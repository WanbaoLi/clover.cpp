"""Clover HMC command line. Each completed trajectory atomically commits the selected gauge and RNG together."""
import argparse
from contextlib import closing
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import subprocess
import tempfile
from time import perf_counter
import numpy as np
from engine import Dynamics,HERE,PROFILE,require,validate_cg_budget,worker_grid
from gauge_io import chroma_restart,read_lime,write_lime
from qdp_stream import QDPNoise,pack_seed

FORMAT=1
def fingerprint(u):
    return hashlib.sha256(np.ascontiguousarray(u,dtype="<c16").view(np.uint8)).hexdigest()
def checked_gauge(u,dims,tol):
    require(u.shape==(math.prod(dims),4,3,3) and u.size>0,"gauge shape does not match geometry")
    err=det=0.0
    identity=np.eye(3)
    det_roundoff=32*np.finfo(np.result_type(u.real.dtype,np.float32)).eps
    for start in range(0,len(u),2048):
        links=u[start:start+2048]
        require(np.isfinite(links).all(),"nonfinite gauge")
        gram=links.conj().swapaxes(-1,-2)@links
        block_err=float(np.max(np.abs(gram-identity)))
        a,b,c=links[...,0,:],links[...,1,:],links[...,2,:]
        determinants=(a[...,0]*(b[...,1]*c[...,2]-b[...,2]*c[...,1])
                      + a[...,1]*(b[...,2]*c[...,0]-b[...,0]*c[...,2])
                      + a[...,2]*(b[...,0]*c[...,1]-b[...,1]*c[...,0]))
        block_det=float(np.max(np.abs(determinants-1)))
        require(math.isfinite(block_err) and math.isfinite(block_det),"nonfinite SU(3) error")
        # Keep the original LU-based decision for fields at the acceptance boundary.
        if abs(block_det-tol)<=det_roundoff:
            block_det=float(np.max(np.abs(np.linalg.det(links)-1)))
        err=max(err,block_err)
        det=max(det,block_det)
    require(max(err,det)<tol,f"SU(3) check failed: unitarity={err:.4g}, det={det:.4g}")
    return dict(max_unitarity_error=err,max_determinant_error=det)
def load_input(path,dims):
    shape=(math.prod(dims),4,3,3)
    if path is None:return np.broadcast_to(np.eye(3,dtype="=c16"),shape).copy()
    if path.suffix.lower() in (".ildg",".lime",".qio"):u=read_lime(path,dims)
    elif path.suffix==".npz":
        with np.load(path,allow_pickle=False) as data:u=data["gauge"].astype("=c16")
    elif path.suffix==".npy":u=np.load(path,allow_pickle=False).astype("=c16")
    elif path.suffix==".bin":u=np.fromfile(path,dtype="<c16").reshape(shape).astype("=c16",copy=False)
    else:raise ValueError("input suffix must be .npy, .npz, .bin, .qio, .ildg or .lime")
    checked_gauge(u,dims,5e-6)
    return np.ascontiguousarray(u)
def statepath(folder,number):return folder/f"trajectory-{number:08d}.npz"
def read_state(path):
    with np.load(path,allow_pickle=False) as f:
        u=f["gauge"].astype("=c16"); state=json.loads(str(f["metadata"]))
    require(state["format"]==FORMAT and state["profile"]==PROFILE["source_xml_sha256"],"incompatible checkpoint")
    require(fingerprint(u)==state["gauge_sha256"],"checkpoint gauge SHA256 mismatch")
    require(path.name==statepath(path.parent,state["number"]).name,"checkpoint number mismatch")
    checked_gauge(u,state["dims"],5e-6)
    return u,state
def latest(folder):
    files=[p for p in folder.glob("trajectory-*.npz") if re.fullmatch(r"trajectory-\d{8,}\.npz",p.name)]
    if not files:raise ValueError("no completed checkpoint")
    return max(files,key=lambda p:int(p.stem.split("-")[-1]))
def sync_directory(folder):
    fd=os.open(folder,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)
def commit(folder,u,state):
    state["gauge_sha256"]=fingerprint(u)
    dest=statepath(folder,state["number"]);require(not dest.exists(),"checkpoint already exists")
    with tempfile.NamedTemporaryFile(dir=folder,prefix=".checkpoint-",delete=False) as f:
        tmp=Path(f.name)
        try:
            np.savez(f,gauge=np.asarray(u,dtype="<c16"),metadata=json.dumps(state,allow_nan=False,separators=(",",":")))
            f.flush();os.fsync(f.fileno())
        except BaseException:
            tmp.unlink(missing_ok=True);raise
    try:os.replace(tmp,dest);sync_directory(folder)
    finally:tmp.unlink(missing_ok=True)
def validate_light_solver(solver,dims):
    grid=worker_grid(dims)
    if solver=="mg":
        dims=[n//p for n,p in zip(dims,grid)]
        result=subprocess.run([str(HERE/"clover-hmc-light"),"--check-mg",*map(str,dims)],
                              stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        require(result.returncode==0,"native QXS MG check failed: "+result.stdout.strip())
def fresh(args):
    if args.chroma_restart is not None:
        require(args.gauge is None,"the Chroma checkpoint already supplies the initial gauge")
        dims,number,words=chroma_restart(args.chroma_restart)
        require(args.dims is None or tuple(args.dims)==dims,"Chroma restart geometry mismatch")
        settings=PROFILE["default_integration"]
        require((args.tau,args.outer_steps,args.inner_steps)==(settings["tau"],settings["outer_steps"],settings["inner_steps"]),
                "Chroma restart requires the embedded trajectory length and integration steps")
        gauge=args.chroma_restart
        pack_seed(words)
        rng_state={"bit_generator":"QDP47","words":words}
    else:
        dims=tuple(args.dims or [4]*4)
        number=0
        gauge=args.gauge
        require(args.seed>=0,"seed must be nonnegative")
        rng_state=np.random.PCG64(args.seed).state
    require(all(n>=2 and n%2==0 for n in dims),"all four extents must be even and at least 2")
    require(math.prod(dims)<=2147483647//72,"volume too large for this Bridge++ configuration")
    require(args.strange_solver in ("shared","independent"),"unsupported strange solver environment")
    validate_light_solver(args.light_solver,dims)
    cfg=dict(tau=args.tau,outer_steps=args.outer_steps,inner_steps=args.inner_steps,
             light_solver=args.light_solver,strange_solver=args.strange_solver)
    if args.chroma_restart is not None:cfg["delta_h_method"]="chroma"
    require(0<cfg["tau"] and math.isfinite(cfg["tau"]) and min(cfg["outer_steps"],cfg["inner_steps"])>0,"invalid integration parameters")
    require(0<=args.chrono_depth<=7,"chrono depth must be 0..7")
    validate_cg_budget(args.light_solver,args.light_cg_max_iterations)
    if args.light_cg_max_iterations is not None:
        cfg["light_cg_max_iterations"]=args.light_cg_max_iterations
    if args.chrono_depth:cfg["chrono_depth"]=args.chrono_depth
    u=load_input(gauge,dims);state=dict(format=FORMAT,profile=PROFILE["source_xml_sha256"],
            dims=dims,config=cfg,number=number,seed=args.seed,rng=rng_state,step=None)
    commit(args.output,u,state)
    return u,state
def chroma_kinetic_energy(momentum):
    # Subtract the kinetic zero point on each link before the lattice reduction.
    # Bound scratch storage independently of the lattice volume.
    totals=[[] for _ in range(4)]
    for begin in range(0,len(momentum),16384):
        block=momentum[begin:begin+16384]
        for mu,partial in enumerate(totals):
            links=block[:,mu]
            link_norms=np.sum(links.real**2+links.imag**2,axis=(1,2))
            partial.append(float(np.sum(link_norms-4.)))
    value=0.
    for partial in totals:
        value+=math.fsum(partial)
    return value

def chroma_action_offsets(volume):
    # Native gauge actions subtract the cold field. The clover blocks use
    # hopping normalization. Restore the reference monomials' constant terms
    # before summing the potential in the order of the production XML.
    ref=PROFILE["chroma_action_conventions"]
    c_rect=-1.0/(20*ref["gauge_u0"]*ref["gauge_u0"])
    return {"lw_tree_gauge":-ref["gauge_beta"]*volume*(6+12*c_rect),
            "logdet_1flav":-12*volume*math.log(4+ref["strange_mass"]),
            "logdet_2flav":-24*volume*math.log(4+PROFILE["masses"]["has4"][0])}

def chroma_energy_difference(start,end,offsets):
    def potential(state):
        value=0.
        for tag in PROFILE["hamiltonian"]:
            value+=state["action"][tag]+offsets.get(tag,0.)
        return value
    initial=potential(start)
    final=potential(end)
    kinetic=end.get("chroma_kinetic",end["kinetic"])-start.get("chroma_kinetic",start["kinetic"])
    potential_change=final-initial
    return dict(initial_potential=initial,final_potential=final,delta_kinetic=kinetic,
                delta_potential=potential_change,delta_h=kinetic+potential_change)

def step(dyn,u,state,rng,timings=False):
    # A checkpoint fixes the arithmetic used for all subsequent accept decisions.
    method=state["config"].get("delta_h_method","total")
    require(method in ("chroma","total"),"unknown Hamiltonian difference method")
    begin=perf_counter()
    p,refresh=dyn.refresh(u,rng)
    refreshed=perf_counter()
    start=dyn.energy(u,p)
    if method=="chroma":start["chroma_kinetic"]=chroma_kinetic_energy(p)
    energy_ready=perf_counter()
    candidate,q=dyn.evolve(u,p,state["config"])
    evolved=perf_counter()
    end=dyn.energy(candidate,q)
    if method=="chroma":end["chroma_kinetic"]=chroma_kinetic_energy(q)
    energy_finished=perf_counter()
    group=checked_gauge(candidate,state["dims"],1e-11)
    reference=None
    if method=="chroma":
        reference=chroma_energy_difference(start,end,chroma_action_offsets(math.prod(state["dims"])))
        delta=reference["delta_h"]
    else:
        delta=end["total"]-start["total"]
    require(math.isfinite(delta),"nonfinite Hamiltonian change")
    probability=1. if delta<=0 else math.exp(-delta)
    # Match Chroma's conditional RNG consumption for a negative Hamiltonian change.
    uniform=None if delta<=0 else float(rng.random())
    accepted=uniform is None or uniform<=probability
    record=dict(accepted=accepted,delta_h=delta,probability=probability,uniform=uniform,
                initial=start,final=end,heatbath=refresh,proposal_sha256=fingerprint(candidate),
                initial_sha256=fingerprint(u),proposal_su3=group)
    if reference is not None:record["chroma_energy_difference"]=reference
    grid=worker_grid(state["dims"])
    if math.prod(grid)>1:record["native_mpi_grid"]=grid
    if dyn.mg_work:record["light_mg"]=dict(dyn.mg_work)
    if dyn.strange_work:record["strange_krylov"]=dict(dyn.strange_work)
    if dyn.chrono_depth:record["md_chronological_predictor"]=dict(dyn.chrono_stats)
    if timings:
        finished=perf_counter()
        record["wall_seconds"]={
            "refresh":refreshed-begin,
            "initial_energy":energy_ready-refreshed,
            "molecular_dynamics":evolved-energy_ready,
            "final_energy":energy_finished-evolved,
            "gauge_check_metropolis_and_metadata":finished-energy_finished,
            "step_before_checkpoint":finished-begin,
        }
    return candidate if accepted else u,record
def main():
    parser=argparse.ArgumentParser(description="Bridge++ cfg3000 Clover/Stout 2+1 HMC with MPI and OpenMP CPU workers")
    sub=parser.add_subparsers(dest="command",required=True)
    start=sub.add_parser("run",help="create a fresh chain in a new or empty directory")
    resume=sub.add_parser("resume",help="continue from the newest committed trajectory")
    info=sub.add_parser("inspect",help="verify and print the newest checkpoint")
    export=sub.add_parser("export",help="write a double-precision SciDAC/ILDG gauge from a checkpoint")
    for p in (start,resume,info,export):p.add_argument("--output",required=True,type=Path,help="chain directory")
    for p in (start,resume):
        p.add_argument("--trajectories",type=int,default=1,help="additional trajectories")
        p.add_argument("--timings",action="store_true",help="save wall times for this invocation's trajectory phases in each new checkpoint")
    for p in (info,export):p.add_argument("--number",type=int,help="specific trajectory index; default is the newest")
    export.add_argument("--lime",type=Path,required=True)
    stream=start.add_mutually_exclusive_group(required=True)
    stream.add_argument("--seed",type=int,help="integer seed for a new PCG64 chain")
    stream.add_argument("--chroma-restart",type=Path,help="continue an archived cfg3000-profile Chroma LIME checkpoint with its QDP++ seed")
    start.add_argument("--dims",type=int,nargs=4)
    start.add_argument("--gauge",type=Path,help="initial gauge; default is the identity")
    default=PROFILE["default_integration"]
    start.add_argument("--tau",type=float,default=default["tau"])
    start.add_argument("--light-solver",choices=("cg","mg"),default="cg",help="light symmetric-Schur normal-equation solver")
    start.add_argument("--light-cg-max-iterations",type=int,
                       help="override the CPU CG iteration cap for the six light factors; residual targets stay at the profile values")
    start.add_argument("--strange-solver",choices=("shared","independent"),
                       default=os.environ.get("CLOVER_HMC_STRANGE_SOLVER","shared"),help="strange shifted-system strategy")
    start.add_argument("--chrono-depth",type=int,default=0,help="MD light normal-equation minimum-residual history (0..7)")
    start.add_argument("--outer-steps",type=int,default=default["outer_steps"])
    start.add_argument("--inner-steps",type=int,default=default["inner_steps"])
    args=parser.parse_args();folder=args.output.resolve();args.output=folder
    if args.command in ("export","inspect"):
        path=statepath(folder,args.number) if args.number is not None else latest(folder)
        u,state=read_state(path)
        if args.command=="inspect":print(json.dumps(state,indent=2));return
        require(not args.lime.exists(),"export destination already exists")
        with tempfile.NamedTemporaryFile(dir=args.lime.parent,prefix=".export-",delete=False) as f:tmp=Path(f.name)
        try:
            write_lime(tmp,u,state["dims"],state["number"])
            require(np.array_equal(read_lime(tmp,state["dims"]),u),"export readback mismatch")
            os.link(tmp,args.lime);sync_directory(args.lime.parent)
        finally:tmp.unlink(missing_ok=True)
        print(state["gauge_sha256"]);return
    require(args.trajectories>0,"trajectory count must be positive")
    for k in ("OMPI_COMM_WORLD_SIZE","PMI_SIZE","MV2_COMM_WORLD_SIZE"):
        require(int(os.getenv(k,"1"))==1,"this executable requires a single MPI rank")
    if args.command=="run":folder.mkdir(parents=True,exist_ok=True)
    require(folder.is_dir(),"chain directory does not exist")
    with (folder/".clover-hmc.lock").open("a") as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc:raise ValueError("another writer holds this chain") from exc
        if args.command=="run":
            require(not any(p.name!=".clover-hmc.lock" for p in folder.iterdir()),"run requires a new or empty chain directory")
            u,state=fresh(args)
        else:
            u,state=read_state(latest(folder))
            validate_light_solver(state["config"].get("light_solver","cg"),state["dims"])
        if state["rng"]["bit_generator"]=="QDP47":
            rng=QDPNoise(state["dims"],state["rng"]["words"])
            def save_rng():
                return {"bit_generator":"QDP47","words":rng.words}
        else:
            rng=np.random.Generator(np.random.PCG64())
            rng.bit_generator.state=state["rng"]
            def save_rng():
                return rng.bit_generator.state
        cfg=state["config"]
        # Checkpoints predating the shared Krylov solver were produced using independent shifts.
        light=cfg.get("light_solver","cg");strange=cfg.get("strange_solver","independent")
        with tempfile.TemporaryDirectory(prefix="clover-hmc-") as scratch, closing(Dynamics(state["dims"],Path(scratch),cfg.get("chrono_depth",0),light,strange,cfg.get("light_cg_max_iterations"))) as dyn:
            for _ in range(args.trajectories):
                u,rec=step(dyn,u,state,rng,timings=args.timings);state=dict(state,number=state["number"]+1,rng=save_rng(),step=rec)
                commit(folder,u,state)
                print(json.dumps(dict(trajectory=state["number"],delta_h=rec["delta_h"],accepted=rec["accepted"],
                                      gauge_sha256=state["gauge_sha256"])),flush=True)
if __name__=="__main__":
    try:main()
    except (OSError,ValueError,RuntimeError,KeyError) as exc:
        print("clover-hmc:",exc,file=sys.stderr);sys.exit(1)
