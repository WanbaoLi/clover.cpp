"""Bridge++ process adapter and the two-level force-gradient update."""
import contextlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
PROFILE=json.loads((HERE/"profile.json").read_text())
GROUPS=PROFILE["groups"]
MASSES=PROFILE["masses"]
SEQUENCE=PROFILE["hamiltonian"]
SOLVERS=PROFILE["solvers"]

def project(p):
    a=(p-p.conj().swapaxes(-1,-2))/2
    return a-np.trace(a,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)/3

def require(ok,message):
    if not ok: raise ValueError(message)

def worker_grid(dims):
    value=os.environ.get("CLOVER_HMC_GRID","1 1 1 1")
    try:
        grid=tuple(int(x) for x in value.split())
    except ValueError as exc:
        raise ValueError("CLOVER_HMC_GRID requires four positive integers in x y z t order") from exc
    require(len(grid)==4 and all(p>0 and n%(2*p)==0 for n,p in zip(dims,grid)),
            "CLOVER_HMC_GRID requires four divisors with even local extents")
    return grid

class Native:
    def __init__(self,kind,dims,scratch,native_env=None):
        self.kind=kind;self.folder=scratch/kind
        self.folder.mkdir()
        self.shape=(math.prod(dims),4,3,3)
        threads="1"
        if kind in ("light","strange"):
            threads=os.environ.get("CLOVER_HMC_FERMION_THREADS","1")
        elif kind=="logdet":
            threads=os.environ.get("CLOVER_HMC_LOGDET_THREADS","1")
        elif kind=="gauge":
            threads=os.environ.get("CLOVER_HMC_GAUGE_THREADS","1")
        env=dict(os.environ,OMP_NUM_THREADS=threads)
        if native_env:env.update(native_env)
        setting=os.environ.get("CLOVER_HMC_MPIEXEC")
        if setting is None:
            # Some MPI runtimes cannot initialize concurrent singleton worlds
            # without their process manager. Honor an explicit empty setting
            # for sites using a working direct singleton launch.
            mpiexec=shutil.which("mpiexec")
            launcher=[mpiexec] if mpiexec else []
        else:
            launcher=shlex.split(setting)
        grid=worker_grid(dims)
        ranks=math.prod(grid)
        if ranks>1 and not launcher:launcher=["mpiexec"]
        wrapper=shlex.split(os.environ.get("CLOVER_HMC_RANK_WRAPPER",""))
        command=[*launcher,"-n",str(ranks)] if launcher else []
        command.extend([*wrapper,str(HERE/("clover-hmc-"+kind)),*map(str,dims)])
        if ranks>1:command.extend(map(str,grid))
        self.proc=subprocess.Popen(command,cwd=self.folder,
                                   stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                   bufsize=1,text=True,env=env)
    def close(self):
        try:
            # MPI launchers may keep rank zero's stdin open after the controller closes its pipe.
            self.proc.stdin.write("-1 0 0 0 0\n")
            self.proc.stdin.flush()
        except BrokenPipeError:
            pass
        try:
            self.proc.stdin.close()
        except BrokenPipeError:
            pass
        diagnostic=self.proc.stdout.read()
        status=self.proc.wait()
        self.proc.stdout.close()
        if status and sys.exc_info()[0] is None:
            raise RuntimeError(f"{self.kind} worker exited with status {status}: {diagnostic[-1200:]}")
    def call(self,u,mode=0,mass=0.,den=9.,tol=1e-12,limit=200):
        np.ascontiguousarray(u,dtype="=c16").tofile(self.folder/"u.bin")
        try:
            self.proc.stdin.write(f"{mode} {mass:.17g} {den:.17g} {tol:.17g} {limit}\n")
            self.proc.stdin.flush()
        except BrokenPipeError as exc:
            raise RuntimeError(self.kind+" worker exited before request") from exc
        tail=[]; errors=[]; rows=[]
        for line in self.proc.stdout:
            if line.strip()=="REQUEST_DONE":return rows
            tail=(tail+[line.rstrip()])[-10:]
            if any(text in line.lower() for text in ("what():","exception","failed","error in","fatal")):
                errors=(errors+[line.rstrip()])[-8:]
            fields=line.split()
            if fields and fields[0] in ("LOGDET","RATIO","SHIFT","ACTION","REFRESH","CHRONO","KRYLOV","MG_NORMAL"):
                rows.append(fields)
        raise RuntimeError(self.kind+" worker failed: "+"\n".join(dict.fromkeys(errors+tail)))
    def gauge(self,name):
        return np.fromfile(self.folder/name,dtype="=c16").reshape(self.shape)

def validate_cg_budget(solver,limit):
    if limit is not None:
        require(solver=="cg","a light CG budget requires the CG light solver")
        require(type(limit) is int and 0<limit<=2147483647,"light CG iteration limit must be a positive native integer")

class Dynamics:
    def __init__(self,dims,scratch,chrono_depth=0,light_solver="cg",strange_solver="shared",light_cg_max_iterations=None):
        require(0<=chrono_depth<=7,"chrono depth must be 0..7")
        require(light_solver in ("cg","mg"),"unsupported light solver")
        require(strange_solver in ("shared","independent"),"unsupported strange solver")
        validate_cg_budget(light_solver,light_cg_max_iterations)
        self.light_cg_max_iterations=light_cg_max_iterations
        self.chrono_depth=chrono_depth;self.chrono_stats={};self.strange_work={};self.mg_work={}
        self.dims=dims;self.vol=math.prod(dims);self.shape=(self.vol,4,3,3)
        self.work={};self.stack=contextlib.ExitStack()
        options={"light":{"CLOVER_HMC_LIGHT_SOLVER":light_solver},
                 "strange":{"CLOVER_HMC_STRANGE_SOLVER":strange_solver}}
        try:
            for kind in ("light","strange","logdet","gauge"):
                w=Native(kind,dims,scratch,options.get(kind));self.stack.callback(w.close);self.work[kind]=w
        except BaseException:
            self.stack.close();raise
        self.phi={};self.cache={}
    def close(self):
        self.stack.close()
    def parameters(self,tag):
        m,den=MASSES.get(tag,(0,None));target,limit=SOLVERS[tag]
        if tag in MASSES and self.light_cg_max_iterations is not None:
            limit=self.light_cg_max_iterations
        return dict(mass=m,den=9. if den is None else den,tol=target,limit=limit)
    def residuals(self,tag,rows,refresh,force_only=False):
        tol,_=SOLVERS[tag]
        if tag in MASSES:
            row=next(r for r in rows if r[0]==("REFRESH" if refresh else "RATIO"))
            values=[float(row[1] if refresh else row[4])]
        else:
            values=[float(r[5]) for r in rows if r[0]=="SHIFT"]
            expected={"heat":15} if refresh else ({"md":12} if force_only else {"action":15,"md":12})
            for system,count in expected.items():
                require(sum(r[:2]==["SHIFT",system] for r in rows)==count,"incomplete strange "+system+" shifted solves")
            require(len(values)==sum(expected.values()),"unexpected strange shifted system")
        for residual in values:
            require(math.isfinite(residual) and 0<=residual<=tol*1.5,
                    f"solver true residual exceeded for {tag}: measured {residual:.17g}, target {tol:.17g}")
        return max(values)
    def record_mg_work(self,tag,rows):
        for r in rows:
            if r[0]!="MG_NORMAL":continue
            st=self.mg_work.setdefault(tag,{}).setdefault(r[1],dict(solves=0,mg_count=0,max_post_mg_schur_residual=0.))
            st["solves"]+=1;st["mg_count"]+=int(r[2])
            st["max_post_mg_schur_residual"]=max(st["max_post_mg_schur_residual"],float(r[3]))
    def record_strange_work(self,rows):
        for r in rows:
            if r[0]!="KRYLOV": continue
            st=self.strange_work.setdefault(r[1],dict(mode=r[2],calls=0,poles_per_call=int(r[3]),
                           iteration_normal_applications=0,polish_normal_applications=0,verified_poles=0))
            st["calls"]+=1;st["iteration_normal_applications"]+=int(r[4])
            st["polish_normal_applications"]+=int(r[5]);st["verified_poles"]+=int(r[6])
    def refresh(self,u,rng):
        self.cache.clear();self.phi.clear();self.chrono_stats.clear();self.strange_work.clear();self.mg_work.clear()
        g=self.work["gauge"]
        raw=rng.standard_normal(self.shape)+1j*rng.standard_normal(self.shape)
        g.call(raw,mode=3);p=g.gauge("moved.bin")
        # Array site order is x-fast. All odd sites have 12 independent complex variates.
        coords=np.indices(tuple(reversed(self.dims)),dtype=np.int32)
        odd=(coords.sum(axis=0).ravel()%2)==1
        result={}
        for tag in SEQUENCE:
            if tag not in MASSES and tag!="rat_3strange": continue
            w=self.work["light" if tag in MASSES else "strange"]
            eta=np.zeros((self.vol,4,3),dtype="=c16"); shape=(self.vol//2,4,3)
            eta[odd]=(rng.standard_normal(shape)+1j*rng.standard_normal(shape))*math.sqrt(.5)
            eta.tofile(w.folder/"phi.bin");rows=w.call(u,mode=2,**self.parameters(tag))
            self.phi[tag]=(w.folder/"refreshed.bin").read_bytes()
            if tag=="rat_3strange":self.record_strange_work(rows)
            else:self.record_mg_work(tag,rows)
            result[tag]=dict(eta_norm2=float(np.vdot(eta,eta).real),max_true_residual=self.residuals(tag,rows,True))
        return p,result
    def field(self,group,u,md=False,force_only=False):
        cache_key=(group,True) if md and self.chrono_depth else group
        hit=self.cache.get(cache_key)
        if hit is not None and hit[0] is u:return hit[1:]
        if force_only:
            # Keep only full evaluations in the cache for future energy queries.
            self.cache.pop(cache_key,None)
        terms={};residuals={}
        total=None
        for tag in GROUPS[group]:
            if tag.startswith("logdet"):
                w=self.work["logdet"]
                if tag not in terms:
                    rows=w.call(u)
                    for n in (1,2):
                        key=f"logdet_{n}flav"; vals=[float(r[5]) for r in rows if r[:2]==["LOGDET",str(n)] and r[3]=="0"]
                        require(len(vals)==2,"incomplete Clover determinant");terms[key]=math.fsum(vals)
                nf={"logdet_1flav":1,"logdet_2flav":2}[tag]
                force=w.gauge(f"nf{nf}_cb0.bin")
                force+=w.gauge(f"nf{nf}_cb1.bin")
            elif tag=="lw_tree_gauge":
                w=self.work["gauge"];rows=w.call(u)
                terms[tag]=float(next(r for r in rows if r[0]=="ACTION")[2]);force=w.gauge("force.bin")
            else:
                w=self.work["light" if tag in MASSES else "strange"]
                (w.folder/"phi.bin").write_bytes(self.phi[tag])
                mode=10+self.chrono_depth if md and self.chrono_depth and tag in MASSES else 0
                strange_force_only=(tag=="rat_3strange" and force_only)
                rows=w.call(u,mode=3 if strange_force_only else mode,**self.parameters(tag))
                residuals[tag]=self.residuals(tag,rows,False,strange_force_only)
                if tag=="rat_3strange":self.record_strange_work(rows)
                else:self.record_mg_work(tag,rows)
                if mode:
                    r=next(x for x in rows if x[0]=="CHRONO")
                    st=self.chrono_stats.setdefault(tag,dict(solves=0,with_history=0,iterations=0,max_true_residual=0.))
                    st["solves"]+=1;st["with_history"]+=int(int(r[2])>0);st["iterations"]+=int(r[5])
                    st["max_true_residual"]=max(st["max_true_residual"],float(r[4]))
                if tag in MASSES:
                    terms[tag]=float(next(r for r in rows if r[0]=="RATIO")[2])
                    if MASSES[tag][1] is None:terms[tag]-=6*self.vol
                elif not strange_force_only:
                    terms[tag]=float(next(r for r in rows if r[0]=="ACTION")[2])
                force=w.gauge("force.bin")
            # Match Python sum's original order and initial zero, retaining one running field.
            if total is None:
                total=force
                np.add(0,total,out=total)
            else:
                total+=force
            del force
        require(np.isfinite(total).all() and all(math.isfinite(x) for x in terms.values()),"nonfinite force or action")
        if not force_only:
            self.cache[cache_key]=(u,total,terms,residuals)
        return total,terms,residuals
    def drift(self,u,p,dt):
        w=self.work["gauge"];p.tofile(w.folder/"p.bin");w.call(u,mode=1,mass=dt)
        return w.gauge("moved.bin")
    def kick(self,u,p,dt,group,force_only=False):
        f,_,_=self.field(group,u,md=True,force_only=force_only)
        return project(p+(.5*dt)*f)
    def middle(self,u,p,h,group):
        force_only=(group=="outer")
        f,_,_=self.field(group,u,md=True,force_only=force_only)
        shift=self.drift(u,project(.5*h*h/24*f),1.)
        return self.kick(shift,p,h*2/3,group,force_only=force_only)
    def inner(self,u,p,length,n):
        h=length/n;p=self.kick(u,p,h/6,"inner")
        for j in range(n):
            u=self.drift(u,p,h/2);p=self.middle(u,p,h,"inner");u=self.drift(u,p,h/2)
            p=self.kick(u,p,h/(6 if j==n-1 else 3),"inner")
        return u,p
    def evolve(self,u,p,config):
        n=config["outer_steps"]; m=config["inner_steps"]; h=config["tau"]/n
        p=self.kick(u,p,h/6,"outer")
        for j in range(n):
            u,p=self.inner(u,p,h/2,m);p=self.middle(u,p,h,"outer");u,p=self.inner(u,p,h/2,m)
            p=self.kick(u,p,h/(6 if j==n-1 else 3),"outer")
        return u,p
    def energy(self,u,p):
        _,outer,ro=self.field("outer",u);_,inner,ri=self.field("inner",u)
        terms=dict(outer,**inner); kinetic=float(np.vdot(p,p).real-16*self.vol)
        return dict(total=kinetic+math.fsum(terms[t] for t in SEQUENCE),kinetic=kinetic,
                    action=terms,true_residuals=dict(ro,**ri))
