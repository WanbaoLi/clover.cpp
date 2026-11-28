"""Exercise full and even/odd native BiCGStab restarts on a nontrivial gauge field."""
import csv
import hashlib
import math
import struct
import subprocess
import sys
import tempfile
from array import array
from pathlib import Path

launcher,count_flag,*driver=sys.argv[1:]
grid='1,1,1,2'

def gauge(path,variant):
    xml=b'<ildgFormat><field>su3gauge</field><precision>64</precision><lx>4</lx><ly>4</ly><lz>4</lz><lt>4</lt></ildgFormat>'
    body=bytearray()
    for site in range(4**4):
        for mu in range(4):
            phase=.014*(1+(site*17+3*mu+variant)%41)
            z=complex(math.cos(phase),math.sin(phase))
            for value in (z,0,0,0,z.conjugate(),0,0,0,1):
                body.extend(struct.pack('>dd',value.real,value.imag))
    with path.open('wb') as f:
        for name,flags,data in ((b'ildg-format',0x8000,xml),(b'ildg-binary-data',0x4000,body)):
            f.write(struct.pack('>IHHQ128s',0x456789ab,1,flags,len(data),name));f.write(data);f.write(bytes((-len(data))%8))

def run(name,extra,code,needle=None,**changed):
    args={'method':'bicgstab','gauge':str(cfg),'lattice':'4,4,4,4','grid':grid,'source':'point','position':'0,0,0,0',
          'mass':'1','csw':'1.15','rho':'.125','threads':'1','rhs':'0,3','fields':'full','bicg-iters':'2','bicg-cycles':'300',
          'tol':'1e-12','accept':'5e-12','output':str(root/name)}
    args.update(changed)
    command=[launcher,count_flag,'2',*driver]
    for k,v in args.items():command.extend(['--'+k,v])
    result=subprocess.run(command+extra,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=40,text=True)
    assert (result.returncode==0)==(code==0) and (not needle or needle in result.stdout),result.stdout
    return root/name

def table(folder,name):
    with (folder/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))

with tempfile.TemporaryDirectory(prefix='clover-solver-checkpoint-') as temp:
    root=Path(temp);cfg=root/'gauge.lime';gauge(cfg,0)
    baseline=run('baseline',[],0)
    first=run('first',['--checkpoint'],1,**{'bicg-cycles':'1'})
    assert (first/'FAILED').is_file() and not (first/'COMPLETE').exists()
    earlier={int(r['rhs']):r for r in table(first,'columns.tsv')}
    assert set(earlier)=={0,3} and all(r['accepted']=='0' for r in earlier.values())
    resume=['--resume-from',str(first)]
    second=run('second',resume,0,rhs='3,0');assert (second/'COMPLETE').is_file()
    assert [r['rhs'] for r in table(second,'columns.tsv')]==['3','0']
    meta={r['key']:r['value'] for r in table(second,'run.tsv')};assert meta['initial_guess']=='checkpoint'
    local_values=4**4//2*24;size=local_values*8
    for slot,rhs in enumerate((0,3)):
        saved=float(earlier[rhs]['independent_relative']);assert abs(float(meta['resume_initial_relative_rhs_'+str(rhs)])-saved)<1e-14
        manifest=(first/('checkpoint-rhs-%d.ready'%rhs)).read_text().splitlines()
        for record in manifest[-2:]:
            rank,where,gauge_hash,field_hash=record.split('\t')
            shard=first/('checkpoint-rhs-%d-%s.f64'%(rhs,where))
            assert shard.stat().st_size==size and hashlib.sha256(shard.read_bytes()).hexdigest()==field_hash
            native=array('d');native.frombytes(shard.read_bytes())
            published=(first/('field_'+where.replace(',','_')+'.f64')).read_bytes()
            assert published[slot*size:(slot+1)*size]==array('d',(v*.2 for v in native)).tobytes()
    for f in baseline.glob('field_*.f64'):
        a=f.read_bytes();b=(second/f.name).read_bytes();assert len(a)==len(b)==2*size
        assert a[:size]==b[size:] and a[size:]==b[:size]
    middle=run('middle',resume,1,rhs='3',fields='none',**{'bicg-cycles':'1'})
    assert not list(middle.glob('field_*.f64')) and (middle/'checkpoint-rhs-3.ready').is_file()
    last=run('last',['--resume-from',str(middle)],0,rhs='3')
    for f in baseline.glob('field_*.f64'):
        assert (last/f.name).read_bytes()==f.read_bytes()[size:]
    run('wrong-rhs',resume,1,'no committed checkpoint',rhs='3,1')
    for key,value in [('position','1,0,0,0'),('source','wall'),('mass','1.01'),('csw','1.14'),('rho','.12'),('grid','2,1,1,1')]:
        run('wrong-'+key,resume,1,'mismatch',**{key:value})
    gauge(cfg,1);run('wrong-gauge',resume,1,'mismatch');gauge(cfg,0)
    shard=next(first.glob('checkpoint-rhs-3-*.f64'));raw=shard.read_bytes()
    shard.write_bytes(raw[:-1]);run('truncated',resume,1,'shard',rhs='3')
    shard.write_bytes(bytes([raw[0]^1])+raw[1:]);run('corrupt',resume,1,'shard',rhs='3');shard.write_bytes(raw)
    m=first/'checkpoint-rhs-3.ready';lines=m.read_text().splitlines();i=next(i for i,l in enumerate(lines) if l.startswith('residual\t'))
    lines[i]='residual\t0.875';m.write_text('\n'.join(lines)+'\n');run('bad-residual',resume,1,'independently reapplying',rhs='3')

    # An odd-site source exercises the inhomogeneous term in postProp.
    eo={'method':'eo-bicgstab','position':'1,0,0,0','threads':'2'}
    direct=run('eo-direct',[],0,**eo)
    capped=run('eo-capped',['--checkpoint'],1,**eo,**{'bicg-cycles':'1'})
    run('eo-wrong-position',['--resume-from',str(capped)],1,'mismatch',method='eo-bicgstab',rhs='0')
    middle=run('eo-middle',['--resume-from',str(capped)],1,**eo,rhs='3,0',**{'bicg-cycles':'1'})
    final=run('eo-final',['--resume-from',str(middle)],0,**eo,rhs='0,3')
    assert {r['key']:r['value'] for r in table(capped,'run.tsv')}['initial_guess']=='zero'
    for old,new in ((capped,middle),(middle,final)):
        meta={r['key']:r['value'] for r in table(new,'run.tsv')};assert meta['initial_guess']=='checkpoint'
        for slot,r in enumerate(table(old,'columns.tsv')):
            rhs=r['rhs'];saved=float(r['independent_relative']);assert r['accepted']=='0' and saved>1e-8
            assert abs(float(meta['resume_initial_relative_rhs_'+rhs])-saved)<1e-13
            assert abs(float(meta['resume_eo_projected_relative_rhs_'+rhs])-saved)<1e-13
            manifest=(old/('checkpoint-rhs-'+rhs+'.ready')).read_text().splitlines()
            for line in manifest[-2:]:
                _,where,_,digest=line.split('\t');raw=(old/('checkpoint-rhs-'+rhs+'-'+where+'.f64')).read_bytes()
                assert len(raw)==size and hashlib.sha256(raw).hexdigest()==digest
                native=array('d');native.frombytes(raw)
                published=(old/('field_'+where.replace(',','_')+'.f64')).read_bytes()
                assert published[slot*size:(slot+1)*size]==array('d',(v*.2 for v in native)).tobytes()
    for f in direct.glob('field_*.f64'):
        assert f.read_bytes()==(final/f.name).read_bytes()
    # Each solver can consume a verified full field produced by the other one.
    full=run('full-at-odd',[],0,position='1,0,0,0',rhs='0')
    full_first=run('full-at-odd-capped',['--checkpoint'],1,position='1,0,0,0',rhs='0',**{'bicg-cycles':'1'})
    to_eo=run('full-to-eo',['--resume-from',str(full_first)],0,**eo,rhs='0')
    to_full=run('eo-to-full',['--resume-from',str(capped)],0,position='1,0,0,0',rhs='0')
    for target in (to_eo,to_full):
        err=norm=0.
        for f in full.glob('field_*.f64'):
            a=array('d');a.frombytes(f.read_bytes());b=array('d');b.frombytes((target/f.name).read_bytes());assert len(a)==len(b)
            err+=sum((x-y)**2 for x,y in zip(a,b));norm+=sum(x*x for x in a)
        assert math.sqrt(err/norm)<1e-10 and float(table(target,'columns.tsv')[0]['independent_relative'])<5e-12

    # Compare all fields after a cold start, two saved CPU iterations, and a
    # resumed source reordering when the two operators cannot coexist in RAM.
    compact_direct=run('compact-direct',['--compact-eo'],0,**eo)
    compact_first=run('compact-first',['--compact-eo','--checkpoint'],1,**eo,**{'bicg-cycles':'1'})
    compact_mid=run('compact-middle',['--compact-eo','--resume-from',str(compact_first)],1,**eo,
                    rhs='3,0',**{'bicg-cycles':'1'})
    compact_last=run('compact-last',['--compact-eo','--resume-from',str(compact_mid)],0,**eo,rhs='0,3')
    for smaller,baseline in ((compact_direct,direct),(compact_first,capped),(compact_mid,middle),(compact_last,final)):
        for f in baseline.glob('field_*.f64'):
            assert (smaller/f.name).read_bytes()==f.read_bytes()
        for f in baseline.glob('checkpoint-rhs-*.ready'):
            assert (smaller/f.name).read_bytes()==f.read_bytes()
        a=table(smaller,'columns.tsv');b=table(baseline,'columns.tsv')
        assert len(a)==len(b) and all(x[k]==y[k] for x,y in zip(a,b) for k in ('independent_relative','nconv','outer_loops','accepted'))
        mode={r['key']:r['value'] for r in table(smaller,'run.tsv')}
        assert mode['compact_eo']=='1'
        for r in a:assert float(mode['compact_eo_setup_seconds_rhs_'+r['rhs']])>0
    cross=run('full-to-compact',['--compact-eo','--resume-from',str(full_first)],0,**eo,rhs='0')
    for f in to_eo.glob('field_*.f64'):
        assert (cross/f.name).read_bytes()==f.read_bytes()
    meta={r['key']:r['value'] for r in table(cross,'run.tsv')}
    for k in ('resume_initial_relative_rhs_0','resume_eo_projected_relative_rhs_0'):
        expected={r['key']:r['value'] for r in table(to_eo,'run.tsv')}[k]
        assert abs(float(meta[k])-float(expected))<1e-14
    run('compact-invalid',['--compact-eo'],1,'requires --method eo-bicgstab',rhs='0')
