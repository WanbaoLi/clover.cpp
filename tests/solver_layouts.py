"""Exercise the installed solver's layout decisions before large field allocation."""
import subprocess
import sys

launcher,count_flag,*solver_command=sys.argv[1:]

def plan(ranks,shape,*extra,expect=None,error=None):
    result=subprocess.run([launcher,count_flag,str(ranks),*solver_command,"--plan-layout","--lattice",shape,*extra],
                          text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=45)
    if error:
        assert result.returncode and error in result.stdout,result.stdout
    else:
        assert result.returncode==0 and expect in result.stdout,result.stdout

plan(2,"4,4,4,4",expect="method=bicgstab")
plan(2,"24,24,24,64",expect="method=mg")
plan(2,"32,32,32,64",expect="method=mg")
plan(2,"48,48,48,96","--method","bicgstab",expect="grid=1,1,1,2")
plan(2,"64,64,64,128",expect="method=mg")
plan(1,"64,64,64,128",error="field indexing")
plan(1,"2,4,4,4","--method","mg",error="no compatible MG grid")
