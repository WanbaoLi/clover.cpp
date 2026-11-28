"""Write an identity SU(3) field in double-precision big-endian ILDG LIME format."""
import argparse
import math
from pathlib import Path
import struct

cli=argparse.ArgumentParser()
cli.add_argument("file",type=Path)
cli.add_argument("shape",help="X,Y,Z,T")
args=cli.parse_args()
shape=[int(n) for n in args.shape.split(",")]
if len(shape)!=4 or min(shape)<1:
    cli.error("shape requires four positive integers")
def header(name,length,flag):
    return struct.pack(">IHHQ128s",0x456789ab,1,flag,length,name.encode())
matrix=struct.pack(">18d",*(part for i in range(3) for j in range(3) for part in (float(i==j),0.0)))
site=matrix*4
xml=("<?xml version=\"1.0\"?><ildgFormat><version>1.0</version><field>su3gauge</field><precision>64</precision>"+
     "".join(f"<{tag}>{n}</{tag}>" for tag,n in zip(("lx","ly","lz","lt"),shape))+"</ildgFormat>").encode()
with args.file.open("xb") as f:
    f.write(header("ildg-format",len(xml),0x8000));f.write(xml);f.write(bytes(-len(xml)%8))
    count=math.prod(shape)
    f.write(header("ildg-binary-data",count*len(site),0x4000))
    block=site*1024
    for _ in range(count//1024): f.write(block)
    f.write(site*(count%1024))
