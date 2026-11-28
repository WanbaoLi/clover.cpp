"""Single-file SciDAC/ILDG SU(3) gauge records, x-fast and big-endian on disk."""
from datetime import datetime, timezone
import hashlib
import math
import struct
import xml.etree.ElementTree as ET
import zlib
import numpy as np

HEADER=struct.Struct(">IHHQ128s")
MAGIC=0x456789ab
def checksum(payload,volume):
    raw=memoryview(payload).cast("B"); size=len(raw)//volume;a=b=0
    for site in range(volume):
        c=zlib.crc32(raw[site*size:(site+1)*size]);k=site%29;l=site%31
        a^=((c<<k)|(c>>(32-k)))&0xffffffff
        b^=((c<<l)|(c>>(32-l)))&0xffffffff
    return a,b
def read_lime(path,dims):
    info={}; data=None; expected=None; vol=math.prod(dims)
    with open(path,"rb") as f:
        while True:
            h=f.read(HEADER.size)
            if not h:break
            magic,ver,flags,n,tag=HEADER.unpack(h);tag=tag.rstrip(b"\0").decode("ascii")
            if (magic,ver)!=(MAGIC,1):raise ValueError("invalid LIME header")
            if tag in ("ildg-binary-data","scidac-binary-data"):
                if data is not None or n not in (vol*288,vol*576): raise ValueError("unexpected gauge record size/count")
                data=f.read(n)
                if len(data)!=n:raise ValueError("truncated LIME gauge data")
            elif tag in ("ildg-format","scidac-private-file-xml","scidac-private-record-xml","scidac-checksum"):
                if n>4*1024*1024:raise ValueError("unexpected LIME metadata size")
                tree=ET.fromstring(f.read(n).rstrip(b"\0"))
                doc={node.tag.split("}")[-1]:node.text for node in tree}
                if tag=="scidac-checksum":expected=(int(doc["suma"],16),int(doc["sumb"],16))
                else:info[tag]=doc
            else:f.seek(n,1)
            f.seek((-n)%8,1)
    if data is None:raise ValueError("LIME file contains no gauge field")
    declared=[]
    if "ildg-format" in info:
        field=info["ildg-format"]
        if field["field"]!="su3gauge":raise ValueError("ILDG input is not an SU(3) gauge")
        declared.append(tuple(int(field[k]) for k in ("lx","ly","lz","lt")))
        if int(field["precision"])!=len(data)*8//(vol*72):raise ValueError("inconsistent ILDG precision")
    if "scidac-private-file-xml" in info:
        declared.append(tuple(map(int,info["scidac-private-file-xml"]["dims"].split())))
    if "scidac-private-record-xml" in info:
        q=info["scidac-private-record-xml"]
        representation="D" if len(data)==vol*576 else "F"
        required={"recordtype":"0","colors":"3","datacount":"4","precision":representation,
                  "datatype":"QDP_"+representation+"3_ColorMatrix"}
        if any(q.get(k)!=v for k,v in required.items()):raise ValueError("unsupported SciDAC data type")
    if not declared or any(tuple(dims)!=d for d in declared):raise ValueError("LIME geometry mismatch or missing")
    if expected is not None and checksum(data,vol)!=expected:raise ValueError("SciDAC checksum mismatch")
    return np.frombuffer(data,dtype=">c16" if len(data)==vol*576 else ">c8").reshape(vol,4,3,3).astype("=c16")
def chroma_restart(path):
    """Read the embedded production HMC profile and its checkpoint seed."""
    with open(path,"rb") as f:
        while True:
            header=f.read(HEADER.size)
            if len(header)!=HEADER.size:
                raise ValueError("Chroma restart requires embedded SciDAC record XML")
            magic,version,flags,size,raw_tag=HEADER.unpack(header)
            if (magic,version)!=(MAGIC,1):
                raise ValueError("invalid Chroma restart LIME header")
            if raw_tag.rstrip(b"\0")==b"scidac-record-xml":
                if size>4*1024*1024:
                    raise ValueError("unexpected Chroma restart metadata size")
                doc=ET.fromstring(f.read(size).rstrip(b"\0"))
                break
            f.seek(size+(-size)%8,1)
    traj=doc.find("HMCTrj")
    if traj is None or traj.find("nrow") is None:
        raise ValueError("missing HMCTrj geometry")
    dims=tuple(map(int,traj.findtext("nrow").split()))
    if len(dims)!=4:
        raise ValueError("expected four HMCTrj dimensions")
    traj.remove(traj.find("nrow"))
    canonical=ET.canonicalize(ET.tostring(traj),strip_text=True)
    # Both archived L24 checkpoints contain this exact action/integrator.
    if hashlib.sha256(canonical.encode()).hexdigest()!="6c5085765b6f9f1c2cd02516de03c5cf04be6e6442e9ae9b92e2b927d4d83dc6":
        raise ValueError("unsupported Chroma HMCTrj profile")
    mc=doc.find("MCControl")
    if mc is None or mc.findtext("NWarmUpUpdates")!="0" or mc.findtext("Cfg/reunit")!="false":
        raise ValueError("unsupported Chroma warm-up or gauge reunitarization setting")
    number=int(mc.findtext("StartUpdateNum","-1"))
    if number<0:
        raise ValueError("invalid Chroma starting update")
    words=[int(node.text) for node in mc.findall("RNG/Seed/elem")]
    return dims,number,words

def write_lime(path,u,dims,number):
    v=math.prod(dims);full=np.ascontiguousarray(u,dtype=">c16"); data=memoryview(full).cast("B")
    a,b=checksum(data,v)
    def xml(root,pairs,namespace=""):
        node=ET.Element(root,{"xmlns":namespace} if namespace else {})
        for k,val in pairs:ET.SubElement(node,k).text=str(val)
        return ET.tostring(node,encoding="utf-8",xml_declaration=True)+b"\0"
    records=[
        (0x8000,"scidac-private-file-xml",xml("scidacFile",[("version","1.1"),("spacetime",4),("dims"," ".join(map(str,dims))),("volfmt",0)])),
        (0x4000,"scidac-file-xml",xml("cloverHMC",[("format",1)])),
        (0x8000,"scidac-private-record-xml",xml("scidacRecord",[("version","1.1"),("date",datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S %Y UTC")),("recordtype",0),
            ("datatype","QDP_D3_ColorMatrix"),("precision","D"),("colors",3),("spins",1),("typesize",144),("datacount",4)])),
        (0,"scidac-record-xml",xml("cloverHMC",[("trajectory",number)])),
        (0,"ildg-format",xml("ildgFormat",[("version","1.0"),("field","su3gauge"),("precision",64),*zip(("lx","ly","lz","lt"),dims)],"http://www.lqcd.org/ildg")),
        (0,"ildg-binary-data",data),
        (0x4000,"scidac-checksum",xml("scidacChecksum",[("version","1.0"),("suma",format(a,"x")),("sumb",format(b,"x"))]))]
    with open(path,"wb") as f:
        for flags,tag,payload in records:
            f.write(HEADER.pack(MAGIC,1,flags,len(payload),tag.encode()));f.write(payload);f.write(b"\0"*((-len(payload))%8))
