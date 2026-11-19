"""Relaxed in-silico PCR on one chromosome: both primers located by their 3' 15 bases (<=1 mismatch there),
5' end mismatches tolerated (BARC primers often carry a 5' GCG tail); product 50-1500 bp."""
import json,sys,gzip,urllib.request,io,re
ch=sys.argv[1];url=sys.argv[2];marks=sys.argv[3].split(',')
P=json.load(open(sys.argv[4]))
comp=str.maketrans('ACGT','TGCA');rc=lambda s:s.translate(comp)[::-1]
g=gzip.GzipFile(fileobj=urllib.request.urlopen(url,timeout=600));seq=[]
for line in io.TextIOWrapper(g):
    if line[0]!='>': seq.append(line.rstrip())
S=''.join(seq).upper();del seq
def find(core):
    out=[];i=S.find(core)
    while i>=0: out.append(i);i=S.find(core,i+1)
    return out
res={}
for m in marks:
    f,r=P[m][:2];hits=[]
    # F on + : 3' core = f[-15:] ; R on - strand appears as rc(r), its 3' end is the left end -> rc(r)[:15]
    for a,b,lab in ((f,r,'F/R'),(r,f,'R/F')):
        for i in find(a[-15:]):
            st=i-(len(a)-15)
            win=S[i:i+1500]
            j=win.find(rc(b)[:15])
            if j>=0:
                end=i+j+len(b);prod=end-st
                mmA=sum(1 for x,y in zip(S[st:st+len(a)],a) if x!=y);mmB=sum(1 for x,y in zip(S[end-len(b):end],rc(b)) if x!=y)
                if 50<=prod<=1500: hits.append((st,end,prod,mmA,mmB,lab))
    res[m]=hits
json.dump(res,open(f'relax2_{ch}.json','w'));print(ch,{k:[(h[0],h[2],h[3],h[4]) for h in v] for k,v in res.items()})
