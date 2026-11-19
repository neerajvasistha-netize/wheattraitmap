"""Genome-wide in-silico PCR (one gz chromosome per run, chunked, low memory).
Each primer is located by its 3' 15 bases (exact); the partner must lie 50-1500 bp downstream in the opposite orientation."""
import json,sys,gzip,urllib.request,io,time,ahocorasick
P=json.load(open(sys.argv[3]))
comp=str.maketrans('ACGTN','TGCAN');rc=lambda s:s.translate(comp)[::-1]
A=ahocorasick.Automaton();W={}
for m,v in P.items():
    f,r=v[0].upper(),v[1].upper()
    for lab,p in (('F',f),('R',r)):
        W.setdefault(p[-15:],[]).append((m,lab,'+',p))      # primer on + strand, 3' end at right
        W.setdefault(rc(p)[:15],[]).append((m,lab,'-',p))   # primer binding the other strand
for k in W: A.add_word(k,k)
A.make_automaton()
def mm(a,b): return sum(1 for x,y in zip(a,b) if x!=y)
def scan(name,seq,base,out):
    hits={}
    for end,key in A.iter(seq):
        st=end-14
        for m,lab,strand,p in W[key]:
            s0=st-(len(p)-15) if strand=='+' else st
            hits.setdefault(m,[]).append((s0,lab,strand,p))
    for m,hs in hits.items():
        plus=[h for h in hs if h[2]=='+'];minus=[h for h in hs if h[2]=='-']
        for s0,lab,_,p in plus:
            for s1,lab2,_,q in minus:
                e=s1+len(q)
                if lab2!=lab and 50<=e-s0<=1500 and s0>=0 and e<=len(seq):
                    out.append((m,name,base+s0,base+e,e-s0,mm(seq[s0:s0+len(p)],p),mm(seq[s1:e],rc(q))))
def run(url):
    for att in range(4):
        try:
            g=gzip.GzipFile(fileobj=urllib.request.urlopen(url,timeout=600));name=None;buf=[];bl=0;base=0;out=[];OV=3000;CH=20_000_000
            for line in io.TextIOWrapper(g):
                if line[0]=='>':
                    if name: scan(name,''.join(buf).upper(),base,out)
                    name=line[1:].split()[0];buf=[];bl=0;base=0
                else:
                    s=line.rstrip();buf.append(s);bl+=len(s)
                    if bl>=CH:
                        sq=''.join(buf).upper();scan(name,sq,base,out);base+=len(sq)-OV;buf=[sq[-OV:]];bl=OV
            if name: scan(name,''.join(buf).upper(),base,out)
            return sorted(set(out))
        except Exception as e: print('retry',e,flush=True);time.sleep(10)
    print('FAILED');return []
o=run(open(sys.argv[1]).read().strip());json.dump(o,open(sys.argv[2],'w'));print('DONE',len(o))
