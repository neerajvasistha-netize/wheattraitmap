"""Scan one gz chromosome FASTA for SSR primer pairs (in-silico PCR, <=1 mismatch per primer outside the last 8 bases)
and for probe sequences (16-mer seeds, clusters within 3 kb)."""
import json,sys,gzip,urllib.request,io,time
import ahocorasick
Q=json.load(open('queries.json'))
comp=str.maketrans('ACGTN','TGCAN');rc=lambda s:s.translate(comp)[::-1]
A=ahocorasick.Automaton();W={}
def add(seed,info):
    W.setdefault(seed,[]).append(info)
for m,(f,r) in Q['primers'].items():
    for lab,p in (('F',f),('R',r)):
        for strand,s in (('+',p),('-',rc(p))):
            # 3' 12-mer seed: on + strand the 3' end is the right end; on - strand it is the left end
            seed=s[-12:] if strand=='+' else s[:12]
            add(seed,('P',m,lab,strand,s))
K=16
for m,seqs in Q['probes'].items():
    for si,s in enumerate(seqs):
        for strand,t in (('+',s),('-',rc(s))):
            for i in range(0,len(t)-K+1,4):
                sd=t[i:i+K]
                if 'N' not in sd: add(sd,('S',m,si,strand,i))
for k in W: A.add_word(k,k)
A.make_automaton()
def scan(name,seq,base,P,S):
    for end,key in A.iter(seq):
        st=end-len(key)+1
        for info in W[key]:
            if info[0]=='P':
                _,m,lab,strand,s=info
                p=st-(len(s)-12) if strand=='+' else st
                if p<0 or p+len(s)>len(seq): continue
                mm=sum(1 for a,b in zip(seq[p:p+len(s)],s) if a!=b)
                if mm<=2: P.append((m,lab,strand,base+p,len(s),mm))
            else:
                S.append((info[1],base+st))
def stream(url):
    for att in range(4):
        try:
            r=urllib.request.urlopen(url,timeout=600);g=gzip.GzipFile(fileobj=r)
            name=None;buf=[];blen=0;base=0;P=[];S=[];OV=2000;CH=20_000_000;res={}
            def flush(final=False):
                pass
            for line in io.TextIOWrapper(g):
                if line[0]=='>':
                    if name: scan(name,''.join(buf).upper(),base,P,S);res[name]=(P,S)
                    name=line[1:].split()[0];buf=[];blen=0;base=0;P=[];S=[]
                else:
                    s=line.rstrip();buf.append(s);blen+=len(s)
                    if blen>=CH:
                        sq=''.join(buf).upper();scan(name,sq,base,P,S);base+=len(sq)-OV;buf=[sq[-OV:]];blen=OV
            if name: scan(name,''.join(buf).upper(),base,P,S);res[name]=(P,S)
            return res
        except Exception as e: print('retry',e,flush=True);time.sleep(10)
    print('FAILED',url);return {}
u=open(sys.argv[1]).read().strip();t0=time.time();res=stream(u)
out={}
for name,(P,S) in res.items():
    out[name]={'P':sorted(set(P)),'S':sorted(set(S))}
json.dump(out,open(sys.argv[2],'w'));print('DONE',u.split('/')[-1],{k:(len(v['P']),len(v['S'])) for k,v in out.items()},'%.0fs'%(time.time()-t0))
