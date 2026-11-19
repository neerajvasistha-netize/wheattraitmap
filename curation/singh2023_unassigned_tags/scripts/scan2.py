"""Chunked version of scan.py (low memory): stream gz FASTA per chromosome, scan in 20 Mb windows with overlap."""
import json,sys,gzip,urllib.request,io,time
import ahocorasick
meta=json.load(open('meta.json'))
comp=str.maketrans('ACGT','TGCA');rc=lambda s:s.translate(comp)[::-1]
A=ahocorasick.Automaton();Q={}
for rid,m in meta.items():
    for strand,seq in (('+',m['trimmed']),('-',rc(m['trimmed']))):
        L=len(seq);k=L//4
        for j in range(4): Q.setdefault(seq[j*k:(j+1)*k],[]).append((rid,strand,seq,j*k))
for key in Q: A.add_word(key,key)
A.make_automaton()
OV=200;CH=20_000_000
def scan(name,seq,base,out,seen):
    for end,key in A.iter(seq):
        st=end-len(key)+1
        for rid,strand,full,off in Q[key]:
            p=st-off
            if p<0 or p+len(full)>len(seq):continue
            gp=base+p
            if (rid,strand,gp) in seen:continue
            w=seq[p:p+len(full)];mm=sum(1 for a,b in zip(w,full) if a!=b)
            if mm<=3: seen.add((rid,strand,gp));out.append({'id':rid,'chr':name,'pos':gp+1,'strand':strand,'mm':mm,'len':len(full)})
def stream(url,out):
    for att in range(4):
        try:
            r=urllib.request.urlopen(url,timeout=600);g=gzip.GzipFile(fileobj=r)
            name=None;buf=[];blen=0;base=0;seen=set();tmp=[]
            for line in io.TextIOWrapper(g):
                if line[0]=='>':
                    if name: scan(name,''.join(buf).upper(),base,tmp,seen)
                    name=line[1:].split()[0];buf=[];blen=0;base=0;seen=set()
                else:
                    s=line.rstrip();buf.append(s);blen+=len(s)
                    if blen>=CH:
                        seq=''.join(buf).upper();scan(name,seq,base,tmp,seen)
                        base+=len(seq)-OV;buf=[seq[-OV:]];blen=OV
            if name: scan(name,''.join(buf).upper(),base,tmp,seen)
            out+=tmp;return
        except Exception as e:
            print('retry',url,e,flush=True);time.sleep(10)
    print('FAILED',url,flush=True)
if __name__=='__main__':
    urls=[l.strip() for l in open(sys.argv[1]) if l.strip()];res=[]
    for u in urls:
        t0=time.time();n=len(res);stream(u,res);print(u.split('/')[-1],len(res)-n,'hits %.0fs'%(time.time()-t0),flush=True)
    json.dump(res,open(sys.argv[2],'w'))
    print('DONE',flush=True)
