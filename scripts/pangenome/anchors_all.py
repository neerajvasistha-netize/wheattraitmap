"""Genome-wide anchors: every Chinese Spring v1.1 HC gene (CDS from RefSeq v2.1 via the v1.1->v2.1 map)
matched to a cultivar's CDS with minimap2 (mappy, asm10). Same rules as pangenome_anchors.py:
identity >= 0.95, coverage of the CS CDS >= 0.6, best-scoring hit, same-chromosome hit preferred."""
import gzip,json,re,sys,os,time,threading,queue
import mappy
MIN_ID,MIN_COV=0.95,0.6
HDR=re.compile(r"(?:chromosome|primary_assembly):[^:]+:([^:]+):(\d+):(\d+)");GENE=re.compile(r"gene:(\S+)")
def chrname(c):
    m=re.search(r"([1-7][ABD])$",c.replace("chr","").replace("Chr",""));return m.group(1) if m else None
def fasta(p):
    op=gzip.open if p.endswith('.gz') else open;n=None;s=[]
    with op(p,'rt') as fh:
        for l in fh:
            if l[0]=='>':
                if n: yield n,''.join(s)
                n=l[1:].strip();s=[]
            else: s.append(l.strip())
    if n: yield n,''.join(s)
ref=json.load(open('/mnt/user-data/uploads/WheatNHR_liftover/output/_ref_cache.json'))['ref']
qlen={n.split()[0]:len(s) for n,s in fasta('cs_query_all.fa')}
idx=mappy.Aligner('cs_query_all.fa',preset='asm10',n_threads=2)
def run(cds,out,limit=None):
    t0=time.time();best={};lock=threading.Lock();Q=queue.Queue(4000);n=[0]
    def worker():
        buf=mappy.ThreadBuffer()
        while True:
            it=Q.get()
            if it is None: break
            name,seq=it;m=HDR.search(name);gm=GENE.search(name);gene=gm.group(1) if gm else name.split()[0]
            pos=(chrname(m.group(1)),(int(m.group(2))+int(m.group(3)))/2e6) if m else None
            for h in idx.map(seq,buf=buf):
                ident=h.mlen/max(h.blen,1);cov=(h.r_en-h.r_st)/max(qlen[h.ctg],1)
                if ident<MIN_ID or cov<MIN_COV: break
                same=pos is not None and pos[0]==ref[h.ctg][0];key=(same,h.mlen)
                with lock:
                    if h.ctg not in best or key>best[h.ctg][0]: best[h.ctg]=(key,gene,pos,round(ident,4))
                break
    th=[threading.Thread(target=worker) for _ in range(2)];[t.start() for t in th]
    for name,seq in fasta(cds):
        n[0]+=1;Q.put((name,seq))
        if limit and n[0]>=limit: break
    for _ in th: Q.put(None)
    [t.join() for t in th]
    k=0
    with open(out,'w') as fh:
        fh.write("v11_gene\tchr\tmid_Mb\tidentity\tcultivar_gene\n")
        for g,(key,gene,pos,ident) in best.items():
            if pos and pos[0]: fh.write("%s\t%s\t%.4f\t%s\t%s\n"%(g,pos[0],pos[1],ident,gene));k+=1
    print(os.path.basename(cds),n[0],'CDS scanned;',k,'CS genes anchored; %.0f s'%(time.time()-t0),flush=True)
def all_cultivars():
    import glob
    os.makedirs('anchors_full',exist_ok=True)
    for cds in sorted(glob.glob('/mnt/user-data/uploads/WheatNHR_liftover/pangenome/*.cds.all.fa.gz')):
        asm=os.path.basename(cds).split('.cds.all')[0]
        run(cds,'anchors_full/anchors_%s.tsv'%asm)

if __name__=="__main__" and len(sys.argv)>1 and sys.argv[1]=="all": all_cultivars()
elif __name__=="__main__":
    run(sys.argv[1],sys.argv[2],int(sys.argv[3]) if len(sys.argv)>3 else None)
