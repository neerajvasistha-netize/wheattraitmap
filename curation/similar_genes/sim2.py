import mappy,sys,json,time,re
part,nparts=int(sys.argv[1]),int(sys.argv[2]);lim=int(sys.argv[3]) if len(sys.argv)>3 else 0
A={g:mappy.Aligner(f'sub_{g}.fa',preset='asm20',best_n=30) for g in 'ABDU'}
out={};t0=time.time();i=-1
for n,s,_ in mappy.fastx_read('../pan/cs_query_all.fa'):
    i+=1
    if i%nparts!=part: continue
    best={}
    for g,a in A.items():
        for h in a.map(s):
            if h.ctg==n: continue
            idn=h.mlen/max(1,h.blen);qc=(h.q_en-h.q_st)/len(s);tc=(h.r_en-h.r_st)/h.ctg_len
            if idn>=0.80 and qc>=0.5 and tc>=0.5:
                sc=h.mlen
                if h.ctg not in best or best[h.ctg][3]<sc: best[h.ctg]=(round(idn*100,1),round(qc*100),round(tc*100),sc)
    out[n]=[[c,v[0],v[1],v[2]] for c,v in sorted(best.items(),key=lambda x:-x[1][3])[:12]]
    if lim and len(out)>=lim: break
json.dump(out,open(f'sim_{part}.json','w'));print('done',len(out),'%.0fs'%(time.time()-t0))
