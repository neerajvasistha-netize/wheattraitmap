import json,csv,glob,os,re,math,sys
sys.path.insert(0,'/home/claude/atlas/liftover')
from collections import Counter,defaultdict
U='/mnt/user-data/uploads'
M=json.load(open('/home/claude/atlas/merged_v13_lifted.json'))
FR='RefSeq v2.1'
ALEN={"1A":598.6,"1B":700.4,"1D":498.6,"2A":787.2,"2B":812.7,"2D":656.4,"3A":754.0,"3B":851.9,"3D":619.5,"4A":754.2,"4B":673.7,"4D":514.9,"5A":713.3,"5B":714.8,"5D":569.9,"6A":622.6,"6B":731.1,"6D":495.3,"7A":744.5,"7B":764.1,"7D":642.8}
num=lambda v:None if v in (None,'') else float(v)
def co(r):
    a=(r.get('asm') or {}).get(FR)
    return a
def pos(c):
    s,e,p=num(c.get('start')),num(c.get('end')),num(c.get('peak'))
    return p if p is not None else ((s+e)/2 if s is not None and e is not None else s)
def iv(c):
    s,e=num(c.get('start')),num(c.get('end'));p=pos(c)
    if p is None: return None
    if s is not None and e is not None: return (s-2.5,e+2.5) if e-s<=20 else (p-10,p+10)
    return (p-2.5,p+2.5)
g=defaultdict(list)
for r in M['reports']:
    c=co(r)
    if not c or c.get('chr') not in ALEN or pos(c) is None: continue
    g[(r['disease'],c['chr'])].append((iv(c),r))
loci=[]
for (d,ch),arr in g.items():
    arr.sort(key=lambda x:x[0][0]);cur=None
    for v,r in arr:
        if not cur or v[0]>cur['hi']: cur={'disease':d,'chr':ch,'lo':max(0,v[0]),'hi':v[1],'reps':[r['id']],'src':{r['source']}};loci.append(cur)
        else: cur['hi']=max(cur['hi'],v[1]);cur['reps'].append(r['id']);cur['src'].add(r['source'])
by=defaultdict(list)
for l in loci: by[l['chr']].append(l)
CH=list(ALEN)
regions=[]
for ch in CH:
    arr=sorted(by[ch],key=lambda l:l['lo']);cur=None
    for l in arr:
        if not cur or l['lo']>cur['hi']: cur={'chr':ch,'lo':l['lo'],'hi':l['hi'],'loci':[l]};regions.append(cur)
        else: cur['hi']=max(cur['hi'],l['hi']);cur['loci'].append(l)
regions=[r for r in regions if len({l['disease'] for l in r['loci']})>1]
for i,R in enumerate(regions,1):
    R['id']='SR%03d'%i;R['diseases']=sorted({l['disease'] for l in R['loci']});R['n_reports']=sum(len(l['reps']) for l in R['loci']);R['n_studies']=len(set().union(*[l['src'] for l in R['loci']]))
print('loci',len(loci),'shared regions',len(regions))
# genes on v2.1 (keyed by v1.1 ID)
cache=json.load(open(U+'/data/WheatNHR_liftover/output/_ref_cache.json'))
g21={k:tuple(v) for k,v in cache['v21'].items() if v[0] in ALEN}
from shared_region_genes import gene_class
cls={}
for row in csv.DictReader(open(U+'/For_database/annotation/iwgsc_refseqv1.0_FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned.TAB',encoding='utf-8',errors='replace'),delimiter='\t'):
    k=row['Gene-ID'].split('.')[0].replace('01G','02G')
    cls[k]=(gene_class(row.get('Pfam-IDs','')+' '+row.get('Pfam-IDs-(Description)',''),row.get('Human-Readable-Description','')),row.get('Human-Readable-Description',''))
zone=lambda ch,p:'distal' if (p/ALEN[ch]<0.1 or p/ALEN[ch]>0.9) else 'interior'
N=Counter();K=defaultdict(Counter)
for gname,(ch,p) in g21.items():
    z=zone(ch,p);N[z]+=1;N['all']+=1
    c=cls.get(gname,('unannotated',''))[0];K[z][c]+=1;K['all'][c]+=1
print('genes',N,'NLR',K['all']['NLR'],'WAK',K['all']['WAK/WAKL'])
def sf(k,n,p):
    if k<=0: return 1.0
    lp,lq=math.log(p),math.log(1-p);s=0.0
    for i in range(k,n+1): s+=math.exp(math.lgamma(n+1)-math.lgamma(i+1)-math.lgamma(n-i+1)+i*lp+(n-i)*lq)
    return min(1.0,s)
# presence in cultivars
cult={}
for f in sorted(glob.glob('/home/claude/pan/anchors_*.tsv')):
    nm=re.sub(r'anchors_Triticum_aestivum_','',os.path.basename(f)).split('.')[0]
    cult[nm]={r['v11_gene']:r['chr'] for r in csv.DictReader(open(f),delimiter='\t')}
tested=set(); [tested.update(d) for d in cult.values()]
qset=set(l[1:].strip() for l in open('/home/claude/pan/cs_all.fa') if l.startswith('>'))
names=list(cult)
summ=[];genes_out=[]
for R in regions:
    gs=sorted((p,gname) for gname,(ch,p) in g21.items() if ch==R['chr'] and R['lo']<=p<=R['hi'])
    n=len(gs);cnt=Counter(cls.get(x,('unannotated',''))[0] for _,x in gs)
    # zone-weighted expectation
    expN={k:sum(K[zone(R['chr'],p)][k]/N[zone(R['chr'],p)] for p,_ in gs) for k in ('NLR','WAK/WAKL')}
    pz={k:(expN[k]/n if n else 0) for k in expN}
    row={'region':R['id'],'chr':R['chr'],'start_Mb_v2.1':round(R['lo'],2),'end_Mb_v2.1':round(R['hi'],2),'diseases':'; '.join(R['diseases']),'n_diseases':len(R['diseases']),'reports':R['n_reports'],'studies':R['n_studies'],'HC_genes':n}
    for k in ('NLR','WAK/WAKL'):
        x=cnt[k];eg=n*K['all'][k]/N['all']
        row[k]=x;row[k+'_expected_genome']=round(eg,1);row[k+'_fold_genome']=round(x/eg,2) if eg else '';row[k+'_p_genome']=sf(x,n,K['all'][k]/N['all']) if n else 1
        row[k+'_expected_zone']=round(expN[k],1);row[k+'_fold_zone']=round(x/expN[k],2) if expN[k] else '';row[k+'_p_zone']=sf(x,n,pz[k]) if n and pz[k]>0 else 1
    for k in ('LRR-RLK','lectin/S-locus RLK','other kinase','ABC transporter (PDR)','PR protein'): row[k]=cnt[k]
    nt=0;pres=Counter();var=[]
    for p,x in gs:
        c,desc=cls.get(x,('unannotated',''))
        rr={'region':R['id'],'chr':R['chr'],'gene_v1.1':x,'pos_Mb_v2.1':round(p,4),'class':c,'description':desc}
        if x in qset:
            nt+=1;miss=0
            for nm in names:
                ok=cult[nm].get(x)==R['chr'];rr[nm]=int(ok);pres[nm]+=ok;miss+=not ok
            if c in ('NLR','WAK/WAKL') and miss: var.append(x)
        genes_out.append(rr)
    row['genes_tested_pangenome']=nt
    row['mean_%present_cultivars']=round(100*sum(pres.values())/(nt*len(names)),1) if nt else ''
    row['NLR/WAK absent in >=1 cultivar']=len(var)
    summ.append(row)
m=len(summ)*2
for r in summ:
    for k in ('NLR','WAK/WAKL'):
        r[k+'_sig_zone_bonf']=r[k+'_p_zone']<0.05/m
        for t in ('genome','zone'): r[k+'_p_'+t]=float('%.2g'%r[k+'_p_'+t])
json.dump({'summary':summ,'genes':genes_out,'names':names,'bonf':0.05/m,'nlr_presence':None},open('sr_v21.json','w'))
sig=[r for r in summ if r['NLR_sig_zone_bonf'] or r['WAK/WAKL_sig_zone_bonf']]
print('tests',m,'threshold',0.05/m,'significant regions (zone baseline):',len(sig))
for r in sorted(sig,key=lambda r:r['NLR_p_zone']):
    print(r['region'],r['chr'],r['start_Mb_v2.1'],r['end_Mb_v2.1'],r['diseases'],'| NLR',r['NLR'],'fold',r['NLR_fold_zone'],'p',r['NLR_p_zone'],'| WAK',r['WAK/WAKL'],r['WAK/WAKL_fold_zone'],r['WAK/WAKL_p_zone'])
# PAV baseline
allt=[x for x in g21 if x in qset]
def rate(gs): return round(100*sum(cult[nm].get(x)==g21[x][0] for x in gs for nm in names)/(len(gs)*len(names)),1)
print('presence: all genes',rate(allt),'NLR',rate([x for x in allt if cls.get(x,('',))[0]=='NLR']),'WAK',rate([x for x in allt if cls.get(x,('',))[0]=='WAK/WAKL']))
