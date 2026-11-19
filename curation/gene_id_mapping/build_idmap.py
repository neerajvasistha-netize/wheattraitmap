import gzip,json,re,os,glob,collections
OUT='/home/claude/wheatnhr-atlas/web/idmap';os.makedirs(OUT,exist_ok=True)
G=json.load(open('/home/claude/wheatnhr-atlas/web/genes_all.json'))
uni=set()
for c,rows in G['g'].items():
    for r in rows: uni.add(f"{c}:{r[1]}")
def key(gid):
    m=re.match(r'TraesCS([1-7][ABD]|U)02G0*(\d+)',gid)
    if not m: return None
    k=f"{m.group(1)}:{int(m.group(2))}"; return k
def dump(name,M):
    D={}
    for x,ks in M.items():
        ks=sorted(ks); D[x]="|".join(ks)
    gzip.open(f'{OUT}/{name}.json.gz','wt',compresslevel=9).write(json.dumps(D,separators=(',',':')))
    print(name,len(D),os.path.getsize(f'{OUT}/{name}.json.gz')//1024,'KB', 'inUni%.1f'%(100*sum(any(k in uni for k in v.split('|')) for v in D.values())/max(1,len(D))))
def tsv(f,dbs,norm):
    M=collections.defaultdict(set)
    for ln in gzip.open(f,'rt'):
        p=ln.rstrip('\n').split('\t')
        if len(p)<5 or p[4] not in dbs: continue
        k=key(p[0]);
        if k: M[norm(p[3])].add(k)
    return M
dl='dl/'
dump('entrez',tsv(dl+'Triticum_aestivum.IWGSC.63.entrez.tsv.gz',{'EntrezGene'},lambda x:x))
dump('refseq',tsv(dl+'Triticum_aestivum.IWGSC.63.refseq.tsv.gz',{'RefSeq_mRNA','RefSeq_mRNA_predicted','RefSeq_peptide','RefSeq_peptide_predicted'},lambda x:x.split('.')[0].upper()))
dump('uniprot',tsv(dl+'Triticum_aestivum.IWGSC.63.uniprot.tsv.gz',{'Uniprot/SPTREMBL','Uniprot/SWISSPROT','UniProtKB_all'},lambda x:x.split('.')[0].split('-')[0].upper()))
# pangenome
for f in sorted(glob.glob('../pan/anchors_full/*.tsv')):
    M=collections.defaultdict(set);pre=None
    for i,ln in enumerate(open(f)):
        if i==0: continue
        p=ln.rstrip('\n').split('\t'); g=p[4]
        if g.startswith('TraesCS'): break
        k=key(p[0])
        if not k: continue
        m=re.match(r'(Traes[A-Z]+)',g); pre=m.group(1)
        M[re.sub(r'\.\d+$','',g).upper()].add(k)
    if M: dump('pan_'+pre,M)
# After building: for f in *.json.gz; do base64 -w0 $f > ${f%.json.gz}.gz.txt; done   (the page reads the .gz.txt files)
