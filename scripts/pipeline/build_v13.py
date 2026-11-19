import json,re
M=json.load(open('/home/claude/atlas/merged_v13_lifted.json'))
old=json.load(open('/home/claude/atlas/db_seed.json'))
seed={'reports':M['reports'],'sources':M['sources'],'papers':old['papers']}
# papers tab: mark papers whose DOI is now in the atlas as extracted
dois=set()
for k,v in M['sources'].items():
    for d in re.findall(r'10\.\d{4,9}/[^\s"\',;]+',json.dumps(v)): dois.add(d.lower().rstrip('.'))
for p in seed['papers']:
    m=re.search(r'10\.\d{4,9}/[^\s"\',;]+',p.get('link') or '')
    if m and m.group(0).lower().rstrip('.') in dois and p.get('status','To search')=='To search': p['status']='Added to atlas'
U=json.load(open('/home/claude/atlas/papers_uploaded.json'))
have={(re.search(r'10\.\d{4,9}/[^\s"\',;]+',p.get('link') or '') or [''])[0].lower().rstrip('.') for p in seed['papers']}
for p in seed['papers']:
    d=(re.search(r'10\.\d{4,9}/[^\s"\',;]+',p.get('link') or '') or [''])[0].lower().rstrip('.')
    u=next((x for x in U if x['link'][16:].lower()==d),None) if d else None
    if u and p.get('status') in ('To search',None,''): p['status']=u['status'];p['notes']=(p.get('notes') or '')+(' ' if p.get('notes') else '')+u['notes']
seed['n_reviewed']=len(U)
seed['papers']+= [u for u in U if u['link'][16:].lower() not in have]
json.dump(seed,open('/home/claude/atlas/db_seed_v13.json','w'))
d='/home/claude/atlas/web/'
sj=json.dumps(seed,ensure_ascii=False,separators=(',',':')).replace('</','<\\/')
html=open(d+'head.html').read()+open(d+'body.html').read()+'<script type="application/json" id="seed">'+sj+'</script>\n<script type="application/json" id="gokegg" data-src="gokegg_lite.json"></script>\n<script type="application/json" id="markerslite" data-src="markers_lite.json"></script>\n<script type="application/json" id="primerslite">'+open(d+'primers_lite.json').read()+'</script>\n<script type="application/json" id="genesall" data-src="genes_all.json"></script>\n<script>\n'+open(d+'app.js').read()+'\n</script>\n</body>\n</html>\n'
open(d+'index.html','w').write(html)
print('reports',len(seed['reports']),'sources',len(seed['sources']),'page MB',round(len(html)/1e6,2))
