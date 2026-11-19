import json,urllib.request,urllib.parse,time,csv
Q={
"Fusarium head blight":'(Fusarium graminearum[tiab] OR "Fusarium head blight"[tiab] OR "head blight"[tiab])',
"Fusarium crown rot":'("crown rot"[tiab] OR "Fusarium pseudograminearum"[tiab])',
"Spot blotch":'("spot blotch"[tiab] OR "Bipolaris sorokiniana"[tiab] OR "Cochliobolus sativus"[tiab])',
"Tan spot":'("tan spot"[tiab] OR "Pyrenophora tritici-repentis"[tiab])',
"Septoria nodorum blotch":'("Parastagonospora nodorum"[tiab] OR "Stagonospora nodorum"[tiab] OR "Septoria nodorum"[tiab])',
"Septoria tritici blotch":'("Zymoseptoria tritici"[tiab] OR "Septoria tritici"[tiab] OR "Mycosphaerella graminicola"[tiab])',
"Wheat blast":'("wheat blast"[tiab] OR "Magnaporthe oryzae"[tiab] OR "Pyricularia oryzae"[tiab])',
"Sharp eyespot":'("sharp eyespot"[tiab] OR "Rhizoctonia cerealis"[tiab])',
"Eyespot":'(eyespot[tiab] OR Oculimacula[tiab] OR Pseudocercosporella[tiab])',
"Take-all":'("take-all"[tiab] OR Gaeumannomyces[tiab])',
"Leaf rust":'("leaf rust"[tiab] OR "Puccinia triticina"[tiab])',
"Stripe rust":'("stripe rust"[tiab] OR "yellow rust"[tiab] OR "Puccinia striiformis"[tiab])',
"Stem rust":'("stem rust"[tiab] OR "Puccinia graminis"[tiab])',
"Powdery mildew":'("powdery mildew"[tiab] OR "Blumeria graminis"[tiab])',
"Common bunt":'("common bunt"[tiab] OR "Tilletia caries"[tiab] OR "Tilletia foetida"[tiab])',
"Dwarf bunt":'("dwarf bunt"[tiab] OR "Tilletia controversa"[tiab])',
"Karnal bunt":'("Karnal bunt"[tiab] OR "Tilletia indica"[tiab])',
"Loose smut":'("loose smut"[tiab] OR "Ustilago tritici"[tiab])',
"Wheat streak mosaic virus":'("wheat streak mosaic"[tiab])',
"Barley yellow dwarf virus":'("barley yellow dwarf"[tiab])',
"Wheat dwarf virus":'("wheat dwarf virus"[tiab])',
"Wheat yellow mosaic virus":'("wheat yellow mosaic"[tiab])',
"Soil-borne cereal mosaic virus":'("soil-borne"[tiab] AND mosaic[tiab])',
"Bacterial leaf streak":'("bacterial leaf streak"[tiab] OR "Xanthomonas translucens"[tiab])',
"Cereal cyst nematode":'("cereal cyst nematode"[tiab] OR "Heterodera avenae"[tiab] OR "Heterodera filipjevi"[tiab])',
"Root-knot nematode":'("root-knot"[tiab] OR Meloidogyne[tiab])',
"Cephalosporium stripe":'("Cephalosporium stripe"[tiab] OR "Cephalosporium gramineum"[tiab])',
"Anthracnose":'(anthracnose[tiab] AND Colletotrichum[tiab])',
"Crown rust":'("crown rust"[tiab] OR "Puccinia coronata"[tiab])',
"Leaf spot complex":'("leaf spot"[tiab] AND (Alternaria[tiab] OR Bipolaris[tiab]))',
}
base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(u):
    for i in range(4):
        try: return json.load(urllib.request.urlopen(u,timeout=60))
        except Exception as e: time.sleep(2)
    raise
rows={}
for dis,q in Q.items():
    term=f'(wheat[tiab] OR Triticum[tiab]) AND (RNA-seq[tiab] OR transcriptom*[tiab] OR "differentially expressed"[tiab] OR "gene expression"[tiab]) AND {q} AND open access[filter]'
    r=get(base+'esearch.fcgi?db=pmc&retmode=json&retmax=200&term='+urllib.parse.quote(term))
    ids=r['esearchresult']['idlist'];print(dis,r['esearchresult']['count']);time.sleep(0.4)
    for i in range(0,len(ids),100):
        s=get(base+'esummary.fcgi?db=pmc&retmode=json&id='+','.join(ids[i:i+100]));time.sleep(0.4)
        for k in s['result'].get('uids',[]):
            d=s['result'][k];doi=next((a['value'] for a in d.get('articleids',[]) if a['idtype']=='doi'),'')
            pmc=next((a['value'] for a in d.get('articleids',[]) if a['idtype']=='pmcid'),'PMC'+k)
            rec=rows.setdefault(pmc,{'pmcid':pmc,'title':d.get('title',''),'journal':d.get('fulljournalname',''),'year':d.get('pubdate','')[:4],'doi':doi,'first_author':(d.get('authors') or [{'name':''}])[0]['name'],'diseases':set()})
            rec['diseases'].add(dis)
with open('candidates.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['pmcid','year','first_author','title','journal','doi','diseases'])
    for r in sorted(rows.values(),key=lambda r:r['year'],reverse=True): w.writerow([r['pmcid'],r['year'],r['first_author'],r['title'],r['journal'],r['doi'],'; '.join(sorted(r['diseases']))])
print('total',len(rows))
