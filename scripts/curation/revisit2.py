import json,re,collections,sys
sys.path.insert(0,'/home/claude/p5run')
from extract import DISQ,MARK,CITE,body,tabley,CHRLEN,CHRLEN21,ENC
d2t=json.load(open('d2t.json'));M=json.load(open('atlas/merged_v13f.json'))
P=json.load(open('markerpos_all.json'))
meta={}
for f in ['p3/zotero_items.json','p5/txt/_items.json']:
    for x in json.load(open(f)):
        if x.get('doi'): meta.setdefault(x['doi'].lower(),x)
for x in json.load(open('p4/match.json'))+json.load(open('p4/sel.json')):
    if x.get('doi'): meta.setdefault(x['doi'].lower(),{'subject':[x.get('dis','')],'title':''})
doi2src={}
for k,v in M['sources'].items():
    u=v.get('url','') if isinstance(v,dict) else ''
    m=re.search(r'10\.\d{4,}/\S+',u)
    if m: doi2src[m.group(0).lower().rstrip('.')]=k
def nn(n): return re.sub(r'[\s(].*','',n).lower().replace('–','-')
allnames=set(nn(r['name']) for r in M['reports'])
srcmark=collections.defaultdict(set)
for r in M['reports']:
    for t in re.split(r'[;,/\s]+',str(r.get('marker') or '')): srcmark[r['source']].add(t.lower())
    srcmark[r['source']].add(nn(r['name']))
def norm(m):
    m=str(m).strip().lower().replace('-','').replace('_','')
    m=re.sub(r'^x(?=(gwm|wms|wmc|barc|cfd|cfa|gdm|psp|wpt|cfe|cwem|mag|ksum|swes|cwm|cnl))','',m).replace('wms','gwm')
    return re.sub(r'^(gwm|wmc|barc|cfd|cfa|gdm|psp)0+(\d)',r'\1\2',m)
def mpos(m,ch):
    L=P.get(norm(m)) or []
    L=[x for x in L if x[0]==ch]
    if not L: return None
    v21=[x for x in L if x[2]=='RefSeq v2.1']
    return (v21 or L)[0]
TC=[('Sb|SB|sb|BB|Bb|Bp|BP|bp|Hb|Lb','Spot blotch'),('Ts|TS|ts|Tan|Tsn|Nec|Chl','Tan spot'),('Snb|SNB|snb|Sng|Snl|Sn|sn|Gb','Septoria nodorum blotch'),('Stb|STB|stb|Std','Septoria tritici blotch'),('Fhb|FHB|fhb|Fhs|fhs|Fhi|Fh|Fdk|FDK|Don|DON|Dis|Isr|Ihb|Ttr|Ns','Fusarium head blight'),('Fcr|FCR|fcr|Crr|Cr|CR','Fusarium crown rot'),('Eye|Pch|pch|Es|ES','Eyespot'),('She|Shes|Sheb|Se|SE|Rs','Sharp eyespot'),('Blast|Bl|BL|Wb|WB|wb|Mot','Wheat blast'),('Ggt|Tak|Ta|TA','Take-all')]
AMB={'lb','hb','dis','ns','nec','chl','ar','gb','rs','ta'}
def disq(q,pd=None):
    code=re.match(r'Q([A-Za-z]+)\.',q)
    if not code: return None
    c=code.group(1).lower()
    if c in AMB: return pd
    for pat,d in TC:
        if re.fullmatch(pat.lower(),c): return d
    return None
SUBJ={'spot blotch':'Spot blotch','tan spot':'Tan spot','tan_spot':'Tan spot','snb':'Septoria nodorum blotch','stb':'Septoria tritici blotch','fhb':'Fusarium head blight','crown_rot':'Fusarium crown rot','crown rot':'Fusarium crown rot','blast':'Wheat blast','wheat blast':'Wheat blast','sharp_eyespot':'Sharp eyespot','sharp eyespot':'Sharp eyespot','eyespot':'Eyespot','take_all':'Take-all','take-all':'Take-all'}
def paperdis(doi,fs,txt):
    s=doi2src.get(doi)
    if s:
        c=collections.Counter(r['disease'] for r in M['reports'] if r['source']==s)
        if c: return c.most_common(1)[0][0]
    m=meta.get(doi)
    if m and m.get('subject'):
        k=str(m['subject'][0]).lower().replace('_',' ')
        for a,b in SUBJ.items():
            if a.replace('_',' ')==k: return b
    for f in fs:
        b=f.split('/')[-1].lower()
        for a,d in SUBJ.items():
            if b.startswith(a.replace(' ','_')): return d
    return None
QR=re.compile(r'\b(Q[A-Za-z][a-zA-Z]{0,6}\.[A-Za-z0-9]{2,12}[-.]?[1-7][ABD][SL]?(?:\.\d{1,2})?[a-z]?)\b')
SNPM=re.compile(r'\b(AX[-_]\d{8,9}|IW[AB]\d{1,6}|BS\d{8}_\d+|BobWhite_c\d+_\d+|Excalibur_(?:rep_)?c\d+_\d+|Kukri_(?:rep_)?c\d+_\d+|RAC875_(?:rep_)?c\d+_\d+|Tdurum_contig\d+_\d+|wsnp_\w+?_\d+|Ku_c\d+_\d+|Ra_c\d+_\d+|tplb\d+\w+_\d+|BS\d{8}|CAP\d+_c\d+_\d+|IACX\d+|IAAV\d+|D_GDS\w+|GENE-\d+_\d+|TA\d{6}-\d{4}|wPt-?\d+|(?:X)?(?:gwm|wmc|barc|cfd|cfa|gdm)\d+|S[1-7][ABD]_\d{5,9}|chr[1-7][ABD]_\d{5,9})\b')
PV=re.compile(r'(\d\.\d+E-\d+|\d\.\d+e-\d+|\d\.\d+\s*×\s*10\s*[-−–]\s*\d+)')
BPR=re.compile(r'(?<![\d.,A-Za-z_\-|:])(\d{1,3}(?:,\d{3}){2}|\d{6,9})(?![\d,]*\.\d)(?![\dA-Za-z_|])')
def rowbp(seg,ch):
    out=[]
    for b in BPR.findall(seg):
        v=int(b.replace(',',''))/1e6
        if 0.05<v<=CHRLEN21.get(ch,900)+3: out.append(round(v,4))
    return out
def designator(q):
    m=re.match(r'Q[A-Za-z]+\.([A-Za-z0-9]+?)[-.]?[1-7][ABD]',q); return m.group(1).lower() if m else ''
new=[];st=collections.Counter();mqseen=set()
for doi,fs in d2t.items():
    doi=doi.lower();txt=''
    for f in fs:
        try: txt+=open(f,errors='ignore').read()+'\n'
        except: pass
    if len(txt)<3000: continue
    B=body(txt);dis=paperdis(doi,fs,txt);src=doi2src.get(doi)
    have=srcmark.get(src,set())
    refv='RefSeq v2.1' if re.search(r'refseq\s*v?\s*2\.1|iwgsc\s*v2\.1',B,re.I) else 'RefSeq v1.0'
    # ---- named QTL ----
    names=QR.findall(B)
    dc=collections.Counter(designator(q) for q in names if disq(q,dis))
    if dc:
        top=dc.most_common(1)[0][1]
        own={d for d,c in dc.items() if c>=3 and c>=0.3*top}
    else: own=set()
    lines=[' '.join(l.split()) for l in B.split('\n')]
    done=set()
    for q in dict.fromkeys(names):
        d=disq(q,dis)
        if not d or designator(q) not in own: continue
        k=nn(q)
        if k in done or k in have or k in allnames: continue
        done.add(k)
        ch=re.search(r'-?([1-7][ABD])',q.split('.',1)[1]).group(1)
        # context: all lines mentioning q (+ next line)
        ctx=[];pos=[];mks=[]
        for i,l in enumerate(lines):
            if q in l and not CITE.search(l):
                seg=l+' '+(lines[i+1] if i+1<len(lines) else '')
                ctx.append(l)
                for m in SNPM.findall(seg):
                    p=mpos(m,ch)
                    if p and m not in mks: mks.append(m);pos.append(p)
        if not ctx: continue
        good=[l for l in ctx if tabley(l) or (re.search(r'detected|identified|mapped|designated|this study|were found|was found|explained|located',l,re.I) and not re.search(r'reported|previous|earlier|et al',l,re.I))]
        if not good: st['qtl_cited_only']+=1;continue
        ctx=good
        rec={'name':q,'cls':'QTL','disease':d,'chr':ch,'doi':doi,'src':src,'marker':' / '.join(mks[:4]),'evidence':ctx[0][:240]}
        rb=[]
        for l in ctx:
            if tabley(l): rb+=rowbp(l,ch)
        if rb:
            v=sorted(set(rb))[:2]
            rec.update(ref=refv,start=v[0] if len(v)>1 else None,end=v[-1] if len(v)>1 else None,peak=v[0] if len(v)==1 else None,grade='B',how='position given in paper table')
            st['qtl_placed_row']+=1
        elif pos:
            # use the frame of the first hit; keep positions of the same ref
            ref=pos[0][2];v=sorted(x[1] for x in pos if x[2]==ref)
            if v[-1]-v[0]>120: v=[v[0]] if len(v)==1 else v[:1]
            rec.update(ref=ref,start=v[0] if len(v)>1 else None,end=v[-1] if len(v)>1 else None,peak=v[0] if len(v)==1 else None,grade='B',how='flanking/peak marker position ('+pos[0][3]+')')
            st['qtl_placed']+=1
        else:
            rec.update(ref=refv,start=None,end=None,peak=None,grade='D',how='named in paper; no physical position found')
            st['qtl_unplaced']+=1
        new.append(rec)
    # ---- meta-QTL rows ----
    if re.search(r'meta-?QTL|meta-analysis|MQTL',B,re.I):
        for l in lines:
            for mq in re.findall(r'\b(MQTL[-_.]?([1-7][ABD])[-_.]?\d{0,2}(?:\.\d{1,2})?)\b',l):
                nm,ch=mq
                k=nn(nm)
                if k in have or (doi,k) in mqseen or CITE.search(l) or not tabley(l): continue
                mqseen.add((doi,k))
                ps=[];mk=[]
                for m in SNPM.findall(l):
                    p=mpos(m,ch)
                    if p: ps.append(p);mk.append(m)
                rec={'name':nm,'cls':'MQTL','disease':dis,'chr':ch,'doi':doi,'src':src,'marker':' / '.join(mk[:2]),'evidence':l[:240]}
                rb=rowbp(l,ch)
                if len(ps)>=1:
                    ref=ps[0][2];v=sorted(x[1] for x in ps if x[2]==ref)
                    if v[-1]-v[0]>150: v=v[:1]
                    rec.update(ref=ref,start=v[0] if len(v)>1 else None,end=v[-1] if len(v)>1 else None,peak=v[0] if len(v)==1 else None,grade='B',how='MQTL flanking markers ('+ps[0][3]+')')
                elif len(rb)>=2:
                    v=sorted(rb[:2]);rec.update(ref=refv,start=v[0],end=v[1],peak=None,grade='B',how='MQTL interval given in paper')
                else:
                    rec.update(ref=refv,start=None,end=None,peak=None,grade='D',how='no physical position found')
                st['mqtl']+=1;new.append(rec)
    # ---- MTAs / QTNs from table rows (GWAS) ----
    gwas=len(re.findall(r'\bGWAS\b|genome-wide association|association mapping|marker[- ]trait association|\bMTAs?\b|\bQTNs?\b',B,re.I))>=3
    if not gwas: continue
    logp=bool(re.search(r'-?\s*log\s*10\s*\(?\s*[pP]|−log10|-log10',B))
    qtn=bool(re.search(r'\bQTNs?\b|mrMLM|FASTmrEMMA|ISIS EM-BLASSO|pLARmEB',B))
    seenm=set()
    for l in lines:
        for seg in re.split(r'\s{6,}',l):
            if CITE.search(seg) or not tabley(seg): continue
            ms=SNPM.findall(seg)
            if not ms or len(ms)>3: continue
            if re.search(r'\b(YR|LR|SR|PM|Yr|Lr|Sr|Pm)[_\d]|rust|mildew|\bKB\b|\bLS\b',seg): continue
            if QR.search(seg) or re.search(r'\b(Lr|Yr|Sr|Pm)\d|\b(PH|HD|DTH|DH|TKW|GY|TGW|SL|PHT|DTA|DTM)\b|yield|height|heading|maturity|grain weight',seg): continue
            if not PV.search(seg):
                if not (logp and re.search(r'(?<![\d.])([3-9]\.\d{1,3}|[1-3]\d\.\d{1,3})(?![\d.%])',seg)): continue
            for m in ms[:1]:
                km=m.lower()
                if km in have or km in seenm: continue
                chm=re.findall(r'(?<![A-Za-z0-9])([1-7][ABD])(?![A-Za-z0-9])',seg)
                enc=re.match(r'(?:S|chr)([1-7][ABD])_(\d+)',m)
                rb=rowbp(seg,chm[0]) if chm and not enc else []
                if enc: ch=enc.group(1);p=(ch,int(enc.group(2))/1e6,refv,'marker-encoded')
                elif rb: ch=chm[0];p=(ch,rb[0],refv,'position given in paper table')
                else:
                    ch=chm[0] if chm else None
                    p=mpos(m,ch) if ch else None
                    if not p and not ch:
                        L=P.get(norm(m)) or []
                        if L: p=L[0];ch=p[0]
                if not ch: continue
                seenm.add(km)
                rec={'name':m,'cls':'QTN' if qtn else 'MTA','disease':dis,'chr':ch,'doi':doi,'src':src,'marker':m,'evidence':seg[:240]}
                if p: rec.update(ref=p[2],start=None,end=None,peak=round(p[1],4),grade='B',how='SNP position ('+p[3]+')');st['mta_placed']+=1
                else: rec.update(ref=refv,start=None,end=None,peak=None,grade='D',how='no physical position found');st['mta_unplaced']+=1
                pvm=PV.search(seg)
                if pvm: rec['stat']='p '+pvm.group(1)
                new.append(rec)
print(st,len(new))
print(collections.Counter((r['src'] is None) for r in new))
print(collections.Counter(r['disease'] for r in new))
json.dump(new,open('revisit_new.json','w'))
