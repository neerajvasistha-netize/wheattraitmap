#!/usr/bin/env python3
"""Add the remaining wheat diseases (bunts, smut, viruses, bacterial leaf streak,
cereal cyst nematode) to the WheatNHR-Atlas master file."""
import json, re, os, sys, shutil

HOME = '/home/claude'
MASTER = os.path.join(HOME, 'atlas/merged_v13f.json')
BACKUP = os.path.join(HOME, 'atlas/merged_v13f_preNewDis.json')

d = json.load(open(MASTER))
reports, sources = d['reports'], d['sources']
have = {r['id'] for r in reports}
nxt = max(int(r['id'][1:]) for r in reports if re.fullmatch(r'V\d+', r['id'])) + 1

MP = json.load(open(os.path.join(HOME, 'markerpos_all.json')))
def mkey(s): return re.sub(r'[^a-z0-9]', '', str(s).lower())

def place(marker, chrom=None, prefer='RefSeq v1.0'):
    """-> (pos_Mb, assembly, note) or (None, None, None)"""
    e = MP.get(mkey(marker))
    if not e:
        return None, None, None
    cands = [x for x in e if chrom is None or x[0] == chrom]
    if not cands:
        return None, None, None
    for pref in (prefer, 'RefSeq v2.1', 'RefSeq v1.0'):
        for x in cands:
            if x[2] == pref:
                return x[1], x[2], x[3]
    x = cands[0]
    return x[1], x[2], x[3]

NEW = []
def add(name, cls, disease, chrom, start, end, peak, ref, marker, pve, stat,
        source, grade, note, evidence=''):
    global nxt
    rid = 'V%04d' % nxt
    while rid in have:
        nxt += 1; rid = 'V%04d' % nxt
    have.add(rid); nxt += 1
    r = {'id': rid, 'name': name, 'cls': cls, 'disease': disease, 'chr': chrom,
         'start': start, 'end': end, 'peak': peak, 'ref': ref, 'marker': marker,
         'pve': '' if pve in (None, '') else str(pve), 'stat': stat,
         'source': source, 'grade': grade, 'curation': 'auto', 'note': note}
    if evidence:
        r['evidence'] = evidence
    NEW.append(r)
    return rid

def mid(a, b):
    if a is None or b is None: return None
    return round((a + b) / 2.0, 4)

# ---------------------------------------------------------------- sources
SRC = {
 'Pundir2023BMC': ("Pundir S, Singh R, Singh VK, Sharma S, Balyan HS, Gupta PK, Sharma S (2023) "
   "Mapping of QTLs and meta-QTLs for Heterodera avenae Woll. resistance in common wheat "
   "(Triticum aestivum L.). BMC Plant Biology 23:529",
   "https://doi.org/10.1186/s12870-023-04526-y", "Pundir et al. (2023)"),
 'Pundir2022SR': ("Pundir S, Sharma R, Kumar D, Singh VK, Chaturvedi D, Kanwar RS, Roder MS, Borner A, "
   "Ganal MW, Gupta PK, Sharma S, Sharma S (2022) QTL mapping for resistance against cereal cyst "
   "nematode (Heterodera avenae Woll.) in wheat (Triticum aestivum L.). Scientific Reports 12:9586",
   "https://doi.org/10.1038/s41598-022-12988-7", "Pundir et al. (2022)"),
 'Chaturvedi2023SR': ("Chaturvedi D, Pundir S, Singh VK, Kumar D, Sharma R, Roder MS, Sharma S, Sharma S "
   "(2023) Identification of genomic regions associated with cereal cyst nematode (Heterodera avenae "
   "Woll.) resistance in spring and winter wheat. Scientific Reports 13:5916",
   "https://doi.org/10.1038/s41598-023-32737-8", "Chaturvedi et al. (2023)"),
 'Muellner2021TAG': ("Muellner AE, Buerstmayr M, Eshonkulov B, Hole D, Michel S, Hagenguth JF, Pachler B, "
   "Pernold R, Buerstmayr H (2021) Comparative mapping and validation of multiple disease resistance QTL "
   "for simultaneously controlling common and dwarf bunt in bread wheat. Theoretical and Applied "
   "Genetics 134:489-503", "https://doi.org/10.1007/s00122-020-03708-8", "Muellner et al. (2021)"),
 'Lunzer2023TAG': ("Lunzer M, Buerstmayr M, Grausgruber H, Mullner AE, Fallbacher I, Buerstmayr H (2023) "
   "Wheat (Triticum aestivum) chromosome 6D harbours the broad spectrum common bunt resistance gene Bt11. "
   "Theoretical and Applied Genetics 136:207", "https://doi.org/10.1007/s00122-023-04452-5",
   "Lunzer et al. (2023)"),
 'Emebiri2021FPS': ("Emebiri L, Hildebrand S, Tan MK, Juliana P, Singh PK, Fuentes-Davila G, Singh RP (2021) "
   "Pre-emptive breeding against Karnal bunt infection in common wheat: combining genomic and agronomic "
   "information to identify suitable parents. Frontiers in Plant Science 12:675859",
   "https://doi.org/10.3389/fpls.2021.675859", "Emebiri et al. (2021)"),
 'Thambugala2020BMC': ("Thambugala D, Menzies JG, Knox RE, Campbell HL, McCartney CA (2020) Genetic "
   "analysis of loose smut (Ustilago tritici) resistance in Sonop spring wheat. BMC Plant Biology 20:314",
   "https://doi.org/10.1186/s12870-020-02525-x", "Thambugala et al. (2020)"),
 'Kumar2018PLOS': ("Kumar S, Knox RE, Singh AK, DePauw RM, Campbell HL, Isidro-Sanchez J, Clarke FR, "
   "Pozniak CJ, N'Daye A, Meyer B, Sharpe A, Ruan Y, Cuthbert RD, Somers D, Fedak G (2018) High-density "
   "genetic mapping of a major QTL for resistance to multiple races of loose smut in a tetraploid wheat "
   "cross. PLoS ONE 13:e0192261", "https://doi.org/10.1371/journal.pone.0192261", "Kumar et al. (2018)"),
 'Manan2025TAG': ("Manan F, Li X, Szabo-Hever A, Shi G, Elias E, Faris J, Xu SS, Liu Z (2025) Genome-wide "
   "association and biparental mapping revealed a major quantitative trait locus associated with seedling "
   "resistance to bacterial leaf streak in durum wheat. Theoretical and Applied Genetics 139:10",
   "https://doi.org/10.1007/s00122-025-05111-7", "Manan et al. (2025)"),
 'Enteshari2026BMC': ("Enteshari M, Massah A, Mehrabi R, Talebi R, Soorni A, Talebi M (2026) Genome-wide "
   "association mapping and characterization of genetic loci associated with resistance against wheat "
   "streak mosaic virus. BMC Plant Biology 26:1727", "https://doi.org/10.1186/s12870-026-09310-2",
   "Enteshari et al. (2026)"),
 'Xie2022FPS': ("Xie Y, Nachappa P, Nalam VJ, Pearce S (2022) Genomic and molecular characterization of "
   "wheat streak mosaic virus resistance locus 2 (Wsm2) in common wheat (Triticum aestivum L.). "
   "Frontiers in Plant Science 13:928949", "https://doi.org/10.3389/fpls.2022.928949", "Xie et al. (2022)"),
 'Bian2025TPG': ("Bian R, Liu N, Xu Y, Su Z, Chai L, Bernardo A, St Amand P, Rupp J, Pumphrey M, Fritz A, "
   "Zhang G, Jordan KW, Bai G (2025) A novel quantitative trait locus for barley yellow dwarf virus "
   "resistance and kernel traits on chromosome 2D of a wheat mutant. The Plant Genome 18:e20548",
   "https://doi.org/10.1002/tpg2.20548", "Bian et al. (2025)"),
 'Buerstmayr2023TAG': ("Buerstmayr M, Buerstmayr H (2023) Two major quantitative trait loci control wheat "
   "dwarf virus resistance in four related winter wheat populations. Theoretical and Applied Genetics "
   "136:103", "https://doi.org/10.1007/s00122-023-04349-3", "Buerstmayr and Buerstmayr (2023)"),
 'Pfrieme2022FPS': ("Pfrieme AK, Ruckwied B, Habekuss A, Will T, Stahl A, Pillen K, Ordon F (2022) "
   "Identification and validation of quantitative trait loci for wheat dwarf virus resistance in wheat "
   "(Triticum spp.). Frontiers in Plant Science 13:828639", "https://doi.org/10.3389/fpls.2022.828639",
   "Pfrieme et al. (2022)"),
 'Kawaguchi2024BS': ("Kawaguchi K, Ohki T, Ishikawa G, Sayama M, Terasawa Y, Oda S, Fujita M, Ito M, "
   "Hatta K (2024) Identification of a major QTL conferring resistance to wheat yellow mosaic virus "
   "derived from the winter wheat 'Hokkai 240'. Breeding Science 74:232-239",
   "https://doi.org/10.1270/jsbbs.23079", "Kawaguchi et al. (2024)"),
 'Chen2025NC': ("Chen Y, Kong D, Wang Z, Liu J, Wang L, Dai K, Ji J, Chen W, Tang X, Wen M, Zhang X, "
   "Zhang H, Jiao C, Sun L, Wang H, Fei X, Guo H, Sun B, Tao X, Wang W, et al. (2025) A wheat CC-NBS-LRR "
   "protein Ym1 confers WYMV resistance by recognizing viral coat protein. Nature Communications 16:3630",
   "https://doi.org/10.1038/s41467-025-58816-0", "Chen et al. (2025)"),
}
for k, (cit, url, short) in SRC.items():
    if k not in sources:
        sources[k] = {'citation': cit, 'url': url, 'short': short}

CCN = 'Cereal cyst nematode'

# ======================================================= 1. CCN, Pundir 2023
# Table 2 (LOD/PVE) + Table S3 (physical positions), IWGSC RefSeq v1.0
for nm, c, s, e, lod, pve, meth in [
    ('QCcnr.ccsu-1B.1', '1B', 49.940952, 49.941015, 2.80, 7, 'CIM'),
    ('QCcnr.ccsu-1B.2', '1B', 49.927166, 49.927207, 2.76, 7, 'CIM'),
    ('QCcnr.ccsu-1B.3', '1B', 59.193613, 70.663028, 2.95, 6, 'ICIM'),
    ('QCcnr.ccsu-2A.1', '2A', 31.045433, 31.045496, 3.62, 14, 'CIM/ICIM'),
    ('QCcnr.ccsu-2A.2', '2A', 50.038871, 50.038902, 2.59, 7, 'CIM/ICIM'),
    ('QCcnr.ccsu-2A.3', '2A', 42.194565, 43.476992, 2.88, 7, 'CIM'),
    ('QCcnr.ccsu-2A.4', '2A', 746.982564, 747.266450, 3.35, 7, 'ICIM'),
    ('QCcnr.ccsu-3A.1', '3A', 432.953612, 498.529974, 2.53, 9, 'CIM'),
]:
    add(nm, 'QTL', CCN, c, s, e, mid(s, e), 'RefSeq v1.0', '', pve,
        'LOD %.2f' % lod, 'Pundir2023BMC', 'B',
        'Main-effect QTL for Heterodera avenae resistance, HUW 468 x C 306 RIL population '
        '(%s); physical interval from Table S3, LOD/PVE from Table 2.' % meth)

# previously reported loci, positions anchored by Pundir et al. (2023) Table S3
for nm, cls, c, s, e, orig in [
    ('QCcn-1B', 'QTL', '1B', 189.338346, 189.338365, 'Al-Ateeq et al.'),
    ('QCre.srd-1B', 'QTL', '1B', 664.520900, 684.861809, 'Williams et al. (2006); Jayatilake et al. (2013)'),
    ('QCre-ma2A', 'QTL', '2A', 2.471848, 2.471868, 'Cui et al.'),
    ('QCre.pau-2A', 'QTL', '2A', 2.320665, 2.320686, 'Singh et al.'),
    ('Cre5', 'DG', '2A', 28.204878, 33.704080, 'Williams et al. (2006); Jayatilake et al. (2013)'),
    ('QCcn-3A', 'QTL', '3A', 712.489366, 721.225071, 'Al-Ateeq et al.'),
]:
    add(nm, cls, CCN, c, s, e, mid(s, e), 'RefSeq v1.0', '', '', '',
        'Pundir2023BMC', 'C',
        'Previously reported locus for Heterodera avenae resistance; physical position anchored '
        'from linked SSR/SNP markers by Pundir et al. (2023) Table S3. Original study: %s.' % orig)

# meta-QTL (Table 4)
for nm, c, peak, ci, pve, nq, fl in [
    ('MQTL1.2A', '2A', 47.2, 3.1, 8.6, 5, 'BS00121553-BS00185908'),
    ('MQTL2.2A', '2A', 29.8, 6.5, 17.5, 2, 'gwm071d-BS00166296'),
]:
    add(nm, 'MQTL', CCN, c, None, None, peak, 'RefSeq v1.0', fl, pve,
        'CI %.1f cM' % ci, 'Pundir2023BMC', 'B',
        'Meta-QTL from %d H. avenae QTL projected on a 76,668-marker consensus map '
        '(Pundir et al. 2023, Table 4); peak given in Mb, confidence interval in cM.' % nq)

# ======================================================= 2. CCN, Pundir 2022
P22 = [
 ('QCcn.ha-2D', 'Excalibur_rep_c67599_2154', '2D', '9.23-11.37', 0.2284, 650.325169),
 ('QCcn.ha-2A', 'Kukri_rep_c68068_95', '2A', '5.30-9.05', 0.03971, 766.400394),
 ('QCcn.ha-5B', 'IACX6116', '5B', '4.39-5.86', 0.04186, 657.711642),
 ('QCcn.ha-7B', 'wsnp_Ex_c26747_35976442', '7B', '3.67-4.42', 0.02799, 174.571382),
 ('QCcn.ha-6B', 'BS00045453_51', '6B', '3.26-4.39', 7.02e-06, None),
 ('QCcn.ha-7D', 'wsnp_Ex_rep_c68671_67525179', '7D', '3.63-4.34', 0.01733, None),
 ('QCcn.ha-5A', 'wsnp_Ku_c21275_31007309', '5A', '3.37-4.18', 0.003085, None),
 ('QCcn.ha-4B', 'Tdurum_contig15737_728', '4B', '3.28-4.17', 0.01052, None),
 ('QCcn.ha-3B', 'Tdurum_contig15050_103', '3B', '3.34-4.09', 0.001591, None),
 ('QCcn.ha-1A', 'Ku_c3523_1959', '1A', '3.31-4.25', 0.002364, None),
 ('QCcn.ha-7A', 'wsnp_Ku_c6065_10682531', '7A', '3.02-3.75', 0.01782, None),
 ('QCcn.ha-5D', 'wms0174', '5D', '4.41', 0.00035, None),
 ('QCcn.ha-2B.1', 'wsnp_RFL_Contig3917_4326857', '2B', '4.25', 0.00045, None),
 ('QCcn.ha-4A', 'Kukri_c29142_1046', '4A', '4.11', 0.01452, None),
 ('QCcn.ha-3A', 'BS00086051_51', '3A', '3.65', 0.00144, None),
 ('QCcn.ha-6A', 'RAC875_c29850_102', '6A', '3.63', 0.00015, None),
 ('QCcn.ha-2B.2', 'Tdurum_contig12159_468', '2B', '3.59', 0.01943, None),
 ('QCcn.ha-1D', 'wsnp_RFL_Contig2036_1264133', '1D', '3.56', 0.000279, None),
 ('QCcn.ha-1B', 'Ku_c1932_1583', '1B', '3.54', 0.004933, None),
]
for nm, mk, c, lod, r2, tab in P22:
    if tab is not None:
        pos, asm, msrc = tab, 'RefSeq v1.0', 'position of linked candidate gene, Table 3'
    else:
        pos, asm, msrc = place(mk, c)
    add(nm, 'QTL', CCN, c, None, None, None if pos is None else round(pos, 4),
        asm or 'RefSeq v1.0', mk, round(r2 * 100, 3), 'LOD %s' % lod,
        'Pundir2022SR', 'B' if pos is not None else 'D',
        'QTL for H. avenae resistance, ITMI doubled-haploid population (CIM); '
        'PVE given as R2 x 100.' + ('' if pos is not None else ' Peak marker could not be placed physically.'),
        'peak marker %s%s' % (mk, '; ' + msrc if msrc else ''))

# ======================================================= 3. CCN, Chaturvedi 2023
import openpyxl
xl = os.path.join(HOME, 'bio/dl/nem/x_PMC10090075/41598_2023_32737_MOESM2_ESM.xlsx')
wb = openpyxl.load_workbook(xl, data_only=True)
# Table S7 gene locations for the 11 significant MTAs
s7 = {}
def sheet(wb, tag):
    for n in wb.sheetnames:
        if re.sub(r'\s+', '', n).lower().startswith(re.sub(r'\s+', '', tag).lower()):
            return wb[n]
    raise KeyError(tag + ' / ' + str(wb.sheetnames))
ws = sheet(wb, 'SupplementaryTableS7')
cur = None
for row in ws.iter_rows(min_row=3, values_only=True):
    mk, chrom, loc = row[1], row[3], row[6]
    if mk: cur = str(mk).strip()
    if cur and chrom and loc and re.match(r'^[\d,]+-[\d,]+$', str(loc).strip()):
        a, b = [float(x.replace(',', '')) / 1e6 for x in str(loc).strip().split('-')]
        if cur not in s7:
            s7[cur] = (str(chrom).strip(), a, b)
# Table S6: unique significant MTAs
ws = sheet(wb, 'SupplementaryTableS6')
seen = {}
for row in ws.iter_rows(min_row=3, values_only=True):
    snp, chrom, cm, p = row[2], row[3], row[4], row[5]
    if not snp or not chrom: continue
    snp = str(snp).strip(); chrom = str(chrom).strip().replace(' ', '')
    try: p = float(p)
    except (TypeError, ValueError): p = None
    k = (snp, chrom)
    if k not in seen or (p is not None and (seen[k] is None or p < seen[k])):
        seen[k] = p
for (snp, chrom), p in sorted(seen.items(), key=lambda x: (x[0][1], x[0][0])):
    if snp in s7 and s7[snp][0].replace(' ', '') == chrom:
        c, a, b = s7[snp]
        pos, asm, ev = mid(a, b), 'RefSeq v1.0', 'SNP in %s:%.4f-%.4f Mb (Table S7)' % (chrom, a, b)
    else:
        pos, asm, msrc = place(snp, chrom)
        ev = 'marker %s%s' % (snp, '; ' + msrc if msrc else '')
    add(snp, 'MTA', CCN, chrom, None, None, None if pos is None else round(pos, 4),
        asm or 'RefSeq v1.0', snp, '', '' if p is None else 'p %.2e' % p,
        'Chaturvedi2023SR', 'B' if pos is not None else 'D',
        'Marker-trait association for H. avenae (cyst count) in a panel of 180 spring and winter '
        'wheat accessions (GWAS, multiple models).' +
        ('' if pos is not None else ' Marker could not be placed physically.'), ev)

# ======================================================= 4. common / dwarf bunt
CB, DB = 'Common bunt', 'Dwarf bunt'
for nm, dis, c, s, e, lod, pve, pop in [
    ('QBt.ifa-1BS', CB, '1B', 7.99, 22.07, 33.0, 41.8, 'MP-CON'),
    ('QBt.ifa-1AL', CB, '1A', 498.54, 513.66, 27.5, 32.5, 'MP-CON'),
    ('QBt.ifa-7AL', CB, '7A', 722.30, 736.69, 17.7, 18.4, 'MP-CON'),
    ('QBt.ifa-1AL', DB, '1A', 498.54, 516.67, 16.2, 28.0, 'MP-CON'),
    ('QBt.ifa-7AL', DB, '7A', 722.30, 736.69, 9.1, 14.4, 'MP-CON'),
    ('QBt.ifa-7DS', DB, '7D', 12.47, 15.32, 5.9, 13.6, 'MP-BLI'),
]:
    add(nm, 'QTL', dis, c, s, e, mid(s, e), 'RefSeq v1.0', '', pve, 'LOD %.1f' % lod,
        'Muellner2021TAG', 'A',
        'QTL for %s incidence, Blizzard/Bonneville x Rainer RIL populations (%s, multiple QTL '
        'mapping across field experiments); 15K SNP array, physical interval on IWGSC RefSeq v1.0.'
        % (dis.lower(), pop))

# Lunzer et al. 2023, IWGSC RefSeq v2.1
for nm, cls, c, s, e, peak, lod, pve, mk, pop in [
    ('Bt11', 'DG', '6D', 482.8, 495.2, 492.6, 22.4, 49.6, 'BobWhite_c13435_700', 'MP-PR1/PR2/PL/MM'),
    ('QTL_4B', 'QTL', '4B', 11.5, 28.6, 13.6, 8.5, 24.1, 'Excalibur_c64418_447', 'MP-PR1/PL'),
    ('QTL_4B', 'QTL', '4B', 662.9, 671.4, 671.4, 4.4, 10.2, 'AX-158542232', 'MP-MM'),
    ('QTL_1B', 'QTL', '1B', 12.2, 46.9, 28.2, 19.4, 42.3, 'AX-158570920', 'MP-PR2'),
    ('QTL_1A', 'QTL', '1A', 355.2, 549.1, 403.9, 4.2, 6.7, 'AX-94404955', 'MP-PR2'),
    ('QTL_2A', 'QTL', '2A', 0.3, 35.1, 31.8, 10.5, 28.3, 'AX-94684111', 'MP-MM'),
    ('QTL_7B', 'QTL', '7B', 7.1, 26.4, 10.1, 12.0, 21.6, 'AX-94810990', 'MP-PR2'),
]:
    add(nm, cls, CB, c, s, e, peak, 'RefSeq v2.1', mk, pve, 'LOD up to %.1f' % lod,
        'Lunzer2023TAG', 'A' if cls == 'DG' else 'B',
        ('Broad-spectrum common bunt resistance gene Bt11, mapped in four bi-parental populations '
         '(%s); interval and peak on IWGSC RefSeq v2.1.' % pop) if cls == 'DG' else
        ('Common bunt resistance QTL, %s; interval and peak marker on IWGSC RefSeq v2.1. '
         'Not formally designated in the source.' % pop))

# ======================================================= 5. Karnal bunt
KB = 'Karnal bunt'
for mk, c, pos, p, fdr, r2 in [
    ('2282741', '1A', 481.52, 2.08e-05, 3.33e-02, 26.91),
    ('1249729', '2A', 723.62, 1.45e-04, 1.16e-01, 7.58),
    ('1037716', '3B', 618.02, 7.64e-07, 6.12e-03, 29.47),
    ('993727', '4A', 719.43, 4.44e-05, 5.08e-02, 26.31),
    ('1128414', '5A', 618.27, 6.07e-06, 2.43e-02, 27.87),
    ('989877', '6B', 683.23, 3.06e-04, 1.88e-01, 8.09),
]:
    add('QTL_' + c, 'MTA', KB, c, None, None, pos, 'RefSeq v1.0', mk, r2,
        'p %.2e; FDR %.2e' % (p, fdr), 'Emebiri2021FPS', 'B',
        'DArTseq marker-trait association for Karnal bunt (Tilletia indica) resistance in the '
        'WH542 x HD29 and WH542 x W485 populations; marker sequence BLASTed to IWGSC RefSeq v1.0. '
        'Not formally designated in the source.',
        'peak DArTseq marker %s' % mk)

# ======================================================= 6. loose smut
LS = 'Loose smut'
for nm, c, mk, lod, pve, note in [
    ('QUt.mrc-3B', '3B', 'BS00030581_51', 4.21, 5.64, 'races T9 and T39, loose smut incidence'),
    ('QUt.mrc-4B', '4B', 'JD_c11606_1380', 3.81, 3.09, 'races T9 and T39 (ICIM), loose smut incidence'),
    ('QUt.mrc-5B', '5B', 'Excalibur_c11605_156', 65.25, 43.31, 'races T2, T9 and T39, loose smut incidence'),
    ('QUt.mrc-7B', '7B', 'Kukri_rep_c71778_644', 86.49, 89.24, 'race T2, loose smut incidence'),
    ('QUtp.mrc-5B', '5B', 'Excalibur_c11605_156', 74.96, 94.11, 'races T2, T9 and T39, fully smutted index'),
]:
    pos, asm, msrc = place(mk, c)
    add(nm, 'QTL', LS, c, None, None, None if pos is None else round(pos, 4),
        asm or 'RefSeq v1.0', mk, pve, 'LOD up to %.2f' % lod, 'Thambugala2020BMC',
        'B' if pos is not None else 'D',
        'QTL for Ustilago tritici resistance, Diamont/TD-14 doubled-haploid population (%s). '
        'Physical position inferred from the closest marker by the atlas.' % note,
        'closest marker %s%s' % (mk, '; ' + msrc if msrc else ''))

for nm, c, mk, lod, pve, races in [
    ('QUt.spa-3A.2', '3A', 'Kukri_c10751_264', 9.9, 21.2, 'T26'),
    ('QUt.spa-6B.2', '6B', 'Tdurum_contig76997_244', 29.2, 74.9, 'T26, T32, T33 and race mixtures'),
    ('QUt.spa-7A.2', '7A', 'Tdurum_contig67992_182', 6.0, 8.0, 'T33'),
]:
    pos, asm, msrc = place(mk, c)
    add(nm, 'QTL', LS, c, None, None, None if pos is None else round(pos, 4),
        asm or 'RefSeq v1.0', mk, pve, 'LOD up to %.1f' % lod, 'Kumar2018PLOS',
        'B' if pos is not None else 'D',
        'QTL for Ustilago tritici resistance in the tetraploid (durum) Strongfield/Blackbird '
        'doubled-haploid population, race(s) %s. Physical position inferred from the associated '
        'marker by the atlas.' % races,
        'associated marker %s%s' % (mk, '; ' + msrc if msrc else ''))

# ======================================================= 7. bacterial leaf streak
BLS = 'Bacterial leaf streak'
for mk, c, pos, p, q, r2, trait in [
    ('RAC875_c13610_2646', '6A', 1.202829, 2.1061e-8, 0.00028, 6.72, 'average infection type'),
    ('BS00021688_51', '5A', 657.084054, 4.8118e-6, 0.03145, 4.42, 'average infection type'),
    ('Excalibur_c16470_623', '5A', 656.981441, 1.275e-5, 0.03334, 4.02, 'average infection type'),
    ('Excalibur_c16470_494', '5A', 656.981570, 1.275e-5, 0.03334, 4.02, 'average infection type'),
    ('Excalibur_c16470_278', '5A', 656.981786, 1.275e-5, 0.03334, 4.02, 'average infection type'),
    ('RAC875_c13610_2646', '6A', 1.202829, 2.5322e-7, 0.00331, 5.65, 'average % water-soaking'),
    ('BS00021688_51', '5A', 657.084054, 6.251e-6, 0.04086, 4.31, 'average % water-soaking'),
]:
    add(mk, 'MTA', BLS, c, None, None, pos, 'RefSeq v1.0', mk, r2,
        'p %.2e; q %.5f' % (p, q), 'Manan2025TAG', 'B',
        'Marker-trait association for seedling resistance to bacterial leaf streak (Xanthomonas '
        'translucens pv. undulosa) in the Global Durum Panel (%s).' % trait)
for nm, c, s, e, lod, r2, pop, fl in [
    ('QBls.zhl-6A.1', '6A', 8.610346, 15.577414, 13.0, 23.0, 'RP336', '6A_8610346-6A_15577414'),
    ('QBls.zhl-6A.2', '6A', 598.189921, 614.872222, 4.1, 11.3, 'RP336', '6A_598189921-6A_614872222'),
    ('QBls.zhl-3A', '3A', 1.214628, 20.232548, 4.0, 12.0, 'RP336', '3A_1214628-3A_20232548'),
]:
    add(nm, 'QTL', BLS, c, s, e, mid(s, e), 'RefSeq v1.0', fl, r2, 'LOD up to %.1f' % lod,
        'Manan2025TAG', 'A',
        'QTL for bacterial leaf streak reaction (infection type and %% water-soaked area) in the '
        'durum %s biparental population; flanking marker names carry the physical position.' % pop)
pos, asm, msrc = place('IWA2743', '5A')
pos2, _, _ = place('IWA7044', '5A')
s, e = (min(x for x in (pos, pos2) if x), max(x for x in (pos, pos2) if x)) if (pos or pos2) else (None, None)
add('QBls.zhl-5A', 'QTL', BLS, '5A', s, e, mid(s, e), asm or 'RefSeq v1.0',
    'IWA2743-IWA7044', 9.5, 'LOD 4.1', 'Manan2025TAG', 'B' if s else 'D',
    'QTL for bacterial leaf streak infection type in the durum BP025 population; physical interval '
    'inferred from the flanking 90K markers by the atlas.',
    'flanking markers IWA2743-IWA7044' + ('; ' + msrc if msrc else ''))

# ======================================================= 8. WSMV
WSMV = 'Wheat streak mosaic virus'
add('Wsm2', 'DG', WSMV, '3B', 17.1, 18.9, 18.0, 'RefSeq v1.0', 'IWA6024/exome SNPs', '',
    '', 'Xie2022FPS', 'A',
    'Wheat streak mosaic virus resistance locus 2, delimited to a 1.8 Mb telomeric interval on 3BS '
    '(17.1-18.9 Mb, IWGSC RefSeq v1.0) containing 142 annotated genes; linkage mapping in a doubled-'
    'haploid population plus exome capture. The Wsm2 haplotype is absent from Chinese Spring and the '
    'other sequenced genomes, so the causal variant is not represented in RefSeq v1.0.')
for mk, c, bp, p, maf, fdr, eff, model, dpi in [
    ('4005455', '1B', 382503775, None, 0.312, 0.000, 0.208, 'FarmCPU', 35),
    ('1133336-57-C/A', '2A', 746982564, None, 0.171, 0.045, None, 'FarmCPU, BLINK', 35),
    ('1070207-37-T/C', '2B', 221486089, None, 0.102, 0.046, None, 'FarmCPU, BLINK', 28),
    ('2254436-32-G/A', '2B', 221486615, None, 0.129, 0.046, None, 'FarmCPU, BLINK, MLMM', 28),
    ('1161805', '2D', 650655958, None, 0.367, 0.000, 0.197, 'FarmCPU', 35),
    ('986698-65-G/T', '3A', 484644137, None, 0.056, 0.046, None, 'FarmCPU, BLINK', 28),
    ('1128927-47-C/T', '3B', 192028551, None, 0.075, 0.006, -0.266, 'FarmCPU', 35),
    ('3064908-13-G/C', '3B', 701828362, None, 0.249, 0.010, -0.159, 'FarmCPU', 35),
    ('2303262', '4B', 653653182, None, 0.335, 0.020, -0.132, 'FarmCPU', 35),
    ('3943527', '4B', 27519570, None, 0.263, 0.020, 0.161, 'FarmCPU', 35),
    ('1184257-26-C/G', '5A', 693326861, None, 0.321, 0.021, 0.150, 'FarmCPU', 35),
    ('3385268-44-A/C', '5B', 435443286, None, 0.445, 0.036, -0.111, 'FarmCPU', 35),
    ('1104744', '6A', 12977901, None, 0.364, 0.014, -0.128, 'FarmCPU', 35),
    ('1092672', '6B', 692728550, None, 0.347, 0.026, None, 'FarmCPU, BLINK', 35),
    ('1095048-53-G/A', '6D', 429945514, None, 0.460, 0.036, -0.117, 'FarmCPU', 35),
    ('1217909-44-G/A', '7A', 94707004, None, 0.101, 0.028, 0.187, 'FarmCPU', 35),
    ('3023224-23-C/T', '7A', 472488999, None, 0.173, 0.010, -0.220, 'FarmCPU', 35),
]:
    add(mk, 'MTA', WSMV, c, None, None, round(bp / 1e6, 4), 'RefSeq v1.0', mk, '',
        'FDR %.3f; MAF %.3f%s' % (fdr, maf, '' if eff is None else '; effect %.3f' % eff),
        'Enteshari2026BMC', 'B',
        'DArTseq marker-trait association for wheat streak mosaic virus resistance in a panel of '
        '175 wheat genotypes (%s; %d dpi).' % (model, dpi))

# ======================================================= 9. BYDV
BYDV = 'Barley yellow dwarf virus'
add('QByd.hwwg-2DL', 'QTL', BYDV, '2D', 642.42, 643.00, 642.71, 'RefSeq v2.1', '', 44.5,
    'LOD up to 21.55', 'Bian2025TPG', 'A',
    'Major QTL for barley yellow dwarf virus resistance on 2DL, JagMut1095 x Jagger RIL population; '
    'detected in three field seasons and on BLUEs (PVE 11.6-44.5%), physical interval on IWGSC '
    'RefSeq v2.1. The same interval also carries QTL for kernel traits in this population.')

# ======================================================= 10. WDV
WDV = 'Wheat dwarf virus'
for nm, c, s, e, lod, pve, pop, fl in [
    ('QTL_6A', '6A', 599.37, 607.34, 48.9, 73.9, 'MIxA40, MUxA40, P1314xA40, P1314xA39',
     'TGWA25K-TG0301 / Tdurum_contig75700_411 - AX-95197581 / AX-89595830'),
    ('QTL_1B', '1B', 1.4217, 458.69, 11.7, 18.5, 'P1314xA40, P1314xA39',
     'BS00050522_51 / TGWA25K-TG0025 - GENE-3626_308 / Tdurum_contig70304_781'),
]:
    add(nm, 'QTL', WDV, c, s, e, mid(s, e), 'RefSeq v1.0', fl, pve, 'LOD up to %.1f' % lod,
        'Buerstmayr2023TAG', 'A' if c == '6A' else 'D',
        'QTL for wheat dwarf virus resistance in four related winter wheat populations (%s), '
        'multiple QTL mapping across three years. Not formally designated in the source.%s' % (
            pop, '' if c == '6A' else ' The physical interval implied by the flanking markers spans '
            'most of the chromosome and should be treated as unreliable; the genetic interval is 0-14 cM.'))

xl2 = os.path.join(HOME, 'bio/dl/rest/x_PMC9047360/Data_Sheet_1.xlsx')
PF = []
from lxml import etree
t = etree.parse(os.path.join(HOME, 'bio/dl/rest/PMC9047360.xml'))
for tw in t.iter('table-wrap'):
    cap = ''.join(tw.find('caption').itertext()) if tw.find('caption') is not None else ''
    if 'Quantitative trait loci' not in cap: continue
    for tr in tw.iter('tr'):
        cells = [re.sub(r'\s+', '', ''.join(x.itertext())).strip() for x in tr]
        if len(cells) < 11 or not re.match(r'^WDV_', cells[0]): continue
        nm, mk, c, cm, ival, lod, r2 = cells[0], cells[1], cells[3], cells[4], cells[5], cells[7], cells[8]
        PF.append((nm, mk, c, cm, ival, lod, r2))
TRAIT = {'PH': 'plant height under WDV infection', 'Yield': 'grain yield under WDV infection',
         'E': 'extinction (ELISA)', 'IR': 'infection rate', 'TKW': 'thousand-kernel weight under WDV infection'}
for nm, mk, c, cm, ival, lod, r2 in PF:
    tr = nm.split('_')[1] if len(nm.split('_')) > 2 else ''
    pos, asm, msrc = place(mk, c)
    add(nm, 'MTA', WDV, c, None, None, None if pos is None else round(pos, 4),
        asm or 'RefSeq v1.0', mk, r2, 'LOD %s; interval %s cM' % (lod, ival),
        'Pfrieme2022FPS', 'B' if pos is not None else 'D',
        'Marker-trait association for wheat dwarf virus resistance/tolerance (%s) in a panel of '
        'winter wheat accessions (GWAS, FDR < 0.05). Physical position inferred from the peak '
        'marker by the atlas; the source reports genetic positions only.'
        % TRAIT.get(tr, 'trait ' + tr),
        'peak marker %s%s' % (mk, '; ' + msrc if msrc else ''))

# ======================================================= 11. WYMV
WYMV = 'Wheat yellow mosaic virus'
add('QTL_2A', 'QTL', WYMV, '2A', None, None, None, 'RefSeq v1.0', 'snp4212-snp4215', 43.8,
    'LOD up to 18.3', 'Kawaguchi2024BS', 'D',
    'Major QTL for wheat yellow mosaic virus resistance from the winter wheat Hokkai 240, '
    'Hokkai 240 x Nanbukomugi RIL population; effective against pathotypes II and III '
    '(PVE 23.2-43.8%). Genetic position 15.7-17.7 cM on 2A; the flanking SNP markers (snp4212, '
    'snp4215) are study-specific and could not be placed physically. Not formally designated '
    'in the source.')
add('Ym1', 'DG', WYMV, '2D', None, None, None, 'RefSeq v1.0', '', '', '',
    'Chen2025NC', 'B',
    'Ym1, a CC-NBS-LRR gene on 2DL conferring wheat yellow mosaic virus resistance by recognising '
    'the viral coat protein. The candidate gene TraesFLD2D01G565300 is present in Fielder but '
    'absent from Chinese Spring, whose orthologous 2DL interval carries the NLR cluster '
    'TraesCS2D01G510000-TraesCS2D01G511700; no Chinese Spring coordinate is therefore given.')

# ---------------------------------------------------------------- write
print('new records:', len(NEW))
from collections import Counter
print(Counter(r['disease'] for r in NEW))
print('unplaced:', sum(1 for r in NEW if r['peak'] is None and r['start'] is None))
if '--write' in sys.argv:
    shutil.copy(MASTER, BACKUP)
    reports.extend(NEW)
    json.dump(d, open(MASTER, 'w'))
    print('written; total reports', len(reports), 'sources', len(sources))
