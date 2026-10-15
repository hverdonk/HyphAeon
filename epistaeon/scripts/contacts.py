import numpy as np, warnings
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
from _fetch import cif
P=MMCIFParser(QUIET=True)
b=_COORD
S={k:P.get_structure(k,cif(k))[0] for k in ['2HHB','1MBO','1U19','2Q1H','3RY9','3GN8','4P6X']}

def res(st,ch,num):
    c=S[st][ch]
    for r in c:
        if r.id[1]==num and r.id[0]==' ': return r
    return None
def hetres(st,ch,name):
    for r in S[st][ch]:
        if r.resname==name: return r
    return None
def mind(r1,r2,heavy=True):
    a1=np.array([a.coord for a in r1 if not heavy or a.element!='H'])
    a2=np.array([a.coord for a in r2 if not heavy or a.element!='H'])
    return float(np.min(np.linalg.norm(a1[:,None,:]-a2[None,:,:],axis=-1)))
def cbd(r1,r2):
    def rep(r):
        if 'CB' in r: return r['CB'].coord
        return r['CA'].coord
    return float(np.linalg.norm(rep(r1)-rep(r2)))
def iface(st,chA,chB,cut=5.0):
    A=[r for r in S[st][chA] if r.id[0]==' ']; B=[r for r in S[st][chB] if r.id[0]==' ']
    n=0; pairs=[]
    for r1 in A:
        c1=np.array([a.coord for a in r1])
        for r2 in B:
            c2=np.array([a.coord for a in r2])
            if np.min(np.linalg.norm(c1[:,None]-c2[None],axis=-1))<=cut:
                n+=1; pairs.append((r1.id[1],r2.id[1]))
    return n,pairs

print("="*78); print("§4 TABLE 2 CLAIM VALIDATION"); print("="*78)
print("\n-- Hemoglobin 2HHB (A/C = alpha, B/D = beta) --")
n11,_=iface('2HHB','A','B'); print(f"  a1b1 packing face, residue pairs <=5.0A : {n11}    [paper claims 55]")
n12,_=iface('2HHB','A','D'); print(f"  a1b2 sliding switch, residue pairs <=5.0A: {n12}    [paper claims 37]")
h87=res('2HHB','A',87); b51=res('2HHB','B',51)
print(f"  alpha-His87 <-> beta-site51 : min-heavy {mind(h87,b51):.1f}A, CB-CB {cbd(h87,b51):.1f}A   [paper claims 18.4A]")
# also to the partner beta chain D
b51d=res('2HHB','D',51)
print(f"  alpha-His87(A) <-> beta51(D): min-heavy {mind(h87,b51d):.1f}A, CB-CB {cbd(h87,b51d):.1f}A")
k40=res('2HHB','A',40); h146d=res('2HHB','D',146); h146b=res('2HHB','B',146)
print(f"  Perutz a40Lys(A)<->b146His(D): min-heavy {mind(k40,h146d):.1f}A   (salt bridge if <4A)")
print(f"  Perutz a40Lys(A)<->b146His(B): min-heavy {mind(k40,h146b):.1f}A")
# epistatic pair named in Table 1
a114=res('2HHB','A',114); b116=res('2HHB','B',116); b116d=res('2HHB','D',116)
print(f"  Table1 pair a114<->b116 (B): min-heavy {mind(a114,b116):.1f}A ; (D): {mind(a114,b116d):.1f}A")
a57=res('2HHB','A',57); a64=res('2HHB','A',64)
print(f"  Case1 drivers a57<->a64: min-heavy {mind(a57,a64):.1f}A, CB-CB {cbd(a57,a64):.1f}A")
for lbl,r in [('a57',a57),('a64',a64)]:
    dmin=min(mind(r,x) for x in S['2HHB']['B'] if x.id[0]==' ')
    print(f"    {lbl} min distance to any beta(B) residue: {dmin:.1f}A")

print("\n-- Steroid receptor: AncGR2 3GN8 (Anc numbering) --")
s106=res('3GN8','A',106); l111=res('3GN8','A',111); y27=res('3GN8','A',27); p36=res('3GN8','A',36); l29=res('3GN8','A',29)
dex=hetres('3GN8','A','DEX')
print(f"  pos106={s106.resname} pos111={l111.resname} pos27={y27.resname} pos29={l29.resname} pos36={p36.resname}")
for lbl,r in [('106',s106),('111',l111),('27',y27),('29',l29),('36',p36)]:
    print(f"    residue {lbl:>3s} ({r.resname}) -> dexamethasone: min-heavy {mind(r,dex):.1f}A")
print(f"  pocket(106) <-> remote 27: CB-CB {cbd(s106,y27):.1f}A ; min-heavy {mind(s106,y27):.1f}A   [paper claims 15.2A / >15A]")
print(f"  pocket(111) <-> remote 27: CB-CB {cbd(l111,y27):.1f}A ; min-heavy {mind(l111,y27):.1f}A")
print(f"  pocket(106) <-> remote 29: CB-CB {cbd(s106,l29):.1f}A")
print(f"  pocket(111) <-> remote 29: CB-CB {cbd(l111,l29):.1f}A")
print(f"  pocket(106) <-> pos 36   : CB-CB {cbd(s106,p36):.1f}A")

# The SWS1 sites use the original paper's conventional bovine numbering.
# These coordinates are homolog-template geometry, not measured SWS1 contacts.
import json
sws1 = json.loads(open(os.path.join(_HERE, '..', 'data', 'numbering_offsets.json')).read())["sws1_benchmark"]
print("\n-- SWS1 sites on the bovine 1U19 homolog template --")
ret = hetres('1U19', 'A', 'RET')
for site in sws1['site_map']:
    r = res(site['template_pdb'], site['template_chain'], site['template_author_position'])
    print(f"  {site['substitution']}: human CDS {site['human_cds_position_1based']}, "
          f"codon column {site['alignment_codon_column_1based']}, "
          f"template {r.resname}{site['template_author_position']}; "
          f"template retinal distance {mind(r, ret):.1f}A")
print("  Template distances do not establish SWS1 functional coupling or dark-noise kinetics.")

print("\n-- Myoglobin 1MBO --")
mb=[r for r in S['1MBO']['A'] if r.id[0]==' ']
hem=hetres('1MBO','A','HEM')
charged={'ARG','LYS','HIS','ASP','GLU'}
net=sum(1 for r in mb if r.resname in ('ARG','LYS'))-sum(1 for r in mb if r.resname in ('ASP','GLU'))
print(f"  {len(mb)} residues; Arg+Lys={sum(1 for r in mb if r.resname in ('ARG','LYS'))}, Asp+Glu={sum(1 for r in mb if r.resname in ('ASP','GLU'))}, crude net charge (no His) = {net:+d}")
print(f"  heme present: {hem.resname if hem else 'NONE'}; proximal His93 -> heme Fe: {mind(res('1MBO','A',93),hem):.1f}A")
