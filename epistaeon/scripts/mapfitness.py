import numpy as np, warnings, itertools; warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
P=MMCIFParser(QUIET=True); CO=_COORD
STD=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
def model(p): return P.get_structure(p,CO+p+'.cif')[0]
def reslist(m,ch):
    return [r for r in m[ch] if r.id[0]==' ' and r.resname in STD]
def coords(rs): return [np.array([a.coord for a in r]) for r in rs]
def dmat(cA,cB):
    n,m=len(cA),len(cB); D=np.zeros((n,m))
    for i,a in enumerate(cA):
        for j,b in enumerate(cB):
            D[i,j]=np.min(np.linalg.norm(a[:,None]-b[None],axis=-1))
    return D
def classes(D, selfpair=False, seqsep=None, nums=None):
    if selfpair:
        iu=np.triu_indices(D.shape[0],k=1); v=D[iu]
        if seqsep is not None and nums is not None:
            sep=np.abs(np.array(nums)[iu[0]]-np.array(nums)[iu[1]]); v=v[sep>=seqsep]
    else: v=D.ravel()
    tot=len(v)
    return tot, (v<=5).sum()/tot, ((v>5)&(v<=10)).sum()/tot, (v>10).sum()/tot, v

print("="*86)
print("A.  DISTANCE-CLASS COMPOSITION  (can a detected pair be separated into contact vs allosteric?)")
print("    direct <=5A | ambiguous 5-10A | allosteric >10A   [|i-j|>=5 to drop trivial backbone neighbours]")
print("="*86)
store={}
for pdb,chs in [('2HHB',['A','B']),('1MBO',['A']),('1U19',['A'])]:
    m=model(pdb)
    for ch in chs:
        rs=reslist(m,ch); nums=[r.id[1] for r in rs]; D=dmat(coords(rs),coords(rs))
        store[(pdb,ch)]=(rs,nums,D)
        tot,a,b,c,v=classes(D,selfpair=True,seqsep=5,nums=nums)
        print(f"  {pdb} chain {ch}: {len(rs):3d} residues, {tot:6d} pairs |  <=5A {a*100:5.1f}%  5-10A {b*100:5.1f}%  >10A {c*100:5.1f}%  | max {v.max():.0f}A")
# hemoglobin inter-chain
m=model('2HHB')
rsA=reslist(m,'A'); rsB=reslist(m,'B'); rsC=reslist(m,'C'); rsD=reslist(m,'D')
Dab=dmat(coords(rsA),coords(rsB)); Dad=dmat(coords(rsA),coords(rsD))
for lbl,D in [('alpha1-beta1 (A-B)',Dab),('alpha1-beta2 (A-D)',Dad)]:
    tot,a,b,c,v=classes(D)
    print(f"  2HHB {lbl}: {tot:6d} pairs |  <=5A {a*100:5.1f}%  5-10A {b*100:5.1f}%  >10A {c*100:5.1f}%  | max {v.max():.0f}A")

print("\n"+"="*86)
print("B.  CRYSTALLOGRAPHIC-COPY STABILITY  (does the class survive choosing a different copy?)")
print("="*86)
def copycmp(pdb,ch1,ch2,label):
    m=model(pdb); r1=reslist(m,ch1); r2=reslist(m,ch2)
    n1={r.id[1] for r in r1}; n2={r.id[1] for r in r2}; common=sorted(n1&n2)
    f1=[r for r in r1 if r.id[1] in common]; f2=[r for r in r2 if r.id[1] in common]
    D1=dmat(coords(f1),coords(f1)); D2=dmat(coords(f2),coords(f2))
    iu=np.triu_indices(len(common),k=1)
    sep=np.abs(np.array(common)[iu[0]]-np.array(common)[iu[1]])
    k=sep>=5
    v1=D1[iu][k]; v2=D2[iu][k]
    def cls(v): return np.where(v<=5,0,np.where(v<=10,1,2))
    c1,c2=cls(v1),cls(v2)
    flip=(c1!=c2).mean()
    hard=((c1==0)&(c2==2))|((c1==2)&(c2==0))
    print(f"  {label}: {k.sum()} pairs | median |Δd| {np.median(np.abs(v1-v2)):.2f}A, 95th pct {np.percentile(np.abs(v1-v2),95):.2f}A")
    print(f"      class disagreement between copies: {flip*100:.2f}%   contact<->allosteric flips: {hard.sum()} ({hard.mean()*100:.3f}%)")
copycmp('2HHB','A','C','2HHB alpha copies A vs C')
copycmp('2HHB','B','D','2HHB beta  copies B vs D')
copycmp('1U19','A','B','1U19 rhodopsin copies A vs B')

print("\n"+"="*86)
print("C.  HEMOGLOBIN INTER-CHAIN AMBIGUITY  (an alpha-beta pair has two inequivalent realisations)")
print("="*86)
cls=lambda v: np.where(v<=5,0,np.where(v<=10,1,2))
c_ab, c_ad = cls(Dab), cls(Dad)
dis=(c_ab!=c_ad).mean(); hard=((c_ab==0)&(c_ad==2))|((c_ab==2)&(c_ad==0))
print(f"  alpha-beta pairs: {Dab.size}")
print(f"  class differs between alpha1beta1 and alpha1beta2: {dis*100:.1f}%")
print(f"  outright contact<->allosteric flips: {hard.sum()} pairs ({hard.mean()*100:.1f}%)")
print(f"  pairs that are a contact in AT LEAST one interface: {((c_ab==0)|(c_ad==0)).sum()}")
print(f"  pairs that are a contact in BOTH interfaces:        {((c_ab==0)&(c_ad==0)).sum()}")
print("  => an alpha-beta pair CANNOT be classified without first fixing the interface convention.")

print("\n"+"="*86)
print("D.  COFACTOR-MEDIATED COUPLING  (pairs far apart in residue space but bridged by haem/retinal)")
print("="*86)
def cofactor(pdb,ch,name):
    m=model(pdb); rs=reslist(m,ch); nums=[r.id[1] for r in rs]
    het=[r for r in m[ch] if r.resname==name]
    if not het: print(f"  {pdb} {ch}: {name} absent"); return
    hc=np.array([a.coord for a in het[0]])
    dh=np.array([np.min(np.linalg.norm(np.array([a.coord for a in r])[:,None]-hc[None],axis=-1)) for r in rs])
    D=dmat(coords(rs),coords(rs)); iu=np.triu_indices(len(rs),k=1)
    sep=np.abs(np.array(nums)[iu[0]]-np.array(nums)[iu[1]]); k=sep>=5
    v=D[iu][k]; both=(dh[iu[0]]<=5)&(dh[iu[1]]<=5); both=both[k]
    lig=(both&(v>10)).sum()
    print(f"  {pdb} {ch} / {name}: {(dh<=5).sum()} residues line the cofactor (<=5A)")
    print(f"      pairs BOTH lining cofactor yet >10A apart (would be mis-called allosteric): {lig}")
    print(f"      ... of {both.sum()} cofactor-lining pairs total ({100*lig/max(1,both.sum()):.0f}%)")
cofactor('2HHB','A','HEM'); cofactor('1MBO','A','HEM'); cofactor('1U19','A','RET')
