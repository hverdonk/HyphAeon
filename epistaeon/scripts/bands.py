import numpy as np, warnings; warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
from scipy.stats import fisher_exact
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
P=MMCIFParser(QUIET=True); CO=_COORD
STD=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
def rl(m,ch): return [r for r in m[ch] if r.id[0]==' ' and r.resname in STD]
def dm(rs):
    C=[np.array([a.coord for a in r]) for r in rs]; n=len(C); D=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            D[i,j]=D[j,i]=np.min(np.linalg.norm(C[i][:,None]-C[j][None],axis=-1))
    return D
print("="*86)
print("E.  COVERAGE OF THE PAPER'S OWN BANDS  (<=5A direct; 10-25A allosteric; everything else undefined)")
print("="*86)
tots={}
for pdb,ch in [('2HHB','A'),('2HHB','B'),('1MBO','A'),('1U19','A')]:
    m=P.get_structure(pdb,CO+pdb+'.cif')[0]; rs=rl(m,ch); nums=np.array([r.id[1] for r in rs])
    D=dm(rs); iu=np.triu_indices(len(rs),k=1)
    sep=np.abs(nums[iu[0]]-nums[iu[1]]); v=D[iu][sep>=5]; tots[(pdb,ch)]=v
    direct=(v<=5).mean(); gap=((v>5)&(v<10)).mean(); allo=((v>=10)&(v<=25)).mean(); far=(v>25).mean()
    print(f"  {pdb} {ch}: direct {direct*100:5.1f}% | UNDEFINED 5-10A {gap*100:5.1f}% | allosteric 10-25A {allo*100:5.1f}% | UNDEFINED >25A {far*100:5.1f}%")
    print(f"        => {100*(gap+far):.1f}% of pairs fall outside both defined bands")
print("\n"+"="*86)
print("F.  POWER FOR THE §5 BENCHMARK: 'Allosteric Interface Specificity, OR>4.0, Fisher p<1e-5'")
print("="*86)
for (pdb,ch),v in tots.items():
    bg=(v<=5).mean()   # background contact rate
    print(f"  {pdb} {ch}: background contact rate (<=5A) = {bg*100:.2f}%")
    for OR in (4.0,):
        # enriched rate implied by OR against background odds
        o=bg/(1-bg); oe=OR*o; pe=oe/(1+oe)
        found=None
        for n in range(5,4001,1):
            k=int(round(pe*n))
            if k<1: continue
            tbl=[[k,n-k],[int(round(bg*20000)),20000-int(round(bg*20000))]]
            try:
                if fisher_exact(tbl,alternative='greater')[1]<1e-5: found=(n,k,pe); break
            except Exception: pass
        if found:
            n,k,pe=found
            print(f"      OR=4 implies {pe*100:.1f}% contact rate among detected pairs;"
                  f" need >= {n} detected pairs ({k} of them true contacts) for Fisher p<1e-5")
