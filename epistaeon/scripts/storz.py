import numpy as np, warnings; warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
P=MMCIFParser(QUIET=True)
m=P.get_structure('x',_COORD + '2HHB.cif')[0]
def res(ch,n):
    for r in m[ch]:
        if r.id[1]==n and r.id[0]==' ': return r
def mind(r1,r2):
    a=np.array([x.coord for x in r1]); b=np.array([x.coord for x in r2])
    return float(np.min(np.linalg.norm(a[:,None]-b[None],axis=-1)))
A=[50,57,60,64,71]; B=[62,72,128,135]
print("Storz 2009 deer-mouse ALPHA sites in human 2HHB chain A:")
for n in A:
    r=res('A',n)
    dB=min(mind(r,x) for x in m['B'] if x.id[0]==' ')
    dD=min(mind(r,x) for x in m['D'] if x.id[0]==' ')
    hem=[x for x in m['A'] if x.resname=='HEM'][0]
    print(f"  a{n:<4d} {r.resname:4s}  ->nearest beta1(B) {dB:5.1f}A  ->nearest beta2(D) {dD:5.1f}A  ->own heme {mind(r,hem):5.1f}A")
print("Storz 2009 deer-mouse BETA sites in human 2HHB chain B:")
for n in B:
    r=res('B',n)
    dA=min(mind(r,x) for x in m['A'] if x.id[0]==' ')
    hem=[x for x in m['B'] if x.resname=='HEM'][0]
    print(f"  b{n:<4d} {r.resname:4s}  ->nearest alpha1(A) {dA:5.1f}A  ->own heme {mind(r,hem):5.1f}A")
print("\nPairwise min-heavy distance matrix among the 9 Storz sites (A-chain alpha, B-chain beta):")
sites=[('A',n) for n in A]+[('B',n) for n in B]
lab=[f"a{n}" for n in A]+[f"b{n}" for n in B]
print("        "+"".join(f"{l:>7s}" for l in lab))
for i,(c1,n1) in enumerate(sites):
    row=f"  {lab[i]:<6s}"
    for j,(c2,n2) in enumerate(sites):
        row += "      ." if i==j else f"{mind(res(c1,n1),res(c2,n2)):7.1f}"
    print(row)
