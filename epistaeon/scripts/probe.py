import sys
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
from _fetch import cif
AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}
def aslist(d,k):
    v=d.get(k); return [] if v is None else (v if isinstance(v,list) else [v])
def chainmap(path):
    d=MMCIF2Dict(path)
    ch=aslist(d,'_pdbx_poly_seq_scheme.pdb_strand_id'); mon=aslist(d,'_pdbx_poly_seq_scheme.mon_id')
    num=aslist(d,'_pdbx_poly_seq_scheme.pdb_seq_num')
    out={}
    for c,m,n in zip(ch,mon,num):
        if n in ('?','.') or not n.lstrip('-').isdigit(): continue
        out.setdefault(c,{})[int(n)]=AA3to1.get(m,m)
    return out
def probe(path, chain, positions, label):
    cm=chainmap(path)
    if chain not in cm:
        print(f"  {label}: chain {chain} ABSENT"); return
    got=[]
    for p in positions:
        got.append(f"{p}={cm[chain].get(p,'--MISSING--')}")
    print(f"  {label} [chain {chain}]: " + "  ".join(got))

base=_COORD
print("CASE 1 HEMOGLOBIN 2HHB -- white paper names a57, a64, a114, b116, a40Lys, b146His, a-His87, b-site51")
probe(base+'2HHB.cif','A',[40,57,64,87,114],'alpha')
probe(base+'2HHB.cif','B',[51,116,146],'beta ')

print("\nCASE 2 MYOGLOBIN 1MBO -- B-helix ~20-35, E-helix ~58-77; surface charge sites")
probe(base+'1MBO.cif','A',[20,27,35,58,66,77,93],'Mb   ')

print("\nCASE 3 STEROID RECEPTOR -- white paper names Ser106, Leu111, Tyr27, Thr36")
for pdb,ch,lab in [('2Q1H','A','AncCR (2Q1H)'),('3RY9','A','AncGR1(3RY9)'),('3GN8','A','AncGR2(3GN8)'),(cif('4P6X'),'A','humGR (4P6X)')]:
    probe(pdb if pdb.endswith('.cif') else cif(pdb),ch,[27,36,106,111],lab)
print("  -- same positions in human-GR numbering (verified offset +531, see align_ladder.py):")
probe(cif('4P6X'),'A',[27+531,36+531,106+531,111+531],'humGR+531')

import json
sws1 = json.loads(open(os.path.join(_HERE, '..', 'data', 'numbering_offsets.json')).read())["sws1_benchmark"]
print("\nCASE 4 SWS1 -- explicitly mapped bovine 1U19 homolog-template sites")
probe(base+'1U19.cif', 'A', [site['template_author_position'] for site in sws1['site_map']], 'SWS1 homolog template')
