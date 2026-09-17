from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Align import PairwiseAligner, substitution_matrices
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
AA3to1={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}
def L(d,k):
    v=d.get(k); return [] if v is None else (v if isinstance(v,list) else [v])
CO=_COORD
UP=os.path.join(_HERE,'..','data','uniprot')+os.sep
def upseq(a):
    s=open(UP+a+'.fasta').read().split('\n'); return ''.join(s[1:])
def pdbchain(pdb,ch):
    d=MMCIF2Dict(CO+pdb+'.cif'); seqres=[];obs=[]
    for c,m,n in zip(L(d,'_pdbx_poly_seq_scheme.pdb_strand_id'),L(d,'_pdbx_poly_seq_scheme.mon_id'),L(d,'_pdbx_poly_seq_scheme.pdb_seq_num')):
        if c!=ch: continue
        seqres.append(AA3to1.get(m,'X'))
        obs.append(None if n in('?','.') else int(n))
    return ''.join(seqres),obs
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score=-11; al.extend_gap_score=-1; al.mode='global'
print("="*84)
print("NUMBERING CROSSWALK: PDB author numbering  vs  UniProt canonical numbering")
print("="*84)
for pdb,ch,acc,label in [('2HHB','A','P69905','Hb alpha'),('2HHB','B','P68871','Hb beta'),
                         ('1MBO','A','P02185','Myoglobin'),('1U19','A','P02699','Rhodopsin')]:
    ps,pn=pdbchain(pdb,ch); us=upseq(acc)
    aln=al.align(ps,us)[0]; iP,iU=aln.indices
    offs={}; mism=[]
    for a,b in zip(iP,iU):
        if a<0 or b<0: continue
        if pn[a] is None: continue
        offs[(b+1)-pn[a]]=offs.get((b+1)-pn[a],0)+1
        if ps[a]!=us[b]: mism.append(f"PDB {pn[a]}{ps[a]}/UP {b+1}{us[b]}")
    off=sorted(offs.items(),key=lambda x:-x[1])
    print(f"\n{pdb} chain {ch}  ({label}, {acc})")
    print(f"   PDB SEQRES {len(ps)} aa | UniProt {len(us)} aa")
    print(f"   offset UniProt = PDB + {off[0][0]}   ({off[0][1]}/{sum(offs.values())} residues consistent)"
          + (f"   OTHER OFFSETS: {off[1:]}" if len(off)>1 else ""))
    print(f"   sequence mismatches vs UniProt: {len(mism)}" + (f"  -> {', '.join(mism[:6])}" if mism else ""))
    unmod=[i+1 for i,n in enumerate(pn) if n is None]
    print(f"   unmodelled SEQRES positions: {len(unmod)}" + (f" -> {unmod[:12]}" if unmod else " (none)"))
