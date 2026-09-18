from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Align import PairwiseAligner, substitution_matrices
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
from _fetch import cif
AA3to1={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}
def aslist(d,k):
    v=d.get(k); return [] if v is None else (v if isinstance(v,list) else [v])
def seqnum(p,chain='A'):
    d=MMCIF2Dict(p); s=[];n=[]
    for c,m,num in zip(aslist(d,'_pdbx_poly_seq_scheme.pdb_strand_id'),aslist(d,'_pdbx_poly_seq_scheme.mon_id'),aslist(d,'_pdbx_poly_seq_scheme.pdb_seq_num')):
        if c!=chain or num in('?','.') or not num.lstrip('-').isdigit(): continue
        s.append(AA3to1.get(m,'X')); n.append(int(num))
    return ''.join(s), n

al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score=-11; al.extend_gap_score=-1; al.mode='global'

def diff(pA,pB,nameA,nameB,chA='A',chB='A'):
    sA,nA=seqnum(pA,chA); sB,nB=seqnum(pB,chB)
    aln=al.align(sA,sB)[0]
    iA,iB=aln.indices
    subs=[];ident=0;algn=0;gaps=0
    for a,b in zip(iA,iB):
        if a<0 or b<0: gaps+=1; continue
        algn+=1
        if sA[a]==sB[b]: ident+=1
        else: subs.append((nA[a],sA[a],nB[b],sB[b]))
    print(f"\n{nameA} -> {nameB}: {algn} aligned cols, {ident} identical ({100*ident/algn:.1f}%), {gaps} gap cols, {len(subs)} substitutions")
    off=set(x[2]-x[0] for x in subs)
    print(f"  numbering offsets seen among substitutions: {sorted(off)[:6]}")
    return subs,nA,nB,sA,sB

b=_COORD
print("="*78);print("STEROID RECEPTOR EPISTATIC LADDER (alignment-based)");print("="*78)
s1,_,_,_,_ = diff(cif('2Q1H'),cif('3RY9'),'AncCR (2Q1H)','AncGR1 (3RY9)')
print('  '+', '.join(f"{a}{na}{bb}" for na,a,nb,bb in s1))
s2,_,_,_,_ = diff(cif('3RY9'),cif('3GN8'),'AncGR1 (3RY9)','AncGR2 (3GN8)')
print('  '+', '.join(f"{a}{na}{bb}" for na,a,nb,bb in s2))
s3,_,_,_,_ = diff(cif('2Q1H'),cif('3GN8'),'AncCR (2Q1H)','AncGR2 (3GN8)')
print('  '+', '.join(f"{a}{na}{bb}" for na,a,nb,bb in s3))

print("\n"+"="*78);print("ANC-NUMBERING  <->  HUMAN GR (NR3C1) CROSSWALK via 4P6X");print("="*78)
s4,nA,nB,sA,sB = diff(cif('3GN8'),cif('4P6X'),'AncGR2 (3GN8)','human GR (4P6X)')
# report offset at the four white-paper positions
sA_,nA_=seqnum(cif('3GN8')); sB_,nB_=seqnum(cif('4P6X'))
aln=al.align(sA_,sB_)[0]; iA,iB=aln.indices
m={}
for a,bx in zip(iA,iB):
    if a>=0 and bx>=0: m[nA_[a]]=(nB_[bx], sA_[a], sB_[bx])
for pos in (27,36,106,111):
    if pos in m:
        h,aa_anc,aa_hum=m[pos]
        print(f"  Anc position {pos:3d} ({aa_anc}) == human GR {h} ({aa_hum})   [offset +{h-pos}]")
    else:
        print(f"  Anc position {pos:3d}: no aligned human GR counterpart")
