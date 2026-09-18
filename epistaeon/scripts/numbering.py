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
    # pdb_seq_num is filled even for residues without coordinates; auth_seq_num is '?' for those
    for c,m,n,a in zip(L(d,'_pdbx_poly_seq_scheme.pdb_strand_id'),L(d,'_pdbx_poly_seq_scheme.mon_id'),
                       L(d,'_pdbx_poly_seq_scheme.pdb_seq_num'),L(d,'_pdbx_poly_seq_scheme.auth_seq_num')):
        if c!=ch or m not in AA3to1: continue
        seqres.append(AA3to1[m])
        obs.append((int(n), a not in ('?','.')))
    return ''.join(seqres),obs
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score=-11; al.extend_gap_score=-1; al.mode='global'
# free end gaps, so a cleaved initiator Met is an end gap rather than forcing a spurious internal gap
al.end_gap_score=0
print("="*84)
print("NUMBERING CROSSWALK: PDB author numbering  vs  UniProt canonical numbering")
print("="*84)
CHECKS=[('2HHB','A','P69905','Hb alpha'),('2HHB','B','P68871','Hb beta'),
        ('1MBO','A','P02185','Myoglobin, sperm whale'),('1U19','A','P02699','Rhodopsin, bovine'),
        ('7PRX','A','P04150','Glucocorticoid receptor LBD'),
        # non-human structures against the HUMAN protein (the TOGA2 reference species)
        ('1MBO','A','P02144','Myoglobin vs HUMAN'),('1U19','A','P08100','Rhodopsin vs HUMAN')]
for pdb,ch,acc,label in CHECKS:
    ps,pn=pdbchain(pdb,ch); us=upseq(acc)
    aln=al.align(ps,us)[0]; iP,iU=aln.indices
    offs={}; mism=[]; unmod_mism=[]
    cols=list(zip(iP,iU))
    paired=[k for k,(a,b) in enumerate(cols) if a>=0 and b>=0]
    # a gap column is internal only if residue pairs exist on both sides of it
    internal_gaps=[f"UP {b+1}" if b>=0 else f"PDB {pn[a][0]}" for k,(a,b) in enumerate(cols)
                   if (a<0 or b<0) and paired and paired[0]<k<paired[-1]]
    for a,b in cols:
        if a<0 or b<0: continue
        n,modelled=pn[a]
        offs[(b+1)-n]=offs.get((b+1)-n,0)+1
        if ps[a]!=us[b]:
            (mism if modelled else unmod_mism).append(f"PDB {n}{ps[a]}/UP {b+1}{us[b]}")
    off=sorted(offs.items(),key=lambda x:-x[1])
    print(f"\n{pdb} chain {ch}  ({label}, {acc})")
    print(f"   PDB SEQRES {len(ps)} aa | UniProt {len(us)} aa")
    print(f"   offset UniProt = PDB + {off[0][0]}   ({off[0][1]}/{sum(offs.values())} residues consistent)"
          + (f"   OTHER OFFSETS: {off[1:]}" if len(off)>1 else ""))
    print(f"   sequence mismatches vs UniProt (modelled residues): {len(mism)}" + (f"  -> {', '.join(mism[:6])}" if mism else ""))
    if unmod_mism: print(f"   (mismatches at residues without coordinates, e.g. construct tags: {', '.join(unmod_mism)})")
    print(f"   internal indels vs UniProt: {len(internal_gaps)}" + (f" -> {internal_gaps[:6]}" if internal_gaps else ""))
    unmod=[n for n,ok in pn if not ok]
    print(f"   unmodelled SEQRES positions: {len(unmod)}" + (f" -> {unmod[:12]}" if unmod else " (none)"))

# The model's positions are TOGA2 alignment columns in human (hg38) coordinates.
# For the non-human structures, confirm the structure species' row has no indels
# relative to hg38, so the human-residue mapping carries over unchanged.
print("\n" + "="*84)
print("TOGA2 ALIGNMENT: structure species vs hg38 (indels would break the mapping)")
print("="*84)
AL=os.path.join(_HERE,'..','data','alignments','toga2')+os.sep
if os.path.exists(AL+'assemblies_and_species.tsv'):
    tab={}
    for line in open(AL+'assemblies_and_species.tsv').read().splitlines()[1:]:
        f=line.split('\t'); tab.setdefault(f[1],[]).append(f[6])
    def read(p):
        s={};k=None
        for line in open(p):
            line=line.strip()
            if line.startswith('>'): k=line[1:]; s[k]=''
            else: s[k]+=line
        return s
    for gene,species in [('MB','Physeter macrocephalus'),('RHO','Bos taurus')]:
        seqs=read(AL+gene+'.codonified.fa'); h=seqs['hg38']
        for asm in [x for x in tab.get(species,[]) if x in seqs]:
            q=seqs[asm]
            hc=[h[i:i+3] for i in range(0,len(h),3)]; qc=[q[i:i+3] for i in range(0,len(q),3)]
            h_only=sum(1 for x,y in zip(hc,qc) if x!='---' and y=='---')
            q_only=sum(1 for x,y in zip(hc,qc) if x=='---' and y!='---')
            print(f"   {gene}: {species} ({asm}): {h_only} human residues without a counterpart, {q_only} residues in human-gap columns")
else:
    print("   (data/alignments/toga2 not present; skipped)")
