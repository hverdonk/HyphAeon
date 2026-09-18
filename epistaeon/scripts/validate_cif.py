"""Inventory each mmCIF: entities, chains, organisms, ligands, numbering, gaps."""
import sys, glob, os
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

AA3to1 = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

def aslist(d, k):
    v = d.get(k)
    if v is None: return []
    return v if isinstance(v, list) else [v]

def report(path):
    d = MMCIF2Dict(path)
    eid = d.get('_entry.id', ['?'])[0]
    title = d.get('_struct.title', ['?'])[0]
    res = d.get('_refine.ls_d_res_high', ['n/a'])[0]
    print(f"\n{'='*78}\n{eid}  ({res} A)  {title}\n{'='*78}")

    # entities
    ent_id   = aslist(d,'_entity.id')
    ent_type = aslist(d,'_entity.type')
    ent_desc = aslist(d,'_entity.pdbx_description')
    ent_cnt  = aslist(d,'_entity.pdbx_number_of_molecules')
    src_ent  = aslist(d,'_entity_src_gen.entity_id')
    src_org  = aslist(d,'_entity_src_gen.pdbx_gene_src_scientific_name')
    src_gene = aslist(d,'_entity_src_gen.pdbx_gene_src_gene')
    nat_ent  = aslist(d,'_entity_src_nat.entity_id')
    nat_org  = aslist(d,'_entity_src_nat.pdbx_organism_scientific')
    syn_ent  = aslist(d,'_pdbx_entity_src_syn.entity_id')
    syn_org  = aslist(d,'_pdbx_entity_src_syn.organism_scientific')
    org = {}
    for e,o in zip(src_ent, src_org): org[e] = o
    for e,o in zip(nat_ent, nat_org): org.setdefault(e, o)
    for e,o in zip(syn_ent, syn_org): org.setdefault(e, o+" (synthetic)")
    gene = {e:g for e,g in zip(src_ent, src_gene)}

    print("ENTITIES:")
    for i,e in enumerate(ent_id):
        t = ent_type[i] if i < len(ent_type) else '?'
        if t == 'water': continue
        desc = ent_desc[i] if i < len(ent_desc) else '?'
        n = ent_cnt[i] if i < len(ent_cnt) else '?'
        extra = []
        if org.get(e): extra.append(f"org={org[e]}")
        if gene.get(e) and gene[e] not in ('?','.'): extra.append(f"gene={gene[e]}")
        print(f"  [{e}] {t:14s} x{n:<3s} {desc}" + (f"   ({'; '.join(extra)})" if extra else ""))

    # polymer chains: author numbering + sequence
    asym   = aslist(d,'_pdbx_poly_seq_scheme.pdb_strand_id')
    mon    = aslist(d,'_pdbx_poly_seq_scheme.mon_id')
    seqnum = aslist(d,'_pdbx_poly_seq_scheme.pdb_seq_num')
    # pdb_seq_num is populated even for residues without coordinates;
    # auth_seq_num is '?' for those, so it is the field that marks modelled residues.
    authnum = aslist(d,'_pdbx_poly_seq_scheme.auth_seq_num')
    entid  = aslist(d,'_pdbx_poly_seq_scheme.entity_id')
    chains = {}
    for ch, m, n, a, e in zip(asym, mon, seqnum, authnum, entid):
        chains.setdefault(ch, {'ent':e, 'res':[]})['res'].append((n if a not in ('?','.') else '?', m))
    print("POLYMER CHAINS (author numbering, '.' = unmodelled):")
    global CHAINS; CHAINS = {}
    for ch, info in chains.items():
        obs = [(n,m) for n,m in info['res'] if n not in ('?','.')]
        nums = [int(n) for n,_ in obs if n.lstrip('-').isdigit()]
        seq = ''.join(AA3to1.get(m,'X') for _,m in obs)
        CHAINS[(os.path.basename(path)[:4], ch)] = {int(n):m for n,m in obs if n.lstrip('-').isdigit()}
        gaps = []
        if nums:
            s = sorted(set(nums))
            for a,b in zip(s, s[1:]):
                if b != a+1: gaps.append(f"{a}->{b}")
        rng = f"{min(nums)}-{max(nums)}" if nums else "n/a"
        print(f"  chain {ch} (entity {info['ent']}): {len(obs)} modelled residues, author range {rng}"
              + (f", GAPS: {', '.join(gaps[:8])}{' ...' if len(gaps)>8 else ''}" if gaps else ", no internal gaps"))
        print(f"      {seq[:70]}{'...' if len(seq)>70 else ''}")

    # ligands
    het = {}
    hn = aslist(d,'_pdbx_nonpoly_scheme.mon_id')
    hc = aslist(d,'_pdbx_nonpoly_scheme.pdb_strand_id')
    for m,c in zip(hn,hc):
        if m == 'HOH': continue
        het.setdefault(m,set()).add(c)
    if het:
        print("LIGANDS/HETATM: " + ", ".join(f"{m}(chains {','.join(sorted(c))})" for m,c in sorted(het.items())))
    else:
        print("LIGANDS/HETATM: none")

for p in sorted(sys.argv[1:]):
    report(p)

