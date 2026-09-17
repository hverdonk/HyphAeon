"""Blocker 2 evidence: hemoglobin interface identity and T/R state dependence.

Shows that (a) alpha-beta contact pairs are disjoint between the a1b1 packing
interface and the a1b2 sliding interface, and (b) the sliding interface turns
over half its contacts between the deoxy (T) and oxy (R) states.

2DN1's asymmetric unit is only an alpha-beta dimer, so its biological assembly
operator must be applied to recover the tetramer -- itself a worked example of
why contacts must be computed on the assembly, not the asymmetric unit.
"""
import os, warnings, numpy as np
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

_HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, _HERE)
from _fetch import cif, COORD

STD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER '
          'THR TRP TYR VAL'.split())
P = MMCIFParser(QUIET=True)


def _l(d, k):
    v = d.get(k)
    return [] if v is None else (v if isinstance(v, list) else [v])


def residues(model, ch):
    return [r for r in model[ch] if r.id[0] == ' ' and r.resname in STD]


def coords(rs):
    return [np.array([a.coord for a in r]) for r in rs]


def contacts(X, rx, Y, ry, cut=5.0):
    out = set()
    for i, a in enumerate(X):
        for j, b in enumerate(Y):
            if np.min(np.linalg.norm(a[:, None] - b[None], axis=-1)) <= cut:
                out.add((rx[i].id[1], ry[j].id[1]))
    return out


def operators(path):
    d = MMCIF2Dict(path)
    ops = []
    for i, oid in enumerate(_l(d, '_pdbx_struct_oper_list.id')):
        R = np.array([[float(_l(d, f'_pdbx_struct_oper_list.matrix[{r}][{c}]')[i])
                       for c in (1, 2, 3)] for r in (1, 2, 3)])
        t = np.array([float(_l(d, f'_pdbx_struct_oper_list.vector[{r}]')[i])
                      for r in (1, 2, 3)])
        ops.append((oid, R, t))
    return ops


def main():
    T = P.get_structure('T', COORD + '2HHB.cif')[0]
    tA, tB, tD = residues(T, 'A'), residues(T, 'B'), residues(T, 'D')
    t11 = contacts(coords(tA), tA, coords(tB), tB)
    t12 = contacts(coords(tA), tA, coords(tD), tD)

    print('2HHB (deoxy, T state)')
    print(f'  alpha1beta1 packing : {len(t11)} pairs <=5 A')
    print(f'  alpha1beta2 sliding : {len(t12)} pairs <=5 A')
    print(f'  pairs contacting at BOTH interfaces: {len(t11 & t12)}')
    print('  => the two contact sets are disjoint; the interface convention '
          'decides every inter-chain call.\n')

    p = cif('2DN1')
    R = P.get_structure('R', p)[0]
    rA, rB = residues(R, 'A'), residues(R, 'B')
    CA, CB = coords(rA), coords(rB)
    got = []
    for oid, rot, tr in operators(p):
        got.append((oid, contacts(CA, rA, [(c @ rot.T) + tr for c in CB], rB)))
    got.sort(key=lambda kv: -len(kv[1]))
    r11, r12 = got[0][1], got[1][1]

    print('2DN1 (oxy, R state; tetramer rebuilt from assembly operators)')
    for label, tset, rset in (('alpha1beta1 packing', t11, r11),
                              ('alpha1beta2 sliding', t12, r12)):
        j = len(tset & rset) / max(1, len(tset | rset))
        print(f'  {label:20s} T={len(tset):3d} R={len(rset):3d} '
              f'shared={len(tset & rset):3d} T-only={len(tset - rset):3d} '
              f'R-only={len(rset - tset):3d} Jaccard={j:.2f}')
    print('  => the packing interface is state-invariant; the sliding '
          'interface turns over ~half its contacts.')


if __name__ == '__main__':
    main()
