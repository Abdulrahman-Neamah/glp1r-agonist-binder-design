"""
GLP-1R groove analysis and design-target preparation (Anthropic x Adaptyv Challenge 3).

Usage:   python 01_glp1r_groove_analysis.py 6X18.pdb [OUT_DIR]
Input:   PDB 6X18 (GLP-1 bound to GLP-1R; chain R = receptor, chain P = GLP-1). Download from files.rcsb.org/download/6X18.pdb
Output (numpy/pandas only):
  glp1r_groove_residues.csv      ECD residues contacted by the peptide (4.5 A), with exposure on the free receptor
  glp1r_patch_map.csv            ECD surface patches (10 A radius) ranked by exposed apolar area
  GLP1R_ECD_29-129.pdb           ECD only (101 residues); NOT recommended as a design target (see docs/lessons.md)
  GLP1R_full_3chains.pdb         ECD (chain R 29-129) + helix bundle as gap-free chains S (136-338) and T (344-423)
Also prints: which peptide residues touch the ECD versus the helix bundle, hotspot-set spans, and the helix-bundle axis.
Numbering = UniProt P43220 (the PDB numbering of chain R matches it). Structural analysis only.
"""
import collections
import itertools
import os
import sys

import numpy as np
import pandas as pd

ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}
MAXASA = {'A': 106, 'R': 248, 'N': 157, 'D': 163, 'C': 135, 'Q': 198, 'E': 194, 'G': 84, 'H': 184, 'I': 169,
          'L': 164, 'K': 205, 'M': 188, 'F': 197, 'P': 136, 'S': 130, 'T': 142, 'W': 227, 'Y': 222, 'V': 142}
RAD = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8}
ECD_END = 129                       # last resolved ECD residue in 6X18 (residues 130-135 are unresolved)
SEGMENTS = {"R": (29, 129), "S": (136, 338), "T": (344, 423)}   # gap-free pieces of the receptor
HOTSPOT_SETS = {"Set 1 (helix side)": [30, 32, 35, 39], "Set 2 (groove centre)": [32, 39, 69, 91]}


def load(path):
    res = collections.OrderedDict()
    for l in open(path):
        if l.startswith("ATOM"):
            k = (l[21], int(l[22:26]))
            r = res.setdefault(k, {"name": l[17:20].strip(), "atoms": {}, "el": {}})
            n = l[12:16].strip()
            r["atoms"][n] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            r["el"][n] = (l[76:78].strip() or n[0])[:1]
    return res


def dmat(A, B):
    return np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))


def xyz(r):
    return np.array(list(r["atoms"].values()))


def sasa(res, probe=1.4, npts=96):
    """Shrake-Rupley solvent-accessible area per residue (residues given as an ordered dict)."""
    keys = list(res)
    at, own, rad = [], [], []
    for i, k in enumerate(keys):
        for n, x in res[k]["atoms"].items():
            at.append(x); own.append(i); rad.append(RAD.get(res[k]["el"][n], 1.7))
    at, rad, own = np.array(at), np.array(rad) + probe, np.array(own)
    idx = np.arange(npts) + 0.5
    phi, th = np.arccos(1 - 2 * idx / npts), np.pi * (1 + 5 ** 0.5) * idx
    sph = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    D = dmat(at, at)
    out = np.zeros(len(keys))
    for i in range(len(at)):
        nb = np.where((D[i] < rad[i] + rad) & (np.arange(len(at)) != i))[0]
        pts = at[i] + rad[i] * sph
        free = np.ones(npts, bool)
        for j in nb:
            free &= (((pts - at[j]) ** 2).sum(1) >= rad[j] ** 2)
        out[own[i]] += 4 * np.pi * rad[i] ** 2 * free.mean()
    return dict(zip(keys, out))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "."
    os.makedirs(out_dir, exist_ok=True)
    st = load(sys.argv[1])
    R = collections.OrderedDict((k[1], v) for k, v in st.items() if k[0] == "R")
    P = collections.OrderedDict((k[1], v) for k, v in st.items() if k[0] == "P")
    pep = np.vstack([xyz(v) for v in P.values()])
    print(f"Receptor chain R: residues {min(R)}-{max(R)}; peptide chain P: "
          f"{''.join(ONE[v['name']] for v in P.values())}")

    # peptide residues: contacts with ECD versus helix bundle
    print("\nGLP-1 residue: contacting ECD residues / contacting helix-bundle residues (4.5 A)")
    line = []
    for pn, pv in P.items():
        ecd = {n for n, v in R.items() if n <= ECD_END and (dmat(xyz(v), xyz(pv)) <= 4.5).any()}
        tm = {n for n, v in R.items() if n > ECD_END and (dmat(xyz(v), xyz(pv)) <= 4.5).any()}
        line.append(f"{ONE[pv['name']]}{pn}:{len(ecd)}/{len(tm)}")
    print("  " + "  ".join(line))

    # groove residues on the ECD, exposure on the free receptor
    rk = collections.OrderedDict((("R", n), v) for n, v in R.items())
    sa = sasa(rk)
    rs = {n: sa[("R", n)] / MAXASA[ONE[R[n]["name"]]] for n in R}
    rows = []
    for n, v in R.items():
        if n <= ECD_END:
            c = int((dmat(xyz(v), pep) <= 4.5).sum())
            if c:
                rows.append(dict(residue=f"{ONE[v['name']]}{n}", number=n, atom_contacts=c, exposed_fraction=round(rs[n], 2)))
    groove = pd.DataFrame(rows).sort_values("atom_contacts", ascending=False)
    groove.to_csv(os.path.join(out_dir, "glp1r_groove_residues.csv"), index=False)
    print(f"\nECD residues touching the peptide: {len(groove)}; top 8: {', '.join(groove.residue.head(8))}")

    # apolar patch map on the ECD
    ca = {n: v["atoms"]["CA"] for n, v in R.items()}
    ns = [n for n in R if n <= ECD_END]
    apolar = set("AVLIMFWYC")
    asa = {n: rs[n] * MAXASA[ONE[R[n]["name"]]] for n in ns}
    gset = set(groove.number)
    prow = []
    for c in ns:
        if rs[c] < 0.2:
            continue
        pat = [n for n in ns if np.linalg.norm(ca[n] - ca[c]) <= 10 and rs[n] >= 0.2]
        ap = sum(asa[n] for n in pat if ONE[R[n]["name"]] in apolar)
        tot = sum(asa[n] for n in pat)
        prow.append(dict(center=f"{ONE[R[c]['name']]}{c}", exposed_residues=len(pat), apolar_area=int(ap),
                         apolar_fraction=round(ap / tot, 2), groove_residues_in_patch=len(set(pat) & gset)))
    pd.DataFrame(prow).sort_values("apolar_area", ascending=False).to_csv(
        os.path.join(out_dir, "glp1r_patch_map.csv"), index=False)

    # hotspot-set spans (alpha-carbon distances)
    print("\nHotspot sets:")
    for name, s in HOTSPOT_SETS.items():
        span = max(np.linalg.norm(ca[a] - ca[b]) for a, b in itertools.combinations(s, 2))
        print(f"  {name}: {', '.join(ONE[R[n]['name']] + str(n) for n in s)}; largest CA-CA gap {span:.0f} A")

    # helix-bundle axis
    tm = np.array([ca[n] for n in ca if 136 <= n <= 420])
    ecdxyz = np.array([ca[n] for n in ca if n <= ECD_END])
    cen = tm.mean(0)
    ax = np.linalg.svd(tm - cen)[2][0]
    if np.dot(ecdxyz.mean(0) - cen, ax) < 0:
        ax = -ax
    print(f"\nHelix-bundle axis: membrane-top estimate (95th percentile of helix CA) = {np.percentile((tm - cen) @ ax, 95):.1f} A "
          f"above the bundle centre; ECD spans {((ecdxyz - cen) @ ax).min():.0f}-{((ecdxyz - cen) @ ax).max():.0f} A")

    # design targets
    lines = [l for l in open(sys.argv[1]) if l.startswith("ATOM") and l[21] == "R"]
    with open(os.path.join(out_dir, "GLP1R_ECD_29-129.pdb"), "w") as fh:
        fh.writelines([l for l in lines if 29 <= int(l[22:26]) <= 129] + ["END\n"])
    out = []
    for new, (lo, hi) in SEGMENTS.items():
        out += [l[:21] + new + l[22:] for l in lines if lo <= int(l[22:26]) <= hi] + ["TER\n"]
    with open(os.path.join(out_dir, "GLP1R_full_3chains.pdb"), "w") as fh:
        fh.writelines(out + ["END\n"])
    print(f"\nWrote groove table, patch map and two target files to {out_dir}")
