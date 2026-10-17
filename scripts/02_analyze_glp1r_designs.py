"""
Check BoltzGen finalists against the FULL GLP-1R (groove coverage, clashes with the helix bundle, membrane-depth guide).

Usage:
    python 02_analyze_glp1r_designs.py RESULTS_DIR 6X18.pdb --hotspots 30,32,35,39 --out finalists.csv

RESULTS_DIR is an unzipped BoltzGen/Tamarind result folder containing
final_ranked_designs/final_designs_metrics_20.csv and final_ranked_designs/final_20_designs/*.pdb.
The design target must have been the three-chain file from 01_glp1r_groove_analysis.py (BoltzGen relabels the chains
A = ECD 101 residues, B = 203, C = 80; the binder is chain D). Output has NO sequences (metrics only).

Columns: ipTM/ipSAE/PAE (BoltzGen), hotspots touched, native-footprint coverage, residues of the helix bundle touched,
atoms closer than 2.5 A to the bundle, and binder atoms below fixed levels on the bundle axis.
Axis geometry: fitted to helix CA atoms 136-420 of 6X18; membrane top estimated as the 95th percentile (26.6 A). This is a rough
guide in the active-state frame of 6X18, not a membrane or micelle model. Predicted structures only.
"""
import argparse
import collections
import glob
import os

import numpy as np
import pandas as pd

ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def load(path):
    res = collections.OrderedDict()
    for l in open(path):
        if l.startswith("ATOM"):
            r = res.setdefault((l[21], int(l[22:26])), {"name": l[17:20].strip(), "atoms": {}})
            r["atoms"][l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return res


def xyz(rs):
    return np.array([a for r in rs for a in r["atoms"].values()])


def dmat(A, B):
    return np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_dir")
    ap.add_argument("receptor_pdb")
    ap.add_argument("--hotspots", default="30,32,35,39")
    ap.add_argument("--binder-chain", default="D")
    ap.add_argument("--out", default="glp1r_finalists_analysis.csv")
    a = ap.parse_args()
    hot = {int(n): None for n in a.hotspots.split(",")}

    ref = load(a.receptor_pdb)
    RR = {k[1]: v for k, v in ref.items() if k[0] == "R"}
    PP = [v for k, v in ref.items() if k[0] == "P"]
    ca = {n: v["atoms"]["CA"] for n, v in RR.items()}
    for n in hot:
        hot[n] = f"{ONE[RR[n]['name']]}{n}"
    tm = np.array([ca[n] for n in ca if 136 <= n <= 420])
    ecd = np.array([ca[n] for n in ca if n <= 129])
    cen = tm.mean(0)
    ax = np.linalg.svd(tm - cen)[2][0]
    if np.dot(ecd.mean(0) - cen, ax) < 0:
        ax = -ax
    ztop = float(np.percentile((tm - cen) @ ax, 95))
    pepxyz = xyz(PP)
    native = {n for n, v in RR.items() if n <= 129 and (dmat(xyz([v]), pepxyz) <= 4.5).any()}

    base = os.path.join(a.results_dir, "final_ranked_designs")
    metrics = pd.read_csv(os.path.join(base, "final_designs_metrics_20.csv"))
    rows = []
    for _, r in metrics.iterrows():
        f = glob.glob(os.path.join(base, "final_20_designs", f"rank*_{r.id}.pdb"))[0]
        st = load(f)
        cnt = collections.Counter(k[0] for k in st)
        if (cnt["A"], cnt["B"], cnt["C"]) != (101, 203, 80):
            raise SystemExit(f"{f}: chains A/B/C have {cnt['A']}/{cnt['B']}/{cnt['C']} residues, expected 101/203/80. "
                             "Was the design target the three-chain file?")
        E = [v for k, v in st.items() if k[0] == "A"]
        Bm = [v for k, v in st.items() if k[0] == "B"] + [v for k, v in st.items() if k[0] == "C"]
        Bd = [v for k, v in st.items() if k[0] == a.binder_chain]
        # put the design into the 6X18 frame by fitting its ECD onto the reference ECD
        P = np.array([v["atoms"]["CA"] for v in E])
        Q = np.array([ca[29 + i] for i in range(101)])
        R, t = kabsch(P, Q)
        bx, tmx = xyz(Bd) @ R.T + t, xyz(Bm) @ R.T + t
        z = (bx - cen) @ ax
        near = dmat(bx, tmx).min(1)
        cont = {29 + i for i, v in enumerate(E) if (dmat(xyz([v]), xyz(Bd)) <= 4.5).any()}
        rows.append(dict(
            design=r.id.replace("design_spec_", "d"), iptm=round(r.design_to_target_iptm, 2),
            ipsae=round(r.design_to_target_ipsae, 2), min_pae=round(r.min_design_to_target_pae, 1),
            passes_boltzgen_filters=bool(r.pass_filters), hotspots_touched=",".join(v for n, v in hot.items() if n in cont) or "-",
            native_footprint_covered=f"{len(cont & native)}/{len(native)}",
            bundle_residues_touched=sum(1 for v in Bm if (dmat(xyz([v]), xyz(Bd)) <= 4.5).any()),
            atoms_closer_than_2p5A_to_bundle=int((near < 2.5).sum()),
            atoms_below_membrane_top=int((z < ztop).sum()), atoms_below_20A_level=int((z < 20).sum()),
            atoms_below_15A_level=int((z < 15).sum()), binder_atoms=len(bx),
            ecd_fit_rmsd_A=round(float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean())), 2)))
    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(f"membrane-top estimate {ztop:.1f} A; native GLP-1 footprint on the ECD: {len(native)} residues")
    print(df.drop(columns=["binder_atoms", "ecd_fit_rmsd_A"]).to_string(index=False))
    print(f"\nWrote {a.out} ({len(df)} designs)")
