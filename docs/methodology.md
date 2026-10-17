> Design sequences are intentionally not included in this repository (they will be published after the competition results). Metrics are from predicted structures only; nothing has been tested in the lab.

# GLP1R (Challenge 3): methodology details for the final 20 designs

## Target and site
- Human GLP-1R (UniProt P43220), PDB 6X18 chain R (GLP-1 bound), numbering = UniProt. ECD = residues 29-129; helix bundle = 136-423.
- The groove: the 17 ECD residues contacted by the C-terminal half of GLP-1 (including E68, W39, S31, L32, W91, Y69, V30, T35). The peptide's N-terminal half contacts only the helix bundle.
- The groove is greasy (apolar surface patches of about 360-440 A^2 within 10 A), so it suits a de novo binder.

## Design (two BoltzGen batches, via Tamarind Bio)
- Protein-binder task, 70-100 residues, 80 designs per batch, binding-site filter on. The target contained the ECD **and** the helix bundle (three gap-free chains).
- Batch 1 hotspots: V30, L32, T35, W39 (9 of 80 passed BoltzGen filters; 20 designs with ipTM 0.5 or more; best ipTM 0.67, ipSAE 0.45).
- Batch 2 hotspots: L32, W39, Y69, W91 (3 of 80 passed; 16 designs with ipTM 0.5 or more; best ipTM 0.67, ipSAE 0.45).
- Two earlier runs against the ECD alone were discarded: every binder overlapped the helix bundle (median 139 and 195 atoms closer than 2.5 A).
- Selection: finalists were checked for hotspot contacts, coverage of the native footprint, clashes with the helix bundle, and depth relative to an estimated membrane top (rough; the detergent micelle is not modelled). The 20 submitted designs were chosen from 22 candidates that all passed the platform's novelty check (3/4 or higher), using a simple score: ipSAE + 0.5 x ipTM + 0.05 if the design passed BoltzGen filters, minus penalties for clashing atoms and for atoms deep below the membrane-top estimate. Weights are heuristic.

## Limits
- Nothing here predicts receptor activation. These are groove binders, and the challenge text notes that a binder filling only the groove may not activate GLP1R. Activation was not engineered.
- Predicted structures only; no affinity prediction. The ECD orientation relative to the helix bundle comes from the active-state structure and may differ in the assay construct (detergent, thermostabilised, tagged).

## Submitted designs
BoltzGen ipTM / ipSAE / min PAE; native GLP-1 footprint residues covered; atoms closer than 2.5 A to the helix bundle; binder atoms below the 15 A level on the receptor axis (a stricter cut about 12 A under the estimated membrane top at 26.6 A).

| Rank | Design | Length | Scores | Footprint | Clash atoms | Atoms below 15 A level |
|---|---|---|---|---|---|---|
| 1 | glp1r_b2_d34 | 80 | 0.66 / 0.37 / 2.8 | 13/17 | 2 | 6 |
| 2 | glp1r_b2_d07 | 92 | 0.58 / 0.33 / 4.8 | 17/17 | 2 | 2 |
| 3 | glp1r_d35 | 91 | 0.64 / 0.44 / 4.0 | 16/17 | 4 | 48 |
| 4 | glp1r_d17 | 75 | 0.61 / 0.32 / 4.3 | 15/17 | 0 | 3 |
| 5 | glp1r_d19 | 77 | 0.63 / 0.39 / 4.7 | 16/17 | 4 | 35 |
| 6 | glp1r_b2_d13 | 81 | 0.58 / 0.31 / 5.0 | 15/17 | 4 | 0 |
| 7 | glp1r_d24 | 90 | 0.66 / 0.39 / 4.1 | 14/17 | 3 | 30 |
| 8 | glp1r_b2_d25 | 82 | 0.67 / 0.45 / 3.7 | 14/17 | 8 | 36 |
| 9 | glp1r_d32 | 86 | 0.58 / 0.38 / 4.3 | 13/17 | 1 | 38 |
| 10 | glp1r_d11 | 74 | 0.58 / 0.27 / 5.1 | 13/17 | 2 | 0 |
| 11 | glp1r_d37 | 83 | 0.66 / 0.41 / 4.0 | 16/17 | 4 | 38 |
| 12 | glp1r_d10 | 87 | 0.60 / 0.41 / 4.8 | 15/17 | 5 | 45 |
| 13 | glp1r_d01 | 74 | 0.58 / 0.36 / 5.0 | 13/17 | 2 | 36 |
| 14 | glp1r_b2_d30 | 81 | 0.61 / 0.37 / 4.8 | 16/17 | 6 | 15 |
| 15 | glp1r_b2_d38 | 90 | 0.60 / 0.40 / 4.2 | 14/17 | 2 | 46 |
| 16 | glp1r_d14 | 77 | 0.56 / 0.22 / 6.3 | 14/17 | 2 | 0 |
| 17 | glp1r_b2_d01 | 90 | 0.57 / 0.32 / 5.3 | 14/17 | 5 | 9 |
| 18 | glp1r_d03 | 78 | 0.61 / 0.31 / 4.9 | 15/17 | 1 | 49 |
| 19 | glp1r_d13 | 81 | 0.64 / 0.36 / 4.1 | 14/17 | 6 | 38 |
| 20 | glp1r_d04 | 87 | 0.53 / 0.30 / 4.9 | 12/17 | 6 | 2 |
