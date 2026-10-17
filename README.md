# GLP-1R agonist design: groove binders against the full receptor

Computational workflow for the third challenge of the
[Anthropic x Adaptyv protein design competition](https://proteinbase.com/competitions/anthropic-adaptyv-2026):
design a de novo protein that **binds human GLP-1R** (the receptor behind semaglutide and tirzepatide, supplied in detergent) and **activates it**
in a cell-based reporter assay.

**Status:** 20 designs submitted; lab results pending. **Design sequences are not included here** and will be added after the results are
published. Everything below is based on *predicted* structures and scores. Nothing has been tested in the lab.

> \*\*Activation was not engineered.\*\* BoltzGen and Boltz-2 cannot predict whether a binder switches a receptor on, and the organizers note that a binder that
> only fills the peptide-binding groove of the extracellular domain (ECD) may not activate GLP-1R. The goal I could work towards was strong, valid groove binders.

## Workflow

1. **Find the site** (`scripts/01\_glp1r\_groove\_analysis.py`, PDB 6X18): the C-terminal half of GLP-1 docks in an ECD groove (17 contacting residues) while its
N-terminal half contacts only the helix bundle. The groove is greasy (about 360-440 A^2 of exposed apolar surface per 10 A patch).
Two compact hotspot sets: V30, L32, T35, W39 (13 A across) and L32, W39, Y69, W91 (16 A across).
2. **Design** (BoltzGen via Tamarind Bio): protein-binder task, 70-100 residues, binding-site filter on, two batches of 80 designs. The target is the **whole receptor**
split into three gap-free chains (ECD 29-129, helix bundle 136-338 and 344-423), with hotspots on the ECD only.
3. **Check against the full receptor** (`scripts/02\_analyze\_glp1r\_designs.py`): hotspot contacts, coverage of the 17 native GLP-1 contact residues, clashes with the
helix bundle, and a rough membrane-depth guide. `scripts/03\_check\_full\_receptor\_clashes.py` is the general version of the clash check.
4. **Novelty check** (competition platform, automatic): designs needed a score of 3/4 or higher. All 22 candidates passed.
5. **Select 20** with a simple heuristic score (`docs/methodology.md`).

## Key results (predicted)

|Run|Target|Designs|Passed BoltzGen filters|ipTM 0.5 or more|Best ipTM / ipSAE|Outcome|
|-|-|-|-|-|-|-|
|1|ECD only|100|47|17|0.62 / 0.23|invalid: all binders overlap the helix bundle (median 139 clashing atoms)|
|2|ECD only (form error)|80|42|about 20|0.58 / about 0.3|invalid: median 195 clashing atoms|
|3|Full receptor|80|9|20|0.67 / 0.45|usable: 15 candidates kept|
|4|Full receptor|80|3|16|0.67 / 0.45|usable: 7 candidates kept|

Per-design tables (no sequences): `results/`. Methods and the selection formula: [`docs/methodology.md`](docs/methodology.md).

### The near-miss

Runs 1 and 2 looked good by score, but every binder had wrapped into the space of the helix bundle I had left out of the target. Run 2 was also a form
mistake: only one chain was entered. Giving the designer the whole receptor fixed it. Details and other lessons: [`docs/lessons.md`](docs/lessons.md).

## Limits

* Predicted structures only; no affinity or activation prediction.
* The ECD orientation relative to the bundle comes from the active-state structure (6X18) and may differ in the assay construct (detergent, thermostabilising
mutations, FLAG and His tags). The micelle is not modelled, so the depth counts are a rough guide.
* The selection score is heuristic.

## Reproducing

```bash
pip install -r requirements.txt
python scripts/01\_glp1r\_groove\_analysis.py 6X18.pdb out/        # groove table, patch map, and the target files (download 6X18 from files.rcsb.org)
# design with BoltzGen (e.g. on Tamarind Bio) against out/GLP1R\_full\_3chains.pdb, chains R,S,T, hotspots on chain R
python scripts/02\_analyze\_glp1r\_designs.py RESULTS\_DIR 6X18.pdb --hotspots 30,32,35,39 --out finalists.csv
```

`cp designs\_template.csv designs.csv` for your own sequences (git-ignored). The scripts need only numpy and pandas; they have not been tested by anyone but the author,
so expect small fixes.

## Related

My earlier repository covers the first two challenges (EGFR and TNF-alpha), including another case where a partial target gave misleading scores:
[egfr-conditional-binder](https://github.com/Abdulrahman-Neamah/egfr-conditional-binder).

## Notes

* Tools: [BoltzGen](https://github.com/HannesStark/boltzgen) and [Tamarind Bio](https://www.tamarind.bio) (free plan), PDB 6X18. Each has its own license; please check them before reuse.
* The workflow and scripts were developed with help from Claude (Anthropic).
* Competition terms (publication license, effect on patentability) apply to the designs, not to this code.

## License

Code: MIT (see `LICENSE`).

