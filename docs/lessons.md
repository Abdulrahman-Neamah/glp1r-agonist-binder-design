# Lessons from GLP-1R (Anthropic x Adaptyv Challenge 3)

Goal: a de novo protein that binds human GLP-1R (supplied in detergent) and activates it in a cell assay. BoltzGen and Boltz-2 cannot predict
activation and I did not try to engineer it; the organizers note that a groove-only binder may not activate the receptor.
Run summaries: `results/glp1r_boltzgen_runs.csv`.

1. **A partial target creates a loophole.** Two runs against the extracellular domain (ECD) alone gave good-looking scores (47 of 100 designs passed
   BoltzGen's filters, best ipTM 0.62), but every binder overlapped the helix bundle (median 139 atoms closer than 2.5 A; none had a clean fit).
   The binders had wrapped into space the omitted bundle occupies. Fix: give the designer the whole receptor, split into gap-free chains
   (ECD, then the two helix-bundle pieces), with hotspots only on the ECD. The same thing happened to me on a TNF-alpha dimer target.
2. **Check `settings.json` after every run.** My second attempt at the full-receptor job only had chain R in `targetChains`, so the helix bundle was missing
   and the results were invalid again. The Target Chains field must list every chain (R,S,T).
3. **A greasy site is much easier than a polar one.** The GLP-1R groove has about 360-440 A^2 of exposed apolar surface per 10 A patch (TNF-alpha's groove core
   had 24-133 A^2) and the tool found binders at once: 20 and 16 designs with ipTM 0.5 or more in two 80-design batches, against none usable in five TNF runs.
   Compact hotspot sets on one chain (13 and 16 A across) worked.
4. **Rough geometry still helps.** I fitted an axis to the helix bundle (`scripts/01_glp1r_groove_analysis.py`) and counted binder atoms below fixed levels
   (`scripts/02_analyze_glp1r_designs.py`). Many binders reach down toward the membrane zone; the micelle is not modelled and the ECD orientation comes from
   the active-state structure, so treat the depth counts as a guide and not a membrane simulation.
5. **Selection is heuristic.** 22 candidates passed the platform's novelty check (3/4 or higher); 20 were chosen by a simple score
   (ipSAE + 0.5 x ipTM + 0.05 if BoltzGen filters passed, minus penalties for clashing atoms and deep placement). The weights are my judgement.
