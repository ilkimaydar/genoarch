# Model card: Genoarch engine 0.1.0

## Summary

A phenotype-similarity ranker with a fitted temperature. Inputs: HPO terms marked present or absent. Outputs:
probabilities over diseases and genes, uncertainty measures, an open-world (abstain) flag, next-phenotype
suggestions and gene hypotheses. No neural network and no patient data are part of the model; the only fitted
parameters are the blend weight w, the absent-term weight beta and the temperature T.

## Intended use

Research prioritisation in undiagnosed disease studies. Not for diagnosis, screening, treatment or any
clinical decision. Not a medical device.

## Parameters

Fitted on the calibration slice (3,619 Phenopacket Store cases from a disease-disjoint split), per mode and
disease universe:

| Mode | Universe | w | beta | T |
| --- | --- | --- | --- | --- |
| blended (default) | all (default) | 0.05 | 0 | 0.173602 |
| blended | omim | 0.10 | 0 | 0.172491 |
| asymmetric | all | 0 | 0 | 0.18215 |
| asymmetric | omim | 0 | 0 | 0.19002 |
| symmetric | all | 0.5 | 0 | 0.203571 |
| symmetric | omim | 0.5 | 0 | 0.196938 |

Selection: w and beta by disease mean reciprocal rank on a grid (w in {0.02, 0.05, 0.1, 0.2, 0.3}; beta in
{0, 0.02, 0.05, 0.1, 0.25, 0.5, 1}); T by minimum negative log-likelihood of the true disease; default mode and
universe by gene mean reciprocal rank on the calibration slice. Fixed constants: alpha 0.05, 100 null patients,
eps 0.05, q_max 0.99, information gain over the top 50 diseases, hypotheses outside the top 10 diseases with at
least 2 matches of IC >= 3.

## Evaluation (test slice: 6,726 cases, 537 diseases)

| Configuration | Gene top-1 | Gene top-10 | Disease top-1 | Disease top-10 | ECE | NLL |
| --- | --- | --- | --- | --- | --- | --- |
| blended / all (default) | 0.658 | 0.818 | 0.613 | 0.764 | 0.119 | 2.772 |
| blended / omim | 0.649 | 0.812 | 0.629 | 0.780 | 0.105 | 2.502 |
| asymmetric / all | 0.638 | 0.809 | 0.601 | 0.763 | 0.121 | 2.813 |
| asymmetric / omim | 0.638 | 0.799 | 0.609 | 0.772 | 0.107 | 2.549 |
| symmetric / all | 0.475 | 0.768 | 0.393 | 0.698 | 0.157 | 3.844 |
| symmetric / omim | 0.466 | 0.758 | 0.404 | 0.713 | 0.136 | 3.620 |

Default configuration, further metrics: gene top-3 0.737, disease top-3 0.695, gene MRR 0.712, disease MRR
0.667, ECE before temperature scaling 0.609, Brier score of the top-1 0.152, AUROC of top-1 confidence for a
correct top-1 0.891; 48.0% of cases have confidence of at least 0.5 and 91.0% of those are correct. The
reliability diagram shows under-confidence between 0.1 and 0.5 (for example, bin 0.1-0.2: mean confidence
0.146, accuracy 0.363).

Open-world test (200 test cases, 100 null patients each): abstain rate 0.245 when the true disease is in the
index, 0.285 after removing every disease linked to the causal gene, and 0.960 for random phenotype sets of the
same size and specificity. AUROC of the p-value for "gene removed" vs "in index": 0.538. The test is calibrated
under the null but is a weak detector of missing diseases when similar diseases exist.

Next-phenotype simulation (192 cases, start from half of the present terms, one question): the top suggestion
could be answered from the full record in 54 cases and improved the true disease's rank in 44.4% of them (mean
log2 rank gain 0.207); a random candidate could be answered in 15 cases and improved the rank in 6.7% of them
(mean gain -0.026).

Logs: `exp/platform_eval/run_20261118T224839Z/` (summary.json, test_cases.jsonl, open_world_cases.jsonl,
next_phenotype_sim.jsonl, run.log, the script). A copy of summary.json is in `docs/evaluation/`.

## Design decision log

* Open-world null: an earlier variant standardised null patients in-sample. On 150 cases it abstained on only
  78.7% of random phenotype sets (expected about 95%). Leave-one-out standardisation with a pooled variance floor
  gave 96.7% on random sets with 21.3% on in-index cases. Log: `exp/platform_eval/null_standardization_prototype/`.
* Absent terms: every beta > 0 lowered calibration MRR (see the grid in summary.json), so beta = 0 by default.

## Caveats

* The evaluation uses solved phenopacket cases (Phenopacket Store 0.1.27) as a retrospective proxy for undiagnosed
  cases; no open cohort of truly undiagnosed patients is available to us, so the reported metrics describe solved,
  published cases. Performance on undiagnosed cohorts is expected to be lower.
* Duplicate OMIM and Orphanet entries split probability in the "all" universe.
* Ties in the asymmetric score are broken by disease position; the blended mode reduces such ties.
* Gene hypotheses reflect phenotype overlap only.
