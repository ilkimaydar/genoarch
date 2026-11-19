# Data card

## Sources

| File | Release | Bytes | sha256 | Public source |
| --- | --- | --- | --- | --- |
| `hp.json` | HPO 2026-09-01 | 21,915,519 | `a7b3a012e7b4007a35a7cf8da35f2373b54cc16907f80d13d62470d31f833501` | github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/hp.json |
| `hp.obo` | HPO 2026-09-01 | 10,863,613 | `93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22` | same release, `hp.obo` |
| `phenotype.hpoa` | annotations 2026-09-02 | 35,816,037 | `e89aa39c8f97bf5a52c5d160f9250a603681d630ade8ec2679d84ac5aece1f72` | same release, `phenotype.hpoa` |
| `genes_to_disease.txt` | 2026-09-01 release | 1,492,039 | `7937a1b92779dc8614e39ace623b443c7c71dd15e4de11d5e52b14e331ab5978` | same release |
| `genes_to_phenotype.txt` | 2026-09-01 release | 20,821,704 | `507a17bff9c49e6329fbd88b1f91734fa7e9fdea5c06304044c0892eb6ab248c` | same release |
| `all_phenopackets.zip` | Phenopacket Store 0.1.27 | 19,431,098 | `d0c70005bb09b87035087516252b6b539fbc0415be7c664b73c9cde7eedf205a` | github.com/monarch-initiative/phenopacket-store/releases/download/0.1.27/all_phenopackets.zip |

The full release files are not redistributed in this repository; `scripts/fetch_data.sh` downloads and verifies them.

## Processing into the index

* Terms: the Phenotypic abnormality subtree (HP:0000118 and descendants): 19,178 of 19,894 current terms.
  Labels, synonyms and `is_a` links are copied unchanged. 588 obsolete terms map to their replacements and
  alternative ids map to current terms.
* Disease annotations: `phenotype.hpoa` has 286,651 rows. Rows with qualifier NOT are excluded from the
  positive profile. 17,736 disease-term pairs fall outside the Phenotypic abnormality subtree (aspects I
  8,978, C 8,557, H 124, M 77: inheritance, clinical course, past history, modifiers) and are not used for
  similarity. Several frequency observations of one
  pair are pooled (all fractions: sum of numerators / sum of denominators; otherwise the maximum; HPO frequency
  terms use class midpoints). 3,018 pairs with pooled frequency 0 are dropped. Unknown frequency is stored as
  unknown and read as 0.5 where a frequency is needed.
* Result: 12,865 diseases (8,463 OMIM, 4,355 Orphanet, 47 DECIPHER) with 264,562 annotations.
* Information content: IC(t) = ln(N / n_t), N = 12,865, n_t = diseases annotated to t or a descendant;
  rounded to 1e-6.
* Genes: `genes_to_disease.txt` links of 5,550 genes; 5,268 genes keep at least one link to an annotated
  disease (3,441 links point to diseases without phenotype annotations). Association types: OMIM links are
  MENDELIAN or POLYGENIC; Orphanet links are UNKNOWN in this release.
* Gene phenotype profiles are derived as the union of each gene's disease annotations. They equal the
  `genes_to_phenotype.txt` profile for 4,381 of 5,268 genes. For the other 887 genes the file has extra terms and
  no gene has union terms missing from it; the extra rows are 1,930 rows whose disease-term pair has pooled
  frequency 0 and 131 rows whose pair is NOT-qualified in `phenotype.hpoa`, both excluded from the index on purpose.
* Genes: 4,972 have at least one OMIM link; 296 are linked only through Orphanet.

Logs: `exp/platform_eval/index_stats_20261118T231000Z/` (stats.json from `scripts/index_stats.py`, index_build_report.json) for the
values above; index file sha256 `9ff8cc6bf2fa9d8bd64a65e56aecf49e15ca48e7d3eb10f1e20bd88c53aade02`.

## Evaluation data

Phenopacket Store 0.1.27: 10,377 phenopackets. Present terms are features without `excluded`; absent terms are
features with `excluded: true`; the truth disease is the interpretation's diagnosis; the truth gene is the
gene symbol of the genomic interpretation. 29 cases have a disease id that is not annotated in this HPO release
and 3 have no usable present term; 10,345 cases remain. Split by disease id (sha256 of the id, first 8 hex digits,
modulo 10 below 3 = calibration): 3,619 calibration cases, 6,726 test cases from 537 diseases. These are solved
phenopacket cases used as a retrospective proxy for undiagnosed cases; no open cohort of truly undiagnosed patients
is available to us.

## Known biases

* Annotation depth varies by disease and source (median 12 annotations for OMIM entries, 22 for Orphanet).
* Phenopacket Store over-represents some genes (the largest cohort has 462 cases) and holds published, solved
  cases.
* Excluded findings in case reports are often pertinent negatives for the true disease; this is why the fitted
  absent-term weight is 0.

## Licenses and attribution

* HPO: free to use under the terms on https://human-phenotype-ontology.github.io/license.html (cite the HPO,
  display file versions, do not alter content or relationships; products should state that they use the HPO).
  Genoarch shows the versions in the page footer, on the Method page and in `GET /health`.
* OMIM and Orphanet: identifiers, names and annotations arrive through the HPO files. Their own terms of use
  were not checked for this build (UNVERIFIED); check them before commercial use.
* Phenopacket Store: used for calibration and evaluation only; the repository stores no case records. The
  parity patients list Phenopacket ids and HPO term ids for 16 test cases. Its license was not checked for
  this build (UNVERIFIED).
