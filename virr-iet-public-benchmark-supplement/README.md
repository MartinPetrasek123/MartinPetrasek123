# VIRR IET Public Benchmark Supplement

Supplementary artefacts for the manuscript:

**Information-Energetic Thermodynamics: A Public-Benchmark Demonstration of the VIRR Framework for Physical Decision Efficiency**

Author: Martin Petrásek

Processing date: 2026-08-01

## Contents

| ID | File | Purpose | SHA-256 |
| --- | --- | --- | --- |
| S0 | `supplement/VLSIE2026_IET_VIRR_Public_Benchmark_Demonstration.pdf` | Rendered manuscript PDF | `5c6ac03577fa3b6d72d3a962c282a854e71f624320c4e03d7678d743d0104592` |
| S0-DOCX | `supplement/VLSIE2026_IET_VIRR_Public_Benchmark_Demonstration.docx` | Editable manuscript in the VLSIE/IET Word template | `4da6de9077074e0c6250c872c14acabd1490e52fc46b5b9dee0084f8a06eb868` |
| S1 | `supplement/mlperf_tiny_v1_3_virr_reproducibility_table.csv` | Derived 27-row MLPerf Tiny v1.3 VIRR_task calculation table | `f91ab7618e2042d714b1faa3022834243ea869177a359b6e28dadac84ac40c86` |
| S2 | `supplement/reproduce_mlperf_tiny_v1_3_virr.py` | Parser used to reproduce S1 from a local MLPerf Tiny v1.3 results checkout | `12b2c24812fd7affa70da536e9416ae5ba17176c7e7fdef261d7873a131e8066` |

## Public input dataset

- Dataset: MLPerf Tiny v1.3 public results repository
- Repository: <https://github.com/mlcommons/tiny_results_v1.3>
- Commit used in the manuscript: `cd605ea0981306693c9e19765a15801c776988dc`
- Number of parsed power-result rows in S1: `27`

## Reproduce S1

```bash
git clone https://github.com/mlcommons/tiny_results_v1.3.git
cd tiny_results_v1.3
git checkout cd605ea0981306693c9e19765a15801c776988dc
cd ..
python3 supplement/reproduce_mlperf_tiny_v1_3_virr.py tiny_results_v1.3 reproduced.csv
```

The parser computes:

```text
VIRR_task [outcomes/J] = 1,000,000 / median_energy_uJ_per_outcome
```

For anomaly detection, AUC is reported as a benchmark quality statistic and is not interpreted as a direct count of correct outcomes.

## Notes

The public MLPerf Tiny calculation demonstrates task-level computability and cross-submission applicability. It is not a claim of full physical VIRR because the public benchmark files do not include microscopic state trajectories, residence-time validation, or entropy-production estimates.
