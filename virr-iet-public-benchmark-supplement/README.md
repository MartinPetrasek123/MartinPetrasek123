# VIRR IET Public Benchmark Supplement

Supplementary artefacts for the manuscript:

**Information-Energetic Thermodynamics: A Public-Benchmark Demonstration of the VIRR Framework for Physical Decision Efficiency**

Author: Martin Petrásek

Processing date: 2026-08-01

## Contents

| ID | File | Purpose | SHA-256 |
| --- | --- | --- | --- |
| S0 | `supplement/VLSIE2026_IET_VIRR_Public_Benchmark_Demonstration.pdf` | Rendered manuscript PDF | `008e21b2eeb7b6f28f0c7fbe89762155f818f172e47f35bae63ba9a377ea2ba0` |
| S0-DOCX | `supplement/VLSIE2026_IET_VIRR_Public_Benchmark_Demonstration.docx` | Editable manuscript in the VLSIE/IET Word template | `84de263856fb5dc3334b457cf2ca701d0e9c8852f7cb28677020fe5e036d6107` |
| S1 | `supplement/mlperf_tiny_v1_3_virr_reproducibility_table.csv` | Derived 27-row MLPerf Tiny v1.3 VIRR_task calculation table | `f91ab7618e2042d714b1faa3022834243ea869177a359b6e28dadac84ac40c86` |
| S2 | `supplement/reproduce_mlperf_tiny_v1_3_virr.py` | Parser used to reproduce S1 from a local MLPerf Tiny v1.3 results checkout | `e72d853918152416f7e94f90a8789e1c68e83b36f3ff3f5878bfac80e006f656` |

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
