# VIRR IET Public Benchmark Supplement

Supplementary artefacts for the manuscript:

**Information-Energetic Thermodynamics: A Public-Benchmark Demonstration of the VIRR Framework for Physical Decision Efficiency**

Author: Martin Petrásek

Processing date: 2026-08-01

## Contents

| ID | File | Purpose | SHA-256 |
| --- | --- | --- | --- |
| S0 | `supplement/VLSIE2026_IET_VIRR_Public_Benchmark_Demonstration.pdf` | Rendered manuscript PDF | `bd65b3c73c4dce02f69b1e4ebc1d44d0868e8a4ece41d565f0e50b97ae4c41b4` |
| S0-DOCX | `supplement/VLSIE2026_IET_VIRR_Public_Benchmark_Demonstration.docx` | Editable manuscript in the VLSIE/IET Word template | `0ec3d8f2783a9735993430702dd797f90c5db4cc76743bbc533115825c7685c9` |
| S1 | `supplement/mlperf_tiny_v1_3_virr_reproducibility_table.csv` | Derived 27-row MLPerf Tiny v1.3 VIRR_task calculation table | `f91ab7618e2042d714b1faa3022834243ea869177a359b6e28dadac84ac40c86` |
| S2 | `supplement/reproduce_mlperf_tiny_v1_3_virr.py` | Parser used to reproduce S1 from a local MLPerf Tiny v1.3 results checkout | `e23a2b682c2c73edaa187fbbbe387a213ec03bcfc57d5a27eb60db994f6c0074` |

## Public input dataset

- Dataset: MLPerf Tiny v1.3 public results repository
- Repository: <https://github.com/mlcommons/tiny_results_v1.3>
- Commit used in the manuscript: `cd605ea0981306693c9e19765a15801c776988dc`
- Number of parsed power-result rows in S1: `27`

## Reproduce S1 byte-for-byte

```bash
python3 -m pip install numpy
mkdir -p work
git clone https://github.com/mlcommons/tiny_results_v1.3.git
cd tiny_results_v1.3
git checkout cd605ea0981306693c9e19765a15801c776988dc
cd ..
mv tiny_results_v1.3 work/tiny_results_v1.3
python3 supplement/reproduce_mlperf_tiny_v1_3_virr.py work/tiny_results_v1.3 reproduced.csv
cmp supplement/mlperf_tiny_v1_3_virr_reproducibility_table.csv reproduced.csv
sha256sum reproduced.csv
```

The final checksum should be `f91ab7618e2042d714b1faa3022834243ea869177a359b6e28dadac84ac40c86`. The `work/tiny_results_v1.3` path is intentional because S1 records each public MLPerf source file in the `source_file` column.

The parser computes:

```text
VIRR_task [outcomes/J] = 1,000,000 / median_energy_uJ_per_outcome
```

For anomaly detection, AUC is reported as a benchmark quality statistic and is not interpreted as a direct count of correct outcomes.

## Notes

The public MLPerf Tiny calculation demonstrates task-level computability and cross-submission applicability. It is not a claim of full physical VIRR because the public benchmark files do not include microscopic state trajectories, residence-time validation, or entropy-production estimates.
