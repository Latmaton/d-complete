# Verification record

All eight cases were computed from scratch and passed on 27 November 2026. Each run constructed its residue families and terminal lower bounds, then checked the moments, tails, scales, and endpoint inequalities using exact arithmetic. The saved records are in [results](../results/README.md).

## 1. Completed computations

The exact rational lower bounds, targets, elapsed times, and run directories are collected in [run_summary.csv](run_summary.csv). Decimals in the table below are for reference only.

| Bases | Status | Paper target | Computed lower bound, decimal for reference | Record |
| --- | --- | --- | --- | --- |
| `(5,2)` | Passed | `9/10` | `0.9001036130335016` | [Summary](../results/5_2_20261127T131118_450218Z/summary.json) |
| `(-2,3)` | Passed | `159/500` | `0.3185187934794411` | [Summary](../results/-2_3_20261127T131118_446921Z/summary.json) |
| `(-2,-3)` | Passed | `31/100` | `0.3197266069585859` | [Summary](../results/-2_-3_20261127T131118_432152Z/summary.json) |
| `(-5,2)` | Passed | `13/30` | `0.4421710832652359` | [Summary](../results/-5_2_20261127T131128_451719Z/summary.json) |
| `(5,-2)` | Passed | `7/20` | `0.3662568930393349` | [Summary](../results/5_-2_20261127T131242_587491Z/summary.json) |
| `(-5,-2)` | Passed | `56/125` | `0.4668356533545396` | [Summary](../results/-5_-2_20261127T131354_171877Z/summary.json) |
| `(7,-2)` | Passed | `1/500000000` | `2.7369352241470965e-9` | [Summary](../results/7_-2_20261127T131119_512352Z/summary.json) |
| `(-7,-2)` | Passed | `1/2100000000` | `6.477584120040208e-10` | [Summary](../results/-7_-2_20261127T131123_636991Z/summary.json) |

The positive case gives the exact bound

$$
\frac{27908486048646169011}{31005859375000000000}
>\frac{9001}{10000}>\frac9{10}.
$$

All eight summaries record `status="passed"`, `current_stage="complete"`, all four stage timings, and the matching hashes of the six core source files. The saved-record check with `--require-all` returned exit status `0`, with eight consistent passed records and no incomplete or inconsistent records. See [saved_results.txt](saved_results.txt).

The console logs are retained for [the positive case](full_positive.txt), [the `(-2,3)` case](full_minus2_3.txt), and [the other six cases](full_other_six.txt).

## 2. Checks and retained records

| Check | Record |
| --- | --- |
| Fixed parameters and inexpensive constant conditions | [parameters.txt](parameters.txt) |
| The same parameter checks under `python -O` | [optimized_parameters.txt](optimized_parameters.txt) |
| Complete supplementary tables against their parameter source | [supplement.txt](supplement.txt) |
| Recurrences against independent small enumerations and boundary cases | [check_small.txt](check_small.txt) |
| The 36 saved-record regression cases included in this package | [Recorded output excerpt](check_saved_results_tests.txt) |
| Complete saved records for all eight cases | [saved_results.txt](saved_results.txt) |
| Saved positive joint arrays against their height distributions and total `5**11` | [positive_joint_check.json](positive_joint_check.json) |

These records come from the completed verification session on 27 November 2026. The regression log contains the recorded output for the 36 included cases. Preparing this publication package did not rerun the computations or tests.

The checks can be invoked from the package root with:

```bash
python main.py --check-parameters
python -O main.py --check-parameters
python export_parameters.py --check
python check_small.py
python check_saved_results_tests.py
python check_saved_results.py --require-all
```

## 3. Scope and reproduction

Parameter checks verify fixed data and inexpensive conditions. Small enumerations test recurrences on explicitly bounded examples. The saved-record checker recomputes moments, tails, scales, and endpoints from the saved counts. It does not reconstruct residue chains or terminal families. The optional joint array file is outside that checker's scope and was checked separately as recorded above.

To reconstruct all eight cases afresh, use a separate output directory:

```bash
python main.py --all --output reproduced_results
python check_saved_results.py --results reproduced_results --require-all
```

A `running`, `failed`, or `interrupted` record is not a completed verification. Interrupted calculations are not resumed automatically. Source hashes identify source file contents; they are not signatures authenticating a calculation.

The mathematical deduction from the finite inequalities to lower asymptotic density uses the proofs in the paper and the supplementary notes. After changing parameters or algorithms, check their hypotheses and repeat the relevant computations. See [the overview](../docs/00_overview.md) and [the saved-result notes](../docs/05_saved_results.md).
