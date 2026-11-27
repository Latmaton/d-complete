# Checking the consistency of saved results

To check results produced by the current package, run the following command from the package root:

```bash
python check_saved_results.py --require-all
```

The script uses only the Python standard library. It reads the current `parameters.py` and `results/*/summary.json`. To select another results directory, use:

```bash
python check_saved_results.py --results /path/to/results
```

By default, the checker audits every run directory it finds. It reports the cases covered by complete, consistent records, but does not require all eight cases. To require all eight cases, use:

```bash
python check_saved_results.py --require-all
```

To check a single case and require at least one complete record for that case, use:

```bash
python check_saved_results.py --pair -2 -3
```

The reversed order `--pair -3 -2` selects the same fixed case, as in `main.py`. `--pair` and `--require-all` cannot be combined. The checker reads the basic summary fields in every run directory so that missing or malformed records remain errors. After that check, records for other pairs are reported as outside the selected pair and their numerical contents are not audited.

The `--source-archive` option explicitly selects the source bytes used for the hash comparison. The ZIP archive must contain exactly one member with each of the six core filenames listed below, at the archive root. Without this option, the script checks the current source files. There is no automatic fallback to the archive. The script prints which source version it is checking and does not execute code from the archive.

For each record with `status="passed"`, the script checks the following items. All rational arithmetic uses `fractions.Fraction`, without floating-point approximations.

1. The recorded `source_sha256` values must exactly match the SHA-256 hashes of the six core source files `main.py`, `parameters.py`, `residues.py`, `terminals.py`, `arithmetic.py`, and `reporting.py`. These files are read from the current package or from the explicitly selected source archive.
2. The bases, paper target, internal verification target, period data, scales, partition intervals, and settings for positive bases must agree entry by entry with the current fixed parameters. The comparison of `parameters.json` also checks JSON types, so an integer cannot be replaced by a Boolean or a floating-point number. The numbers of scales, intervals, endpoints, and terminal requests must be complete. The `chunk` parameter may take the positive integer value explicitly recorded for that run.
3. The height counts for each size bound must be nonnegative integers within the permitted height range. The count at height zero must be 1, the total must equal the modulus $P$, and the height-distribution CSV and JSON must agree.
4. The script recomputes the exponential moments and finite tail bounds from the saved height distributions, then checks `moments.json` and `tails.csv`. It also recomputes the scale coefficients from the saved terminal lower bounds and checks `scales.csv`.
5. Using the fixed choices on each partition interval and the resulting scale coefficients, the script recomputes every individual contribution, numerator, and normalized value at every endpoint. It checks `endpoints.json` and its agreement with `endpoints.csv`. Finally, it verifies that the lower bound in `summary.json` equals the minimum of all endpoint values and strictly exceeds the internal target.
6. A passed record must have a `stages_seconds` object with exactly the four stages `residues`, `tails`, `terminals`, and `density`. Each time must be a finite, nonnegative JSON number, not a Boolean. This checks the completeness of the recorded stage metadata. It does not independently establish that a stage was executed.

The optional `current_stage` field is informational and is not checked by this script. The current main program writes `complete` in that field when a run passes; the checker's completion requirements use `status` and the four-stage metadata described above.

The only accepted statuses are `passed`, `running`, `failed`, and `interrupted`. Records with one of the last three statuses are listed as incomplete and are not checked as passed records. Their presence in the selected scope causes exit status 3, while the checker continues auditing completed records. A missing or unknown status is an error, as is a malformed pair, an unparseable `summary.json`, or a run subdirectory without `summary.json`. Any direct subdirectory of the selected results directory is treated as a run directory. Store unrelated directories elsewhere.

| Exit status | Meaning |
| --- | --- |
| `0` | At least one complete record was consistent, no incomplete or erroneous record was found in the selected scope, and any requested coverage requirement was met. |
| `1` | A record is malformed or inconsistent, a run subdirectory has no summary, the results directory is missing, or the selected source archive is invalid. This takes precedence over status 3. |
| `2` | A command-line usage error, such as an unsupported pair or combining `--pair` with `--require-all`. |
| `3` | No complete record passed, a legitimate incomplete record exists in the selected scope, or a case required by `--pair` or `--require-all` lacks a complete consistent record. |

An empty results directory returns 3. In default mode, one complete case may return 0, but the printed coverage list identifies that case and explicitly says that full coverage was not required. Multiple runs of the same pair count as one covered case. Use `--require-all` when checking whether the records cover all eight cases. A previous incomplete run in the selected scope still causes status 3 even if a later complete run of the same pair exists. The checker never changes or removes old records.

A directory containing one complete case returns 0 in default mode and 3 with `--require-all`. A directory containing complete, consistent records for all eight cases and no incomplete or malformed records returns 0 with `--require-all`.

The regression checks can be run with:

```bash
python check_saved_results_tests.py
```

The regression script selects a passed `(-2,-3)` run whose source hashes match the current six core files. It copies that run into a temporary directory before testing partial coverage, so the supplied results directory may contain all eight cases. Use `--results /path/to/results` to select another results parent. Tests create their own incomplete records and damaged copies in temporary storage. The supplied records and source files are left unchanged.

**This check verifies only the consistency of saved records. It does not replace constructing the residue families and terminal lower bounds from scratch.** A height-count sum of $P$ checks only the recorded coverage count, without rechecking an actual representation for every residue. The terminal lower bounds are used here as saved data, without reconstructing their representation families or interval families. Source hashes identify the source version declared by a run record. They are not signatures authenticating the record. The optional `joint_counts.npz` is outside the scope of this script. Full reproduction still requires rerunning `main.py` as described in the main instructions.

For the mathematical roles and derivation of the parameters, see [01_parameters.md](01_parameters.md). The checker itself is [check_saved_results.py](../check_saved_results.py).
