# Checking the consistency of saved results

To check results produced by the current package, run the following command from the package root:

```bash
python check_saved_results.py
```

The script uses only the Python standard library. It reads the current `parameters.py` and `results/*/summary.json`. To select another results directory, use:

```bash
python check_saved_results.py --results /path/to/results
```

The original recorded runs are kept in `validation/recorded_results`. Their numerical data and recorded source hashes are unchanged. The exact six core source files used for those runs are preserved in `validation/recorded_source.zip`. Check those historical records with:

```bash
python check_saved_results.py --results validation/recorded_results --source-archive validation/recorded_source.zip
```

The `--source-archive` option explicitly selects the source bytes used for the hash comparison. The ZIP archive must contain exactly one member with each of the six core filenames listed below, at the archive root. Without this option, the script checks the current source files. There is no automatic fallback to the archive. The script prints which source version it is checking and does not execute code from the archive.

For each record with `status="passed"`, the script checks the following items. All rational arithmetic uses `fractions.Fraction`, without floating-point approximations.

1. The recorded `source_sha256` values must exactly match the SHA-256 hashes of the six core source files `main.py`, `parameters.py`, `residues.py`, `terminals.py`, `arithmetic.py`, and `reporting.py`. These files are read from the current package or from the explicitly selected source archive.
2. The bases, paper target, internal verification target, period data, scales, partition intervals, and settings for positive bases must agree entry by entry with the current fixed parameters. The numbers of scales, intervals, endpoints, and terminal requests must be complete. The `chunk` parameter may take the positive integer value explicitly recorded for that run.
3. The height counts for each size bound must be nonnegative integers within the permitted height range. The count at height zero must be 1, the total must equal the modulus \(P\), and the height-distribution CSV and JSON must agree.
4. The script recomputes the exponential moments and finite tail bounds from the saved height distributions, then checks `moments.json` and `tails.csv`. It also recomputes the scale coefficients from the saved terminal lower bounds and checks `scales.csv`.
5. Using the fixed choices on each partition interval and the resulting scale coefficients, the script recomputes every individual contribution, numerator, and normalized value at every endpoint. It checks `endpoints.json` and its agreement with `endpoints.csv`. Finally, it verifies that the lower bound in `summary.json` equals the minimum of all endpoint values and strictly exceeds the internal target.

Records with status `running`, `failed`, `interrupted`, or any other status except `passed` are listed individually and are not checked as passed records. Other cases still running do not prevent checks of completed cases. If any record marked `passed` is missing, damaged, or inconsistent, the script reports its directory and exits with a nonzero status. An unparseable `summary.json` is also an error. The script continues checking the other records.

Exit status 0 means that no record inconsistency was found. Status 1 indicates an inconsistency, a missing results directory, or an invalid selected source archive. Command-line usage errors usually return 2. If there are no results, or all runs are incomplete, the script explicitly states that there are no completed records to check and may still return 0. Thus, exit status 0 alone does not establish that all eight cases passed. Read the individual statuses and the final counts as well.

**This check verifies only the consistency of saved records. It does not replace constructing the residue families and terminal lower bounds from scratch.** A height-count sum of \(P\) checks only the recorded coverage count, without rechecking an actual representation for every residue. The terminal lower bounds are used here as saved data, without reconstructing their representation families or interval families. Source hashes identify the source version declared by a run record. They are not signatures authenticating the record. The optional `joint_counts.npz` is outside the scope of this script. Full reproduction still requires rerunning `main.py` as described in the main instructions.

For the mathematical roles and derivation of the parameters, see [01_parameters.md](01_parameters.md). The checker itself is [check_saved_results.py](../check_saved_results.py).
