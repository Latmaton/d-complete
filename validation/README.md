# Validation record

This page distinguishes checks performed on the English edition from earlier full computations. Translating comments and messages changes source hashes, so the earlier summaries have not been relabelled with the new hashes.

## 1. Checks of the English edition

| Check | Result and record |
| --- | --- |
| `python main.py --check-parameters` | All eight cases passed the preliminary conditions. [Log](parameters.txt) |
| `python -O main.py --check-parameters` | The same checks passed with optimization enabled. [Log](optimized_parameters.txt) |
| `python export_parameters.py --check` | All generated tables match the parameter source. [Log](supplement.txt) |
| `python check_small.py` | Independent small enumerations and explicit failure cases passed. [Log](check_small.txt) |
| `python main.py -2 -3 --quiet` | One complete run of the English source passed. [Log](english_smoke.txt) |
| `python check_saved_results.py` | The new run's saved records are consistent with the current source. [Log](saved_results.txt) |
| Fixed parameter comparison | All eight parameter collections are exactly unchanged. [Record](parameter_equivalence.json) |
| English text and document structure | No Chinese characters remain in the documentation or current code. Links and blocks are valid. [Record](english_edition_check.json) |
| Core source comparison | Only explanatory text changed in the six core files. [Record](translation_check.json) |

The source comparison removes docstrings and allows translation of existing non-ASCII explanatory strings. It retains every identifier, operation, numeric literal, numerical string constant, and formatted expression. This checks the scope of the translation. It is not a new proof of the mathematics or a substitute for a full run.

The small tests cover residue searches, Boolean recurrences, positive height and size distributions, minimum odd heights, signed window counts, initial intervals with nonintegral endpoints, interval packing, direct double sums, complete-period tail bounds, and invalid arguments. A separate `python -O` subprocess confirms that the 26 explicit rejection cases remain active.

## 2. Earlier complete computations

The directory [recorded_results](recorded_results) contains the results produced before the English translation. Their exact numerical values and original source hashes are preserved. The six source files that produced them are stored at the root of [recorded_source.zip](recorded_source.zip). Its `original_logs/` directory contains the verbatim original console logs and comparison records.

The archived calculations used `python main.py 5 2 --output results --quiet` for the positive case. The other cases were run through `main.run_case(pair, verbose=False)` in the order `(-2,-3), (7,-2), (-7,-2), (-5,2), (5,-2), (-5,-2), (-2,3)`. Each successful run constructed the required residue and terminal counts afresh.

| Bases | Archived status | Paper target | Computed lower bound, decimal for reference |
| --- | --- | --- | --- |
| `(5,2)` | Incomplete | `9/10` | Not available |
| `(-2,3)` | Passed | `159/500` | `0.318518793479` |
| `(-2,-3)` | Passed | `31/100` | `0.319726606959` |
| `(-5,2)` | Passed | `13/30` | `0.442171083265` |
| `(5,-2)` | Passed | `7/20` | `0.366256893039` |
| `(-5,-2)` | Passed | `56/125` | `0.466835653355` |
| `(7,-2)` | Passed | `1/500000000` | `2.73693522415e-09` |
| `(-7,-2)` | Passed | `1/2100000000` | `6.47758412004e-10` |

Each original `summary.json` contains the exact rational lower bound, stage times, environment, and source hashes. Seven cases completed all finite checks and strictly exceeded their internal targets. The positive record only completed the residue and tail stages. Its preserved `running` status is an incomplete historical snapshot, not evidence that a process is currently active, and is not counted as a successful verification.

To check the archived records against their actual source version, run:

```bash
python check_saved_results.py --results validation/recorded_results --source-archive validation/recorded_source.zip
```

This check passed for all seven successful archived cases, with no inconsistencies. The incomplete positive case is explicitly skipped. See the [English audit log](recorded_results_check.txt). The checker recomputes the moments, tails, scales, and endpoint expressions from the saved counts. It does not reconstruct the residue chains or terminal families.

## 3. Interpretation

Passing the preliminary checks establishes only that the fixed parameters and inexpensive conditions are valid. Passing the small tests establishes agreement on the stated finite tests. A successful complete run verifies all finite conditions for that instance. The deduction of lower asymptotic density still uses the proofs in the paper and the mathematical notes.

The English edition has not been presented as a new full recomputation of all eight cases. It includes the unchanged earlier evidence, checks that the translation preserves the computational formulas, and records the tests and complete example actually run with the English source.

After changing mathematical parameters or algorithms, regenerate the tables and repeat the relevant verification. When checking saved records, explicitly select the source that generated them. A source hash identifies a version but is not a signature certifying the authenticity of a result.

If an argument or directory check fails before initialization, no summary may be created. Ordinary exceptions or keyboard interruptions write a failed or interrupted status. A process killed directly by the operating system may leave a `running` summary. Any status other than `passed` is not a successful verification. Interrupted computations are not resumed automatically.
