# Verification results

This directory contains the eight completed verification runs. Each timestamped subdirectory contains the output for one pair of bases, with the bases and UTC start time in its directory name. All eight supplied summaries have `status="passed"` and record all four completed stages. The exact lower bounds and verification details are listed in [validation/README.md](../validation/README.md).

From the package root, check the records and require coverage of all eight cases with:

```bash
python check_saved_results.py --require-all
```

The checker reports which cases have complete, consistent records. A run is complete only when `summary.json` has `status="passed"`. Its `current_stage` is then `complete`, and `stages_seconds` contains all four stages. A `running`, `failed`, or `interrupted` record is not a completed verification. The current verification status is listed in [validation/README.md](../validation/README.md).

| File | Meaning |
| --- | --- |
| `summary.json` | Completion status, exact final lower bound, targets, stage times, environment, and source hashes |
| `parameters.json` | The full parameters used for this run |
| `constants.json` | Exact values of the constants in the density criterion |
| `residue_heights.json`, `residue_heights.csv` | Counts of residue representatives at each height |
| `moments.json`, `tails.csv` | Exponential moments and tail upper bounds |
| `terminal_lower_bounds.csv` | The constructed lower bound for each requested terminal count |
| `scales.csv` | The scale coefficients used in the endpoint estimates |
| `endpoints.json`, `endpoints.csv` | Exact endpoint calculations and their minimum |
| `joint_counts.npz`, when saved | The positive case's joint height and size counts |

Fractions in the JSON and CSV files are exact. Decimal displays are only aids to reading. The optional array file does not replace the other records and is not checked by `check_saved_results.py`.

Exit status `0` from the command above means that all eight cases have complete, consistent records and no incomplete or erroneous record was found. Status `3` means that completion or coverage is missing. Status `1` indicates malformed or inconsistent records. A check without `--require-all` may succeed for only the cases present, so always read the coverage list.

The saved-record checker recomputes arithmetic from recorded counts. To reconstruct the residue families and terminal counts afresh, use a separate output directory:

```bash
python main.py --all --output reproduced_results
python check_saved_results.py --results reproduced_results --require-all
```

See [the main README](../README.md) for memory requirements and [the saved-result notes](../docs/05_saved_results.md) for the precise scope of the checks.
