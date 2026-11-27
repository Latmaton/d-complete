# d-complete: verification programs and supplementary material

This package constructs the residue families and terminal counting bounds required for eight density estimates in the paper. It then checks the exponential moments, tail bounds, scale coefficients, and all endpoint inequalities using exact integer and rational arithmetic.

The mathematical notes establish each recurrence before explaining its implementation. The positive terminal construction is presented through lemmas and propositions, with supplementary proofs. Complete parameter tables, numerical records, and independent checks are included.

**Verification status: all eight cases were computed from scratch and passed on 27 November 2026.** The final `--require-all` check also passed. Each supplied run has `status="passed"` and records all four completed stages. See the case table in [validation/README.md](validation/README.md).

These notes were prepared with the assistance of ChatGPT, and the authors reviewed the mathematical logic underlying the programs.

## 1. Getting started

Use Python 3.10 or later and NumPy. No other third-party package is required. The programs use functions and ordinary loops, without classes. From the directory containing this README, run:

```bash
python -m pip install -r requirements.txt
python check_saved_results.py --require-all
```

The second command checks the saved records against the current source and requires complete records for all eight cases. Exit status `0` means that this check passed. Status `3` means that complete records are missing or a run is incomplete. Status `1` indicates a damaged or inconsistent record. The checker prints the covered and missing cases. It checks the arithmetic based on saved counts, rather than reconstructing the representations.

For the preliminary and small independent checks, run:

```bash
python main.py --list
python main.py --check-parameters
python export_parameters.py --check
python check_small.py
python check_saved_results_tests.py
```

These commands list the cases, check the fixed parameters, compare the supplementary tables with their source, test the recurrences against independent small enumerations, and test the record checker's handling of incomplete or damaged records. Their precise scope is recorded in [validation/README.md](validation/README.md).

The small checks remain enabled under `python -O`. A normal invocation of `check_small.py` additionally repeats only the invalid-input checks in an optimized subprocess.

To recompute one complete case, for example:

```bash
python main.py -2 -3 --output my_results
python check_saved_results.py --results my_results --pair -2 -3
```

To recompute all eight cases:

```bash
python main.py --all --output reproduced_results
python check_saved_results.py --results reproduced_results --require-all
```

The `--all` option follows the order in the parameter file, starting with the large case `(5,2)`. Running `(-2,-3)` first provides a short example of the complete workflow. The main program constructs its counts afresh and does not read previously saved results. Using a separate output directory keeps a new calculation distinct from the records supplied with the package.

Other useful commands are:

```bash
python main.py 5 2 --plan
python main.py 5 2 --chunk 100000 --output my_results
python main.py -2 3 --quiet
```

The `--plan` option checks parameters and reports the main array sizes. Changing `--chunk` changes the length of temporary batches, not the mathematical parameters or results. It does not reduce the size of the height tables that must remain in memory. The `--quiet` option retains only the four main stage messages. Use `--save-arrays` to save the positive case's joint counts as a compressed array file. This large intermediate file is omitted by default. Running `python main.py` without arguments displays help.

## 2. Mathematical notes and reading order

| File | Contents |
| --- | --- |
| [docs/00_overview.md](docs/00_overview.md) | Proof dependencies, program structure, and the relation between parameters and results |
| [docs/01_parameters.md](docs/01_parameters.md) | Meaning and use of each parameter, with the corresponding checks |
| [docs/02_residues.md](docs/02_residues.md) | Residue chains, completeness of the search, Boolean recurrences, and the positive construction |
| [docs/03_positive_terminals.md](docs/03_positive_terminals.md) | Positive terminal construction, including proofs of the minimum-height and grouped-shift recurrences |
| [docs/03_signed_terminals.md](docs/03_signed_terminals.md) | Signed terminal construction, multiplication closure, and compressed storage |
| [docs/04_density.md](docs/04_density.md) | Exponential moments, tail phases and rounding, scale coefficients, and endpoint inequalities |
| [docs/05_saved_results.md](docs/05_saved_results.md) | Consistency checks for saved results and their limitations |
| [supplement/README.md](supplement/README.md) | Complete parameter tables for all eight cases |
| [results/README.md](results/README.md) | Result files and how to check them |
| [validation/README.md](validation/README.md) | Checks performed, recorded results, and source hashes |

The source files have the following roles.

| File | Role |
| --- | --- |
| `parameters.py` | Single source of fixed mathematical parameters and computational bounds |
| `residues.py` | Residue families, height distributions, and positive joint counts |
| `terminals.py` | Terminal counting lower bounds |
| `arithmetic.py` | Exact checks of constants, tails, scales, and endpoints |
| `reporting.py` | Exact output files, environment details, and source hashes |
| `main.py` | Command-line entry point and coordination of the four stages |
| `export_parameters.py` | Generation and checking of complete Markdown and CSV parameter tables |
| `check_small.py` | Independent enumerations and boundary checks for the implemented recurrences |
| `check_saved_results.py` | Consistency checks for saved parameters, counts, scales, endpoints, and source hashes |
| `check_saved_results_tests.py` | Regression checks for record validation, completion status, and case coverage |

## 3. Cases and targets

| Program order `(p,q)` | Target in the paper | Target that the program must strictly exceed |
| --- | --- | --- |
| `(5,2)` | `9/10` | `9001/10000` |
| `(-2,3)` | `159/500` | `159/500` |
| `(-2,-3)` | `31/100` | `31/100` |
| `(-5,2)` | `13/30` | `13/30` |
| `(5,-2)` | `7/20` | `7/20` |
| `(-5,-2)` | `56/125` | `56/125` |
| `(7,-2)` | `1/500000000` | `1/500000000` |
| `(-7,-2)` | `1/2100000000` | `1/2100000000` |

Interchanging the bases leaves the represented set unchanged, but changes the direction of the construction. The command line accepts the reversed pair and uses the orientation listed above. These targets are thresholds to be checked, not preloaded conclusions.

The supplement contains all 71 periodic entries, 191 scales, 466 intervals, and 2426 choices within those intervals. No long table is shortened by an ellipsis. Each case includes a readable Markdown page and CSV files for the basic parameters, period, scales, intervals, and expanded choices. Rational values remain exact fractions.

## 4. Output of a complete run

Each run creates a new directory named after the bases and a UTC timestamp. Existing records are not overwritten.

| Output | Contents |
| --- | --- |
| `summary.json` | Status, current stage, targets, exact lower bound when complete, stage times, environment, and core source SHA-256 hashes |
| `parameters.json` | Parameters, scales, and interval choices actually used |
| `constants.json` | Exact values of `c,A,kappa,theta` |
| `residue_heights.csv/json` | Residue counts classified by size bound and height |
| `moments.json` | Moment product `g` and the end-of-period constant `Lambda` |
| `tails.csv` | Tail upper bounds at integer heights |
| `terminal_lower_bounds.csv` | Lower bound for every requested pair `(h,X)` |
| `scales.csv` | Each scale and its final safe coefficients `a,b` |
| `endpoints.csv` | Normalized density expression at both ends of each interval |
| `endpoints.json` | Exact endpoint values and the full contribution of each `j` |

A run has completed all four stages only when its summary has `status="passed"`. While it runs, `current_stage` identifies `residues`, `tails`, `terminals`, or `density`. A successful run ends with `current_stage="complete"`. If an exception or keyboard interruption is caught during one of the four stages, the program attempts to record the completed stages and the failure location. An initialization or output-write failure, or termination by the operating system, may leave an incomplete record. A status of `running` is not a successful verification. The program does not resume an interrupted calculation automatically.

All decisive comparisons use integers or `Fraction`. Decimal output is for reference only. Explicit exceptions perform the mathematical checks, so `python -O` does not disable them. The main program compares core source hashes before and after each full run to detect changes during execution.

To check all current records and require all eight cases, use:

```bash
python check_saved_results.py --require-all
```

To check one case, use:

```bash
python check_saved_results.py --pair -2 -3
```

Without a coverage option, the checker reports whichever cases it finds. A successful check of one case does not establish coverage of the other seven. Any incomplete record in the selected scope causes exit status `3`, even if a later run of the same case passed. Full details are in [docs/05_saved_results.md](docs/05_saved_results.md).

The saved-result checker recomputes the arithmetic from recorded height distributions and terminal lower bounds. It does not reconstruct the representations themselves. Complete reproduction uses `main.py`. Source hashes identify the source bytes associated with a record. They do not authenticate the record or replace a mathematical verification.

## 5. Computational size

The positive case uses an initial odd-integer table with `5**13 = 1220703125` positions. Its two single-byte height arrays alone occupy `2441406250` bytes, in addition to the joint counts and temporary arrays. Reducing `--chunk` reduces temporary memory but leaves these two main arrays unchanged. A complete run therefore requires several GB of memory and substantial time.

The positive height recurrence uses exact integer lookups to extract powers of two. It initializes its first step directly and copies only entries that the next step does not overwrite. These operations implement the recurrence proved in the mathematical notes.

For `(-2,3)`, the terminal window is `[-300000000,300000000]`. Only positions congruent to `1 mod 3` are stored, so the two main Boolean arrays occupy about 400 million bytes. The other four window constructions use radius `10000000`. The terminal bound for each case involving `7` follows directly from two single-term representations. All bounds are listed in the supplement.

These programs verify eight fixed numerical instances of the paper's criterion. After changing parameters, check that the mathematical constructions still apply, regenerate the tables, and repeat the relevant full verification. The [overview](docs/00_overview.md) explains the respective roles of proof and computation.
