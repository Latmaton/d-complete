# Mathematical results and their implementation

This package supports eight lower-density estimates in the paper. Its overall basis is the quantitative criterion `thm:density-criterion`, together with the scale proposition `prop:ds-scales` and the endpoint proposition `pro-density-improve`. The programs construct and check the finite data required by those results.

Each note defines its mathematical objects and proves the recurrence before describing the code. Successful execution does not replace a proof of the recurrence. Conversely, a proof of the recurrence does not imply that a particular numerical choice exceeds the proposed target.

## 1. Correspondence with the paper

| Required object or condition | Mathematical justification | Implementation | Output |
| --- | --- | --- | --- |
| One bounded valid chain for each residue | Completeness of the search and set recurrences in the [residue notes](02_residues.md) | `residue_counts` | `residue_heights.csv` |
| Joint classification by height and size in the positive case | Section 5 of the same notes | `positive_families` | `joint[B][a,b]` in memory, optionally saved |
| `c<1` and the size constant `A` | `def-c`, `def-A` in the paper and the [exact arithmetic notes](04_density.md) | `constants` | `constants.json` |
| `g<1` and a tail bound at the end of a period | `lem:ds-tail` in the paper and the moment and period recurrences in the exact arithmetic notes | `tails` | `moments.json`, `tails.csv` |
| Positive terminal lower bounds | Lemmas, propositions, and Appendices A and B of the [positive terminal notes](03_positive_terminals.md) | `seed_intervals`, `next_intervals`, `pack_intervals` | `terminal_lower_bounds.csv` |
| Signed terminal lower bounds | Window and closure arguments in the [signed terminal notes](03_signed_terminals.md) | `signed_terminal_counts` | The same terminal output file |
| Safe scales `(a_i,b_i)` | `prop:ds-scales` in the paper and the rounding argument in the exact arithmetic notes | `scale_coefficients` | `scales.csv` |
| Endpoint inequalities on each interval | `pro-density-improve` in the paper and the concavity argument in the exact arithmetic notes | `density` | `endpoints.csv/json`, `summary.json` |

Every function in this table appears in the corresponding module in the package root. Appendix C of the positive terminal notes also identifies the arrays, indices, integer types, and rounding directions.

## 2. Why these data suffice

The residue families reduce residue choices to finitely many selected chains. Their height distributions give exponential moments for uniformly chosen residues and hence uniform tail bounds. The terminal programs give proved lower bounds for `T_h(X)`. Combining these two kinds of bounds with nonnegative weights, the scale proposition in the paper gives counting bounds at selected scales. The endpoint proposition then combines those scales into a lower asymptotic density estimate.

The directions of the estimates are essential. Tail bounds and scale abscissae may be rounded upwards. Terminal counts and scale ordinates must be rounded downwards. The implementation preserves these directions. In particular, the final positive terminal bound may be rational and need not be rounded again.

The deduction of asymptotic density from these finite data uses the representation concatenation, injectivity, size estimates, and division by powers of `q` proved in the paper. The program does not identify a finite experimental proportion with an asymptotic density, nor does it treat an unsuccessful finite search as a proof of nonrepresentability.

## 3. Flow of data

`main.py` reads `parameters.py`, checks the inexpensive conditions, and then performs four stages: residue construction, moments and tails, terminal bounds, and scales and endpoints. Each completed stage writes its results. The positive joint counts pass from the residue module to the terminal module. The other cases do not require these counts.

The fixed data specify construction ranges, search bounds, periods, scales, and endpoint choices. They contain no successful height distributions, terminal counts, or final density conclusions. A full run computes these quantities afresh.

`export_parameters.py` generates the supplement from the same parameter source, avoiding a separately maintained handwritten copy of the long tables. Running `python export_parameters.py --check` compares the exported files with their expected contents. Table consistency and mathematical verification are different checks, and both are recorded.

## 4. Distinctions used throughout

- The height of a rational chain, the minimum d-height of an integer representation, and the terminal scale height `H` are different objects. Sections 4 and 5 of the residue notes justify the transformations between them.
- In the positive construction, `cap` bounds the scanned integer numerators. The program only needs a valid chain for every residue within that range. It does not claim to find the optimal family among unrestricted numerators.
- A signed window construction need only produce a valid terminal lower bound. Intermediate integers outside the window may be discarded, so an unmarked integer need not be unrepresentable.
- Representations coming from different size bounds `B` may overlap. The positive terminal recurrence takes their maximum. For a fixed `B`, injectivity separates the images belonging to distinct residues and permits their addition.
- The final value must strictly exceed the target. The internal target `9001/10000` for `(5,2)` is slightly stronger than the paper's `9/10`.

## 5. Relation to the supplied material

The original `PARAMETERS`, `SCALES`, and `INTERVALS` from `data.py` are retained entry by entry in `parameters.py`. Computational ranges, including chunk size, window radii, positive grids, and height bounds, are collected in the same file. The package adds a common entry point, exact output files, explicit checks, complete parameter tables, and independent small enumerations.

The positive terminal notes follow the revised Lemma–Proposition exposition selected by the author. Appendix A proves the minimum odd-height recurrence and the initial interval count. Appendix B proves the equivalent summation grouped by equal shifts. Appendix C connects the results with functions, arrays, and integer ranges. The other notes likewise present proofs before implementation details.

The English edition changes the language of the exposition and the programs' explanatory text. Earlier run records remain attached to their original source snapshot, rather than being relabelled as runs of a new source hash. The [validation record](../validation/README.md) distinguishes those complete computations from checks performed on the English edition.
