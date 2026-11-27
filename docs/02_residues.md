# Residue families: mathematical constructions, finite recurrences, and implementation

This note describes the objects computed by `residues.py`. It first proves the results underlying the searches and recurrences, then explains their finite implementation. Subsequent computations use counts obtained by selecting one chain for each residue, rather than counts of all chains.

## 1. The required construction

Let $p,q$ be coprime integers with $|p|,|q|\geq2$. Fix a width $d\geq1$ and write $P=|p|^d$. A rational chain is a finite sequence

$$
\mathcal C=((\alpha_1,\gamma_1),\ldots,(\alpha_t,\gamma_t)),
\qquad
0\leq\alpha_1<\cdots<\alpha_t<d,
\quad
0\leq\gamma_1<\cdots<\gamma_t.
$$

Define its value and height by

$$
\sigma(\mathcal C)=\sum_{j=1}^t p^{\alpha_j}q^{-\gamma_j},
\qquad
a(\mathcal C)=1+\gamma_t\quad(t>0).
$$

The empty chain has value and height $0$. Since $q$ and $P$ are coprime, the residue of the rational number $q^{-\gamma}$ modulo $P$ is the multiplicative inverse of $q^\gamma$.

Given a size bound $B\geq0$, for each $0\leq r<P$ we seek a chosen chain $\mathcal C^{(B)}(r)$ satisfying

$$
\sigma^{(B)}(r):=\sigma(\mathcal C^{(B)}(r))\equiv r\pmod P,
\qquad |\sigma^{(B)}(r)|\leq B.
\tag{1}
$$

For $r=0$, choose the empty chain. Set $h^{(B)}(r)=a(\mathcal C^{(B)}(r))$. The program returns

$$
c_B(a)=\#\{0\leq r<P:h^{(B)}(r)=a\}.
\tag{2}
$$

Thus $c_B(0)=1$ and $\sum_a c_B(a)=P$ must hold. For $(p,q)=(5,2)$, size information is also needed; see Section 5.

**Independent verification of a chain.** `verify_family(p,q,d,r,terms,B)` checks the range and strict ordering of both exponent sequences term by term. It computes the exact rational sum using `Fraction` and separately computes its residue modulo $P$. It returns the height only if (1) and the empty-chain convention for residue zero both hold. Here `pow(q,-gamma,P)` computes a modular inverse, not a floating-point negative power.

## 2. Complete search with a size constraint

This section applies to all admissible bases. In this project, `signed_families` handles the cases other than the positive case and the two Boolean-recurrence cases; its name does not restrict its mathematical applicability.

### 2.1 Clearing denominators at a fixed height

Fix $H\geq1$ and consider only chains with $\gamma<H$. Multiplication by $q^{H-1}$ gives the integer

$$
T(\mathcal C)=q^{H-1}\sigma(\mathcal C)
=\sum_j p^{\alpha_j}q^{H-1-\gamma_j}.
$$

The size condition is therefore equivalent to

$$
-N\leq T(\mathcal C)\leq N,
\qquad N=\left\lfloor B|q|^{H-1}\right\rfloor.
\tag{3}
$$

This equivalence also holds when $q<0$, since the absolute value of the multiplier is $|q|^{H-1}$. Clearing denominators is used only for the size comparison; the congruence condition is still computed for $\sigma(\mathcal C)$.

### 2.2 Recurrences for the extrema of the integer sum

Keep $H$ fixed and set $w_g=q^{H-1-g}$. To define $\mathcal T_{e,g}$, restrict the exponents by

$$
0\leq\alpha_1<\cdots<\alpha_t<e,
\qquad g\leq\gamma_1<\cdots<\gamma_t<H
$$

Take all corresponding integer sums $\sum_jp^{\alpha_j}w_{\gamma_j}$, including the empty sum $0$. Put

$$
L_e(g)=\min\mathcal T_{e,g},\qquad
U_e(g)=\max\mathcal T_{e,g}.
$$

**Lemma 1.** If $e=0$ or $g=H$, then $L_e(g)=U_e(g)=0$. Otherwise, the following set identity holds:

$$
\mathcal T_{e,g}
=p\mathcal T_{e-1,g}
\ \cup\ \bigl(w_g+p\mathcal T_{e-1,g+1}\bigr)
\ \cup\ \mathcal T_{e,g+1}.
\tag{4}
$$

**Proof.** A chain with no $\alpha=0$ term gives a sum in the first set, as seen by decreasing every $\alpha$ by $1$. A chain containing $(0,g)$ gives a sum in the second set: remove that term and decrease every remaining $\alpha$ by $1$. If a chain has an $\alpha=0$ term but does not contain $(0,g)$, its least $\gamma$ is greater than $g$, so its sum belongs to the third set. Each inverse operation restores a valid chain, proving the identity. The three sets may overlap; this does not affect their union or its extrema. This proves the lemma.

Taking the minimum and maximum of the three sets in (4) gives the recurrences for `low[e][g]` and `high[e][g]`. In particular, when $p<0$, multiplication by $p$ reverses the endpoints, so the program sorts them before use. Computing in increasing $e$ and decreasing $g$ ensures that all required values on the right are already available.

These are the exact extrema, but not every integer in $[L_e(g),U_e(g)]$ need be representable. The interval is used only to discard impossible states; existence is decided by the search below.

### 2.3 Search states and branches

`dfs(e,residue,g,lo,hi)` searches for a normalized remaining chain with $p$-exponents below $e$, denominator exponents at least $g$, rational sum congruent to `residue` modulo $|p|^e$, and integer sum in $[\mathrm{lo},\mathrm{hi}]$. The number of factors of $p$ already removed is $d-e$, so a selected term with current exponent $\alpha=0$ has original exponent $d-e$.

First intersect the target interval with $[L_e(g),U_e(g)]$. An empty intersection means that no solution exists. If the residue is zero and the interval contains zero, the empty chain is a solution. The remaining states split into two cases according to whether $\alpha=0$ occurs.

* If $r\equiv0\pmod{|p|}$, the chain cannot contain an $\alpha=0$ term: such a term would be a unit modulo $|p|$, while all other terms would vanish, contradicting $r\equiv0$. Divide all terms by $p$, reducing the width from $e$ to $e-1$.
* If $r\not\equiv0\pmod{|p|}$, every solution must contain exactly one $\alpha=0$ term. The program tries every possible denominator exponent $g\leq\gamma<H$. It continues only when $r-q^{-\gamma}$ is divisible by $|p|$. Removing this term and dividing by $p$ reduces the width to $e-1$, and all subsequent denominator exponents must exceed $\gamma$.

Division when reducing the modulus from $|p|^e$ to $|p|^{e-1}$ has an exact meaning: if $u$ is divisible by $p$, then

$$
px\equiv u\pmod{|p|^e}
\iff x\equiv u/p\pmod{|p|^{e-1}}.
\tag{5}
$$

This is also valid for negative $p$. The code uses `// p` only after divisibility has been checked, so the direction of rounding does not affect the quotient. Interval division is handled separately. If $p>0$, the condition $\ell\leq px\leq u$ is equivalent to

$$
\left\lceil\ell/p\right\rceil\leq x\leq\left\lfloor u/p\right\rfloor;
$$

if $p<0$, first replace the interval by $[-u,-\ell]$ and then divide by $|p|$. `divide_interval` performs these ceiling and floor operations entirely with integer division.

**Proposition 1.** For fixed $r,B,H$, `search(r,B,H)` returns a chain satisfying the conditions if and only if a valid chain of height at most $H$ exists.

**Proof.** Induct on the remaining width $e$. The empty-chain case and the case $e=0$ are handled correctly. The congruence modulo $|p|$ determines which of the two branches is possible, and the second branch tries every possible $\gamma$. Each congruence and interval transformation is an equivalence, and the remaining width strictly decreases, so the induction hypothesis applies. Lemma 1 excludes only intervals that contain no valid sum, hence no solution is lost. Reversing the transformations preserves the strict ordering of both exponent sequences and satisfies the original congruence and size constraints. This proves the proposition.

States proved infeasible are stored in `failures`. This set is local to one search with fixed $r,B,H$ and records the full state, so failures from other parameter choices are not reused. The code distinguishes the empty chain `()` from the failure value `None`.

### 2.4 Minimum height and output

For fixed $r,B$, the set of chains of height at most $H$ increases with $H$. The least feasible height can therefore be found by bisection. The program first verifies feasibility at `H_MAX`, then maintains an infeasible left endpoint and a feasible right endpoint until they are adjacent. Bounds `B` are processed in increasing order, so a chain found for the preceding bound gives a feasible starting point for the next one.

By Proposition 1, the result is the least height satisfying (1). Once a feasible height lies within `H_MAX`, larger heights cannot improve this minimum. The selected chain is then checked again from its exponent pairs by `verify_family` before being counted in $c_B(a)$. If any residue has no chain within `H_MAX`, the function raises `RuntimeError` rather than returning an incomplete count table.

## 3. Boolean recurrence when the size bound is automatic

### 3.1 Why the size state can be omitted

**Lemma 2.** Every rational chain of width $d$ satisfies

$$
|\sigma(\mathcal C)|\leq |p|^{d-1}.
\tag{6}
$$

**Proof.** The empty chain satisfies the bound. For a chain with $t\geq1$ terms, strict ordering gives $t\leq d$, $\alpha_j\leq d-t+j-1$, and $\gamma_j\geq j-1$. Hence

$$
\begin{aligned}
|\sigma(\mathcal C)|
&\leq\sum_{j=1}^t |p|^{d-t+j-1}|q|^{-(j-1)}\\
&=|p|^{d-1}\sum_{j=1}^t
\frac1{|p|^{t-j}|q|^{j-1}}\\
&\leq |p|^{d-1}\frac{t}{2^{t-1}}
\leq |p|^{d-1}.
\end{aligned}
$$

The last inequality uses $t\leq2^{t-1}$, which follows by induction on $t$. This proves the lemma.

Consequently, if every requested $B$ is at least $|p|^{d-1}$, the existence recurrence already enforces the size condition. The project cases $(7,-2)$ and $(-7,-2)$ both use $d=8$ and $B=823543=7^7$, so they satisfy this condition exactly. `boolean_families` checks it before computation.

### 3.2 The exact recurrence from three cases

Let $E_{e,H}(r)$ mean that there exists a chain of width $e$ and height at most $H$ whose rational sum, reduced modulo $|p|^e$, equals $r$. No additional size constraint is imposed here. The initial values are

$$
E_{e,0}(r)=[r=0],\qquad E_{0,H}(0)=\mathrm{true}.
$$

Square brackets indicate the truth value of a statement. All arguments below are reduced modulo the appropriate modulus, and a quotient subject to a divisibility condition is used only when that condition holds.

**Proposition 2.** For $e,H\geq1$,

$$
\begin{aligned}
E_{e,H}(r)={}&E_{e,H-1}(qr)\\
&\lor\bigl([p\mid r]\land E_{e-1,H}(r/p)\bigr)\\
&\lor\bigl([p\mid r-1]\land E_{e-1,H-1}(q(r-1)/p)\bigr).
\end{aligned}
\tag{7}
$$

**Proof.** If every $\gamma\geq1$, multiplying the rational sum by $q$ lowers the height bound by $1$, giving the first term. If every $\alpha\geq1$, division by $p$ lowers the width by $1$, giving the second term. If neither case holds, both $\alpha=0$ and $\gamma=0$ occur. Since both exponent sequences increase strictly, they must occur together in the first term, $(0,0)$. After removing this term of value $1$, both exponents of every remaining term are at least $1$. Division by $p$ followed by multiplication by $q$ gives the third term. Each inverse operation restores a valid chain, so the disjunction of these three existence conditions is equivalent to the original condition. This proves the proposition.

The program computes in increasing $H$ and, for each $H$, in increasing $e$. The first and third terms read `previous`; the second reads the completed `current[e-1]`. The arrays are exchanged only after all $e$ have been computed, preserving every required value from the preceding height.

Write $A_H=\#\{r:E_{d,H}(r)=\mathrm{true}\}$. By definition, $A_H$ is nondecreasing. The number of residues first represented at this height is $A_H-A_{H-1}$, giving the minimum-height distribution. The program returns when $A_H=P$ and supplies this same distribution for every $B$. It raises an error if complete coverage has not been reached by `H_MAX`.

Here the exact existence recurrence proves that each residue has a chain; individual chains are not stored. Chunking the arrays changes only how many residues are processed together, not recurrence (7).

## 4. The minimum-height table for positive integers

The positive case uses the auxiliary function `odd_heights(limit,width,chunk,verbose)`. For a positive integer $n$, let $h_{\min,k}(n)$ denote the minimum d-height among representations of the form

$$
n=\sum_j5^{\alpha_j}2^{\beta_j},\qquad
0\leq\alpha_1<\cdots<\alpha_t<k,
\quad \beta_1>\cdots>\beta_t\geq0
$$

where the d-height is $1+\beta_1$; set the minimum to $+\infty$ if no representation exists. For the integer $0$, assign the auxiliary height $0$. This d-height of an integer representation and the height of a rational chain in Section 1 are distinct definitions and must not be interchanged.

A complete proof of the minimum-height recurrence appears in [Appendix A of the positive terminal-count note](03_positive_terminals.md). For odd $n$, its central formula is

$$
h_{\min,k}(n)=\min\left\{
h_{\min,k-1}(n),\
1+h_{\min,k-1}\!\left(\frac{n-5^{k-1}}2\right)
\right\},
\tag{8}
$$

where the second term is used only when $n\geq5^{k-1}$. Nonzero even integers are reduced to odd integers by

$$
h_{\min,k}(2^v u)=v+h_{\min,k}(u),\qquad u\text{ is a positive odd integer}
\tag{9}
$$

If the argument in the second term of (8) is zero, that term is assigned the value $1$ directly.

Index $i$ in the returned array corresponds to the positive odd integer $n=2i+1\leq\mathrm{limit}$. The two buffers remain separate: `previous` stores width $k-1$, and `current` stores width $k$. Only `previous` may be read when computing the second term; otherwise a height newly obtained in the same round could be reused, violating (8). Swapping the buffers after a full round preserves the recurrence without retaining a separate table for every $k$.

The implementation initializes width $1$ directly: among positive odd integers, only $1$ has such a representation, of height $1$. For each later width, write `offset=(5**(k-1)-1)//2`. At array index $i$, the nonnegative parent is exactly $i-\mathrm{offset}$. Entries before this offset retain the preceding width's value; every later entry is overwritten from the separate preceding buffer.

To factor a parent into a power of $2$ and an odd integer, the program uses an integer table of $v_2(u)$ for $1\le u<2^{16}$. If the low 16 bits are zero, it adds $16$ to the table value for the high 16 bits. The supported limit is less than $2^{32}$, so every parent is at most $2^{31}-1$ and these two lookups suffice. The zero parent is treated separately as the singleton case. Intermediate heights use `uint16`, which holds even `255+32+1`; they are capped at the marker `255` before being written to the single-byte height buffer. This changes neither the recurrence nor the finite range and uses no floating-point arithmetic.


The implementation requires `limit < 2**32` and `width <= 20`. It uses one-byte `uint8` entries, with $255$ denoting no representation. In any representation of a positive integer in this range, every term is at most $n<2^{32}$, so its $2$-exponent is at most $31$ and its finite d-height is at most $32$. Thus $255$ cannot be confused with a genuine finite height. Addition is performed in a wider integer type before values are capped at this marker.

## 5. $(5,2)$: constructing residue families from odd numerators

### 5.1 Correspondence between integer representations and rational chains

**Lemma 3.** Fix $a\geq1$. Nonempty rational chains of height exactly $a$,

$$
\sigma=\sum_j5^{\alpha_j}2^{-\gamma_j},\qquad
\alpha_j<d,
$$

are in bijection with representations of positive odd integers $S$ of width $d$ and d-height at most $a$, through

$$
\sigma=\frac S{2^{a-1}},\qquad
\beta_j=a-1-\gamma_j.
\tag{10}
$$

**Proof.** The largest $\gamma$ in the rational chain is $a-1$. Reflection gives strictly decreasing nonnegative $\beta$, with least $\beta$ equal to zero. Exactly one term is odd, so the integer sum $S$ is positive and odd, and the largest $\beta$ is at most $a-1$. Conversely, a valid representation of an odd integer $S$ must have least $\beta=0$, since otherwise every term would be even. A d-height at most $a$ means that every $\beta\leq a-1$; reflection therefore gives nonnegative, strictly increasing $\gamma$ with maximum exactly $a-1$. The two operations are inverse to each other. This proves the lemma.

Thus, to decide whether a candidate $(a,S)$ yields a valid chain, it suffices to check

$$
S>0\text{ and odd},\quad h_{\min,d}(S)\leq a,
\quad S\leq B2^{a-1},
\quad r\equiv S(2^{a-1})^{-1}\pmod{5^d}.
\tag{11}
$$

The condition $h_{\min,d}(S)=a$ is not required. Even if $S$ has a smaller minimum d-height, the rational chain obtained from (10) still has height $a$, since its least $\beta=0$.

### 5.2 Finite search range and selection rule

`positive_families` first computes the height table for odd $1\leq S\leq\mathrm{cap}$. It then processes $a=1,2,\ldots,\mathrm{hmax}$ in order; at each fixed $a$, it scans odd $S$ in increasing order. A separate Boolean array `seen` is maintained for every $B$. Initially only residue zero is selected, using the empty chain. When (11) holds and the residue has not yet been selected, the pair $(a,S)$ is fixed for that residue.

Each residue is thus assigned its first occurrence in lexicographic order of $(a,S)$ within the candidate range. At fixed $a$, ordering by $S$ is the same as ordering by $\sigma$. The same `seen` array is used across all heights and chunks, so every residue is counted once. The odd $S$ in one chunk are distinct modulo $5^d$: the program ensures `2*chunk < P`, so the absolute difference of distinct numerators is less than $P$. Multiplication by $(2^{a-1})^{-1}$ preserves this distinction. A batch update therefore cannot add the same residue twice.

The centralized parameter values are

$$
d=11,\quad \mathrm{hmax}=31,\quad
\mathrm{cap}=369098752,\quad
\mathcal B=\{10,11,12,14,22\}.
$$

`cap` and the value $R$ below are passed explicitly from `POSITIVE_SETTINGS` in `parameters.py`. The parameter `cap` bounds the odd numerators in this finite search. A residue might have a chain of smaller height whose odd numerator exceeds `cap`, and the current search might miss it. The claim here is that the selected family is valid and covers every residue; a first occurrence within the truncation is not asserted to have globally minimum height without truncation. The later density proof requires only validity and coverage.

If any residue remains uncovered within `hmax`, the program raises `RuntimeError`. It returns only when `covered == P` and every entry of `seen` is true.

### 5.3 Joint height and size counts

Let $R=32768$. In addition to the height distribution, retain

$$
M_B(a,b)=\#\left\{0\leq r<P:
h^{(B)}(r)=a,\ 
\left\lfloor R\sigma^{(B)}(r)\right\rfloor=b\right\}.
\tag{12}
$$

The variable `counts[a,b]` represents $M_B(a,b)$. By (10), its bin index is computed entirely with integers:

$$
b=\left\lfloor\frac{RS}{2^{a-1}}\right\rfloor.
$$

The endpoint $\sigma=B$ gives $b=RB$ and is included, so the second array dimension has length $RB+1$. The empty chain contributes $M_B(0,0)=1$ separately. Row sums give the height distribution, and the counts must satisfy

$$
c_B(a)=\sum_bM_B(a,b),\qquad
\sum_{a,b}M_B(a,b)=P.
$$

Equation (12) is the joint count used in Proposition 1 of the [positive terminal-count note](03_positive_terminals.md). There, the bin intervals provide a proved interval transfer, so the full chain and exact $\sigma$ for each residue need not be retained.

## 6. Mathematical objects and their implementation

| Mathematical content | Code location | Returned or stored quantity |
|---|---|---|
| Verification of one chain in (1) | `verify_family` | Height after verification |
| Integer extrema in Lemma 1 | `geometry` inside `signed_families` | `low[e][g]`, `high[e][g]` |
| Complete search in Proposition 1 | Internal `search` and `dfs` | An exponent chain, or `None` |
| Height monotonicity and bisection | Residue loop in `signed_families` | `hist[B][a] = c_B(a)` |
| Lemma 2 and Proposition 2 | `boolean_families` | Existence arrays and minimum-height distribution |
| Equations (8) and (9), proved in Appendix A | `odd_heights` | Index $i$ corresponds to $h_{\min,k}(2i+1)$ |
| Lemma 3 and the finite selection rule | `positive_families` | Height distribution and `joint[B]` |
| Equation (12) | `counts[a,b]` | $M_B(a,b)$ |
| Selection of an applicable algorithm | `residue_counts` | `(hist, joint)`; `joint=None` when no joint counts are used |

The mathematical parameters used by `residue_counts(params,verbose=True)` are `p,q,d,B,H_MAX`. The optional computational parameter `chunk` defaults to the centralized `CHUNK`. Setting `verbose=False` only suppresses progress messages. Chunking the Boolean arrays and height tables changes neither the recurrence nor the order in which numerators are scanned.

## 7. Exact arithmetic and implementation limits

The search uses Python arbitrary-precision integers and integer or `Fraction` size bounds, without floating-point comparisons. The NumPy calculations use fixed-width types, so their ranges are checked before allocation and arithmetic:

* The Boolean recurrence checks that the bounds required for multiplying residues by $p,q$ fit in signed $64$-bit integers. Division by a negative base is guarded by divisibility conditions.
* The positive case checks $P\leq2^{32}-1$, so the upper bounds for each joint-count cell and the total coverage fit in unsigned $32$-bit integers.
* The residue product $S(2^{a-1})^{-1}$, binning product $RS$, and size comparison $B2^{a-1}$ are checked against the signed $64$-bit range before evaluation. Joint-count row sums use `uint64`.
* Section 4 justifies the finite height-table range and the marker $255$. Both `previous` and `current` buffers are retained throughout.

All essential conditions are checked with explicit exceptions, without relying on `assert` statements that `python -O` can disable. Invalid input raises `ValueError`; incomplete coverage within the specified range or inconsistent internal counts raise `RuntimeError`.

A total count of $P$ is a necessary coverage check, but does not by itself prove that the chains are valid. Validity follows, respectively, from exact verification of individual chains, the existence recurrence of Proposition 2, and Lemma 3 together with the proved minimum-height recurrence. These arguments and the total-count checks together verify the residue module.
