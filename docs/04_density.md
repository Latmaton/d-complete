# Tail bounds, scales, and density endpoints: proofs and exact verification

This note continues the construction of residue families and terminal counts. The program uses the resulting finite counts to compute tail bounds and scales, then checks the finitely many endpoints in the density criterion. The lemmas and propositions below establish the recurrences; `arithmetic.py` evaluates them exactly.

The relevant results in the manuscript are `lem:ds-tail`, `prop:ds-scales`, `lem:ds-division`, and `pro-density-improve`. This note supplies the finite-length and rounding arguments and relates each step to the implementation. The program checks the finite hypotheses required by the manuscript. The passage from these hypotheses to a lower density bound still uses the manuscript's representation constructions, injectivity arguments, size estimates, and density proposition.

## 1. Inputs, notation, and preliminary checks

Fix coprime integers $p,q$ with $|p|,|q|>1$. Put

$$
P=|p|^d,
$$

and fix a period of length $s$, denoted $(B_i,\ell_i)_{i=1}^s$. For each $i$, the residue program has constructed families $\mathcal C_i(r)$ satisfying the manuscript's conditions, with size bound $B_i$ and height $h_i(r)$. The zero residue uses the empty family, of height $0$. Write

$$
n_i(h)=\#\{0\le r<P:h_i(r)=h\}.
$$

The same selected family is used for every occurrence of a given $B$. The code therefore stores `hist[B][h]` and takes $n_i(h)=\texttt{hist}[B_i][h]$. Before using the histograms, the program checks

$$
n_i(h)\in\mathbb Z_{\ge0},\qquad n_i(0)=1,
\qquad\sum_hn_i(h)=P.
\tag{1}
$$

Equation (1) checks the consistency of the counts. The validity and full residue coverage of the families follow from the preceding construction and its verification; they do not follow merely from a numerical table satisfying (1).

The contraction and size constants in the manuscript are

$$
c=\frac{|q|^{\ell_1+\cdots+\ell_s}}{P^s}<1,
\qquad
A=\frac1{1-c}\sum_{i=1}^s
\frac{B_i|q|^{\ell_i+\cdots+\ell_s-1}}{P^{s-i+1}}.
\tag{2}
$$

Also put

$$
\kappa=\begin{cases}1,&p>0,\ q>0,\\2,&\text{otherwise},\end{cases}
\qquad
\theta=\frac{\kappa\operatorname{ord}_{|q|}(p)}{|q|}.
\tag{3}
$$

`constants(params)` evaluates (2) and (3) exactly. These quantities require neither residue histograms nor terminal counts and should be checked before large computations. `read_inputs(p,q)` calls this check after verifying the following conditions:

- Integer parameters have the required types and ranges, including $d,\ell_i,B_i,m_0,t_{\max},J,H,n,e$. Rational inputs must be integers or `Fraction` objects; floating-point values are not accepted.
- The scale table is nonempty and $X>0$. In the positive case, $e$ selects an available level and $X$ lies on that level's valid grid. In signed cases, $X$ is an integer within the computation window.
- Each interval has exactly $J$ choices $(i,k)$, where $i$ is a valid scale index starting at $1$, and $k$ is an integer.
- The interval endpoints are positive and strictly increasing, adjacent intervals meet, and the final endpoint is exactly $|p|$ times the first.

Here `HEIGHT_LIMITS` bounds the residue search height, not the terminal scale height $H$. For example, some positive-case scales have $H=123$, which is consistent with a residue search height limit of $31$.

## 2. Lemma 1: a uniform tail bound from exponential moments

Extend the heights and $\ell_i$ with period $s$. For a residue tuple of length $j$, define

$$
E_0=0,\qquad E_j=\max\{E_{j-1}+h_j(r_j)-\ell_j,0\}.
$$

Define the normalized tail count by

$$
\tau_j(t)=P^{-j}\#\{(r_1,\ldots,r_j):E_j\ge t\}.
$$

Thus $\tau_j(t)=1$ for all $t\le0$. Fix a rational number $z>1$ and define

$$
G_i=\frac1P\sum_hn_i(h)z^{h-\ell_i},
\qquad g=\prod_{i=1}^sG_i<1,
$$

$$
\Lambda_i=\max\left\{1,\max_{1\le a<s}
\prod_{b=0}^{a-1}G_{i-b}\right\}.
\tag{4}
$$

When $s=1$, omit the inner maximum and set $\Lambda_i=1$. Then, for every $j\ge0,t\ge1$, we have

$$
\tau_j(t)\le\min\{1,\Lambda_jz^{-t}\}.
\tag{5}
$$

**Proof.** Periodicity and $g<1$ give

$$
\Lambda_{i+s}=\Lambda_i,
\qquad \Lambda_i=\max\{1,G_i\Lambda_{i-1}\}.
$$

For $t\ge1$, fixing the last residue and grouping by its height gives the exact recurrence

$$
\tau_j(t)=\frac1P\sum_hn_j(h)
\tau_{j-1}(t+\ell_j-h).
\tag{6}
$$

Initially, $\tau_0(t)=0$. Suppose that $\tau_{j-1}(u)\le\Lambda_{j-1}z^{-u}$ holds for positive integers $u$. For $u\le0$, the same bound holds because $\tau_{j-1}(u)=1$, $\Lambda_{j-1}\ge1$, and $z^{-u}\ge1$. Substitution into (6) gives

$$
\tau_j(t)\le G_j\Lambda_{j-1}z^{-t}
\le\Lambda_jz^{-t}.
$$

Combining this with $\tau_j(t)\le1$ proves (5). ∎

`tails` computes all $G_i$ from the histograms and checks $g<1$ and $\Lambda_s\le C$. The constant $\Lambda_s$ uses the suffix products ending at the $s$th term:

$$
1,\quad G_s,\quad G_sG_{s-1},\quad\ldots,\quad G_s\cdots G_2.
$$

The code forms these products in reverse order. The subsequent tail recurrence uses the forward order $1,2,\ldots,s$. These orders serve different formulas and must not be interchanged. The contraction condition $c<1$ and the exponential-moment condition $g<1$ also serve different purposes and must be checked separately.

## 3. Proposition 2: a finite number of periods with upward rounding

Fix an integer $m_0\ge1$, a rational number $C\ge\Lambda_s$, and a positive integer $D_0$. The fixed inputs in this version use $D_0=10^9$. For $t\le0$, set $w^{(a)}(t)=1$. For $t\ge1$, put

$$
w^{(0)}(t)=\frac1{D_0}
\left\lceil D_0\min\{1,Cz^{-t}\}\right\rceil,
\tag{7}
$$

and, for $1\le a\le m_0s$, define recursively

$$
w^{(a)}(t)=\frac1{D_0}
\left\lceil\frac{D_0}{P}
\sum_hn_a(h)w^{(a-1)}(t+\ell_a-h)\right\rceil.
\tag{8}
$$

Then $w_t=w^{(m_0s)}(t)$ satisfies

$$
1=w_0\ge w_1\ge w_2\ge\cdots\ge0,
\qquad \tau_{Ns}(t)\le w_t\quad(N\ge m_0).
\tag{9}
$$

**Proof.** Fix any $N\ge m_0$. Lemma 1 and periodicity give

$$
\tau_{(N-m_0)s}(t)\le\min\{1,\Lambda_sz^{-t}\}
\le w^{(0)}(t).
$$

This starts at the end-of-period phase; it does not assume that the tail distribution at length $(N-m_0)s$ is zero. By (6) and periodicity, the next step uses the height distribution $n_1$, followed by $n_2,\ldots,n_s$. Since every $n_a(h)$ is nonnegative, substituting an upper bound into (6) preserves the upper-bound property, as does upward rounding. Induction yields

$$
\tau_{(N-m_0)s+a}(t)\le w^{(a)}(t)
\qquad(0\le a\le m_0s).
$$

Taking $a=m_0s$ gives (9). Updating through $m_0$ complete periods returns to the $Ns$ phase required by the manuscript's scale proposition.

The sequence in (7) is nonincreasing and lies in $[0,1]$. If the previous sequence has these properties, each summand in (8) is nonincreasing in $t$, and so is their nonnegative weighted average. Equation (1) puts this average in $[0,1]$. Rounding upward to the $D_0^{-1}$ grid keeps it in that range and preserves monotonicity. Together with the prescribed value for $t\le0$, namely $1$, induction proves all the stated properties. ∎

The program stores the integers $W_t=D_0w_t$. Each step of (8) therefore becomes

$$
W'_t=\left\lceil\frac{\sum_hn_i(h)W_{t+\ell_i-h}}P\right\rceil,
\qquad W_u=D_0\quad(u\le0).
\tag{10}
$$

The integer expression `-(-total // P)` implements upward rounding. No floating-point approximation is used. Replacing `ceil` by rounding to the nearest integer or rounding downward would not preserve the argument.

### Why a finite table suffices

If output is required only for $0\le t\le t_{\max}$, retain the initial table through

$$
L_0=t_{\max}+m_0\sum_{i=1}^s\ell_i.
\tag{11}
$$

If the current table ends at $L$, compute the next one only through $L-\ell_i$. Since $h\ge0$, every accessed index satisfies

$$
t+\ell_i-h\le t+\ell_i\le L.
$$

Nonpositive indices use the value $1$ directly. Thus no unavailable entry is needed. After $m_0s$ steps, the final index is exactly $t_{\max}$, which explains the extra entries in (11).

### Scales whose heights exceed the computed tail range

If a scale has $H>t_{\max}$, extend the computed bound by a constant:

$$
w_t=w_{t_{\max}}\qquad(t>t_{\max}).
\tag{12}
$$

The true tail probabilities are nonincreasing in $t$, so (12) remains a valid upper bound. The extended differences $w_v-w_{v+1}$ vanish for $v\ge t_{\max}$. Hence scale computation needs only `tail_terms=min(H,t_max)` terms. This truncation neither assigns zero to unknown tail probabilities nor discards any positive weight from the chosen bound.

## 4. Proposition 3: valid scales from terminal lower bounds

Following the manuscript, $T_h(X)$ counts distinct integers with $|n|\le X$, $q\nmid n$, and height at most $h$ in some $d$-representation. The terminal construction supplies exact nonnegative rational lower bounds

$$
0\le L_h(X)\le T_h(X).
$$

These may be integers or rational lower bounds stored after downward rounding. Put $M=\min\{H,t_{\max}\}$ and define

$$
\widehat\beta(H,X)=\sum_{v=0}^{M-1}
(w_v-w_{v+1})L_{H-v}(X).
\tag{13}
$$

By (9) and (12), every weight is nonnegative, so the true terminal count may be replaced term by term with $L$. The manuscript's `prop:ds-scales` then gives, for every $\varepsilon>0$ and all sufficiently large integers $K$,

$$
\mathcal R_0\big((A|q|^H+X+\varepsilon)|p|^K\big)
\ge\widehat\beta(H,X)|p|^K.
\tag{14}
$$

For a scale input $(H,X,n,e)$, define

$$
\alpha=\frac{A|q|^H+X}{|p|^n},
\qquad
 a=\frac{\lfloor D_0\alpha\rfloor+1}{D_0},
\qquad
 b=\frac1{D_0^2}
 \left\lfloor\frac{D_0^2\widehat\beta(H,X)}{|p|^n}\right\rfloor.
\tag{15}
$$

Use only scales that pass the check $b>0$. Then

$$
\mathcal R_0(a|p|^K)\ge b|p|^K
\quad\text{for all sufficiently large integers }K.
\tag{16}
$$

**Proof.** Equation (15) ensures $a>\alpha$, even when $\alpha$ itself lies exactly on the grid. Take

$$
\varepsilon=a|p|^n-A|q|^H-X>0.
$$

Substituting $K-n$ for $K$ in (14) makes the left side exactly $\mathcal R_0(a|p|^K)$ and the right side $\widehat\beta(H,X)|p|^{K-n}$. Since $b$ in (15) is rounded downward, the latter is at least $b|p|^K$. ∎

`scale_coefficients` evaluates (13) and (15) directly. The exponent `n` normalizes the scale; `e=count_exponent` selects the positive-case terminal table. They have different meanings: an input may have $n=56,e=57$, so they must not be identified. Nonnegative weights make weaker terminal lower bounds valid. If the tail sequence were not nonincreasing, some weights could be negative and the substitution argument would fail. The program therefore checks monotonicity explicitly before using the sequence.

## 5. Proposition 4: why interval endpoints suffice

Choose $J\ge1$ and a full multiplicative period

$$
0<x_0<x_1<\cdots<x_V=|p|x_0.
$$

For each interval $[x_{v-1},x_v]$ and each $0\le j<J$, fix a scale index $i(v,j)$ and an integer $k(v,j)$. Define

$$
F_v(y)=\sum_{j=0}^{J-1}
\left[|p|^{k(v,j)}b_{i(v,j)}-
\theta\left(|p|^{k(v,j)}a_{i(v,j)}-
\frac y{|q|^j}\right)_+\right],
\tag{17}
$$

where $(u)_+=\max\{u,0\}$. Put

$$
D=\min_{1\le v\le V}\min_{y\in\{x_{v-1},x_v\}}
\frac{F_v(y)}{\kappa(1-|q|^{-J})y}.
\tag{18}
$$

Then the manuscript's `pro-density-improve` gives

$$
\underline d(\mathcal R_{p,q})\ge D.
\tag{19}
$$

**Why the endpoints suffice.** For fixed $v,j$, the function

$$
y\longmapsto-\theta\left(|p|^{k(v,j)}a_{i(v,j)}-y/|q|^j\right)_+
$$

is concave, so $F_v$ is concave. Put $c_J=\kappa(1-|q|^{-J})$. By (18), both endpoints satisfy $F_v(x)\ge c_JDx$. For any point $y=(1-\lambda)x_{v-1}+\lambda x_v$ in the interval, concavity gives

$$
F_v(y)\ge(1-\lambda)F_v(x_{v-1})+\lambda F_v(x_v)
\ge c_JDy.
\tag{20}
$$

This uses the concavity of $F_v$ and the linear function $c_JDy$. It does not require or assert that the ratio $F_v(y)/y$ is concave.

**From interval estimates to lower density.** The manuscript's capacity estimate

$$
\mathcal R_0(Y)-\mathcal R_0(X)\le\theta(Y-X)+O(1)
\quad(0\le X\le Y)
$$

together with (16), gives the termwise lower bounds in (17) for $\mathcal R_0(y|p|^K/|q|^j)$. There are only finitely many choices, so the threshold for sufficiently large $K$ can be chosen uniformly. Coverage of a full multiplicative period ensures that every sufficiently large $x$ can be written as $x=y|p|^K$, with $x_0\le y<|p|x_0$. Combining this with (20) gives

$$
\sum_{j=0}^{J-1}\mathcal R_0(x/|q|^j)
\ge \kappa(1-|q|^{-J})Dx-O(1).
$$

Now use the manuscript's division identity

$$
\mathcal R(x)=\sum_{j=0}^{J-1}\mathcal R_0(x/|q|^j)
+\mathcal R(x/|q|^J)
$$

and take lower limits. Writing $\delta=\underline d(\mathcal R_{p,q})$, we obtain

$$
\delta\ge(1-|q|^{-J})D+|q|^{-J}\delta,
$$

Hence $\delta\ge D$. This also explains the factor $1-|q|^{-J}$ in the denominator of (18). ∎

For each interval, `density` uses the supplied choices $(i,k)$ and evaluates its two endpoints; it does not search for or optimize new choices. A shared endpoint is evaluated separately for each adjacent interval because the choices may differ. Finally, the exact fraction $D$ must **strictly exceed** `params['target']`. For $(5,2)$, the verification target is $9001/10000$ and the paper target is $9/10$; the two targets coincide for the other seven cases. Decimal displays are not used in this comparison.

## 6. Correspondence between mathematics, code, and certificates

| Proved mathematical relation | Function or input | Output from the current run |
|---|---|---|
| Period, parameter ranges, and full multiplicative period | `read_inputs`, `parameters.py` | `parameters.json` |
| Contraction, size, and capacity constants (2), (3) | `constants` | `constants.json` |
| Height distributions of valid residue families and (1) | Constructed by `residues.py`; counts checked again by `tails` | `residue_heights.csv`, `residue_heights.json` |
| Exponential moments and end-of-period $\Lambda_s$ (4) | `tails` | `moments.json` |
| Upward-rounded tail bounds (7)–(12) | `tails` | `tails.csv` |
| $L_h(X)\le T_h(X)$ | `terminals.py`; positive-case proof in `03_positive_terminals.md` | `terminal_lower_bounds.csv` |
| Valid scales (13)–(16) | `scale_coefficients` | `scales.csv` |
| Interval endpoint values and strict target comparison (17), (18) | `density(..., records=records)` | `endpoints.json`, `endpoints.csv` |
| Completion status, source files, and environment | `main.py`, `reporting.py` | `summary.json` |

In the endpoint JSON, `interval_index` and the scale index `i` start at $1$, while `j` starts at $0$. Each record stores `left,right,y,value` and the sum `numerator`. For every $j$, its `contributions` list records the selected `i,k`, the factor $|p|^k$ (`factor`), scale coefficients `a,b`, the deduction `penalty`, and the resulting `contribution`, allowing (17) to be checked term by term. Rational values are stored as exact integer or fraction strings.

The parameter-check mode verifies only formats, ranges, and inexpensive constants. It does not construct residue families, compute terminal lower bounds, or check (18). A `status: passed` record written after a successful full run means that all finite checks were completed in that run. Although the final endpoint formula can be recomputed quickly from existing numerical tables, the full verification also includes the preceding constructions, proofs, and the origin of the counts used in the run.
