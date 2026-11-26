# Terminal counting lower bounds for signed bases

This section explains the signed terminal construction in `terminals.py`. We first prove two transformations of representations, then define a recurrence in a finite window and describe its array implementation. Intermediate integers outside the window may occur in representations of smaller integers. The program therefore counts only the integers it constructs. It does not classify unmarked integers in the window as unrepresentable.

Let $p,q$ be coprime integers with $|p|,|q|>1$. A nonempty $d$-representation has the form

$$
n=\sum_{i=1}^{t}p^{a_i}q^{b_i},\qquad
0\leq a_1<\cdots<a_t,\qquad
b_1>\cdots>b_t\geq0.
$$

Since $p,q$ are coprime, these exponent conditions ensure that no summand divides another. Negative signs do not change divisibility. We call $1+b_1$ the d-height of the representation. For an integer $h\geq0$ and a real number $X\geq0$, define

$$
T_h(X)=\#\{n\in\mathbb Z:|n|\leq X,\ q\nmid n,
\ n\text{ has a nonempty }d\text{-representation of d-height at most }h\}.
$$

In particular, $T_0(X)=0$. The construction allows the empty sum $0$, but this is excluded from $T_h(X)$ because $q\mid0$. Throughout this section, the window radius $L$ and the requested endpoints $X$ are nonnegative integers with $X\leq L$.

## Lemma 1

If $u$ has a nonempty $d$-representation of d-height at most $h$, then so does $pu$. If $h\geq1$ and $u$ has a nonempty $d$-representation of d-height at most $h-1$, then

$$
pu+q^{h-1}
$$

has a nonempty $d$-representation of d-height at most $h$. The second assertion also holds for the empty sum $u=0$.

**Proof.** Multiplication by $p$ increases every $p$-exponent by $1$ and preserves the $q$-exponents and their order. It therefore preserves the representation height and the antichain condition.

For the second assertion, every summand of $pu$ has $p$-exponent at least $1$ and $q$-exponent at most $h-2$. The new summand $q^{h-1}$ has $p$-exponent $0$ and $q$-exponent $h-1$. Placing it first preserves the strict increase of the first exponent sequence and the strict decrease of the second. When $u=0$, use the singleton representation $q^{h-1}$. This proves the lemma.

## Proposition 1

Fix the window $[-L,L]$ and put $\mathcal S_0=\{0\}$. For $h\geq1$, first set

$$
\mathcal E_h=\mathcal S_{h-1}\cup
\{pu+q^{h-1}:u\in\mathcal S_{h-1},\ |pu+q^{h-1}|\leq L\},
$$

and then define

$$
\mathcal S_h=
\{p^jv:v\in\mathcal E_h,\ j\geq0,\ |p^jv|\leq L\}.
\tag{1}
$$

For every $h\geq0$, each nonzero element of $\mathcal S_h$ has a nonempty $d$-representation of d-height at most $h$. Consequently,

$$
T_h(X)\geq
\#\{n\in\mathcal S_h:|n|\leq X,\ q\nmid n\}
\qquad(0\leq X\leq L).
\tag{2}
$$

**Proof.** We use induction on $h$. The assertion holds for $\mathcal S_0$, which has no nonzero elements. Suppose that it holds for $h-1$. Every nonzero element already in $\mathcal S_{h-1}$ retains its valid representation when the allowed height increases. By Lemma 1, each added integer $pu+q^{h-1}$ also has a valid representation of height at most $h$. Repeated application of the first assertion of Lemma 1 shows that each nonzero element in (1) has such a representation. Discarding images outside the window does not invalidate any retained representation. This proves the induction step, and (2) follows.

We do not claim that $\mathcal S_h$ contains every representable integer in the window of height at most $h$. Reversing a representation may pass through an integer outside the window, in which case that construction is not included in this recurrence. The numerical lower bound uses only marked elements, so it does not require this possibility to be excluded.

## Lemma 2

The set in (1) can be obtained by the following finite procedure. First mark $\mathcal E_h$. Then examine the marked nonzero integers $u$ in the successive ranges

$$
|p|^j\leq |u|\leq |p|^{j+1}-1
\qquad(j=0,1,2,\ldots)
$$

Mark $pu$ whenever $|pu|\leq L$. Stop when $|p|^j>\lfloor L/|p|\rfloor$.

**Proof.** Each newly marked integer is obtained by multiplying a marked integer by $p$, so it belongs to the set in (1).

Conversely, suppose that $v\in\mathcal E_h$, $v\ne0$, and $p^kv$ belongs to the window. Then

$$
|v|<|pv|<\cdots<|p^kv|\leq L.
$$

Thus all intermediate terms also belong to the window. If a source satisfies $|p|^j\leq|u|<|p|^{j+1}$, its image satisfies $|pu|\geq|p|^{j+1}$. Hence the image lies in a later absolute-value range, rather than the current or an earlier range. The stated order therefore marks $pv,\ldots,p^kv$ in turn and finds every nonzero element of (1). Multiplication by $p$ leaves zero unchanged, and zero was already marked at the start. This proves the lemma.

The same argument shows that the order in which positive and negative sources are examined does not affect the result. At each step, the code only reads marks in the current absolute-value range and writes their images into later ranges.

## Lemma 3

If $q\nmid n$ and $n$ has a $d$-representation, then the smallest $q$-exponent in that representation is $0$. In particular, when $(p,q)=(-2,3)$, we have

$$
n\equiv1\pmod3.
$$

**Proof.** The smallest $q$-exponent $b_t$ occurs only once. After dividing the representation by $q^{b_t}$, the resulting integer is congruent to $p^{a_t}$ modulo $|q|$ and is therefore coprime to $q$. If $b_t\geq1$, the original integer is divisible by $q$, contrary to the hypothesis. Thus $b_t=0$ and $n\equiv p^{a_t}\pmod{|q|}$. For $p=-2,q=3$, we have $p\equiv1\pmod3$, which gives the assertion.

## Proposition 2

When $(p,q)=(-2,3)$, it suffices to store the integers in the window congruent to $1\pmod3$, without losing any integer constructed by Proposition 1 that is not divisible by $3$. More precisely, put $\mathcal P_0=\varnothing$. For $h\geq1$, first retain $\mathcal P_{h-1}$ and add the integer $1$ if it lies in the window. When $h\geq2$, also add the integers

$$
-2u+3^{h-1}\qquad(u\in\mathcal P_{h-1})
$$

that lie in the window. Finally, take the closure under multiplication by $-2$ as in Lemma 2. Then

$$
\mathcal P_h=\{n\in\mathcal S_h:3\nmid n\}
\qquad(h\geq0).
\tag{3}
$$

**Proof.** Both sides are empty when $h=0$. When $h=1$, Proposition 1 produces $1$ from $0$ and then takes the closure under multiplication by $-2$. Its elements not divisible by $3$ are therefore exactly the integers $(-2)^j$ in the window, in agreement with the construction above.

Suppose that $h\geq2$ and that the assertion holds for $h-1$. Since $3\mid3^{h-1}$,

$$
3\nmid(-2u+3^{h-1})\iff3\nmid u.
$$

Thus the multiples of $3$ in the previous set cannot produce new integers not divisible by $3$. Multiplication by $-2$ also preserves whether an integer is divisible by $3$. Retaining only integers not divisible by $3$ in each operation of Proposition 1 therefore gives precisely the present recurrence. If $L\geq1$, the integer $1$ has already been retained, and if $L=0$, it is not added. Repeatedly adding it does not change the equality. This completes the induction.

By Lemma 3, every integer in $\mathcal P_h$ is congruent to $1\pmod3$. This also follows directly from the recurrence. We have $1\equiv1\pmod3$, multiplication by $-2$ preserves this residue class, and $-2u+3^{h-1}\equiv u\pmod3$ for $h\geq2$. Thus the compressed array needs only this one residue class. Equality (3) compares two constructions in a finite window. It is not a completeness assertion about all representable integers.

## Mathematical objects and code

| Mathematical object or operation | Code | Explanation |
| --- | --- | --- |
| Window radius $L$ | `TERMINAL_LIMITS[(p,q)]` | Configured in `parameters.py` |
| $\mathcal S_{h-1}$ or $\mathcal P_{h-1}$ | `previous` | The previous Boolean array, read without modification during step $h$ |
| $\mathcal E_h$ and its closure under multiplication by $p$ | `current` | Copy the previous marks, add affine images, and complete the closure |
| $u\mapsto pu+q^{h-1}$ | `add_image` | Mark only images in the window |
| The absolute-value order in Lemma 2 | The `bound, upper` loop | Handle positive and negative sources separately, with new images in later absolute-value ranges |
| The right side of (2) | Return value of `signed_terminal_counts` | The dictionary `result[h,X]`, whose values are nonnegative integers |

The general case uses a Boolean array of length $2L+1$, with index $i$ representing the integer $-L+i$. The array records whether an integer has been constructed. If several representations reach the same position, that integer is still counted only once.

For $(-2,3)$, the program sets

$$
\texttt{first}=-L+((1+L)\bmod3),
$$

which is the smallest integer in the window congruent to $1\pmod3$. Index $i$ represents `first + 3*i`. Since this array contains no multiples of $3$, the program directly counts true entries in the requested subinterval. For the other signed cases, it first counts all true entries and then subtracts those at positions representing multiples of $q$.

`add_image` reads true positions in batches of size `chunk`. Batching only limits the length of temporary arrays. It changes neither the source interval nor its images. All necessary checks use explicit exceptions. The program refuses to continue if the window is smaller than a requested endpoint, a requested height or endpoint is invalid, the bases are not coprime, or an intermediate integer exceeds the range of the NumPy integer type being used.

The power $q^{h-1}$ is first computed as a Python integer. The program then checks whether the image of the entire source interval intersects the window. If the two are disjoint, it skips that operation. If images must be computed in NumPy arrays, the program checks the integer range before multiplication. It does not make decisions using values that have overflowed.

### The two fixed instances involving $7$

For $(p,q)=(7,-2)$ or $(-7,-2)$, the scales in this program request only $(h,X)=(1,7)$. The singleton representations $1$ and $p$ are two distinct integers of absolute value at most $7$. Neither is divisible by $q=-2$, and both have d-height $1$. Hence

$$
T_1(7)\geq2.
$$

`terminal_counts` checks that the requested pair is exactly this one and directly returns the lower bound $2$. Neither case requires a window array.

## Scope of verification

`check_small.py` uses set constructions in small windows to check Boolean array indices and closure operations. It also checks height zero, empty requests, invalid requests, and integer boundaries. The window construction provides a terminal counting lower bound. The density conclusion additionally requires the separate scale and endpoint inequality checks in `arithmetic.py`.

The actual commands, windows, requested ranges, and outcomes are recorded in the package's `validation/` directory. The proof in this section does not depend on whether a small example happens to cover the window, and it does not treat unmarked integers as evidence of nonrepresentability.
