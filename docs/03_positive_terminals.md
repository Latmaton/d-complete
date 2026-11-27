# Terminal count lower bounds for $(5,2)$ and their implementation

This document first establishes finite recurrences for terminal count lower bounds, then explains how the program evaluates them. The main text follows the Lemma–Proposition order of the argument. Appendix A proves the recurrence for the initial height table, Appendix B proves the formula for grouping terms with the same shift, and Appendix C relates the mathematical results to the implementation.

The construction has three steps. First, minimum representation heights are computed over a finite range to obtain counts for the initial intervals. Next, representations are concatenated using residue families to transfer these counts to a larger range. Finally, pairwise disjoint intervals are selected and their counts are added to obtain terminal count lower bounds. The proofs below justify each recurrence, while the program evaluates the finite integer quantities involved.

We consider only the positive bases $(5,2)$. A nonempty $d$-representation of an integer is an expression

$$
n=\sum_{j=1}^t5^{\alpha_j}2^{\beta_j},\qquad
0\leq\alpha_1<\cdots<\alpha_t,\quad
\beta_1>\cdots>\beta_t\geq0.
$$

The **d-height** of this representation is $1+\beta_1$. Let $T_h(I)$ denote the number of positive odd integers in the interval $I$ that have a nonempty $d$-representation of d-height at most $h$, and write $T_h(X)=T_h((0,X])$. In particular, $T_0(I)=0$, and $T_h(I)$ is nondecreasing in $h$.

We use the notation from the residue construction. For $B\in\mathcal B=\{10,11,12,14,22\}$ and $0\leq r<P$, let

$$
\sigma^{(B)}(r)
=\sum_{(\alpha,\gamma)\in\mathcal C^{(B)}(r)}
5^\alpha2^{-\gamma},\qquad P=5^{11}.
$$

The exponent pairs in the chosen family satisfy

$$
0\leq\alpha_1<\cdots<\alpha_t<11,\qquad
0\leq\gamma_1<\cdots<\gamma_t,
$$

and

$$
\sigma^{(B)}(r)\equiv r\pmod P,\qquad
0\leq\sigma^{(B)}(r)\leq B.
$$

The **height** of a nonempty family is $h^{(B)}(r)=1+\max\gamma$. For $r=0$, set $\mathcal C^{(B)}(0)=\varnothing$, $\sigma^{(B)}(0)=0$, and $h^{(B)}(0)=0$. Below, $a$ denotes the height of a residue family, and $h$ denotes the upper bound on the d-height of an integer representation. These two heights have different definitions.

Fix

$$
Q=8192,\qquad R=32768,\qquad D_0=10^9,
$$

and set

$$
\nu_m=13+11(m-1),\qquad
W_m=400+2(m-1),\qquad 1\leq m\leq5,
$$

$$
(K_1,K_2,K_3,K_4,K_5)=(40,57,83,108,134).
$$

Then, for $m\geq2$,

$$
5^{\nu_m}=P5^{\nu_{m-1}},\qquad W_m=W_{m-1}+2.
$$

For $0\leq i\leq2Q$, define

$$
I_{m,i}=\left(\frac{(i-W_m)5^{\nu_m}}Q,\frac{i5^{\nu_m}}Q\right].
$$

We will construct nonnegative integers $V_m(h,i)$ for $0\leq h\leq K_m$ and prove that

$$
T_h(I_{m,i})\geq\frac{5^{\nu_m}}{D_0}V_m(h,i).
\tag{1}
$$

## Lemma 1.

For each positive odd integer $n<2\cdot5^{13}$, let

$$
c(n)=\left\lceil\frac{Qn}{5^{13}}\right\rceil.
$$

Let $h_{\min,14}(n)$ be the minimum d-height among all $d$-representations whose $5$-exponents are less than $14$. Set this value to $+\infty$ if no such representation exists. Appendix A gives and proves the recurrence for this minimum height. Since all terms are positive and $5^{14}>2\cdot5^{13}$, restricting the $5$-exponents to be less than $14$ excludes no representation in the present range. For $0\leq h\leq K_1$ and $0\leq i\leq2Q$, let

$$
V_1(h,i)=\left\lfloor\frac{D_0}{5^{13}}\#\{n:1\leq n<2\cdot5^{13},\ 2\nmid n,\ h_{\min,14}(n)\leq h,\ i-W_1<c(n)\leq i\}\right\rfloor.
$$

Then

$$
T_h(I_{1,i})\geq\frac{5^{13}}{D_0}V_1(h,i).
$$

**Proof.** Since $i-W_1$ and $i$ are integers, the definition of the ceiling gives

$$
i-W_1<c(n)\leq i\iff i-W_1<\frac{Qn}{5^{13}}\leq i\iff n\in I_{1,i}.
$$

Since $0\leq i\leq2Q$, every positive odd integer in $I_{1,i}$ satisfies $n\leq2\cdot5^{13}$. This upper bound is even, so in fact $n<2\cdot5^{13}$.

The condition $h_{\min,14}(n)\leq h$ ensures that every integer counted above has a nonempty $d$-representation of d-height at most $h$. Combining this with the interval equivalence just proved gives

$$
\frac{5^{13}}{D_0}V_1(h,i)\leq\#\{n\in I_{1,i}:n>0,\ 2\nmid n,\ h_{\min,14}(n)\leq h\}\leq T_h(I_{1,i}).
$$

This completes the proof.

## Lemma 2.

Fix $B\in\mathcal B$, $0\leq r<P$, and $h\geq1$, and let $a=h^{(B)}(r)$. If $r\ne0$, assume that $1\leq a<h$. Let $u$ be a positive odd integer with a representation of d-height at most $h-a$. Then

$$
\boxed{n=Pu+2^{h-1}\sigma^{(B)}(r)}
\tag{2}
$$

is a positive odd integer with a nonempty $d$-representation of d-height at most $h$. The conclusion also holds for $r=0$, in which case $a=0$ and $n=Pu$.

**Proof.** First suppose that $r\ne0$. In a representation of $u$ satisfying the stated conditions, every $5$-exponent is nonnegative and every $2$-exponent is at most $h-a-1$. Thus the terms of $Pu$ have $5$-exponents at least $11$ and $2$-exponents at most $h-a-1$.

On the other hand,

$$
2^{h-1}\sigma^{(B)}(r)
=\sum_{(\alpha,\gamma)\in\mathcal C^{(B)}(r)}
5^\alpha2^{h-1-\gamma}.
$$

Since $\alpha<11$ and $\gamma\leq a-1$, these terms have $5$-exponents less than $11$ and $2$-exponents at least $h-a$. The exponents within each of the two parts already satisfy the required strict order. After the parts are combined, the $5$-exponents are still strictly increasing and the $2$-exponents strictly decreasing, so (2) gives a valid $d$-representation. All its $2$-exponents are at most $h-1$, so its d-height is at most $h$.

Also, $a<h$ implies that all $2$-exponents in the second part are at least $1$. Since $Pu$ is a positive odd integer and $\sigma^{(B)}(r)>0$, the integer $n$ is positive and odd.

If $r=0$, then $n=Pu$. Multiplication by $P=5^{11}$ increases all $5$-exponents of $u$ by $11$ without changing the $2$-exponents. It therefore preserves both the validity and the d-height of the $d$-representation. This completes the proof.

## Lemma 3.

Fix $B\in\mathcal B$ and $h\geq1$. On the set of pairs $(r,u)$ allowed by Lemma 2, the map

$$
(r,u)\longmapsto Pu+2^{h-1}\sigma^{(B)}(r)
$$

is injective.

**Proof.** If two ordered pairs $(r_1,u_1)$ and $(r_2,u_2)$ have the same image, reduction modulo $P$ and the congruence $\sigma^{(B)}(r)\equiv r\pmod P$ give

$$
2^{h-1}r_1\equiv2^{h-1}r_2\pmod P.
$$

Since $(2,P)=1$ and $0\leq r_1,r_2<P$, we have $r_1=r_2$. The same residue class uses the same fixed value $\sigma^{(B)}(r)$. Substituting back into the original equality gives $Pu_1=Pu_2$, and hence $u_1=u_2$. This completes the proof.

## Lemma 4.

Fix $2\leq m\leq5$ and $1\leq h\leq K_m$. For each nonnegative integer $b$, let

$$
\eta=\frac{Q2^{h-1}}{R5^{\nu_m}},\qquad j_b=\lfloor\eta b\rfloor.
$$

For the parameters above, $\eta\leq1$. Indeed, since $Q/R=1/4$, it suffices to check that $2^{K_m-1}\leq4\cdot5^{\nu_m}$. For $m=2,3,4,5$, these inequalities are respectively

$$
2^{56}\leq4\cdot5^{24},\quad
2^{82}\leq4\cdot5^{35},\quad
2^{107}\leq4\cdot5^{46},\quad
2^{133}\leq4\cdot5^{57}.
$$

Suppose that $0\leq r<P$, $b=\lfloor R\sigma^{(B)}(r)\rfloor$, and $2\leq i\leq2Q$, $i-j_b-2\geq0$. Then

$$
\boxed{u\in I_{m-1,i-j_b-2}\ \Longrightarrow\
Pu+2^{h-1}\sigma^{(B)}(r)\in I_{m,i}.}
\tag{3}
$$

**Proof.** From $b=\lfloor R\sigma^{(B)}(r)\rfloor$ and $0<\eta\leq1$, we have

$$
\eta b\leq
\frac{Q2^{h-1}\sigma^{(B)}(r)}{5^{\nu_m}}
<\eta(b+1)\leq\eta b+1.
$$

Together with $j_b\leq\eta b<j_b+1$, this gives

$$
j_b\leq
\frac{Q2^{h-1}\sigma^{(B)}(r)}{5^{\nu_m}}
<j_b+2.
\tag{4}
$$

If $u\in I_{m-1,i-j_b-2}$, then

$$
i-j_b-2-W_{m-1}<\frac{Qu}{5^{\nu_{m-1}}}\leq i-j_b-2.
\tag{5}
$$

Let $n=Pu+2^{h-1}\sigma^{(B)}(r)$. Since $5^{\nu_m}=P5^{\nu_{m-1}}$, (4) and (5) give

$$
i-W_m=i-j_b-2-W_{m-1}+j_b<\frac{Qn}{5^{\nu_m}}<i-j_b-2+j_b+2=i.
$$

Thus $n\in I_{m,i}$. This completes the proof.

## Proposition 1.

**Definitions and conventions.** Lemma 1 gives $V_1(h,i)$ for the first step. For $2\leq m\leq5$, we compute $V_m(h,i)$ only for $0\leq h\leq K_m$ and $0\leq i\leq2Q$, and set

$$
V_m(0,i)=0,\qquad V_m(h,0)=V_m(h,1)=0.
$$

All subsequent expressions involving $V_{m-1}(t,j)$ from the previous step use the same convention: the value is zero when $j<0$, and when $0\leq j\leq2Q$ and $t>K_{m-1}$, we use $V_{m-1}(K_{m-1},j)$ instead. The recurrence uses only $t\geq0$, and $V_{m-1}(0,j)=0$ when $t=0$. Every nonnegative $j$ that occurs is at most $2Q$.

For each $B\in\mathcal B$, use the joint counts of the residue families defined above,

$$
M_B(a,b)=\#\{0\leq r<P:
 h^{(B)}(r)=a,\ \lfloor R\sigma^{(B)}(r)\rfloor=b\}.
$$

For fixed $B$, only one chain is chosen for each residue class. Thus $M_B(a,b)$ counts residue classes, rather than all possible chains. In particular, $M_B(0,0)=1$ and $\sum_{a,b}M_B(a,b)=P$.

Below we use only nonzero residue classes of heights $15\leq a\leq28$. This range is part of the chosen recurrence and does not mean that no residue classes occur at other heights. Omitting their contributions can only decrease the lower bound. Appendix C explains how the residue program constructs $M_B(a,b)$.

For $1\leq h\leq K_m$ and $2\leq i\leq2Q$, let

$$
\begin{aligned}
V_m(h,i)=\max\Biggl\{&V_m(h-1,i),\\
&\max_{B\in\mathcal B}\Biggl\lfloor\frac1P\Biggl(
V_{m-1}(h,i-2)\\
&\qquad+\sum_{a=15}^{\min(h,28)}\sum_{b\geq0}
M_B(a,b)V_{m-1}(h-a,i-j_b-2)
\Biggr)\Biggr\rfloor\Biggr\}.
\end{aligned}
\tag{6}
$$

Here $j_b$ is defined by Lemma 4 for the current $m,h$, and terms with $M_B(a,b)=0$ do not affect the sum. For fixed $B$, contributions from different residue classes can be added because Lemma 3 guarantees that the resulting integers are distinct. Constructions for different $B$ may overlap, so we take the maximum over $B$ rather than a sum. Taking the maximum with $V_m(h-1,i)$ in (6) uses the fact that the count cannot decrease when the permitted height increases.

Then, for $1\leq m\leq5$, $0\leq h\leq K_m$, and $0\leq i\leq2Q$, this recurrence satisfies

$$
T_h(I_{m,i})\geq\frac{5^{\nu_m}}{D_0}V_m(h,i).
\tag{7}
$$

**Proof.** Lemma 1 proves the case $m=1$. Now let $2\leq m\leq5$, and suppose that, for all $0\leq t\leq K_{m-1}$ and $0\leq j\leq2Q$,

$$
T_t(I_{m-1,j})\geq
\frac{5^{\nu_{m-1}}}{D_0}V_{m-1}(t,j).
\tag{8}
$$

Since $T_t(I)$ is nondecreasing in $t$, this still gives, for every $t\geq0$,

$$
T_t(I_{m-1,j})\geq
\frac{5^{\nu_{m-1}}}{D_0}
V_{m-1}(\min(t,K_{m-1}),j),
\qquad 0\leq j\leq2Q.
\tag{9}
$$

For fixed $m$, we now induct on $h$. When $h=0$, we have $V_m(0,i)=0$, so the claim holds. Let $1\leq h\leq K_m$, and suppose the claim holds for $h-1$ and all $0\leq i\leq2Q$. Since $T_h(I)\geq T_{h-1}(I)$, for all these $i$ we have

$$
T_h(I_{m,i})\geq
\frac{5^{\nu_m}}{D_0}V_m(h-1,i).
\tag{10}
$$

For $i=0,1$, (7) follows directly from the boundary convention. Now fix $2\leq i\leq2Q$ and $B\in\mathcal B$. We will show that the value corresponding to this $B$ in (6) also gives a valid lower bound.

For $r=0$, we have $\sigma^{(B)}(0)=0$. If a positive odd integer $u\in I_{m-1,i-2}$ has a representation of d-height at most $h$, then $n=Pu$ also has a valid representation with the same d-height bound. Moreover,

$$
i-W_m=i-2-W_{m-1}<\frac{Qu}{5^{\nu_{m-1}}}=\frac{QPu}{5^{\nu_m}}\leq i-2\leq i.
$$

Thus $n\in I_{m,i}$. The map $u\mapsto Pu$ is injective. Thus, by (9) and the boundary convention above, the number of new integers produced by the zero residue class is at least

$$
\frac{5^{\nu_{m-1}}}{D_0}V_{m-1}(h,i-2).
\tag{11}
$$

For the nonzero residue classes, fix $15\leq a\leq\min(h,28)$ and $b\geq0$. When $a=h$, the d-height bound for the previous integer is $0$, whereas a nonempty d-representation has d-height at least $1$. Thus no positive odd integer qualifies, and the corresponding term in (6) is also zero because $V_{m-1}(0,j)=0$. Henceforth suppose that $a<h$.

There are $M_B(a,b)$ nonzero residues with $h^{(B)}(r)=a$ and $\lfloor R\sigma^{(B)}(r)\rfloor=b$. If $i-j_b-2<0$, this term is zero by definition. Otherwise, for each such $r$, take all positive odd integers $u\in I_{m-1,i-j_b-2}$ of d-height at most $h-a$. By (9), there are at least

$$
\frac{5^{\nu_{m-1}}}{D_0}
V_{m-1}(h-a,i-j_b-2)
$$

such integers $u$. By Lemmas 2 and 4, each $n=Pu+2^{h-1}\sigma^{(B)}(r)$ is a positive odd integer in $I_{m,i}$ of d-height at most $h$. Lemma 3 shows that these pairs $(r,u)$ produce distinct integers. Thus the fixed pair $(a,b)$ contributes at least

$$
\frac{5^{\nu_{m-1}}}{D_0}
M_B(a,b)V_{m-1}(h-a,i-j_b-2)
$$

integers.

For fixed $B$, each nonzero residue class has a unique height $a$ and size bin $b$, so the contributions from different $(a,b)$ are also disjoint. The images of the zero residue class satisfy $n\equiv0\pmod P$, whereas those of nonzero residue classes satisfy $n\equiv2^{h-1}r\not\equiv0\pmod P$. Hence these two parts are disjoint as well. Combining (11) with the contributions of the nonzero residue classes gives

$$
T_h(I_{m,i})\geq\frac{5^{\nu_{m-1}}}{D_0}\left(V_{m-1}(h,i-2)+\sum_{a=15}^{\min(h,28)}\sum_{b\geq0}M_B(a,b)V_{m-1}(h-a,i-j_b-2)\right).
\tag{12}
$$

Since $5^{\nu_m}=P5^{\nu_{m-1}}$, the right-hand side of (12) equals

$$
\frac{5^{\nu_m}}{D_0P}
\left(V_{m-1}(h,i-2)+
\sum_{a=15}^{\min(h,28)}\sum_{b\geq0}
M_B(a,b)V_{m-1}(h-a,i-j_b-2)\right).
$$

The quantity in parentheses is a nonnegative integer. Dividing it by $P$ and taking the floor still gives a valid lower bound as required in (7). This argument holds for each $B\in\mathcal B$, so the maximum of these lower bounds is also valid. Taking the maximum with (10) shows that $V_m(h,i)$ as defined in (6) satisfies (7). This completes the induction on $h$, and then the induction on $m$.

## Lemma 5.

Fix $1\leq m\leq5$ and $0\leq h\leq K_m$. Set $U_m(h,i)=0$ for $i\leq0$, and for $1\leq i\leq2Q$ define

$$
\boxed{
U_m(h,i)=\max\{U_m(h,i-1),\ V_m(h,i)+U_m(h,i-W_m)\}.
}
\tag{13}
$$

Then, for all $0\leq i\leq2Q$,

$$
\boxed{
T_h\left(\frac{i5^{\nu_m}}Q\right)
\geq\frac{5^{\nu_m}}{D_0}U_m(h,i).
}
\tag{14}
$$

**Proof.** Take any two integers $1\leq j<k\leq2Q$. By the definitions of $I_{m,j}$ and $I_{m,k}$,

$$
I_{m,j}\cap I_{m,k}=\varnothing\iff\frac{j5^{\nu_m}}Q\leq\frac{(k-W_m)5^{\nu_m}}Q\iff k-j\geq W_m.
\tag{15}
$$

Here $j,k$ are arbitrary distinct interval indices, not parameters fixed in advance. The intervals remain disjoint when the relevant endpoints coincide, since $I_{m,k}$ does not contain its left endpoint.

For each integer $0\leq i\leq2Q$, consider all index sets $J$ satisfying

$$
J\subseteq\{1,\ldots,i\},\qquad
k-j\geq W_m\quad\text{for all }j,k\in J\text{ with }j<k
\tag{16}
$$

with the convention $\{1,\ldots,0\}=\varnothing$. We will prove that

$$
U_m(h,i)=\max_{J\text{ satisfying }(16)}\sum_{j\in J}V_m(h,j)
\qquad(0\leq i\leq2Q).
\tag{17}
$$

We induct on $i$. When $i=0$, the only set allowed by (16) is $J=\varnothing$, so both sides of (17) are zero. Let $i\geq1$, and suppose that (17) holds for all nonnegative integers smaller than $i$. For any $J$ satisfying (16), there are two cases: $i\notin J$ and $i\in J$.

If $i\notin J$, then $J\subseteq\{1,\ldots,i-1\}$. By the induction hypothesis, the maximum sum over such sets is $U_m(h,i-1)$.

If $i\in J$, then for each $j\in J\setminus\{i\}$, (16) gives $i-j\geq W_m$, or equivalently $j\leq i-W_m$. Conversely, adjoining $\{i\}$ to any subset of $\{1,\ldots,i-W_m\}$ satisfying (16) produces another set satisfying (16). When $i-W_m\leq0$, that subset must be empty, and $U_m(h,i-W_m)=0$. When $i-W_m>0$, the induction hypothesis gives the maximum sum over such subsets as $U_m(h,i-W_m)$. Therefore, in the case $i\in J$, the maximum sum is exactly

$$
V_m(h,i)+U_m(h,i-W_m).
$$

Taking the maximum of the two cases and using (13) proves (17).

Choose a set $J$ attaining the maximum in (17). By (15) and (16), the intervals $I_{m,j}$ for $j\in J$ are pairwise disjoint, and their right endpoints are all at most $i5^{\nu_m}/Q$. Proposition 1 gives

$$
T_h\left(\frac{i5^{\nu_m}}Q\right)\geq\sum_{j\in J}T_h(I_{m,j})\geq\frac{5^{\nu_m}}{D_0}\sum_{j\in J}V_m(h,j)=\frac{5^{\nu_m}}{D_0}U_m(h,i).
$$

This completes the proof.

## Proposition 2.

For $1\leq m\leq5$, $0\leq i\leq2Q$, and every positive integer $h$,

$$
\boxed{
T_h\left(\frac{i5^{\nu_m}}Q\right)
\geq\frac{5^{\nu_m}}{D_0}
U_m(\min\{h,K_m\},i).
}
\tag{18}
$$

**Proof.** If $h\leq K_m$, the result follows from Lemma 5. If $h>K_m$, then $T_h(X)$ is nondecreasing in $h$, so

$$
T_h\left(\frac{i5^{\nu_m}}Q\right)
\geq T_{K_m}\left(\frac{i5^{\nu_m}}Q\right)
\geq\frac{5^{\nu_m}}{D_0}U_m(K_m,i).
$$

This completes the proof.

Proposition 2 is the result that the terminal program supplies to the subsequent density calculation. It gives a lower bound on the count. It does not require the program to find all representations in the interval, nor does it require the right-hand side to be an integer. We next prove the additional recurrences used to compute these bounds and explain how the implementation corresponds to the results above.

## Appendix A. The recurrence for the initial height table

Lemma 1 uses $h_{\min,14}(n)$. This appendix proves how to compute this function over a finite range and then obtain $V_1(h,i)$ from it.

For $k\geq0$ and each positive integer $n$, let $h_{\min,k}(n)$ denote the minimum d-height among all nonempty $d$-representations whose $5$-exponents are less than $k$. Set this value to $+\infty$ if no such representation exists. In particular, $h_{\min,0}(n)=+\infty$ for every positive integer $n$. We use the convention $e+(+\infty)=+\infty$. There is no need here to define a height for a nonempty representation of $0$.

### Lemma A.1

If $n=2^e u$, where $e\geq0$ and $u$ is a positive odd integer, then

$$
h_{\min,k}(n)=e+h_{\min,k}(u).
\tag{A.1}
$$

**Proof.** In a nonempty $d$-representation, the smallest $2$-exponent occurs exactly once. Denote it by $b$. After dividing every term by $2^b$, exactly one term is odd and all the others are even, so the resulting sum is odd. Hence $b=v_2(n)=e$.

Dividing the entire representation of $n$ by $2^e$ gives a representation of $u$. All $5$-exponents remain unchanged, and all $2$-exponents decrease by $e$, so the d-height also decreases by $e$. Conversely, multiplying any representation of $u$ by $2^e$ gives a representation of $n$ with d-height increased by $e$. These operations are inverse to each other and preserve the restriction that the $5$-exponents are less than $k$. Taking minima proves (A.1). If one of the two integers has no such representation, neither does the other, so the conclusion still holds. This completes the proof.

### Lemma A.2

Let $k\geq1$, and let $n$ be a positive odd integer. Then

$$
h_{\min,k}(n)=
\begin{cases}
h_{\min,k-1}(n),&n<5^{k-1},\\
1,&n=5^{k-1},\\
\displaystyle\min\left\{h_{\min,k-1}(n),\,
1+h_{\min,k-1}\left(\frac{n-5^{k-1}}2\right)\right\},&n>5^{k-1}.
\end{cases}
\tag{A.2}
$$

**Proof.** Divide the valid representations of $n$ into two classes, according to whether they omit or use the $5$-exponent $k-1$.

The first class consists of representations with all $5$-exponents less than $k-1$, and its minimum height is $h_{\min,k-1}(n)$.

In the second class, the term with $5$-exponent $k-1$ has the largest $5$-exponent, and hence also the smallest $2$-exponent. Since $n$ is odd, the proof of Lemma A.1 shows that this smallest exponent is $0$, so this term must be $5^{k-1}$. All terms are positive, so the second class is empty when $n<5^{k-1}$. If $n=5^{k-1}$, the singleton representation has height $1$, which is the smallest possible height of a nonempty representation.

Now suppose that $n>5^{k-1}$. Delete the final term $5^{k-1}$ and divide the remaining terms by $2$. This gives a nonempty $d$-representation of the positive integer $(n-5^{k-1})/2$, with all $5$-exponents less than $k-1$ and d-height $1$ less than that of the original representation.

Conversely, multiply any representation of this positive integer satisfying the exponent restriction by $2$, and then adjoin $5^{k-1}$. This gives a valid representation of $n$. The new term has the largest $5$-exponent and the unique zero $2$-exponent, so it preserves the antichain condition. Thus the minimum height in the second class is exactly the second candidate in (A.2). Taking the minimum over the two classes proves the result.

Equations (A.1) and (A.2) give a complete finite order of computation. Start with height $+\infty$ for every positive integer when $k=0$, and then compute successively for $k=1,\ldots,14$. Only heights of odd integers need to be stored, since the positive integer appearing in (A.2) can be reduced to an odd integer using (A.1). That integer is strictly less than $n$, and its odd part is also strictly less than $n$. Thus, when computing for $n$ up to a given bound, no values beyond that bound are needed.

The function `odd_heights` in `residues.py` evaluates precisely these two recurrences. It uses two arrays. At the start of step $k$, `previous[t]` represents $h_{\min,k-1}(2t+1)$, while the entry written to `current[t]` represents $h_{\min,k}(2t+1)$. Every height needed by the recurrence is read from `previous`, and the singleton case is assigned $1$ directly. The arrays are exchanged only after all entries have been computed, so the current step cannot inadvertently use a height already updated to $k$. Width $1$ is initialized directly, since its only represented positive odd integer is $1$. At later widths the integer factorization in Lemma A.2 is evaluated by two exact 16-bit lookup operations. Only the prefix not overwritten by the recurrence is copied. The implementation details and fixed-width bounds are given in Section 4 of `02_residues.md`.

### Lemma A.3

Let $N(a,j)$ count the integers satisfying the following conditions:

$$
1\leq n<2\cdot5^{13},\qquad 2\nmid n,\qquad
h_{\min,14}(n)=a,\qquad c(n)=j
$$

For $0\leq h\leq K_1$ and $0\leq i\leq2Q$, set

$$
C(h,i)=\sum_{a=0}^{h}\sum_{j=0}^{i}N(a,j).
$$

Then the initial values in Lemma 1 can be written as

$$
V_1(h,i)=\left\lfloor\frac{D_0}{5^{13}}
\left(C(h,i)-C\bigl(h,\max\{i-W_1,0\}\bigr)\right)\right\rfloor.
\tag{A.3}
$$

**Proof.** Since $n>0$, we have $c(n)\geq1$. Thus $N(a,0)=0$, and consequently $C(h,0)=0$.

If $i-W_1\geq0$, subtracting the two cumulative counts retains exactly those integers of height at most $h$ with $i-W_1<c(n)\leq i$. If $i-W_1<0$, every integer with $c(n)\geq1$ automatically satisfies the left endpoint condition, while the subtracted term $C(h,0)$ is zero. Thus both cases give the count in Lemma 1. Substitution into its definition proves (A.3).

The function `seed_intervals` first uses the minimum height table to compute $N(a,j)$, then takes cumulative sums along the height and interval index axes to obtain $C(h,i)$. It then performs the subtraction, multiplication, and integer division in (A.3). Summing first over heights and then over interval indices merely specifies the order of a finite double sum and introduces no further mathematical assumption.

## Appendix B. Grouping terms with the same shift

The recurrence in Proposition 1 sums over size bins $b$. In the computation, many different $b$ give the same shift $j_b$. Grouping these terms reduces the amount of computation. We now prove the grouping formula used by the program.

Fix $m,h$ and $B$, and retain $\eta>0$ from Lemma 4. Since $0\leq\sigma^{(B)}(r)\leq B$, only $0\leq b\leq BR$ can occur. Let

$$
F_B(a,t)=\sum_{0\leq b<t}M_B(a,b),\qquad 0\leq t\leq BR+1.
$$

For each nonnegative integer $j$, let

$$
e_j=\min\left\{BR+1,\left\lceil\frac j\eta\right\rceil\right\},
\qquad
\omega_B(a,j)=F_B(a,e_{j+1})-F_B(a,e_j).
$$

### Lemma B.1

For every $a,j\geq0$, we have

$$
\omega_B(a,j)=\sum_{\substack{0\leq b\leq BR\\\lfloor\eta b\rfloor=j}}M_B(a,b).
\tag{B.1}
$$

Consequently, for $2\leq i\leq2Q$, the sum in Proposition 1 satisfies

$$
\begin{aligned}
&\sum_{a=15}^{\min(h,28)}\sum_{b\geq0}
M_B(a,b)V_{m-1}(h-a,i-j_b-2)\\
&\qquad=
\sum_{a=15}^{\min(h,28)}\sum_{j=0}^{i-2}
\omega_B(a,j)V_{m-1}(h-a,i-j-2).
\end{aligned}
\tag{B.2}
$$

Here heights exceeding the stored range at the previous step are still handled by the convention in Proposition 1.

**Proof.** For an integer $b$,

$$
\lfloor\eta b\rfloor=j
\iff j\leq\eta b<j+1
\iff
\left\lceil\frac j\eta\right\rceil\leq b<
\left\lceil\frac{j+1}\eta\right\rceil.
$$

Intersecting this integer interval with $0\leq b<BR+1$ gives $e_j\leq b<e_{j+1}$. Summing the corresponding $M_B(a,b)$ is exactly the difference of the two prefix sums, which proves (B.1).

For a fixed $j$, the factor $V_{m-1}(h-a,i-j-2)$ is independent of $b$, so these terms can be grouped. If $j>i-2$, the previous interval index is negative and the term is zero by convention. Thus it suffices to sum over $0\leq j\leq i-2$, giving (B.2). This completes the proof.

The program does not compute $\eta$ or $e_j$ in floating point. It uses directly

$$
e_j=\min\left\{BR+1,
\left\lceil\frac{jR5^{\nu_m}}{Q2^{h-1}}\right\rceil\right\}.
$$

The code sets `numerator = Q * 2**(h-1)` and `divisor = R * 5**exponent`, where `exponent` equals $\nu_m$. Thus `ceil_div(j * divisor, numerator)` gives exactly the integer ceiling required above.

### Correspondence with shifted array addition

In `next_intervals`, `prefix[B][a,t]` represents $F_B(a,t)$, `edges[j]` represents $e_j$, and `weights[j]` represents $\omega_B(a,j)$.

For fixed $h,a$, the array `source` stores

$$
\texttt{source}[t]=V_{m-1}(\min\{h-a,K_{m-1}\},t),
\qquad 0\leq t\leq2Q-2.
$$

Index $t$ in the temporary array `total` corresponds to interval index $i=t+2$ in the new table. Thus

```python
total[j:] += weights[j] * source[:len(total) - j]
```

at each index $t\geq j$, adds exactly

$$
\omega_B(a,j)V_{m-1}(h-a,t-j)
=\omega_B(a,j)V_{m-1}(h-a,i-j-2).
$$

This agrees term by term with (B.2). The implementation also splits these array slices into segments of length `chunk` for addition. Each segment uses the same index correspondence, so the sum is unchanged. The initial value of `total` includes the contribution $V_{m-1}(h,i-2)$ from the zero residue class. After summing, the program divides by $P$ and takes the floor, then takes the maximum over $B$ and the previous height. This gives exactly (6). The code also performs integer type conversions. Appendix C gives the bounds that justify them.

## Appendix C. Correspondence between the results and the program

The following discussion refers to `residues.py`, `terminals.py`, and `parameters.py` in the root directory of this package. The function `odd_heights` uses the two arrays described above. The mathematical interpretation in this appendix explains the code, so the validity of the recurrences need not be inferred anew from the array operations.

### C.1 Construction of the joint residue counts

The main text assumes that a valid chain has been chosen for each $B$ and each residue class $r$, and the terminal program uses the joint counts $M_B(a,b)$. The function `positive_families` makes these choices using representations of positive odd integers.

Suppose that a positive odd integer $S$ has a representation with $5$-exponents less than $11$ and d-height at most $a$,

$$
S=\sum_j5^{\alpha_j}2^{\beta_j}.
$$

Set $\gamma_j=a-1-\beta_j$. Since $\alpha_j$ is strictly increasing and $\beta_j$ strictly decreasing, $\gamma_j$ is strictly increasing and nonnegative. Also, since $S$ is odd, the smallest $\beta_j$ is $0$, so the largest $\gamma_j$ is exactly $a-1$. The resulting residue chain therefore has height exactly $a$, and

$$
\sigma=\sum_j5^{\alpha_j}2^{-\gamma_j}
=\frac{S}{2^{a-1}}.
$$

If $S\leq B2^{a-1}$, then $\sigma\leq B$. Its residue class and size bin are respectively

$$
r\equiv S(2^{a-1})^{-1}\pmod P,
\qquad
b=\left\lfloor\frac{RS}{2^{a-1}}\right\rfloor.
$$

Thus the program only needs to know that $h_{\min,11}(S)\leq a$ to ensure the existence of a suitable chain with the stated value of $\sigma$ and height $a$. Here $a$ need not equal the minimum d-height of $S$. These quantities must be distinguished.

For fixed $B$, the program scans in increasing order of $a$, and then in increasing order of $S$. Each residue is recorded only at its first occurrence, with the empty chain selected in advance for the zero residue. All $P$ residues must have been covered by the end of the scan before verification can proceed. Failure to cover every residue within the finite search range means only that this construction is incomplete.

Thus `joint[B][a,b]` is precisely $M_B(a,b)$ in the main text. The main program need not store every chain: Appendix A and the transformation above prove that each record corresponds to a valid chain, and the subsequent calculation depends only on height and size bin. Nor does the finite scan bound need to ensure that every possible chain is found. It need only ensure that a valid choice has been made for every residue.

The implementation in chunks also ensures that distinct $S$ within a chunk give distinct residues. It limits the number of odd integers in each chunk so that `2 * chunk < P`. The absolute difference between two distinct odd integers in a chunk is strictly less than $P$, so they cannot be congruent modulo $P$. The same `seen` array is used across chunks and heights to exclude residues already recorded.

### C.2 Results, functions, and return values

| Mathematical object or result | Function | Return value and use |
| --- | --- | --- |
| Minimum height recurrence in Lemmas A.1 and A.2 | `residues.odd_heights` | A one-dimensional array whose index $t$ corresponds to the odd integer $2t+1$ |
| The residue chain selection above and $M_B(a,b)$ | `residues.positive_families` | `hist, joint`, with `joint[B][a,b]=M_B(a,b)` |
| $V_1(h,i)$ in Lemmas 1 and A.3 | `terminals.seed_intervals` | A two-dimensional integer array `table[h,i]` |
| $V_m(h,i)$ in Lemmas 2–4 and Proposition 1 | `terminals.next_intervals` | A new table formed from the previous table `previous` and `joint` |
| $U_m(h,i)$ in Lemma 5 | `terminals.pack_intervals` | An integer table obtained by optimally selecting pairwise disjoint intervals at each height |
| The lower bound for $T_h(X)$ in Proposition 2 | `terminals.positive_terminals` | A dictionary `result[h,X]` with exact rational values |
| Selection of the terminal algorithm by bases | `terminals.terminal_counts` | Calls `positive_terminals` for $(5,2)$ |

The functions only need to pass these mathematical quantities between them. No classes or further abstractions are needed. The function `odd_heights` is used both for the initial interval calculation and for the residue family construction. The former computes up to $2\cdot5^{13}$ with $5$-exponents less than $14$. The latter uses the finite range specified by `positive_families`, with $5$-exponents less than $11$. These are two applications of the same proved recurrence.

The entry `hist[B][a]` equals $\sum_bM_B(a,b)$ and is used in the subsequent exponential moment and tail bound calculations. The terminal recurrence also needs size bin information, so it uses `joint`. This information must not be discarded by reducing `joint` to `hist` before computing terminal lower bounds.

### C.3 Fixed parameters and array indices

| Notation in the text | Name in the code | Fixed value or meaning |
| --- | --- | --- |
| $P$ | `P` | $5^{11}$ |
| $\mathcal B$ | Keys of `joint` | $\{10,11,12,14,22\}$ |
| $Q$ | `grid` | $8192$ |
| $R$ | `size_grid` | $32768$ |
| $D_0$ | `denominator` | $10^9$ |
| $\nu_1$ | `seed_exponent` | $13$ |
| Initial upper bound on the $5$-exponents | `seed_width` | $14$ |
| $(K_1,\ldots,K_5)$ | `height_caps` | $(40,57,83,108,134)$ |
| $W_1$ | `span` in `seed_intervals` | $400$ |
| $W_m$ | `span` in `pack_intervals` | $400+2(m-1)$ |
| $m$ | `level` | $1,\ldots,5$ |
| $\nu_m$ | `exponent` | $13+11(m-1)$ |

These quantities are parameters of the chosen mathematical construction. The parameter `chunk` only determines how many array entries are processed at a time and does not change the recurrence or the final result. The fixed values are collected in `POSITIVE_SETTINGS` in `parameters.py`, and `positive_terminals` passes them explicitly to the computational functions. The retained height range is specified by `retained_heights=(15,28)`. The increments in the exponent width and interval span are specified by `width_increment=11` and `span_increment=2`, respectively. Any parameter change requires the corresponding assumptions to be checked again. The numerical conclusions for the original parameters cannot simply be reused after changing an individual value.

The following distinctions between array entries are useful.

| Array or entry | Mathematical meaning |
| --- | --- |
| `previous[t]` in `odd_heights` | The minimum height of the odd integer $2t+1$ at the preceding value of $k$ |
| `raw[a,j]` in `seed_intervals`, before both cumulative sums | $N(a,j)$ in Appendix A |
| `raw[h,i]` after both cumulative sums | $C(h,i)$ in Appendix A |
| `previous[h,i]` in `next_intervals` | $V_{m-1}(h,i)$ |
| `table[h,i]` in `next_intervals` | $V_m(h,i)$ |
| `total[t]` in `next_intervals` | The numerator in (6) for target index $i=t+2$ |
| `U[h,i]` in `pack_intervals` | $U_m(h,i)$ |

### C.4 Order of computation and boundary conventions

The residue program first constructs `joint`. The terminal program for positive bases then performs the following steps.

1. Compute the initial minimum height table using Lemmas A.1 and A.2, then obtain $V_1$ using Lemma A.3.
2. For $m=2,3,4,5$, compute $V_m$ from $V_{m-1}$ using Proposition 1. Within each step, compute in increasing order of $h$ so that $V_m(h-1,i)$ is available.
3. At each step where terminal lower bounds are needed, apply Lemma 5 to $V_m$ to obtain $U_m$.
4. Read the requested $(h,X)$ according to Proposition 2 to obtain the lower bounds needed for the subsequent density calculation.

Step 2 always passes on $V_m$. It cannot be replaced by $U_m$ from Step 3. The former counts integers in one short interval, whereas the latter counts integers in a union of pairwise disjoint intervals. Their interval meanings are different.

The following boundary conventions all come from the main text. They are not additional rules introduced for convenient array access.

| Case | Program behavior | Mathematical justification |
| --- | --- | --- |
| $h=0$ | The corresponding row is zero | A nonempty representation cannot have d-height zero |
| Negative previous interval index | Omit the contribution | Omitting part of the construction still gives a lower bound |
| Previous height exceeds $K_{m-1}$ | Use the last row of the previous table | $T_h$ is nondecreasing in $h$ |
| $a=h$ | Read row zero of the previous table, giving zero contribution | The remaining integer requires a nonempty representation |
| $i=0,1$ when $m\geq2$ | The first two columns of the new table are zero | Proposition 1 uses only the trivial lower bound for these columns |
| `pack_intervals` starts at $i=0$ | Gives $U_m(h,0)=0$ | $V_m(h,0)=0$, consistent with the initial value in Lemma 5 |

In particular, the initial entry $V_1(h,1)$ is not set to zero. It is determined by the count in Lemma 1. The convention that the first two columns are zero applies only to the recurrence tables for $m\geq2$.

When reading the results, the scale record specifies $\nu_m$ and

$$
X=\frac{i5^{\nu_m}}Q,\qquad 0\leq i\leq2Q.
$$

The function `positive_terminals` computes $i=QX/5^{\nu_m}$ using exact rational arithmetic and requires it to be an integer in the stated range. It then returns

$$
\texttt{result}[h,X]
=\frac{5^{\nu_m}}{D_0}U_m(\min\{h,K_m\},i).
$$

No further floor is taken here. The integer floors have already been taken in the definition of $V_1$ and in (6), and the return value is retained as an exact rational number using `Fraction`. Proposition 2 already guarantees that this rational number does not exceed the actual count.

### C.5 Bounds for integer arithmetic

The proofs above use mathematical integers. The implementation must also ensure that the intermediate results fit within the fixed-width NumPy integer types. For the fixed parameters in this document, the counting bounds verify this directly.

The initial range contains $5^{13}$ positive odd integers, so each cumulative count is at most $5^{13}$. Hence

$$
5^{13}D_0=1220703125000000000<2^{63}.
$$

The intermediate product of the initial count and $D_0$ fits in `int64`, and $0\leq V_1(h,i)\leq D_0$.

If every value at the previous step is at most $D_0$, then for fixed $B$ the numerator in (6) is at most

$$
\left(1+\sum_{a=15}^{28}\sum_bM_B(a,b)\right)D_0\leq PD_0.
$$

Here $1$ accounts for the zero residue class, while the other terms include only some of the nonzero residue classes. After dividing by $P$, taking the floor, and taking maxima, we still have $0\leq V_m(h,i)\leq D_0$. By induction, this bound holds at every step. Also,

$$
D_0<2^{32},\qquad
PD_0=48828125000000000<2^{64}.
$$

Thus $V_m$ can be stored in `uint32`, and the numerator can be accumulated in `uint64`. Prefix counts are at most $P<2^{32}$, which also fits in `uint32`.

Packing selects at most $2Q$ indices, each with weight at most $D_0$, so

$$
U_m(h,i)\leq2QD_0=16384000000000<2^{64}.
$$

The type `uint64` used by `pack_intervals` is sufficient for these values. Before the final multiplication by $5^{\nu_m}$, the code converts the value to a Python integer using `int(...)` and then constructs a `Fraction`. Thus this large integer multiplication is not performed in a fixed-width NumPy integer type.

The function `odd_heights` uses `255` to represent $+\infty$. The implementation requires $n<2^{32}$. Since all terms are positive, every $2^\beta$ in a representation is at most $n$, so every representation that exists has height at most $32$. Thus `255` cannot be confused with a finite height in the current computational range. The initial bound $2\cdot5^{13}=2441406250<2^{32}$ used here lies within this range.

### C.6 Scope of verification and program interfaces

The necessary checks use `if ...: raise ...`, so they remain active under `python -O`. The function `seed_intervals` checks the integer range of the initial counts. The function `next_intervals` checks the dimensions of the previous array, its zero-height row and zero-endpoint column, the total count, the integer range, and the condition $\eta\leq1$ required for interval transfer. The function `pack_intervals` checks that the span is positive, that the zero-endpoint column is zero, and that the packing sum fits in the integer type. The function `positive_terminals` checks the exponents, terminal heights, and endpoint grid of all requests before the large initial computation.

The lower-level functions retain parameterized interfaces so that the finite formulas can be checked on small arrays. Artificial joint counts can be used to check the implementation of (B.2), but this does not make them valid residue families. For the full computation, `residues.positive_families` constructs the joint counts from scratch and checks coverage of all residues before the terminal program is called.

The scale records supplied by `read_inputs` use Python integers for `H` and `tail_terms`. Direct calls to `positive_terminals` must preserve these types, rather than replace them with equal-valued fractions or floating-point numbers.

The script `check_small.py` checks minimum heights and initial counts by direct enumeration, checks Lemma 5 by enumerating all subsets of intervals, and evaluates the double sum in (6) directly to compare it with the grouped shifted additions. These checks address implementation errors. The proofs in the main text and Appendices A and B establish why the recurrences give valid lower bounds at the stated scales.

The actual commands, parameters, coverage, and pass results for this package are recorded in `validation/`. The full computations establish the finite inequalities required by the density criterion; the small-scale checks test the implementation on independently enumerable examples.
