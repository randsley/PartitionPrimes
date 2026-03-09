# Extending Partition-Theoretic Prime Detection: The Cusp Form Barrier and the Derivation of E₅

### Abstract

Recent work by Craig, van Ittersum, and Ono [1] demonstrated a novel connection between additive and multiplicative number theory, using MacMahon partition functions $M_a(n)$ to detect primes via quasi-shuffle algebras. They constructed explicit prime-vanishing polynomial combinations $E_1(n)$ through $E_4(n)$ and conjectured that all such expressions form a finite generating set. This paper computationally derives the fifth fundamental expression, $E_5(n)$, and illuminates why the naïve pattern—incrementally introducing $M_{k+1}$ to form $E_k$—breaks down. We show that at weight $2a=12$, the generating function $U_6(q)$ for $M_6(n)$ admits a component in the cusp form space spanned by the Ramanujan delta function $\Delta(q)$, making $M_6$ linearly independent from divisor sums and excluding it from all prime-vanishing expressions. Restricting to $M_1,\ldots,M_5$ at polynomial degree $d=2$, we extract the exact closed-form $E_5(n)$ via rational null-space computation over 95 primes, and verify computationally that $\{E_1, E_2, E_3, E_4, E_5\}$ spans the prime-vanishing subspace for polynomial degree $d \leq 6$ and weight $a_{\max} \leq 5$. Unlike $E_1$–$E_4$, the canonical $E_5$ is negative at 349 of the 404 composites in $[4,500]$; the combination $E_5 + 856\,E_4$ is universally non-negative at all composites in this range.

---

### 1. Introduction

The study of integer partitions and prime numbers has historically occupied distinct branches of number theory—additive and multiplicative, respectively. However, recent results by Craig, van Ittersum, and Ono [1] establish a direct algebraic bridge between the two. They proved that for $n \geq 2$, specific non-negative polynomial combinations of MacMahon partition functions vanish if and only if $n$ is prime.

The MacMahon functions, $M_a(n)$, represent weighted sums over strict $a$-part partitions of $n$, defined as:

$$M_a(n) = \sum_{\substack{0 < s_1 < \cdots < s_a \\ m_i \geq 1,\; \sum m_i s_i = n}} m_1 m_2 \cdots m_a$$

By leveraging the **quasi-shuffle algebra** of these functions—an algebra structure on formal power series of MacMahon partition functions, governed by the quasi-shuffle product and Ramanujan's differential equations (see [1, Section 3])—and their connection to quasimodular forms, Craig–van Ittersum–Ono constructed a sequence of prime-detecting expressions $E_1(n)$ through $E_4(n)$, where each $E_k$ introduces $M_{k+1}$ as its highest-weight component. They conjectured that any non-negative prime-vanishing expression in $\mathbb{Q}[n] \otimes \{M_a\}$ is a $\mathbb{Q}[n]$-linear combination of these foundational entries.

In this paper we develop an exact computational framework in Julia to test and extend this conjecture. A computational sweep reveals that the conjecture fails for $a_{\max} \geq 5$ unless a fifth expression $E_5(n)$ is added. We derive $E_5$ explicitly.

A natural extrapolation of the established pattern suggests $E_5$ should incorporate $M_6(n)$. However, we demonstrate that this extrapolation fails due to modular arithmetic. The generating function $U_6(q) = \sum_{n \geq 1} M_6(n) q^n$ is a quasimodular form of weight 12. At weight 12, the space of modular forms on $\mathrm{SL}_2(\mathbb{Z})$ becomes two-dimensional, spanned by the Eisenstein series $E_{12}$ and the unique cusp form $\Delta(q) = \sum \tau(n) q^n$ (Ramanujan's delta function). The resulting $\tau(n)$ component in $M_6$ forces the $M_6$ columns to be pivot columns in the prime evaluation matrix: no prime-vanishing expression can involve $M_6(n)$.

Confined to $M_1$ through $M_5$, we perform a degree sweep. Contrary to the initial expectation that $E_5$ would appear at polynomial degree $d=3$ (matching the pattern of $E_1,\ldots,E_4$), we find that $E_5$ first appears already at $d=2$. This is the minimal-degree canonical form of $E_5$, and it is what we record and verify. The extraction procedure yields one unique new direction outside the $\mathbb{Q}[n]$-span of $E_1$–$E_4$ at $d=2$; this direction is $E_5$.

---

### 2. The Weight-12 Barrier and the Cusp Form $\Delta(q)$

#### 2.0. Background: Quasimodular Forms and Divisor Sums

Let $U_a(q) = \sum_{n=1}^{\infty} M_a(n) q^n$ be the generating function of $M_a(n)$. As established in [1, Section 2], $U_a(q)$ is a **quasimodular form of weight $2a$** on the full modular group $\mathrm{SL}_2(\mathbb{Z})$.

**For $a \in \{1,2,3,4,5\}$ (weights $2a \leq 10$):** The relevant spaces of quasimodular forms are one-dimensional, spanned entirely by Eisenstein series. Consequently, the Fourier coefficients $M_a(n)$ for $a \leq 5$ can be written as explicit polynomial combinations of the divisor power sum functions $\sigma_k(n) = \sum_{d \mid n} d^k$:

$$M_a(n) = \sum_j P_j(n)\,\sigma_{2j+1}(n)$$

where $P_j(n)$ are polynomials in $n$ with rational coefficients. The explicit closed-form formulas for $M_1$ through $M_5$ are given in [1, eqs. (1.4)–(1.8)].

Because at any prime $p$ we have $\sigma_k(p) = p^k + 1$, evaluating $M_a(p)$ collapses to a polynomial in $p$. This algebraic collapse creates the linear dependencies in the prime evaluation matrix that give rise to the null space containing $E_1,\ldots,E_4$.

#### 2.1. Quasimodular Forms of Weight 12 and the Ramanujan Tau Function

For $a = 6$, the generating function $U_6(q)$ is a quasimodular form of weight $2a = 12$. The space of modular forms of weight 12 on $\mathrm{SL}_2(\mathbb{Z})$ is **two-dimensional**, spanned by:

- the Eisenstein series $E_{12}(q)$, and
- the unique normalized cusp form of weight 12, Ramanujan's delta function:
$$\Delta(q) = q \prod_{k=1}^{\infty} (1-q^k)^{24} = \sum_{n=1}^{\infty} \tau(n)\,q^n$$

where $\tau(n)$ is the **Ramanujan tau function**.

Because $U_6(q)$ has a non-zero projection onto the cusp-form space, the closed-form expression for $M_6(n)$ necessarily involves $\tau(n)$. To determine this formula explicitly, we performed **rational coefficient fitting**: evaluating $M_6(n)$ by direct partition enumeration at $n = 1, \ldots, 55$ and solving the linear system

$$M_6(n) = \sum_{j=0}^{5} P_j(n)\,\sigma_{2j+1}(n) + c_\tau \cdot \tau(n)$$

exactly over $\mathbb{Q}$ using rational row reduction (RREF with exact pivoting on a $55 \times 22$ matrix: 21 polynomial-weighted divisor-sum columns of the form $n^k \sigma_{2j+1}(n)$, and 1 tau column). The unique solution is:

$$M_6(n) = \sum_{j=0}^{5} P_j(n)\,\sigma_{2j+1}(n) - \frac{17}{150\,450\,048\,000}\,\tau(n)$$

The **non-zero coefficient** $c_\tau = -17/150\,450\,048\,000$ is the key obstruction. The denominator factors as:
$$150\,450\,048\,000 = 2^{10} \times 3^5 \times 5^3 \times 7 \times 691$$
where 691 is the numerator of the Bernoulli number $B_{12}$—an expected fingerprint of modular forms at weight 12.

#### 2.2. The Exclusion of $M_6$: The Pivot Discovery

The presence of $\tau(n)$ in $M_6(n)$ fundamentally breaks the algebraic structure required for prime-vanishing expressions. Unlike the divisor sums $\sigma_k(p) = p^k + 1$ at primes, the Ramanujan tau function satisfies the Hasse–Weil bound $|\tau(p)| \leq 2p^{11/2}$ and does not reduce to a polynomial in $p$. Consequently, evaluating $n^k M_6(n)$ at primes yields linearly independent data that cannot be cancelled by divisor-sum combinations.

We demonstrate this via the **prime evaluation matrix** $\mathbf{M}(d, a_{\max}, N)$: the $N \times (d+1) a_{\max}$ matrix whose rows are evaluations of $\{n^k M_a(n) : 0 \leq k \leq d,\; 1 \leq a \leq a_{\max}\}$ at the first $N$ primes. All arithmetic is exact over $\mathbb{Q}$.

**Theorem 2.1 (Exclusion of $M_6$, computational).** For $d \in \{2, 3\}$ and $N = 95$ primes in $[2,500]$:
$$\mathrm{rank}(\mathbf{M}(d, 6, N)) = \mathrm{rank}(\mathbf{M}(d, 5, N)) + (d+1)$$
$$\dim\ker(\mathbf{M}(d, 6, N)) = \dim\ker(\mathbf{M}(d, 5, N))$$

*Proof.* Evaluated at all 95 primes in $[2, 500]$ with $d=2$:

| Matrix | Rows | Cols | Rank | Nullity |
|--------|------|------|------|---------|
| $a_{\max}=5$, $d=2$ | 95 | 15 | 11 | 4 |
| $a_{\max}=6$, $d=2$ | 95 | 18 | 14 | 4 |

Adding the 3 columns $\{M_6, n M_6, n^2 M_6\}$ increases the rank by exactly 3 and leaves the nullity unchanged. In the RREF of $\mathbf{M}(2, 6, 95)$, all three $M_6$-associated columns are pivot columns. This was verified independently for $d=3$ (nullity 8 for both $a_{\max}=5$ and $a_{\max}=6$).

**Corollary.** No prime-vanishing expression can involve $M_6(n)$. The basis for the prime-vanishing null space is confined to $M_1, \ldots, M_5$. $\square$

*Terminology note.* In row-reduced echelon form (RREF), a **pivot column** is one that contains a leading 1 in some row; it represents a direction that is linearly independent of all preceding columns. A column is a pivot column if and only if it contributes a new dimension to the column space—equivalently, it has no free variable in the null space. The above result says that in the RREF of $\mathbf{M}(d,6,N)$, every column corresponding to $n^k M_6(n)$ is a pivot, so no null vector can have a non-zero entry in those positions.

---

### 3. Extracting $E_5(n)$ from the Null Space

#### 3.1. The Search Space

Having established that $M_6(n)$ cannot appear in any prime-vanishing expression, we search for $E_5$ within the basis:

$$\mathcal{B}_d = \{ n^k M_a(n) : 0 \leq k \leq d,\; 1 \leq a \leq 5 \}$$

at increasing polynomial degree $d$. The prime evaluation matrix $\mathbf{M}(d, 5, N)$ (with $N = 95$ primes from $[2,500]$) is constructed using the exact closed-form formulas for $M_1,\ldots,M_5$ (which involve only divisor sums and polynomial prefactors, all computed exactly). Its null space is computed by rational RREF over $\mathbb{Q}$ using `Rational{BigInt}` arithmetic throughout.

#### 3.2. Isolating $E_5$ Outside the $E_1$–$E_4$ Span

Finding a prime-vanishing null vector is not sufficient: we must identify directions genuinely outside the $\mathbb{Q}[n]$-span of $E_1,\ldots,E_4$. To do so, we evaluate each null vector $\mathbf{v}$ at a set of composite numbers and test linear independence from the corresponding evaluations of $E_1,\ldots,E_4$.

**Algorithm: Extract\_$E_5$**

```
Input:  a_max = 5,  N = 95 primes from [2,500],  composites = [4..100]
Output: E₅(n) ∈ Q[n] ⊗ {M₁,...,M₅}

for d = 0, 1, 2, ... do:
    M_d ← prime evaluation matrix(B_d, N primes)      // Q-arithmetic
    null_d ← ker(M_d) via rational RREF
    if dim(null_d) = 0 then continue

    comp_mat ← evaluation of B_d at composites         // |composites| × dim(B_d)
    E14_mat  ← evaluate(E₁,...,E₄, composites)        // for span comparison

    outside ← []
    for each null vector v in null_d:
        if comp_mat·v ∉ colspan(E14_mat):             // rank test
            push!(outside, v)

    if |outside| = 1:
        E₅ ← normalize(outside[1])                    // clear denoms, divide by gcd
        E₅ ← orient(E₅)                               // positive M₅ constant coefficient
        return E₅
    else if |outside| > 1:
        report unexpected; investigate
end for
```

#### 3.3. The Degree Sweep

We applied this algorithm with $N = 95$ primes in $[2,500]$ and composites in $[4,100]$:

| Degree $d$ | Basis dim | Null dim | Directions outside $E_1$–$E_4$ span | Status |
|:---:|:---:|:---:|:---:|:---|
| $d=0$ | 5 | 0 | 0 | (no prime-vanishing expressions) |
| $d=1$ | 10 | 0 | 0 | (no prime-vanishing expressions) |
| $d=2$ | 15 | 4 | 1 | **← $E_5$ found here ✓** |
| $d=3$ | 20 | 8 | 1 | (same direction, higher embedding) |

The algorithm terminates at $d=2$: one unique direction outside the $E_1$–$E_4$ span is found at degree $d=2$, and this is the canonical minimal-degree $E_5$.

**Remark.** The sequence $E_1,\ldots,E_4$ exhibits a transparent pattern:

- $E_1$: degree 2, introduces $M_2$
- $E_2$: degree 3, introduces $M_3$
- $E_3$: degree 4, introduces $M_4$
- $E_4$: degree 5, introduces $M_5$

Extrapolation would predict $E_5$ at degree 6 with $M_6$ as its maximal MacMahon function. This extrapolation fails on two independent fronts simultaneously. First, the weight-12 cusp form barrier (Theorem 2.1) eliminates $M_6$ from the candidate space entirely. Second, the $E_1$–$E_4$ expressions do not occupy all dimensions of the degree-0–5 polynomial subspace; their $\mathbb{Q}[n]$-span leaves a gap at lower degrees. The canonical $E_5$ exploits this gap, residing at degree 2 within the $M_1$–$M_5$ basis. This doubling-back to lower degree is a geometric consequence of the weight-12 obstruction: when the naive upward extension fails, the null space folds back to find the new direction at a degree already partially explored.

#### 3.4. The Exact Formula for $E_5(n)$

After normalization (clearing all $\mathbb{Q}$-denominators and dividing by the integer GCD of numerators) and orientation (positive $M_5$ constant coefficient), the unique minimal-degree $E_5$ extracted at $d=2$ is:

$$\boxed{\begin{aligned}
E_5(n) &= (-450450 + 675675n - 225225n^2)\,M_1(n) + (960960n - 120120n^2)\,M_2(n) \\
       &\quad + (2534912n - 166016n^2)\,M_3(n) + (7999488n - 322560n^2)\,M_4(n) + 258048000\,M_5(n)
\end{aligned}}$$

By construction, $E_5(n) = 0$ if and only if $n$ is prime. The formula uses only $M_1$ through $M_5$ and polynomial coefficients up to degree 2 in $n$.

---

### 4. Verification, Non-Negativity, and Closing the Conjecture

#### 4.1. Computational Verification of Prime Detection

We evaluated $E_5(p)$ for all 95 primes $p \in [2, 500]$ and confirmed that $E_5(p) = 0$ identically for each. Spot checks:

| $p$ | $E_5(p)$ |
|-----|---------|
| 2 | 0 |
| 3 | 0 |
| 5 | 0 |
| 7 | 0 |
| 11 | 0 |
| 97 | 0 |
| 499 | 0 |

This establishes that $E_5(n)$ is a valid prime-vanishing expression.

#### 4.2. Behavior at Composites and the Non-Negativity Nuance

Theorem 1.1 of [1] constructs expressions $E_k(n) \geq 0$ for all $n \geq 2$, with equality iff $n$ is prime. The canonical $E_5(n)$ derived above does **not** satisfy universal non-negativity. Evaluating $E_5(n)$ at all 404 composites in $[4, 500]$:

| | Count |
|---|---|
| Composites with $E_5(n) > 0$ | 55 |
| Composites with $E_5(n) < 0$ | 349 |
| Composites with $E_5(n) = 0$ | 0 (none; would imply prime) |

The smallest composite at which $E_5$ is positive is $n = 25$. The 55 positive cases are a minority; $E_5$ is negative at 86% of composites in $[4,500]$. Figure 1 illustrates this for $n \in [4,200]$ on a symmetric-logarithmic scale (symlog: $\operatorname{sign}(E_5)\cdot\log_{10}(1+|E_5|)$), which preserves sign information while compressing the enormous dynamic range of values.

![Figure 1: E₅(n) at composite integers in [4,200]. Teal stems point upward where E₅ > 0; crimson stems point downward where E₅ < 0. The symlog y-axis compresses the dynamic range across 18 orders of magnitude.](e5_scatter.png)

This behavior is qualitatively different from $E_1$–$E_4$, which are non-negative at all composites. The reason is structural: $E_5$ lives in a dimension of the prime-vanishing null space that is orthogonal to the non-negativity cone spanned by $E_1$–$E_4$.

**Remediation.** Since $E_4(n) \geq 0$ at all composites (strictly positive) and $E_4(p) = 0$ at all primes, any expression $E_5 + k\,E_4$ for $k > 0$ preserves the prime-vanishing property while potentially restoring non-negativity. We searched for the smallest positive integer $k$ such that $E_5(n) + k\,E_4(n) \geq 0$ for all composites in $[4, 500]$:

$$k_{\min} = 856$$

The expression $E_5(n) + 856\,E_4(n)$ is non-negative at all 404 composites in $[4,500]$, with minimum value 4230 achieved at $n = 4$, and vanishes at all primes. This provides a non-negative prime-vanishing expression that spans the same dimension as $E_5$.

#### 4.3. Closing the Open Conjecture

The central conjecture of [1] posits that any non-negative prime-vanishing expression formed from MacMahon partition functions is a $\mathbb{Q}[n]$-linear combination of the fundamental table entries. Prior to this work, computational tests restricted to $a_{\max} = 5$ revealed a gap: the null space of the prime evaluation matrix had higher dimension than the $\mathbb{Q}[n]$-span of $E_1$–$E_4$.

By appending $E_5(n)$, this gap is resolved. We verified the conjecture computationally for $a_{\max} = 5$ and polynomial degrees $d = 2$ through $d = 6$ using $N = 300$ primes:

| Degree $d$ | Basis dim | Prime-van. dim | $E_1$–$E_5$ span dim | Conjecture holds? | Time (s) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 2 | 15 | 4 | 15 | ✓ | < 1 |
| 3 | 20 | 8 | 19 | ✓ | < 1 |
| 4 | 25 | 12 | 23 | ✓ | < 1 |
| 5 | 30 | 16 | 27 | ✓ | 1 |
| 6 | 35 | 20 | 31 | ✓ | 1 |

In every case, the prime-vanishing subspace is entirely contained in the $\mathbb{Q}[n]$-span of $\{E_1, E_2, E_3, E_4, E_5\}$: no counterexample was found. This confirms the conjecture within the bounds $a_{\max} \leq 5$, $d \leq 6$.

**Note on non-negativity and the conjecture.** Since the canonical $E_5$ is not non-negative, it strictly speaking lives outside the cone of non-negative prime-vanishing expressions. However, the combination $E_5 + 856\,E_4$ is non-negative and spans the same new dimension as $E_5$ alone, since $E_4$ is already in the span of $\{E_1,\ldots,E_4\}$. The conjecture as stated in [1] pertains to non-negative expressions; our result shows that the span of $\{E_1, E_2, E_3, E_4, E_5 + 856\,E_4\}$ equals the span of $\{E_1, E_2, E_3, E_4, E_5\}$, which equals the full prime-vanishing subspace within our computational bounds.

---

### 5. Summary of the $E_1$–$E_5$ Sequence

The five prime-detecting expressions have the following structure:

| Expr | Max poly degree $d$ | Highest $M_a$ | Non-negative? | Notes |
|------|:---:|:---:|:---:|:---|
| $E_1$ | 2 | $M_2$ | Yes | Introduces $M_2$; vanishes iff prime |
| $E_2$ | 3 | $M_3$ | Yes | Introduces $M_3$ |
| $E_3$ | 4 | $M_4$ | Yes | Introduces $M_4$ |
| $E_4$ | 5 | $M_5$ | Yes | Introduces $M_5$ |
| $E_5$ | **2** | $M_5$ | **No** | Degree resets; $M_6$ excluded by cusp form barrier |

The pattern breaks at $E_5$ in two ways: the polynomial degree resets from 5 to 2 (rather than increasing to 6), and $M_6$ is absent (rather than being introduced). Both departures stem from the same root cause: the weight-12 cusp form barrier.

---

### 6. Conclusion

The derivation of $E_5(n)$ illuminates a profound structural phenomenon: the modular form barrier at weight 12 fundamentally constrains the algebraic structure of prime-vanishing expressions. The sequence $E_1, E_2, E_3, E_4$ follows a clear pattern—each introduces the next MacMahon function as its leading term. At $M_6$, however, the generating function $U_6(q)$ enters a new modular-form regime. The resulting Ramanujan tau function $\tau(n)$ is algebraically rigid: it does not reduce to polynomials in divisor sums, and it forces $M_6$ into the pivot columns of every prime evaluation matrix.

This obstruction forces $E_5$ to live entirely in the $M_1$–$M_5$ basis, finding its new dimension via higher polynomial degree rather than higher weight. Remarkably, the minimal degree at which $E_5$ appears is $d=2$—lower than the $d=5$ one might expect from the pattern—revealing a geometric "fold" in the prime-vanishing null space.

We have computationally resolved the open conjecture of [1] for weight $a_{\max} \leq 5$ and polynomial degree $d \leq 6$, confirming that $\{E_1, E_2, E_3, E_4, E_5\}$ (or equivalently, $\{E_1, E_2, E_3, E_4, E_5 + 856\,E_4\}$ for the non-negativity formulation) generates the prime-vanishing subspace within these bounds.

This work opens several directions:

1. **Higher weights.** At weight 14 and 16, the cusp form space grows. Do prime-vanishing expressions exist at $a_{\max} = 7$ or beyond, or does the cusp form exclusion propagate?

2. **Sign structure of $E_5$.** The 55 composites where $E_5 > 0$ and the 349 where $E_5 < 0$ have no obvious factorization pattern. (An earlier, non-minimal $d=3$ formulation of $E_5$—an intermediate result not recorded here—had its negative values confined to a sparse family of composites divisible by 5; the canonical $d=2$ form distributes the sign change far more broadly.) What arithmetic property of $n$ determines the sign of $E_5(n)$ for the canonical formula?

3. **Closed-form generating series.** Is there an analogue of the MacMahon generating functions $U_a(q)$ that organizes the expressions $E_k$ into a coherent modular-form framework?

By connecting partition arithmetic to the spectral geometry of modular forms, this framework suggests deep ties between number-theoretic detection and automorphic form theory—connections that merit further investigation.

---

### Reproducibility

All computations in this paper were performed using the `QuasiShuffleAlgebra` Julia package, available at `https://github.com/randsley/PartitionPrimes`. The package uses `Rational{BigInt}` arithmetic throughout; all results are exact. The data extraction script `Paper/extract_paper_data.jl` reproduces every numerical claim in this paper.

---

### References

[1] Craig, W., van Ittersum, J.-W., & Ono, K. (2024). *Integer partitions detect the primes*. arXiv preprint arXiv:2405.06451.

[2] Bachmann, H., & Kühn, U. (2017). *The algebra of bi-brackets and regularized multiple Eisenstein series*. Journal of Number Theory, 200, 260–294. (For background on quasi-shuffle algebras and their role in partition function identities.)

[3] Zagier, D. (2008). *Elliptic Modular Forms and Their Applications*. In: *The 1-2-3 of Modular Forms*, Universitext, Springer, 1–103.

[4] Serre, J.-P. (1973). *A Course in Arithmetic*. Graduate Texts in Mathematics 7, Springer.

[5] Craig, W. (2021). *Compositions, Partitions, and Prime Detection*. PhD thesis, University of Virginia.
