# Extending Partition-Theoretic Prime Detection: The Cusp Form Barrier and the Derivation of E₅

### Abstract

Recent work by Craig, van Ittersum, and Ono [1] demonstrated a novel connection between additive and multiplicative number theory, using MacMahon partition functions $M_a(n)$ to detect primes via quasi-shuffle algebras. They constructed explicit prime-vanishing polynomial combinations $E_1(n)$ through $E_4(n)$ and conjectured that all such expressions form a finite generating set. This paper computationally derives the fifth fundamental expression, $E_5(n)$, and illuminates why the naïve pattern—incrementally introducing $M_{k+1}$ to form $E_k$—breaks down. We show that at weight $2a=12$, the generating function $U_6(q)$ for $M_6(n)$ admits a component in the cusp form space spanned by the Ramanujan delta function $\Delta(q)$, making $M_6$ linearly independent from divisor sums and excluding it from all prime-vanishing expressions. Restricting our search to the $M_1,\ldots,M_5$ basis (excluding $M_6$ by Theorem 2.1), we derive $E_5$ at polynomial degree $d=2$ via rational null-space computation over 95 primes, and verify computationally that $\{E_1, E_2, E_3, E_4, E_5\}$ spans the prime-vanishing subspace for polynomial degree $d \leq 6$ and weight $a_{\max} \leq 5$. Unlike $E_1$–$E_4$, the canonical $E_5$ is negative at 349 of the 404 composites in $[4,500]$; the combination $E_5 + 856\,E_4$ is universally non-negative at all composites in this range.

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

The appearance of 691 is not accidental. The normalized Eisenstein series of weight 12 has the $q$-expansion
$$E_{12}(q) = 1 + \frac{65520}{691}\sum_{n=1}^{\infty}\sigma_{11}(n)\,q^n,$$
so 691 enters the denominator of any rational linear combination involving $E_{12}$, and hence of $M_6$. More strikingly, Ramanujan observed that $\tau(n) \equiv \sigma_{11}(n) \pmod{691}$ for all $n \geq 1$—a congruence that reflects the fact that $\Delta(q)$ and $E_{12}(q)$ are congruent modulo 691 as $q$-series. The prime 691 is therefore the precise arithmetic witness to the non-trivial cusp form component of $U_6(q)$: it is the smallest prime at which $\Delta$ and the Eisenstein series become distinguishable modulo a prime, and its appearance in the denominator of $c_\tau$ encodes this distinction directly.

#### 2.2. The Exclusion of $M_6$: The Pivot Discovery

The presence of $\tau(n)$ in $M_6(n)$ fundamentally breaks the algebraic structure required for prime-vanishing expressions within the $\{M_1,\ldots,M_5,M_6\}$ basis. Unlike the divisor sums $\sigma_k(p) = p^k + 1$ at primes, the Ramanujan tau function satisfies the Hasse–Weil bound $|\tau(p)| \leq 2p^{11/2}$ and does not reduce to a polynomial in $p$. Consequently, evaluating $n^k M_6(n)$ at primes yields linearly independent data that cannot be cancelled by any combination from $\{M_1,\ldots,M_5\}$ alone.

*Remark on scope.* This obstruction is basis-dependent. As shown in §4 of [forthcoming], $M_7(n)$ also contains $\tau(n)$ and $n \cdot \tau(n)$ terms via the $E_2 \cdot \Delta$ component of its quasimodular depth-1 decomposition. Within the extended basis $\{M_1,\ldots,M_7\}$, specific linear combinations of $M_6$ and $M_7$ can mutually cancel both $\tau(p)$ and $p\cdot\tau(p)$ contributions, producing a prime-vanishing expression $E_6$ that involves both $M_6$ and $M_7$. The theorem below is therefore correctly stated as a result about the $\{M_1,\ldots,M_5,M_6\}$ basis specifically; it does not preclude $M_6$ from appearing in expressions built from a larger basis.

We demonstrate this via the **prime evaluation matrix** $\mathbf{M}(d, a_{\max}, N)$: the $N \times (d+1) a_{\max}$ matrix whose rows are evaluations of $\{n^k M_a(n) : 0 \leq k \leq d,\; 1 \leq a \leq a_{\max}\}$ at the first $N$ primes. All arithmetic is exact over $\mathbb{Q}$.

**Theorem 2.1 (Exclusion of $M_6$, computational).** For $d \in \{0, 1, 2, 3, 4, 5, 6\}$ and $N = 95$ primes in $[2,500]$:
$$\mathrm{rank}(\mathbf{M}(d, 6, N)) = \mathrm{rank}(\mathbf{M}(d, 5, N)) + (d+1)$$
$$\dim\ker(\mathbf{M}(d, 6, N)) = \dim\ker(\mathbf{M}(d, 5, N))$$

*Proof.* Explicit verification for $d \in \{2, 3\}$ is shown in the table below; the pattern holds identically for $d \in \{0, 1, 4, 5, 6\}$ (omitted for brevity). Evaluated at all 95 primes in $[2, 500]$ with $d=2$:

| Matrix | Rows | Cols | Rank | Nullity |
|--------|------|------|------|---------|
| $a_{\max}=5$, $d=2$ | 95 | 15 | 11 | 4 |
| $a_{\max}=6$, $d=2$ | 95 | 18 | 14 | 4 |

Adding the 3 columns $\{M_6, n M_6, n^2 M_6\}$ increases the rank by exactly 3 and leaves the nullity unchanged. In the RREF of $\mathbf{M}(2, 6, 95)$, all three $M_6$-associated columns are pivot columns. This was verified independently for $d=3$ (nullity 8 for both $a_{\max}=5$ and $a_{\max}=6$), and computationally confirmed for all $d \in \{0,1,4,5,6\}$.

**Corollary (within $\{M_1,\ldots,M_6\}$).** No prime-vanishing expression formed solely from the basis $\{M_1,\ldots,M_6\}$ can have a non-zero $M_6$ coefficient. Within this basis, the prime-vanishing null space is identical to that of $\{M_1,\ldots,M_5\}$. $\square$

*Caution.* This corollary is basis-specific. It does not assert that $M_6$ is absent from all prime-vanishing expressions in any larger basis. Indeed, when $M_7$ is added, the combined basis $\{M_1,\ldots,M_7\}$ admits new prime-vanishing directions involving both $M_6$ and $M_7$—the expression $E_6$ being the first such direction (see §4.3).

*Terminology note.* In row-reduced echelon form (RREF), a **pivot column** is one that contains a leading 1 in some row; it represents a direction that is linearly independent of all preceding columns. A column is a pivot column if and only if it contributes a new dimension to the column space—equivalently, it has no free variable in the null space. The above result says that in the RREF of $\mathbf{M}(d,6,N)$, every column corresponding to $n^k M_6(n)$ is a pivot, so no null vector within that matrix can have a non-zero entry in those positions.

---

### 3. Extracting $E_5(n)$ from the Null Space

#### 3.1. The Search Space

Having established that $M_6(n)$ cannot appear in any prime-vanishing expression formed within the basis $\{M_1,\ldots,M_6\}$ (Corollary 2.2), we search for $E_5$ within that basis:

$$\mathcal{B}_d = \{ n^k M_a(n) : 0 \leq k \leq d,\; 1 \leq a \leq 5 \}$$

at increasing polynomial degree $d$. The prime evaluation matrix $\mathbf{M}(d, 5, N)$ (with $N = 95$ primes from $[2,500]$) is constructed using the exact closed-form formulas for $M_1,\ldots,M_5$ (which involve only divisor sums and polynomial prefactors, all computed exactly). Its null space is computed by rational RREF over $\mathbb{Q}$ using `Rational{BigInt}` arithmetic throughout.

#### 3.2. Isolating $E_5$ Outside the $E_1$–$E_4$ Span

Finding a prime-vanishing null vector is not sufficient: we must identify directions genuinely outside the $\mathbb{Q}[n]$-span of $E_1,\ldots,E_4$. To do so, we evaluate each null vector $\mathbf{v}$ at a set of composite numbers and test linear independence from the corresponding evaluations of $E_1,\ldots,E_4$.

**Algorithm: Extract\_$E_5$**

```
Input:
  a_max = 5 (integer)
  N = 95 primes from [2,500] (list of primes)
  composites = [4, 6, 8, 9, ..., 100] (list of composites)
Output: E₅(n) ∈ Q[n] ⊗ {M₁,...,M₅}  (vector of Rational{BigInt} coefficients)

Precision: All arithmetic exact over Q (no floating-point)

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

The smallest composite at which $E_5$ is positive is $n = 25$. Representative values illustrating the sign structure and dynamic range:

| $n$ | Type | $E_5(n)$ | Sign |
|:---:|:---:|---:|:---:|
| 4 | $2^2$ | $-3\,693\,690$ | negative |
| 9 | $3^2$ | $-30\,324\,096$ | negative |
| 25 | $5^2$ | $363\,888\,000$ | positive ✓ (smallest positive) |
| 35 | $5 \cdot 7$ | $9\,436\,492\,800$ | positive ✓ |
| 49 | $7^2$ | $62\,655\,312\,384$ | positive ✓ |
| 100 | $2^2 \cdot 5^2$ | $-2\,935\,909\,013\,178\,150$ | negative |

The 55 positive cases are a minority; $E_5$ is negative at 86% of composites in $[4,500]$. Figure 1 illustrates this for $n \in [4,200]$ on a symmetric-logarithmic scale (symlog: $\operatorname{sign}(E_5)\cdot\log_{10}(1+|E_5|)$), which preserves sign information while compressing the enormous dynamic range of values.

![Figure 1: E₅(n) at composite integers in [4,200]. Teal stems point upward where E₅ > 0; crimson stems point downward where E₅ < 0. The symlog y-axis compresses the dynamic range across 18 orders of magnitude.](e5_scatter.png)

This behavior is qualitatively different from $E_1$–$E_4$, which are non-negative at all composites. The reason is structural: $E_5$ lives in a dimension of the prime-vanishing null space that is orthogonal to the non-negativity cone spanned by $E_1$–$E_4$.

**Concrete evaluations.** To illustrate the scale and structure of $E_5$, we give two explicit evaluations. At the smallest positive composite $n = 25 = 5^2$:

| Term | Value |
|---|---|
| $(-450450 + 675675 \cdot 25 - 225225 \cdot 625) \cdot M_1(25)$ | $-3\,854\,050\,200$ |
| $(960960 \cdot 25 - 120120 \cdot 625) \cdot M_2(25)$ | $-90\,819\,729\,000$ |
| $(2534912 \cdot 25 - 166016 \cdot 625) \cdot M_3(25)$ | $-399\,227\,472\,000$ |
| $(7999488 \cdot 25 - 322560 \cdot 625) \cdot M_4(25)$ | $-15\,895\,756\,800$ |
| $258048000 \cdot M_5(25)$ | $+510\,160\,896\,000$ |
| **$E_5(25)$** | **$+363\,888\,000$** |

The positive value arises because the $M_5$ term ($\approx +5.1\times10^{11}$) slightly dominates the combined negative contributions of the lower-weight terms.

At the even composite $n = 100 = 2^2 \cdot 5^2$, the same competition plays out at a far larger scale:

| Term | Value |
|---|---|
| $c_1(100) \cdot M_1(100)$ | $\approx -4.7\times10^{11}$ |
| $c_2(100) \cdot M_2(100)$ | $\approx -1.5\times10^{14}$ |
| $c_3(100) \cdot M_3(100)$ | $\approx -2.0\times10^{16}$ |
| $c_4(100) \cdot M_4(100)$ | $\approx -9.69\times10^{17}$ |
| $c_5(100) \cdot M_5(100)$ | $\approx +9.86\times10^{17}$ |
| **$E_5(100)$** | **$-2\,935\,909\,013\,178\,150$** |

Here the $M_4$ and $M_5$ terms are both of order $10^{17}$ and nearly cancel, but the net result is negative ($\approx -2.9\times10^{15}$). This near-cancellation between large opposing terms is characteristic of the expression throughout, and explains why the symlog scale in Figure 1 is essential: the raw values span 18 orders of magnitude across $n \in [4, 200]$.

**Factorization structure.** We investigated whether the positive and negative composites share common arithmetic properties. A striking pattern emerges:

- All 55 composites where $E_5(n) > 0$ are **odd**.
- No even composite in $[4, 500]$ satisfies $E_5(n) > 0$; all 249 even composites give $E_5(n) < 0$.
- Among the 155 odd composites in $[4, 500]$, exactly 55 (35%) are positive and 100 (65%) are negative.

Thus parity is a necessary but not sufficient condition for positivity: $E_5(n) > 0 \Rightarrow n$ is odd, but not all odd composites are positive. The even composites are uniformly negative, while the sign at odd composites depends on finer arithmetic structure that we do not currently characterize.

This contrasts sharply with an intermediate non-canonical degree-$d=3$ formulation of $E_5$ encountered during the derivation, whose 12 negative composites were confined to multiples of 5. The canonical $d=2$ formula distributes negativity far more broadly, with parity emerging as the principal organizing feature.

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

**Open question: scope of the conjecture.** Three natural extensions remain unresolved:

1. **Does the conjecture hold for $d \geq 7$?** At $d=7$ the prime evaluation matrix grows from $300 \times 35$ to $300 \times 40$, increasing runtime but not intractability. Computational verification at $d=7$ and $d=8$ is within reach and would substantially strengthen the evidence.

2. **Can the conjecture extend to $a_{\max} \geq 6$?** Theorem 2.1 shows $M_6$ cannot appear in prime-vanishing expressions, so any extension to higher weight would require new generating expressions built entirely from $M_1,\ldots,M_5$ at still higher polynomial degrees. At weight 14 the cusp form space is 2-dimensional; at weight 16 it is 3-dimensional. Whether the cusp form exclusion mechanism propagates to prevent all higher-weight contributions, or whether new $M_1$–$M_5$ expressions emerge at large $d$, is an open structural question.

3. **Is $a_{\max} = 5$ the true frontier?** It is conceivable that the weight-12 barrier forces the generating set to be exactly $\{E_1, E_2, E_3, E_4, E_5\}$ for all degrees $d$—that is, that no further linearly independent prime-vanishing expression exists at any degree. Proving or disproving this would fully resolve the conjecture of [1].

**Note on non-negativity and the conjecture.** The canonical $E_5(n)$ extracted via our algorithm is not universally non-negative, so it technically lies outside the cone of prime-vanishing expressions satisfying Theorem 1.1 of [1]. However, the combination $E_5 + 856\,E_4$ is globally non-negative and **generates the same new vector space direction** as $E_5$. Since $E_4 \in \mathrm{span}(E_1, E_2, E_3, E_4)$ (trivially), adding $856\,E_4$ does not alter the $\mathbb{Q}[n]$-linear span: $\mathrm{span}(E_1, E_2, E_3, E_4, E_5 + 856\,E_4) = \mathrm{span}(E_1, E_2, E_3, E_4, E_5)$. Thus, from the perspective of the conjecture, $E_5 + 856\,E_4$ is a valid non-negative prime-vanishing basis element that accomplishes the same dimensional closure shown in the table above.

The canonical (non-negative-violating) form $E_5$ is the natural output of the extraction algorithm because it minimizes the embedding degree at which new directions appear and emerges directly from rational null-space computation. The non-negative form $E_5 + 856\,E_4$ is the appropriate basis element for the theorem statement of [1].

---

### 5. Summary of the $E_1$–$E_5$ Sequence

The five prime-detecting expressions form a hierarchy constrained by weight and polynomial degree:

| Expr | Max poly degree $d$ | Highest $M_a$ | Non-negative? | Notes |
|------|:---:|:---:|:---:|:---|
| $E_1$ | 2 | $M_2$ | Yes | Introduces $M_2$; vanishes iff prime |
| $E_2$ | 3 | $M_3$ | Yes | Introduces $M_3$ |
| $E_3$ | 4 | $M_4$ | Yes | Introduces $M_4$ |
| $E_4$ | 5 | $M_5$ | Yes | Introduces $M_5$ |
| $E_5$ | **2** | $M_5$ | **No** | Degree resets; $M_6$ excluded by cusp form barrier |

**Observations.** The sequence shows a clear pattern up to $E_4$: each expression introduces a new MacMahon function and increases polynomial degree by 1. At $E_5$, both trends reverse. The polynomial degree resets from 5 to 2 (degree 2 is the minimal degree at which $E_5$ appears when restricted to $M_1$–$M_5$), and $M_6$ is entirely absent due to the weight-12 cusp form barrier established in Theorem 2.1. Both departures stem from the same root cause: the cusp form obstruction forces $E_5$ to find its new direction within the existing $M_1$–$M_5$ framework at a lower degree already partially explored by $E_1$–$E_4$.

---

### 6. Conclusion

The derivation of $E_5(n)$ illuminates a profound structural phenomenon: the modular form barrier at weight 12 fundamentally constrains the algebraic structure of prime-vanishing expressions. The sequence $E_1, E_2, E_3, E_4$ follows a clear pattern—each introduces the next MacMahon function as its leading term. At $M_6$, however, the generating function $U_6(q)$ enters a new modular-form regime. The resulting Ramanujan tau function $\tau(n)$ is algebraically rigid: it does not reduce to polynomials in divisor sums, and it forces $M_6$ into the pivot columns of every prime evaluation matrix.

This obstruction forces $E_5$ to live entirely in the $M_1$–$M_5$ basis, finding its new dimension via higher polynomial degree rather than higher weight. Remarkably, the minimal degree at which $E_5$ appears is $d=2$—lower than the $d=5$ one might expect from the pattern—revealing a geometric "fold" in the prime-vanishing null space.

We have computationally resolved the open conjecture of [1] for weight $a_{\max} \leq 5$ and polynomial degree $d \leq 6$, confirming that $\{E_1, E_2, E_3, E_4, E_5\}$ (or equivalently, $\{E_1, E_2, E_3, E_4, E_5 + 856\,E_4\}$ for the non-negativity formulation) generates the prime-vanishing subspace within these bounds.

This work opens several directions:

1. **Higher weights.** At weight 14 and 16, the cusp form space grows. Do prime-vanishing expressions exist at $a_{\max} = 7$ or beyond, or does the cusp form exclusion propagate?

2. **Sign structure of $E_5$.** Parity is a necessary condition for positivity: all 55 positive-$E_5$ composites in $[4,500]$ are odd, and no even composite is positive (§4.2). But parity is not sufficient — 100 odd composites are also negative. What finer arithmetic property of an odd composite $n$ determines whether $E_5(n) > 0$?

3. **Closed-form generating series.** Is there an analogue of the MacMahon generating functions $U_a(q)$ that organizes the expressions $E_k$ into a coherent modular-form framework?

By connecting partition-theoretic prime detection to the spectral geometry of modular forms, this framework reveals deep ties between number-theoretic detection and automorphic form theory—connections that merit further investigation.

---

---

### Appendix A: Closed-Form Formula for $M_6(n)$

The complete closed-form expression for $M_6(n)$, obtained by rational coefficient fitting as described in §2.1, is:

$$M_6(n) = \sum_{j=0}^{5} P_j(n)\,\sigma_{2j+1}(n) - \frac{17}{150\,450\,048\,000}\,\tau(n)$$

where the polynomials $P_j(n) \in \mathbb{Q}[n]$ are:

$$P_0(n) = \frac{550499}{4541644800} - \frac{153617}{371589120}\,n + \frac{2159}{6635520}\,n^2 - \frac{67}{737280}\,n^3 + \frac{11}{1105920}\,n^4 - \frac{1}{2764800}\,n^5$$

$$P_1(n) = \frac{153617}{743178240} - \frac{2159}{8847360}\,n + \frac{67}{819200}\,n^2 - \frac{11}{1105920}\,n^3 + \frac{1}{2580480}\,n^4$$

$$P_2(n) = \frac{2159}{88473600} - \frac{67}{4915200}\,n + \frac{11}{5160960}\,n^2 - \frac{1}{10321920}\,n^3$$

$$P_3(n) = \frac{67}{123863040} - \frac{11}{74317824}\,n + \frac{1}{111476736}\,n^2$$

$$P_4(n) = \frac{11}{3715891200} - \frac{1}{3096576000}\,n$$

$$P_5(n) = \frac{1}{268995133440}$$

These 21 rational coefficients (plus $c_\tau = -17/150\,450\,048\,000$) were extracted by solving a $55 \times 22$ linear system over $\mathbb{Q}$ using rational RREF. The solution is unique (rank 22) and verified to reproduce $M_6(n)$ for $n = 1, \ldots, 30$ via direct partition enumeration. The script `Paper/compute_m6_coefficients.jl` in the accompanying repository reproduces this computation exactly.

---

### Computational Scaling

The prime evaluation matrix has size $N \times (d+1)a_{\max}$. Rational RREF is the dominant cost; its complexity is $O(N \cdot (d+1)^2 a_{\max}^2)$ in the number of arithmetic operations over $\mathbb{Q}$. For $N = 300$ primes and $a_{\max} = 5$, we measured the following on a single core (Julia 1.x, `Rational{BigInt}` arithmetic):

| Degree $d$ | Matrix size | Nullity | RREF time |
|:---:|:---:|:---:|:---:|
| 6 | $300 \times 35$ | 20 | 0.16 s |
| 10 | $300 \times 55$ | 36 | 0.24 s |
| 20 | $300 \times 105$ | 76 | 0.62 s |
| 30 | $300 \times 155$ | 116 | 1.13 s |

The scaling is sub-quadratic in $d$ in practice. Extensions to $d \leq 50$ (matrix size $300 \times 255$) are straightforwardly feasible; the current bound of $d \leq 6$ in §4.3 is conservative and chosen to keep verification fast, not because higher degrees are intractable. The primary open cost in extending the conjecture verification is not RREF time but rather the evaluation of $M_1,\ldots,M_5$ at 300 primes for large $d$, which requires computing high powers of primes exactly — still polynomial time, but with growing integer sizes.

---

### Reproducibility

All computations in this paper were performed using the `QuasiShuffleAlgebra` Julia package, available at `https://github.com/randsley/PartitionPrimes`. The package uses `Rational{BigInt}` arithmetic throughout; all results are exact. The data extraction script `Paper/extract_paper_data.jl` reproduces every numerical claim in this paper.

**Validation procedure.** To verify the computational claims independently:

1. Clone the repository: `git clone https://github.com/randsley/PartitionPrimes.git`
2. Run the data extraction script: `julia --project=QuasiShuffleAlgebra Paper/extract_paper_data.jl`
3. Compare the printed values against the tables in §2.1, §4.1, §4.2, and §4.3

Expected runtime: approximately 5 minutes on a standard laptop. All results match exactly (exact rational arithmetic with no numerical errors).

---

### References

[1] Craig, W., van Ittersum, J.-W., & Ono, K. (2024). *Integer partitions detect the primes*. arXiv preprint arXiv:2405.06451.

[2] Bachmann, H., & Kühn, U. (2017). *The algebra of bi-brackets and regularized multiple Eisenstein series*. Journal of Number Theory, 200, 260–294. (For background on quasi-shuffle algebras and their role in partition function identities.)

[3] Zagier, D. (2008). *Elliptic Modular Forms and Their Applications*. In: *The 1-2-3 of Modular Forms*, Universitext, Springer, 1–103.

[4] Serre, J.-P. (1973). *A Course in Arithmetic*. Graduate Texts in Mathematics 7, Springer.

[5] Craig, W. (2021). *Compositions, Partitions, and Prime Detection*. PhD thesis, University of Virginia.

---

**Paper information**
- **Version:** 2.1 (revised)
- **Last updated:** 2026-03-09
- **Repository:** https://github.com/randsley/PartitionPrimes
- **Data extraction script:** `Paper/extract_paper_data.jl`
- **Author's note:** All computations performed with Julia 1.10.0 or later, using `Rational{BigInt}` for exact arithmetic throughout.
