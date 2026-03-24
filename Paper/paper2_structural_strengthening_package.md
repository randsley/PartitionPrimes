# Paper 2 structural strengthening package

This package contains two insertion-ready components for the recovery paper.

## Included files

1. `paper2_section4_full_spanning_statement.tex`  
   Contains:
   - the definition of `O_d`,
   - a lemma showing `E_d` is polynomial on primes,
   - a non-cancellation lemma for `Q(p) tau(p)`,
   - a corollary proving linear independence of polynomial multiples of `tau`,
   - and a full spanning proposition showing
     `O_d = Span([tau], [p tau], ..., [p^(d+1) tau])`
     with dimension `d+2`.

2. `paper2_section5_minimality_rewrite.tex`  
   Contains a rewrite of Section 5 that:
   - uses the stronger Section 4 result,
   - separates structural and computational content,
   - makes the cusp cancellation constraint central,
   - treats degrees 0 and 1 structurally,
   - treats degrees 2 and 3 by exact computation,
   - and identifies degree 4 as the first recovery.

## Usage note

The Section 4 proposition is written conditionally on the prime-level decompositions
for `M_6(p)` and `M_7(p)`. If those remain computationally established in the paper,
this is still appropriate, but the surrounding prose should say that the proposition
uses the decompositions established in Section 3.

If those decompositions are later derived structurally, the proposition can stand
without further qualification.
