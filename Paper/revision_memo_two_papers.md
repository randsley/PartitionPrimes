# Revision Memo for the Two Draft Papers

## Paper 1: *The Weight-12 Cusp Obstruction in Partition-Theoretic Prime Detection*

Overall: this is the stronger of the two drafts. It already has the shape of a theorem paper, but it needs a proof-tightening pass and citation cleanup before submission.

### Abstract
The abstract is good in structure and states a real theorem. Two edits:
- Replace “linearly independent of all polynomials on the primes” with “not equal to any polynomial on the primes” unless you explicitly prove the stronger linear-independence formulation in the body. What the paper actually uses is non-polynomiality and the stronger “polynomial multiple cannot be polynomial” lemma.
- The phrase “for any polynomial degree d≥0” is fine because the theorem is stated globally, but then the proof must remain fully structural and not lean on bounded computation at any critical point. That is almost true now, but Section 4 still needs tightening.

### Section 1.1 Context
This section is clear and motivates the theorem well. Two issues:
- Check the citation to the original Craig–van Ittersum–Ono paper. The reference listed later as [3] does not match the prime-detection paper you are discussing. That needs correction.
- The sentence “This closes the story at weights ≤10” is slightly broad. It would be safer to say it closes the Eisenstein-only level-1 prime-detection characterization in that range.

### Section 1.2 Main Results
This section is strong and well focused. One concrete correction:
- In Remark 1.3 you write “Theorem 1.2 rules out M6...” but Theorem 1.2 is not the obstruction theorem there; this is a numbering mistake and should be fixed. It should refer to Theorem 1.1 or Corollary 1.2, depending on what you mean.

### Section 2 Background
This section is fine as mathematical setup. I would make one stylistic change:
- In Section 2.3, when you say U6 has non-trivial projection onto the cusp-form space, add either a citation or a forward pointer to Lemma 3.2 where the nonzero coefficient c_tau is established. Right now the logic is slightly split across sections.

### Section 3 Closed Form of M6
This section is useful, but it is also where the paper is most vulnerable if a referee asks “what is proved and what is fit?”
- Lemma 3.1 is fine.
- Lemma 3.2 is acceptable if you are explicit that the exact closed form is obtained by exact rational fitting and then used as an exact identity in the paper. But if you want this to feel more theorem-driven, add one sentence justifying why the fit is uniquely determined by the ansatz space and the number of sample points.
- The phrase “Theorem 4.1 below” inside Lemma 3.2 should be checked for numbering consistency after final edits.

### Section 4 Linear Independence of τ over Polynomials at Primes
This is the heart of the paper, and it needs the most tightening.
- Lemma 4.1 is basically sound. The argument using Deligne plus Sato–Tate is good. I would sharpen one sentence: instead of “g(p) is eventually positive or eventually negative,” say “a nonzero polynomial has finitely many real zeros and therefore constant sign for all sufficiently large real arguments.” That reads more like a proof.
- Remark 4.2 is good. Keep it.
- Lemma 4.3 is where I would spend revision time. The second paragraph is too compressed. You should explicitly say:
  1. since Q is nonzero, Q(p)≠0 for all sufficiently large primes outside a finite exceptional set;
  2. hence τ(p)=g(p)/Q(p) on an infinite cofinal set of primes;
  3. a nonzero rational function over Q has constant sign for all sufficiently large real inputs avoiding poles;
  4. this contradicts Sato–Tate sign oscillation.

  That would make the lemma fully referee-proof.

### Section 5 Proof of the Obstruction Theorem
This section is conceptually clean and nearly ready.
- The decomposition into A(p) and Q(p) is exactly the right proof format.
- The phrase “Theorem 4.3 says” should be corrected to “Lemma 4.3 says” if that is the final numbering. There are theorem/lemma mismatches here.
- Once the numbering is fixed and Lemma 4.3 is tightened, this section should work well.

### Section 6 Arithmetic of the Obstruction
This section is a nice bonus. It is not essential to the theorem, but it gives the paper personality.
- Keep it, but make sure every arithmetic remark remains clearly secondary to the proof.
- The reference to “Paper 2, [10]” should be checked against the actual bibliography numbering.

### Section 7 Computational Verification
This section is good as supporting evidence.
- Label it explicitly as computational confirmation of a theorem proved independently, which you already largely do.
- The proposition here is fine. It strengthens trust without carrying proof burden.

### Section 8 Open Questions / Section 9 Conclusion
Both are fine.
- The open questions are good research-facing questions.
- The conclusion is acceptable, but keep it short; the theorem is the center of gravity.

### References
This needs cleanup.
- The original CIO paper citation appears wrong.
- Cross-references to the companion paper and the computational paper must be checked carefully.
- Make sure every cited preprint actually exists in the form stated.

### Bottom line for Paper 1
Revise, do not redesign. The main tasks are:
- fix theorem/lemma numbering,
- correct the bibliography,
- tighten Lemma 4.3,
- add one sentence clarifying the exact-fitting uniqueness in Lemma 3.2.

---

## Paper 2: *Cusp Cancellation and the First Recovery of Prime-Vanishing Relations Beyond Weight 12*

Overall: this draft has a strong idea, but it is not yet a theorem paper in its current form. The main issue is mismatch between rhetoric and proof level.

### Abstract
This is the first place I would revise.
- The abstract says “we prove” several things that are not yet structurally proved in the body: two-dimensionality of O_d, minimality of degree 4, and structural derivation of E6. Right now these are partly structural and partly exact computational.
- Decide now whether the paper is:
  - a theorem paper, or
  - a structural-computational paper.

  In its current state, it is the second, not the first.

### Section 1.1 The Two-Paper Arc
This section is good.
- It explains the motivation cleanly.
- Keep it concise; it helps anchor the companion-paper relationship.

### Section 1.2 Main Results
This section overstates the current level of proof.
- Theorem 1.1 says the smallest degree is 4 and the first such direction is E6. But Section 5 establishes that primarily by exact rational computation, not by a structural theorem.
- Theorem 1.2 says O_d is two-dimensional and E6 is unique up to scalar. The uniqueness via nullity gain is computationally fine; the exact structural two-dimensionality argument is not yet fully rigorous in the quotient-space sense.

My recommendation: either relabel these as computational theorems/propositions, or weaken the verbs in the statements.

### Section 2 Background
This section is useful, but Remark 2.1 is too ambitious.
- The sentence that the quasimodular completion of U7 contributes both τ(p) and pτ(p) is exactly what the paper later still treats through fitting. So here it reads as if already established.
- Rephrase this remark as heuristic motivation, not as settled structural fact.

### Section 3 The Cusp Structure of M7 at Primes
This is the main unresolved gap.
- Lemma 3.1 is not yet a proof of the claimed form
  M7(p)=f7(p)+ατ(p)+β pτ(p), with β≠0.
- The step using E2Δ to get τ(p)+p is fine as a motivational calculation, but it does not get you pτ(p).
- The sentence “to obtain pτ(p), one uses the depth-1 quasimodular decomposition of U7” is not a proof. Then the argument falls back to fitting a 30-column rational system. That means the result is computationally established, not structurally proved.

This section needs one of two fixes:
1. Honest computational framing: state that the prime-level form is established by exact rational fitting, with modular structure as motivation.
2. Actual theorem upgrade: give a real derivation from the quasimodular decomposition of U7, including why β≠0.

### Section 4 The Two-Dimensional Obstruction Space
This section has the right concept but the proof needs sharpening.
- Lemma 4.1 is plausible, though the phrase “any rational function is eventually monotone” is not the best justification for sign behavior. Use “eventually constant sign away from poles” instead.
- Proposition 4.2 is not yet rigorous enough. The sentence saying higher p^jτ(p) still span a two-dimensional space because they lie in Q[p]·[τ] ∪ Q[p]·[pτ] is heuristic, not a proof of dimension in the quotient ev(V_d^(7))/E_d.

You need a precise quotient-space argument. At minimum:
- define the quotient map carefully,
- show the images of [τ] and [pτ] are nonzero and independent,
- show every cusp contribution from the relevant basis lies in their span.

### Section 5 Minimality of Degree 4
This section is the most important one rhetorically, and right now it is still computational at the decisive steps.
- The cusp constraint
  c_τ Q6(p)+(α+βp)Q7(p)=0
  is good and should stay.
- But Lemma 5.1 and Lemma 5.2 prove minimality and existence mainly by exact nullity computations. That is fine for a computational paper, but not for a paper whose abstract says the result is proved structurally.
- Also, the proof references “Theorem 5.1” and “Theorem 5.2” where the section appears to define lemmas. The numbering and labels are inconsistent.

This section should either:
- become explicitly computational, or
- be substantially strengthened to derive degree 4 minimality from dimension counting plus structural constraints.

### Section 6 Explicit Formula for E6
This section is good as a computational result section.
- The explicit formula belongs here.
- Verification and arithmetic content are useful.
- Non-negativity discussion is fine.

But it is a consequence of exact null-space extraction, not yet of the structural theory alone. Keep that distinction explicit.

### Section 7 Bounded Completeness and Weight-24 Barrier
This is a strong section conceptually.
- The caution around weight 24 is good.
- The conjectures are well framed.
- Keep the wording conditional and exploratory, as you already do.

### Section 8 Conclusion
The conclusion currently says the two papers establish a “complete two-paper arc.” That is attractive, but slightly too strong for Paper 2 as drafted.
- If you keep Paper 2 in theorem mode, this sentence invites challenge.
- If you reframe Paper 2 as structural-computational, the conclusion becomes accurate.

### References
This section needs attention.
- The citation [1] again appears to be the wrong CIO paper.
- Make sure the companion paper and computational paper references are consistent with final titles and dates.
- Check all preprint references carefully.

### Bottom line for Paper 2
You need to choose between two revision strategies.

#### Strategy A: make it an honest structural-computational paper
This is the fastest and best current option.
- change “we prove” to “we show computationally” or “we establish in exact computation” where appropriate,
- relabel theorems as computational propositions where needed,
- present Section 3 as modular motivation plus exact rational fit,
- present Section 5 as exact minimality evidence rather than structural theorem.

#### Strategy B: make it a true theorem paper
This is possible, but it is much more work.
You would need to:
- prove the exact prime-level form of M7 structurally,
- prove β≠0 without using fit output as the decisive argument,
- prove Proposition 4.2 rigorously in the quotient-space language,
- prove degree-4 minimality without leaning on nullity computation as the final step.

Right now, Strategy A is the better move.

## Recommendation order
1. **Finish Paper 1 first.**
2. **Reframe Paper 2 as structural-computational unless you are prepared to prove the missing structural steps.**

## Next best follow-up
The next useful deliverable would be:
- a **concrete revision checklist** for Paper 1, and
- a **rewritten abstract + theorem section** for Paper 2 under Strategy A.
