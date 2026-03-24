# Final Submission Checklist — Three-Paper Package

## Package-level consistency
- Confirm the trilogy order is stable everywhere:
  - computational paper = discovery paper
  - Paper 1 = obstruction theorem
  - Paper 2 = recovery paper
- Check that each paper refers to the others with the same titles, years, and numbering.
- Make sure the CIO reference is identical across all three papers.
- Use one consistent notation set across the trilogy:
  - `E_d`, `O_d`, `V_d^(6)`, `V_d^(7)`
  - `c_tau`, `beta`, `E_5`, `E_6`

## Acknowledgements
- Replace or delete the placeholder `[colleagues to be named]` in **Paper 1** and **Paper 2**.
- If you keep acknowledgements, use the same style in all three papers.

## Bibliography and citations
- Verify every bibliography entry is complete:
  - authors
  - title
  - journal or preprint source
  - volume/issue/pages or identifier
  - year
- Check that cross-citations between the three papers are final and consistent.
- Make sure the computational paper cites the two companion papers in the same style those papers cite each other.
- Compile and inspect for:
  - undefined citations
  - duplicate bibliography entries
  - stale labels from old titles

## Cross-references and numbering
- Recompile each paper twice and check for:
  - undefined refs
  - wrong theorem/lemma names
  - stale cleveref/autoref labels
- Specifically confirm:
  - Paper 1: theorem/corollary labels remain correct in text and proof headings
  - Paper 2: proposition/theorem references in Sections 4–5 are all correct after the latest restructuring
- Check table, figure, appendix, and equation references manually in the PDFs.

## Mathematical claim hygiene
- Make sure every statement is clearly marked as one of:
  - structural theorem
  - exact computational proposition
  - conjecture/interpretation
- In Paper 2 especially, keep the hybrid nature explicit:
  - structural in Sections 3–4
  - exact computation in the degree-threshold verification
- In the computational paper, keep bounded claims bounded and avoid drifting back into global language.

## Front matter and abstracts
- Confirm each abstract matches the actual proof level:
  - computational paper = discovery and bounded exact computation
  - Paper 1 = structural obstruction theorem
  - Paper 2 = structural/computational recovery
- Check that each introduction states exactly what that paper contributes, without overclaiming the other two.

## Figures, tables, appendices
- Verify every table and appendix is cited in the main text.
- Check that captions are informative and consistent in style.
- Confirm the appendices contain exactly what the main text promises.
- Check image resolution and print readability in the final PDFs.

## LaTeX hygiene
- Remove unused files or stale insertion snippets from the build tree if they are no longer needed.
- Check for:
  - overfull hboxes
  - underfull warnings that affect visible output
  - duplicate labels
  - unused labels if you want a cleaner source
- Confirm the shared preamble does not create hidden differences in theorem styling or reference naming across papers.

## Journal-facing polish
- Decide whether the computational paper is a preprint-style companion or a submission-ready paper. If submission-ready, match its bibliography/polish level to Papers 1 and 2.
- Standardize title capitalization, theorem naming style, and bibliography style across the set.
- Add author affiliations, emails, ORCID, keywords, MSC classes, or acknowledgements only if the target venue expects them.

## Final read-through order
1. `paper1-obstruction`
2. `paper2-recovery`
3. `paper`
4. Read all three abstracts back-to-back
5. Read all three introductions back-to-back

## Minimal final action list
- Finalize acknowledgements in Papers 1 and 2.
- Verify all bibliography entries and cross-citations across the trilogy.
- Do a clean double compile of all three.
- Read the three abstracts and introductions together.
- Submit.
