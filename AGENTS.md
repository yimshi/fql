# AGENTS.md — Repository rules for Codex (fql R package)

These instructions guide Codex when working in this repository.
They are meant to keep PRs reviewable, preserve statistical behavior, and ensure R package quality.

## Project context
- This is an R package implementing a Flexible Quasi-Likelihood (FQL) model for microbiome count data.
- Primary goals:
  1) Correctness & numerical stability
  2) Reproducibility (deterministic, well-tested)
  3) Engineering quality (R CMD check clean, docs, examples, CI)

## Operating mode (how to work)
- Prefer **small, reviewable PRs**. One PR = one theme.
  - Target: < ~300 net LOC change unless explicitly requested.
- For multi-step work: propose a short plan first (milestones + acceptance checks), then execute milestone-by-milestone.
- Do not broaden scope. Follow the task’s “Scope” and “Non-goals”.

## Non-negotiable quality gates
Every PR must:
- Pass `R CMD check` (0 errors, 0 warnings; minimize notes).
- Pass unit tests (`devtools::test()`).
- Keep the package installable from a clean environment.

If you add or change core fitting logic, you MUST add/extend regression tests that would fail without the change.

## Commands to run (and report in PR)
Run and summarize results in the PR description:
- `R CMD check`
- `R -q -e "devtools::test()"`
Optional but recommended when touching docs/site:
- `R -q -e "devtools::document()"`
- `R -q -e "pkgdown::build_site()"` (only when pkgdown is configured)

## Dependency policy
- Avoid adding new dependencies unless clearly justified.
- Prefer base R / recommended packages.
- If adding a dependency:
  - Put runtime deps in `Imports`.
  - Put dev-only deps in `Suggests` (testthat, roxygen2, knitr, rmarkdown, pkgdown, covr, lintr, etc).
  - Update DESCRIPTION accordingly.

## R package conventions
- Use **roxygen2** for documentation; do not hand-edit `.Rd` files.
- Avoid `exportPattern` in NAMESPACE. Export only intended public functions.
- Keep helpers internal:
  - Use `.`-prefixed function names (e.g., `.update_beta()`) or place in `R/utils-*.R`.
- Never set global options (e.g., `options(warn = -1)`); handle warnings locally.

## API & backward compatibility
- Do not change function signatures or return structures unless the task explicitly requests it.
- If an API change is required:
  - Provide a compatibility path (deprecated alias, clear warning, transition notes).
  - Update README/examples/tests accordingly.

## Numerical stability & safety
- Ensure variance estimates are strictly positive:
  - Use `eps` guards (e.g., `pmax(v, eps)`).
  - Avoid `log(0)` by adding `eps` before log transforms.
- Guard against:
  - Overflow/underflow in `exp(eta)` (consider clamping eta if needed).
  - Singular matrices (use stable solvers; provide clear error messages).
  - Missing/invalid inputs (validate early; informative errors).
- No silent failure:
  - If convergence fails, return a clear `converged = FALSE` with diagnostics.

## Testing guidelines (testthat)
Minimum test coverage expectations:
- Smoke test: fit runs end-to-end on a small synthetic dataset.
- Regression tests for historical bugs (e.g., inner-loop beta usage, log(0) issues).
- Numerical invariants:
  - variance > 0
  - fitted values finite
  - solver returns finite coefficients or clear errors
- If changing math: compare against a baseline (Poisson GLM in Var=mu special case), within tolerances.

## Documentation guidelines
- Public functions must have:
  - @param, @return, @examples (examples must be runnable and fast).
- Keep examples lightweight (avoid long simulations in @examples).
- Prefer vignettes for longer narratives or simulations.

## Repo hygiene
- Do not commit IDE/user artifacts:
  - `.Rhistory`, `.Rproj.user/`, `*.Rproj`, etc.
- Maintain `.Rbuildignore` / `.gitignore` to exclude these.

## PR description template (include in PR body)
- Summary: what changed and why
- Scope: files/areas touched
- Tests: output summary of R CMD check + devtools::test
- Risks/notes: anything reviewers should watch
- Next steps: follow-ups if any (kept out of this PR)

## Review guidelines (optional but preferred)
- Preserve statistical intent; do not “improve” the method without explicit request.
- Keep code readable and well-commented at algorithm boundaries (not every line).
- Prefer vectorized linear algebra (`crossprod`, `tcrossprod`) over slow loops where appropriate.