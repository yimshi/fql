## Release checklist (draft)

- [ ] Run `R CMD check --as-cran .`
- [ ] Run `R -q -e "devtools::test()"`
- [ ] Run `R -q -e "pkgdown::build_site()"` (if website is configured)
- [ ] Confirm `NEWS.md` updated for release version
- [ ] Confirm `inst/CITATION` is current
- [ ] Tag release and publish notes
