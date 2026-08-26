# Private data inputs

Patient-level SEER data are not included in this repository and must not be committed.

After obtaining authorised SEER access, create a private directory with this layout:

```text
private-data/
  raw/
    Breast.dta
    Colon.dta
    Corpus_and_uterus.dta
    Lung_and_bronchus.dta
    Oral_cavity_and_pharynx.dta
    Ovary.dta
    Pancreas.dta
    Rectum.dta
    popmort.dta
  tidy/
  clean/
```

The eight site files are SEER*Stat case-list exports with fields `v1`-`v12` in the
order documented in `docs/variable_dictionary.csv`. The mortality file is a SEER
population mortality export with fields `v1`-`v5`.

Set `PRIVATE_DATA_ROOT` in `code/00_config.do` to the absolute path of `private-data`.
The cleaning scripts write only to its `tidy/` and `clean/` subfolders.

The public `results/` directory contains only aggregate model estimates and aggregate
counts. It contains no direct identifiers, row-level records, or SEER case data.

