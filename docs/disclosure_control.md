# Data access and disclosure control

The analysis uses SEER Research Data. Patient-level data are not redistributed in this
repository. The following safeguards are applied:

- Raw, tidy, and cleaned SEER files remain in a user-specified private directory.
- `.dta`, log, history, workspace, and deployment-account files are excluded by
  `.gitignore`.
- Public CSVs contain model estimates or aggregate cancer-site/sex/age-group counts.
- No patient identifiers, record identifiers, diagnosis dates, survival times, or
  individual covariate combinations are released.
- Shiny deployment metadata and account tokens are excluded.

The code is provided to make cohort construction, modelling, and reporting transparent.
Authorised researchers can rerun it against their own licensed SEER extract.

