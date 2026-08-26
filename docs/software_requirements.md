# Software requirements

## Stata

- Stata 19
- `stpm3`
- `standsurv`
- `stpp`
- `rcsgen`
- `erepost`

The analysis also uses built-in Stata frames, `mlogit`, `drawnorm`, `putdocx`, and
graph commands. Package versions installed for the final analysis should be recorded
alongside any independent rerun.

## R and Shiny

- R 4.5.1
- shiny
- bslib
- dplyr
- tidyr
- stringr
- ggplot2
- scales
- haven
- plotly
- shinyjs

Run `Rscript shiny/install_packages.R` to install missing R packages.

