packages <- c(
  "shiny", "bslib", "dplyr", "tidyr", "stringr", "ggplot2",
  "scales", "haven", "plotly", "shinyjs"
)

missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) {
  install.packages(missing)
}

