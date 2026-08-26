# Shiny application

The application reads `data/Shiny.csv`, a disclosure-safe aggregate result file. It
does not read patient-level SEER data.

## Run locally

From the repository root:

```r
shiny::runApp("shiny")
```

## Deploy

GitHub cannot execute a Shiny server. Deploy the `shiny/` folder to shinyapps.io,
Posit Connect, or another Shiny host. For shinyapps.io:

```r
rsconnect::deployApp("shiny")
```

Account-specific `rsconnect` metadata is intentionally excluded from this repository.
The existing deployment is at
<https://g9v1qs-aisha-vayani.shinyapps.io/Cancer/>.

