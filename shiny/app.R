# Cancer outcomes explorer
#
# This app reads the project results table and provides interactive plots,
# and methods notes.
# The app looks for the disclosure-safe aggregate file in shiny/data/.

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(haven)
library(plotly)
library(shinyjs)

# Helper functions -------------------------------------------------------

read_shiny_data <- function() {
  data_path <- file.path("data", "Shiny.csv")
  if (!file.exists(data_path)) {
    stop("Could not find the aggregate file at shiny/data/Shiny.csv.")
  }
  df <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Convert labelled Stata variables before applying common string cleaning.
  as_chr <- function(x) {
    if (inherits(x, "labelled")) {
      return(as.character(haven::as_factor(x, levels = "labels")))
    }
    if (is.factor(x)) return(as.character(x))
    if (is.character(x)) return(x)
    as.character(x)
  }
  
  df <- df %>%
    mutate(
      site     = as_chr(site),
      estimate = as_chr(estimate),
      agegroup = as_chr(agegroup),
      stage    = as_chr(stage),
      sex      = as_chr(sex),
      race     = as_chr(race),
      tt       = as.numeric(tt),
      PE       = as.numeric(PE),
      lci      = as.numeric(lci),
      uci      = as.numeric(uci),
      n_black  = as.numeric(n_black)
    ) %>%
    mutate(
      site = str_squish(site),
      estimate = str_squish(estimate),
      agegroup = str_squish(agegroup),
      stage = str_squish(stage),
      sex = str_squish(sex),
      race = str_squish(race)
    ) %>%
    mutate(
      sex  = na_if(sex,  ""),
      race = na_if(race, "")
    )
  
  df <- df %>%
    mutate(
      race = case_when(
        str_detect(tolower(race), "^white") ~ "White",
        str_detect(tolower(race), "^black") ~ "Black",
        TRUE ~ race
      )
    )
  
  # Normalise estimate names so the plotting logic can use stable codes.
  df <- df %>%
    mutate(
      estimate_code = tolower(estimate),
      estimate_code = str_replace_all(estimate_code, "\\s+", " "),
      estimate_code = case_when(
        estimate_code %in% c("acs", "ac s", "standardised acs", "standardized acs") ~ "acs",
        estimate_code %in% c("model fit", "model_fit", "modelfit") ~ "model fit",
        estimate_code %in% c("pp", "p p") ~ "pp",
        estimate_code %in% c("ad") ~ "ad",
        estimate_code %in% c("ada", "ad a") ~ "ada",
        estimate_code %in% c("adb", "ad b") ~ "adb",
        estimate_code %in% c("pm", "proportion mediated", "proportion mediated by stage") ~ "pm",
        estimate_code %in% c("tce") ~ "tce",
        estimate_code %in% c("tde") ~ "tde",
        estimate_code %in% c("tie") ~ "tie",
        TRUE ~ estimate_code
      )
    )
  
  df
}

pretty_site <- function(x) {
  x <- str_replace_all(x, "_", " ")
  str_to_title(x)
}

# Values shown to users mapped to the estimate codes in the data.
estimate_category_map <- c(
  "Model fit" = "model fit",
  "Standardised all-cause survival" = "acs",
  "Mediation analysis" = "mediation",
  "Avoidable deaths" = "ad"
)

# Effect labels for mediation and avoidable-death plots.
mediation_effect_map <- c(
  "Total causal effect" = "tce",
  "Direct effects" = "tde",
  "Indirect effects (via stage)" = "tie",
  "All" = "all"
)

ad_effect_map <- c(
  "Total avoidable deaths" = "ad",
  "Avoidable deaths due to other all-cause survival differences" = "ada",
  "Avoidable deaths due to stage differences" = "adb",
  "All" = "all"
)

NIHR_COLOURS <- c(
  navy   = "#193E72",
  orange = "#F29330",
  coral  = "#EA5D4E",
  aqua   = "#2EA9B0",
  purple = "#6667AD",
  green  = "#46A86C",
  grey   = "#ACBCC3",
  pale_grey = "#EFF1F3"
)

series_colours <- function(series) {
  series <- unique(as.character(series))
  base <- c(
    "White" = NIHR_COLOURS[["navy"]],
    "Black" = NIHR_COLOURS[["orange"]],
    "Male" = NIHR_COLOURS[["navy"]],
    "Female" = NIHR_COLOURS[["orange"]]
  )
  
  fallback <- c(
    NIHR_COLOURS[["navy"]],
    NIHR_COLOURS[["orange"]],
    NIHR_COLOURS[["coral"]],
    NIHR_COLOURS[["aqua"]],
    NIHR_COLOURS[["purple"]],
    NIHR_COLOURS[["green"]]
  )
  
  vals <- base[series]
  missing_idx <- which(is.na(vals))
  if (length(missing_idx) > 0) {
    vals[missing_idx] <- fallback[((seq_along(missing_idx) - 1) %% length(fallback)) + 1]
  }
  
  setNames(unname(vals), series)
}

plot_y_axis_title <- function(cat, effect_code = NULL) {
  if (cat == "ad") return("Avoidable deaths")
  if (cat == "mediation" && identical(effect_code, "pm")) return("Proportion mediated by stage (%)")
  if (cat == "mediation") return("Diff in all-cause probabilities of death")
  return("Standardised all-cause survival")
}

y_axis_title <- function(cat) {
  plot_y_axis_title(cat)
}

nice_y_range <- function(values, include_zero = TRUE, pad = 0.08) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  
  if (!length(values)) return(c(0, 1))
  
  lo <- min(values)
  hi <- max(values)
  
  if (include_zero) {
    lo <- min(lo, 0)
    hi <- max(hi, 0)
  }
  
  if (identical(lo, hi)) {
    bump <- ifelse(abs(lo) < 1, 0.1, abs(lo) * 0.1)
    lo <- lo - bump
    hi <- hi + bump
  }
  
  span <- hi - lo
  c(lo - span * pad, hi + span * pad)
}

robust_y_range <- function(values, include_zero = TRUE, pad = 0.10) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  
  if (!length(values)) return(c(0, 1))
  
  full <- nice_y_range(values, include_zero = include_zero, pad = pad)
  
  if (length(values) < 10) return(full)
  
  core <- as.numeric(quantile(values, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE))
  core <- nice_y_range(core, include_zero = include_zero, pad = pad)
  
  if (diff(full) > 2.5 * diff(core)) return(core)
  full
}

standard_plot_theme <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = "Lato") +
    theme(
      text = element_text(family = "Lato"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = NIHR_COLOURS[["pale_grey"]], linewidth = 0.35),
      axis.line = element_line(colour = NIHR_COLOURS[["navy"]], linewidth = 0.45),
      axis.ticks = element_line(colour = NIHR_COLOURS[["navy"]], linewidth = 0.45),
      axis.title = element_text(colour = "#1f2933", size = base_size),
      axis.text = element_text(colour = "#1f2933", size = base_size - 1),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(size = base_size - 1, colour = "#1f2933"),
      legend.text = element_text(size = base_size - 2, colour = "#1f2933"),
      plot.title = element_text(face = "bold", colour = NIHR_COLOURS[["navy"]]),
      plot.subtitle = element_text(colour = "#455a64")
    )
}

build_subtitle <- function(compare_by, sexes, races, age, stage, cat, extra = NULL) {
  bits <- c()
  
  if (cat %in% c("acs", "model fit")) {
    if (compare_by == "Race/Ethnicity") {
      bits <- c(bits, "By Race", if (!is.null(sexes)) paste0(sexes, collapse = ", "))
      bits <- c(bits, paste0("Ages ", age))
      bits <- c(bits, paste0(str_to_title(stage), " stage"))
    } else {
      bits <- c(bits, "By Sex", if (!is.null(races)) paste0(races, collapse = ", "))
      bits <- c(bits, paste0("Ages ", age))
      bits <- c(bits, paste0(str_to_title(stage), " stage"))
    }
  } else {
    bits <- c(bits, paste0("Ages ", age))
  }
  
  if (!is.null(extra) && nzchar(extra)) bits <- c(bits, extra)
  
  paste(bits, collapse = ", ")
}

female_only_site_keys <- c(
  "breast",
  "ovary",
  "uterine corpus",
  "corpus uteri",
  "corpus uteri and uterus"
)

normalise_site_key <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_lower()
}

is_female_only_site <- function(x) {
  normalise_site_key(x %||% "") %in% female_only_site_keys
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

clean_plotly_legend <- function(fig) {
  if (is.null(fig$x$data) || !length(fig$x$data)) return(fig)
  
  extract_series <- function(name) {
    if (is.null(name) || is.na(name) || !nzchar(name)) return(name)
    s <- gsub("^\\(|\\)$", "", name)
    parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(name)
    
    if (tolower(parts[1]) %in% c("model fit", "pp estimator", "pp", "acs")) {
      if (length(parts) >= 2) return(parts[2])
    }
    
    parts[1]
  }
  
  for (i in seq_along(fig$x$data)) {
    tr <- fig$x$data[[i]]
    nm <- extract_series(tr$name)
    
    tr$name <- nm
    tr$legendgroup <- nm
    
    mode <- tr$mode %||% ""
    if (grepl("markers", mode, fixed = TRUE)) {
      tr$showlegend <- FALSE
    }
    
    fig$x$data[[i]] <- tr
  }
  
  fig
}

apply_plotly_style <- function(fig, axis_lwd = 1.15) {
  
  ax_col <- NIHR_COLOURS[["navy"]]
  grid_col <- NIHR_COLOURS[["pale_grey"]]
  
  fig %>% layout(
    legend = list(
      orientation = "h",
      x = 0, xanchor = "left",
      y = 1.18, yanchor = "top",
      bgcolor = "rgba(255,255,255,0.88)",
      bordercolor = grid_col,
      borderwidth = 1,
      font = list(size = 12, color = "#1f2933")
    ),
    margin = list(l = 72, r = 24, t = 62, b = 104),
    hoverlabel = list(
      align = "left",
      bgcolor = "white",
      bordercolor = grid_col,
      font = list(color = "#1f2933", size = 12)
    ),
    
    xaxis = list(
      showline = TRUE,
      mirror = FALSE,
      linecolor = ax_col,
      linewidth = axis_lwd,
      ticks = "outside",
      tickwidth = axis_lwd,
      tickcolor = ax_col,
      showgrid = FALSE,
      zeroline = FALSE,
      automargin = TRUE,
      title = list(text = "Time since diagnosis (years)", standoff = 8)
    ),
    yaxis = list(
      showline = TRUE,
      mirror = FALSE,
      linecolor = ax_col,
      linewidth = axis_lwd,
      ticks = "outside",
      tickwidth = axis_lwd,
      tickcolor = ax_col,
      showgrid = TRUE,
      gridcolor = grid_col,
      gridwidth = 1,
      zeroline = FALSE,
      automargin = TRUE
    ),
    
    plot_bgcolor = "white",
    paper_bgcolor = "white"
  )
}

# Footnote text shown in the app.
FOOTNOTE_LINES <- c(
  "Model fitted",
  "\u2022 Flexible parametric relative survival model fitted on the log cumulative excess hazard scale (stratified by sex) with a period window of 01/01/2021 to 31/12/2022",
  "\u2022 Baseline cumulative excess hazard modelled using restricted cubic splines with 3 degrees of freedom.",
  "\u2022 Main effects included: race, stage, age, race \u00d7 stage, stage \u00d7 age, race \u00d7 age.",
  "\u2022 Age modelled non-linearly using restricted cubic splines with 3 degrees of freedom.",
  "\u2022 Age was modelled continuously for the central 95% of the age distribution, while patients younger than 2.5th and older than 97.5th age percentile were assigned the same expected survival as those at these respective cut-off ages.",
  "\u2022 Time-dependent effects were included for race, stage, age, race \u00d7 stage and race \u00d7 age, each with 2 degrees of freedom.",
  "\u2022 Expected survival was incorporated using SEER population life tables stratified by race, year, age and sex.",
  "",
  "Data sources",
  "\u2022 Incidence data: SEER Program (seer.cancer.gov) SEER*Stat Database: Incidence - SEER Research Limited-Field Data, 21 Registries (excl IL), Nov 2024 Sub (2000\u20132022) - Linked To County Attributes - Time Dependent (1990-2023) Income/Rurality, 1969-2023 Counties, released April 2025, based on the Nov 2024 submission.",
  "\u2022 Expected survival (life tables): SEER Program (seer.cancer.gov) SEER*Stat Database: Expected Survival - U.S. 1970-2021 by individual year (White, Black, Other (AI/API), Ages 0-99, etc.).",
  "",
  "Suggested citation",
  "\u2022 [To be updated when paper is published]"
)

FOOTNOTE_TEXT <- paste(FOOTNOTE_LINES, collapse = "\n")

# Small browser helper for copying bookmarked app state.
COPY_JS <- "
shinyjs.copyToClipboard = function(text) {
  if (navigator.clipboard && window.isSecureContext) {
    navigator.clipboard.writeText(text);
  } else {
    const textArea = document.createElement('textarea');
    textArea.value = text;
    textArea.style.position = 'fixed';
    textArea.style.left = '-9999px';
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    document.execCommand('copy');
    document.body.removeChild(textArea);
  }
};
"

enableBookmarking(store = "url")

# User interface ---------------------------------------------------------

ui <- page_sidebar(
  
  title = tagList(
    div(
      class = "app-title-wrap",
      div("Stage-at-diagnosis and racial disparities in cancer survival", class = "app-title"),
      div("Interactive companion to the population-based SEER mediation analysis", class = "app-subtitle")
    )
  ),
  fillable = TRUE,
  sidebar = sidebar(
    width = 360,
    # The controls are rebuilt when the selected estimate type changes.
    uiOutput("sidebar_controls"),
    hr(),
    uiOutput("sidebar_notes")
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = NIHR_COLOURS[["navy"]],
    secondary = NIHR_COLOURS[["orange"]]
  ),
  
  tags$head(
    tags$style(HTML("
  html, body {
    height: 100%;
    overflow: hidden;
  }
  
  body {
    font-family: Lato, Arial, sans-serif;
    background: #f7f9fb;
  }
  
  .navbar, .navbar.navbar-default, .bslib-page-title {
    background: #193E72 !important;
    color: #fff !important;
  }
  
  .navbar .navbar-brand,
  .navbar .navbar-brand *,
  .bslib-page-title,
  .bslib-page-title * {
    color: #fff !important;
  }
  
  .app-title-wrap {
    display: flex;
    flex-direction: column;
    gap: 2px;
    color: #fff;
  }
  
  .app-title {
    font-weight: 800;
    font-size: 1.05rem;
    line-height: 1.1;
  }
  
  .app-subtitle {
    font-size: 0.78rem;
    color: rgba(255,255,255,0.82);
    line-height: 1.1;
  }
  
  .selectize-dropdown, .selectize-dropdown.form-control {
    z-index: 3000 !important;
  }
  
  .bslib-sidebar-layout > .sidebar {
    max-height: calc(100vh - 72px);
    overflow-y: auto;
  }
  
  .bslib-sidebar-layout > .main {
    min-height: 0;
    overflow: hidden;
  }
  
  .bslib-sidebar-layout [aria-label='Resize sidebar'] {
    display: none !important;
    pointer-events: none !important;
  }
  
  .card-header { padding-top: .5rem; padding-bottom: .5rem; }
  
/* Compact plot header */
.card-header { padding: .35rem .75rem !important; }

.card-header h4 {
  font-size: 1.15rem !important;
  line-height: 1.1 !important;
  margin-bottom: .15rem !important;
  color: #193E72 !important;
  font-weight: 800 !important;
}

.card-header .plot-subtitle {
  font-size: 0.9rem !important;
  line-height: 1.1 !important;
  margin-top: 0 !important;
  white-space: nowrap !important;
  overflow: hidden !important;
  text-overflow: ellipsis !important;
}

.card-header .btn {
  padding: .25rem .45rem !important;
}

  
  /* Fill the available card height with the selected tab. */
.bslib-navs {
  display: flex;
  flex-direction: column;
  flex: 1 1 auto;
  min-height: 0;
}

.bslib-navs .tab-content {
  flex: 1 1 auto;
  min-height: 0;
  overflow: hidden !important;
}

.bslib-navs .tab-pane {
  height: 100%;
  min-height: 0;
}

  




  #plot {
    height: clamp(330px, calc(100vh - 330px), 620px) !important;
  }
  
  #plot .plotly, #plot .svg-container {
    height: 100% !important;
  }
  
  .sidebar-notes {
    font-size: 12px;
    line-height: 1.35;
    color: #344054;
  }
  
  .sidebar-notes ul {
    margin: 6px 0 0 1rem;
    padding: 0;
  }
  
  .sidebar-notes li {
    margin-bottom: 5px;
  }
"))
  ),
  
    div(
      style = "display:flex; flex-direction:column; gap:10px; height:calc(100vh - 96px); min-height:0; overflow:hidden;",
      
      # Fixed-height selection row.
      card(
        fill = FALSE,
        style = "padding: 2px; flex: 0 0 auto;",
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            selectizeInput(
              "site", "Cancer site",
              choices = NULL,
              selected = NULL,
              options = list(
                dropdownParent = "body"
              )
            ),
            
            selectizeInput(
              "estimate_category", "Estimate",
              choices = names(estimate_category_map),
              selected = "Model fit",
              options = list(
                dropdownParent = "body"
              )
            )
          )
        )
      ),
      
      # Main card fills the remaining page height.
      card(
        full_screen = TRUE,
        fill = TRUE,
        style = "flex: 1 1 auto; min-height: 0;",
        card_header(
          div(
            style = "display:flex; align-items:flex-start; justify-content:space-between; gap:12px;",
            div(
              tags$h4(textOutput("plot_title"), style = "margin-bottom: 2px;"),
              tags$div(
                textOutput("plot_subtitle"),
                class = "plot-subtitle",
                style = "color:#555; margin-top:0;"
              )
              
            ),
            div(
              style = "display:flex; gap:10px; align-items:center;",
              tooltip(
                actionButton("btn_copy_link", label = NULL, icon = icon("link"),
                             class = "btn btn-outline-secondary"),
                "Copy link", placement = "bottom"
              )
            )
          )
        ),
        card_body(
          style = "min-height: 0; display:flex; flex-direction:column; overflow:hidden;",
          
          div(
            style = "flex: 1 1 auto; min-height: 0; display:flex; flex-direction:column;",
            navset_tab(
              nav_panel(
                "Graph",
                div(
                  style = "height: 100%; min-height: 0; display:flex; flex-direction:column; gap:6px;",
                  div(
                    style = "flex: 0 0 auto; display:flex; justify-content:flex-end; align-items:center; gap:8px; margin-bottom:-2px;",
                    icon("mouse-pointer"),
                    tags$span("Tap/hover on points for more details", style = "color:#555; font-size: 12px;")
                  ),
                  div(
                    style = "flex: 1 1 auto; min-height: 0; overflow:hidden;",
                    plotlyOutput("plot", height = "100%")
                  )
                )
              ),
              nav_panel(
                "Footnote / Methods & data sources",
                tags$pre(FOOTNOTE_TEXT, style = "font-size: 12px; white-space: pre-wrap;")
              )
            )
          )
          
        )
      )
    ),
  useShinyjs(),
  extendShinyjs(text = COPY_JS, functions = "copyToClipboard")
)

# Server ----------------------------------------------------------------

server <- function(input, output, session) {
  
  # Load data once, then update the site selector from the file contents.
  df <- reactiveVal(NULL)
  
  observeEvent(TRUE, {
    dat <- read_shiny_data()
    df(dat)
    
    observeEvent(df(), {
      message("Unique race labels: ", paste(sort(unique(df()$race)), collapse = " | "))
    })
    
    sites <- sort(unique(dat$site))
    default_site <- sites[1]
    lung_idx <- which(str_detect(tolower(sites), "lung"))
    if (length(lung_idx) > 0) default_site <- sites[lung_idx[1]]
    
    updateSelectizeInput(session, "site", choices = sites, selected = default_site, server = TRUE)
    
  }, once = TRUE)
  
  # Race and sex comparisons are only meaningful for survival estimates.
  compare_allowed <- reactive({
    cat_code <- estimate_category_map[[input$estimate_category]]
    cat_code %in% c("acs", "model fit")
  })
  
  observeEvent(list(input$site, input$estimate_category), {
    req(input$site, input$estimate_category)
    
    cat_code <- estimate_category_map[[input$estimate_category]]
    is_compare <- cat_code %in% c("acs", "model fit")
    is_female_only <- is_female_only_site(input$site)
    
    if (is_compare && is_female_only) {
      updateRadioButtons(session, "compare_by", selected = "Race/Ethnicity")
    }
  }, ignoreInit = TRUE)
  
  # Default to race comparisons for ACS and model-fit views.
  observeEvent(input$estimate_category, {
    if (compare_allowed() && is.null(input$compare_by)) {
      updateRadioButtons(session, "compare_by", selected = "Race/Ethnicity")
    }
  }, ignoreInit = TRUE)
  
  
  # Rebuild sidebar inputs from the data available for the selected view.
  output$sidebar_controls <- renderUI({
    req(df())
    dat <- df()
    
    cat_code <- estimate_category_map[[input$estimate_category]]
    
    sex_choices  <- sort(unique(dat$sex))
    race_choices <- sort(unique(dat$race))
    
    sex_choices  <- sex_choices[!is.na(sex_choices) & nzchar(sex_choices)]
    race_choices <- race_choices[!is.na(race_choices) & nzchar(race_choices)]
    
    age_choices <- sort(unique(na.omit(dat$agegroup)))
    stage_levels <- c("Localised", "Regional", "Distant")
    
    stage_choices <- stage_levels[stage_levels %in% unique(na.omit(dat$stage))]

    
    default_age <- age_choices[min(3, length(age_choices))]
    idx_5564 <- which(str_detect(age_choices, "55"))
    if (length(idx_5564) > 0) default_age <- age_choices[idx_5564[1]]
    
    preferred_default <- "Regional"
    default_stage <- if (preferred_default %in% stage_choices) preferred_default else stage_choices[1]
    
    default_sex <- "Male"
    if (!("Male" %in% sex_choices) && length(sex_choices) > 0) default_sex <- sex_choices[1]
    
    # Sex-specific cancers only have female estimates.
    site_now <- input$site %||% ""
    if (is_female_only_site(site_now)) {
      sex_choices <- intersect(sex_choices, c("Female"))
      default_sex <- "Female"
    }
    
    compare_ui <- NULL
    site_now <- input$site %||% ""
    is_female_only <- is_female_only_site(site_now)
    
    if (compare_allowed()) {
      
      if (is_female_only) {
        compare_ui <- div(
          tags$div(tags$strong("Compare by:"), style = "margin-bottom:6px;"),
          tags$div("Race/Ethnicity", style = "font-weight: 600;"),
          tags$div("Sex comparison is not available for sex-specific cancers.", style = "color:#666; font-size:12px;"),
          tags$hr()
        )
      } else {
        selected_compare <- if (!is.null(input$compare_by) &&
                                input$compare_by %in% c("Sex", "Race/Ethnicity")) {
          input$compare_by
          
        } else {
          "Race/Ethnicity"
        }
        
        compare_ui <- div(
          tags$div(tags$strong("Compare by:"), style = "margin-bottom:6px;"),
          radioButtons(
            "compare_by", NULL,
            choices = c("Sex", "Race/Ethnicity"),
            selected = selected_compare,
            inline = TRUE
          ),
          tags$hr()
        )
      }
        
    } else {
      compare_ui <- div(
        tags$div(tags$strong("Compare by:"), style = "margin-bottom:6px;"),
        tags$div("Not applicable for this estimate type.", style = "color:#666; font-size: 12px;"),
        tags$hr()
      )
    }
      
    is_female_only <- is_female_only_site(input$site)
    cat_code <- estimate_category_map[[input$estimate_category]]
    
    compare_by <- if (cat_code %in% c("acs", "model fit") && is_female_only) {
      "Race/Ethnicity"
    } else {
      input$compare_by %||% "Race/Ethnicity"
    }
    
    
    # Mediation and avoidable-death estimates are compared by sex only.
    compare_controls <- NULL
    
    if (cat_code %in% c("acs", "model fit")) {
      
      if (compare_by == "Sex") {
        compare_controls <- tagList(
          checkboxGroupInput(
            "sex_multi", "Sex",
            choices = sex_choices,
            selected = intersect(c("Male", "Female"), sex_choices)
          ),
          selectInput(
            "race_single", "Race/Ethnicity",
            choices = race_choices,
            selected = if ("White" %in% race_choices) "White" else race_choices[1]
          )
        )
      } else {
        compare_controls <- tagList(
          checkboxGroupInput(
            "race_multi", "Race/Ethnicity",
            choices = race_choices,
            selected = intersect(c("White", "Black"), race_choices)
          ),
          selectInput(
            "sex_single", "Sex",
            choices = sex_choices,
            selected = default_sex
          )
        )
      }
      
      stage_ui <- selectInput(
        "stage_single", "Stage at diagnosis",
        choices = stage_choices,
        selected = default_stage
      )
      
      age_ui <- selectInput(
        "age_single", "Age",
        choices = age_choices,
        selected = default_age
      )
      
      more_ui <- tagList(
        tags$hr(),
        tags$strong("More options"),
        if (cat_code == "model fit") {
          checkboxInput("show_pp", "Show Pohar-Perme estimator", value = FALSE)
        },
        checkboxInput("show_ci", "Show confidence intervals", value = TRUE)
      )
      
      tagList(compare_ui, compare_controls, age_ui, stage_ui, more_ui)
      
    } else if (cat_code == "mediation") {
      
      age_ui <- selectInput("age_single", "Age", choices = age_choices, selected = default_age)
      
      effects_ui <- selectInput(
        "med_effect", "Causal effects",
        choices = names(mediation_effect_map),
        selected = "Total causal effect"
      )
      
      if (is_female_only) {
        compare_ui2 <- div(
          tags$div(tags$strong("Compare by:"), style="margin-bottom:6px;"),
          tags$div("Sex", style="font-weight:600;"),
          tags$div("Sex comparison is not available for sex-specific cancers.", style="color:#666; font-size:12px;"),
          tags$hr()
        )
        sex_ui <- NULL
      } else {
        compare_ui2 <- div(
          tags$div(tags$strong("Compare by:"), style="margin-bottom:6px;"),
          tags$div("Sex", style="font-weight:600;"),
          tags$hr()
        )
        sex_ui <- checkboxGroupInput(
          "sex_med_multi", "Sex",
          choices = sex_choices,
          selected = intersect(c("Male","Female"), sex_choices)
        )
      }
      
      more_ui <- tagList(
        tags$hr(),
        tags$strong("More options"),
        checkboxInput("show_ci", "Show confidence intervals", value = TRUE)
      )
      
      tagList(compare_ui2, sex_ui, age_ui, effects_ui, more_ui)
      
    } else if (cat_code == "ad") {
      
      age_ui <- selectInput("age_single", "Age", choices = age_choices, selected = default_age)
      
      effects_ui <- selectInput(
        "ad_effect", "Avoidable deaths",
        choices = names(ad_effect_map),
        selected = "Total avoidable deaths"
      )
      
      if (is_female_only) {
        compare_ui2 <- div(
          tags$div(tags$strong("Compare by:"), style="margin-bottom:6px;"),
          tags$div("Sex", style="font-weight:600;"),
          tags$div("Sex comparison is not available for sex-specific cancers.", style="color:#666; font-size:12px;"),
          tags$hr()
        )
        sex_ui <- NULL
      } else {
        compare_ui2 <- div(
          tags$div(tags$strong("Compare by:"), style="margin-bottom:6px;"),
          tags$div("Sex", style="font-weight:600;"),
          tags$hr()
        )
        sex_ui <- checkboxGroupInput(
          "sex_ad_multi", "Sex",
          choices = sex_choices,
          selected = intersect(c("Male","Female"), sex_choices)
        )
      }
      
      more_ui <- tagList(
        tags$hr(),
        tags$strong("More options"),
        checkboxInput("show_ci", "Show confidence intervals", value = TRUE)
      )
      
      tagList(compare_ui2, sex_ui, age_ui, effects_ui, more_ui)
      
    } else {
      tagList(compare_ui)
    }
  })
  
  output$sidebar_notes <- renderUI({
    req(input$estimate_category)
    
    cat_code <- estimate_category_map[[input$estimate_category]]
    notes <- character(0)
    
    if (isTRUE(input$show_ci)) {
      notes <- c(notes, "Shaded bands show 95% confidence intervals for the selected estimate.")
    }
    
    if (cat_code == "model fit") {
      if (isTRUE(input$show_pp)) {
        notes <- c(notes, "Diamond markers show Pohar-Perme estimates at their observed values; vertical bars show their 95% confidence intervals.")
      } else {
        notes <- c(notes, "Pohar-Perme estimates are hidden for this view.")
      }
    }
    
    if (cat_code == "mediation") {
      eff <- mediation_effect_map[[input$med_effect %||% "Total causal effect"]]
      if (eff %in% c("tie", "all")) {
        notes <- c(notes, "Hover over indirect-effect points to see the proportion mediated by stage.")
      }
      notes <- c(notes, "Mediation estimates can be negative; the y-axis is scaled around zero when needed.")
    }
    
    if (cat_code == "ad") {
      notes <- c(notes, "N* is the number of Black diagnoses in 2022 for the selected age group, inclusive of those with missing stage at diagnosis.")
    }
    
    if (!length(notes)) {
      notes <- "Use the controls above to update the plotted site, group, stage and age."
    }
    
    div(
      class = "sidebar-notes",
      tags$strong("Notes for this view"),
      tags$ul(lapply(notes, tags$li))
    )
  })
  
  # Filter the results table to the current sidebar choices.
  filtered_data <- reactive({
    req(df(), input$site, input$estimate_category, input$age_single)
    
    dat <- df()
    cat_code <- estimate_category_map[[input$estimate_category]]
    
    est_codes <- character(0)
    
    if (cat_code %in% c("acs", "model fit")) {
      est_codes <- cat_code
      if (cat_code == "model fit" && isTRUE(input$show_pp)) {
        est_codes <- c(est_codes, "pp")
      }
    } else if (cat_code == "mediation") {
      eff <- mediation_effect_map[[input$med_effect %||% "Total causal effect"]]
      est_codes <- if (eff == "all") c("tce", "tde", "tie") else eff
    } else if (cat_code == "ad") {
      eff <- ad_effect_map[[input$ad_effect %||% "Total avoidable deaths"]]
      est_codes <- if (eff == "all") c("ad", "ada", "adb") else eff
    }
    
    out <- dat %>%
      filter(site == input$site) %>%
      filter(estimate_code %in% est_codes) %>%
      filter(agegroup == input$age_single)
    
    # Survival estimates can be shown by race or by sex.
    if (cat_code %in% c("acs", "model fit")) {
      req(input$stage_single)
      out <- out %>% filter(stage == input$stage_single)
      
      is_female_only <- is_female_only_site(input$site)
      compare_by <- if (is_female_only) "Race/Ethnicity" else (input$compare_by %||% "Race/Ethnicity")
      
      if (compare_by == "Sex") {
        req(input$sex_multi, input$race_single)
        out <- out %>%
          filter(race == input$race_single) %>%
          filter(sex %in% input$sex_multi) %>%
          mutate(series = sex)
      } else {
        req(input$race_multi, input$sex_single)
        out <- out %>%
          filter(sex == input$sex_single) %>%
          filter(race %in% input$race_multi) %>%
          mutate(series = race)
      }
    } else {
      out <- out %>% mutate(series = estimate)
    }
    
    # Mediation and avoidable-death rows can otherwise duplicate by race/stage.
    if (cat_code %in% c("mediation", "ad")) {
      is_female_only <- is_female_only_site(input$site)
      
      if (is_female_only) {
        out <- out %>% filter(sex == "Female")
      } else {
        if (cat_code == "mediation") {
          req(input$sex_med_multi)
          out <- out %>% filter(sex %in% input$sex_med_multi)
        } else {
          req(input$sex_ad_multi)
          out <- out %>% filter(sex %in% input$sex_ad_multi)
        }
      }
      
      out <- out %>%
        group_by(site, agegroup, sex, estimate_code, tt) %>%
        summarise(
          PE      = first(PE),
          lci     = first(lci),
          uci     = first(uci),
          n_black = safe_max(n_black),
          .groups = "drop"
        )
      
      if (cat_code == "mediation") {
        eff <- mediation_effect_map[[input$med_effect %||% "Total causal effect"]]
      } else {
        eff <- ad_effect_map[[input$ad_effect %||% "Total avoidable deaths"]]
      }
      
      out <- out %>%
        mutate(
          stage = NA_character_,
          race  = NA_character_,
          series = if (eff == "all") paste0(sex, " - ", toupper(estimate_code)) else sex
        )
    }
    
    # Keep estimates even if a confidence interval is unavailable.
    out <- out %>% filter(!is.na(PE))
    
    out <- out %>% filter(!is.na(tt))
    
    out %>% arrange(series, tt)
  })
  
  # Point-level tooltips for the integer-year markers.
  tooltip_points <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    cat_code <- estimate_category_map[[input$estimate_category]]
    show_ci  <- isTRUE(input$show_ci)
    
    ints <- c(0, 1, 2, 3, 4, 5)
    
    pts <- dat %>%
      filter(tt %in% ints | (estimate_code == "pp" & tt %in% 1:5)) %>%
      mutate(
        compare_label = series,
        point_type = ifelse(estimate_code == "pp", "Pohar-Perme estimator", "Model fit")
      )
    
    # Attach sex-specific N* for avoidable deaths.
    if (cat_code == "ad") {
      
      nstar_by_sex <- dat %>%
        group_by(sex) %>%
        summarise(N_star = suppressWarnings(max(n_black, na.rm = TRUE)), .groups = "drop") %>%
        mutate(N_star = ifelse(is.finite(N_star), N_star, NA_real_))
      
      pts <- pts %>%
        left_join(nstar_by_sex, by = "sex")
      
    } else {
      pts <- pts %>% mutate(N_star = NA_real_)
    }
    
    # Attach saved proportion-mediated estimates to indirect-effect tooltips.
    if (cat_code == "mediation") {
      pm_lookup <- df() %>%
        filter(
          site == input$site,
          agegroup == input$age_single,
          estimate_code == "pm",
          tt %in% ints
        )
      
      if (is_female_only_site(input$site)) {
        pm_lookup <- pm_lookup %>% filter(sex == "Female")
      } else {
        req(input$sex_med_multi)
        pm_lookup <- pm_lookup %>% filter(sex %in% input$sex_med_multi)
      }
      
      pm_lookup <- pm_lookup %>%
        group_by(sex, tt) %>%
        summarise(
          pm_PE  = first(PE),
          pm_lci = first(lci),
          pm_uci = first(uci),
          .groups = "drop"
        )
      
      pts <- pts %>% left_join(pm_lookup, by = c("sex", "tt"))
    } else {
      pts <- pts %>% mutate(pm_PE = NA_real_, pm_lci = NA_real_, pm_uci = NA_real_)
    }
    
    fmt_num <- function(x, digits = 2) {
      ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
    }
    
    # Build the displayed hover text from the plotted row.
    pts <- pts %>%
      rowwise() %>%
      mutate(
        tooltip = {
          time_label <- formatC(tt, format = "f", digits = 0)
          pe_label   <- fmt_num(PE, 2)
          lci_label  <- fmt_num(lci, 2)
          uci_label  <- fmt_num(uci, 2)
          estimate_label <- case_when(
            estimate_code == "pp" ~ "Pohar-Perme estimator",
            estimate_code == "model fit" ~ "Model estimate",
            estimate_code == "acs" ~ "Standardised ACS",
            estimate_code == "pm" ~ "Proportion mediated",
            TRUE ~ toupper(estimate_code)
          )
          
          lines <- c(
            sprintf("<b>%s</b>", compare_label),
            sprintf("<b>Estimate type:</b> %s", estimate_label),
            sprintf("<b>Time since diagnosis:</b> %s years", time_label),
            sprintf("<b>Estimate:</b> %s", pe_label)
          )
          
          
          show_row_ci <- show_ci || (cat_code == "model fit" && estimate_code == "pp")
          
          if (show_row_ci) {
            lines <- c(lines, sprintf("<b>95%% CI:</b> %s to %s", lci_label, uci_label))
          }
          
          if (cat_code == "mediation" && estimate_code == "tie" && is.finite(pm_PE)) {
            pm_label <- sprintf("<b>Proportion mediated by stage:</b> %s%%", fmt_num(pm_PE, 1))
            if (show_ci && is.finite(pm_lci) && is.finite(pm_uci)) {
              pm_label <- sprintf(
                "%s (95%% CI %s to %s)",
                pm_label,
                fmt_num(pm_lci, 1),
                fmt_num(pm_uci, 1)
              )
            }
            lines <- c(lines, pm_label)
          }
          
          if (cat_code == "ad" && is.finite(N_star)) {
            lines <- c(lines, sprintf("<b>N*</b> (Black diagnoses, 2022): %s", formatC(N_star, format = "f", digits = 0)))
          }
          
          paste(lines, collapse = "<br>")
        }
      ) %>%
      ungroup()
    
    pts
  })
  
  
  # Plot title/subtitle
  output$plot_title <- renderText({
    req(input$site, input$estimate_category)
    pretty_site(input$site)
  })
  
  plot_subtitle_text <- reactive({
    req(input$estimate_category, input$age_single)
    cat_code <- estimate_category_map[[input$estimate_category]]
    
    if (cat_code %in% c("acs", "model fit")) {
      compare_by <- input$compare_by %||% "Race/Ethnicity"
      stage <- input$stage_single %||% ""
      age <- input$age_single %||% ""
      sexes <- if (compare_by == "Race/Ethnicity") (input$sex_single %||% "") else (input$sex_multi %||% character(0))
      races <- if (compare_by == "Race/Ethnicity") (input$race_multi %||% character(0)) else (input$race_single %||% "")
      extra <- if (cat_code == "model fit" && isTRUE(input$show_pp)) "Includes Pohar-Perme estimator" else ""
      
      paste0(
        input$estimate_category, " by time since diagnosis\n",
        build_subtitle(compare_by, sexes = as.character(sexes), races = as.character(races),
                       age = age, stage = stage, cat = cat_code, extra = extra)
      )
    } else if (cat_code == "mediation") {
      eff <- input$med_effect %||% "Total causal effect"
      paste0("Mediation analysis (", eff, ") by time since diagnosis\nBy Ages ", input$age_single)
    } else if (cat_code == "ad") {
      eff <- input$ad_effect %||% "Total avoidable deaths"
      paste0("Avoidable deaths (", eff, ") by time since diagnosis\nBy Ages ", input$age_single)
    } else ""
  })
  
  output$plot_subtitle <- renderText(plot_subtitle_text())
  
  # Interactive plot.
  output$plot <- renderPlotly({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    cat_code <- estimate_category_map[[input$estimate_category]]
    show_ci <- isTRUE(input$show_ci)
    
    if (cat_code %in% c("acs", "model fit")) {
      
      is_female_only <- is_female_only_site(input$site)
      compare_by <- if (is_female_only) "Race/Ethnicity" else (input$compare_by %||% "Race/Ethnicity")
      
      # PP rows are plotted as separate points and error bars, not as a line.
      dat_fit <- dat %>% filter(estimate_code != "pp")
      dat_pp  <- dat %>% filter(estimate_code == "pp" & tt %in% 1:5)
      
      pts_all <- tooltip_points()
      pts_fit <- pts_all %>% filter(estimate_code != "pp")
      pts_pp  <- pts_all %>% filter(estimate_code == "pp")
      cols <- series_colours(dat$series)
      
      p <- ggplot(dat_fit, aes(x = tt, y = PE, colour = series, group = series)) +
        geom_line(linewidth = 1.05, alpha = 0.95)
      
      if (show_ci) {
        p <- p + geom_ribbon(
          data = dat_fit,
          aes(ymin = lci, ymax = uci, fill = series, group = series),
          alpha = 0.18, colour = NA, show.legend = FALSE
        )
      }
      
      p <- p + geom_point(
          data = pts_fit,
          aes(x = tt, y = PE, text = tooltip, shape = point_type),
          size = 2.8,
          show.legend = FALSE
        )
      
      if (nrow(dat_pp) > 0) {
        p <- p + geom_errorbar(
          data = dat_pp,
          aes(x = tt, ymin = lci, ymax = uci, colour = series),
          inherit.aes = FALSE,
          width = 0.08, linewidth = 0.8, alpha = 0.95,
          show.legend = FALSE
        )
      }
      
      if (nrow(pts_pp) > 0) {
        p <- p + geom_point(
          data = pts_pp,
          aes(x = tt, y = PE, colour = series, text = tooltip, shape = point_type),
          inherit.aes = FALSE,
          size = 3.4,
          show.legend = FALSE
        )
      }
      
      p <- p + scale_shape_manual(
        values = c("Model fit" = 16, "Pohar-Perme estimator" = 18),
        guide = "none"
      )
      
      p <- p +
        scale_colour_manual(values = cols, drop = FALSE) +
        scale_fill_manual(values = cols, drop = FALSE) +
        scale_x_continuous(limits = c(0, 5), breaks = 0:5, expand = c(0, 0)) +
        scale_y_continuous(
          limits = c(0, 1.03),
          breaks = seq(0, 1, by = 0.1),
          labels = label_number(accuracy = 0.1), expand = c(0, 0)
        ) +
        labs(
          x = "Time since diagnosis (years)",
          y = y_axis_title(cat_code),
          colour = ifelse(compare_by == "Race/Ethnicity", "Race/Ethnicity", "Sex")
        ) +
        standard_plot_theme(base_size = 13)
      
      fig <- suppressWarnings(ggplotly(p, tooltip = "text"))
      fig <- clean_plotly_legend(fig)
      fig <- apply_plotly_style(fig, axis_lwd = 0.8)
      
      fig <- fig %>% layout(
        hovermode = "closest",  
        xaxis = list(showspikes = FALSE, title = list(text = "Time since diagnosis (years)", standoff = 8)),
        yaxis = list(showspikes = FALSE)
      )
      
      fig <- fig %>% config(
        displayModeBar = FALSE,
        displaylogo = FALSE
      )
      
      fig
      
    } else if (cat_code %in% c("mediation", "ad")) {
      
      dat <- dat %>% arrange(series, tt)
      pts <- tooltip_points()
      effect_code <- if (cat_code == "mediation") {
        mediation_effect_map[[input$med_effect %||% "Total causal effect"]]
      } else {
        NULL
      }
      
      if (identical(effect_code, "pm")) {
        dat <- dat %>% filter(tt %in% 0:5)
        pts <- pts %>% filter(tt %in% 0:5)
      }
      
      cols <- series_colours(dat$series)
      
      y_values <- if (show_ci) c(dat$lci, dat$uci, dat$PE) else dat$PE
      y_range <- if (identical(effect_code, "pm")) {
        robust_y_range(y_values, include_zero = TRUE)
      } else {
        nice_y_range(y_values, include_zero = TRUE)
      }
      if (cat_code == "ad") y_range[1] <- min(0, y_range[1])
      oob_fun <- if (identical(effect_code, "pm")) scales::squish else scales::censor
      
      p <- ggplot(dat, aes(x = tt, y = PE, colour = series, group = series)) +
        geom_hline(yintercept = 0, colour = NIHR_COLOURS[["grey"]], linewidth = 0.35) +
        geom_line(linewidth = 1.05, alpha = 0.95)
      
      if (show_ci) {
        p <- p + geom_ribbon(
          aes(ymin = lci, ymax = uci, fill = series),
          alpha = 0.18, colour = NA, show.legend = FALSE
        )
      }
      
      p <- p +
        geom_point(
          data = pts,
          aes(x = tt, y = PE, colour = series, text = tooltip),
          size = 2.8,
          show.legend = FALSE
        ) +
        scale_colour_manual(values = cols, drop = FALSE) +
        scale_fill_manual(values = cols, drop = FALSE) +
        scale_x_continuous(limits = c(0, 5), breaks = 0:5, expand = c(0, 0)) +
        scale_y_continuous(
          limits = y_range,
          labels = label_number(accuracy = ifelse(identical(effect_code, "pm"), 1, 0.01)),
          expand = c(0, 0),
          oob = oob_fun
        ) +
        labs(
          x = "Time since diagnosis (years)",
          y = plot_y_axis_title(cat_code, effect_code),
          colour = if (cat_code == "ad") "Avoidable deaths" else "Causal effect"
        ) +
        standard_plot_theme(base_size = 13)
      
      fig <- suppressWarnings(ggplotly(p, tooltip = "text"))
      fig <- clean_plotly_legend(fig)
      fig <- apply_plotly_style(fig, axis_lwd = 0.8)
      
      fig <- fig %>% layout(
        hovermode = "closest",        
        xaxis = list(showspikes = FALSE, title = list(text = "Time since diagnosis (years)", standoff = 8)),
        yaxis = list(showspikes = FALSE, range = y_range)
      )
      
      fig <- fig %>% config(
        displayModeBar = FALSE,
        displaylogo = FALSE
      )
      
      fig
      
    }
  })
  
  # Copy bookmarked app state to the clipboard.
  observeEvent(input$btn_copy_link, {
    session$doBookmark()
    runjs("shinyjs.copyToClipboard(window.location.href);")
    showNotification("Link copied to clipboard.", type = "message", duration = 2)
  })
  
}

shinyApp(ui, server)
