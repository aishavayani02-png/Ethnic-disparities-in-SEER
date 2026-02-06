# app.R
# ------------------------------------------------------------
# Cancer outcomes explorer (Relative survival / Mediation / Avoidable deaths)
# Data: long format with columns: site, sex, agegroup, estimate, tt, PE, lci, uci, n_black
# ------------------------------------------------------------

library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(plotly)
library(DT)
library(shinyjs)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

# ---- GitHub raw CSV URL ----
DATA_URL <- "https://raw.githubusercontent.com/aishavayani02-png/Ethnic-disparities-in-SEER/main/data/Shiny.csv"
# If your file is in the repo ROOT (not /data), use:
# DATA_URL <- "https://raw.githubusercontent.com/aishavayani02-png/Ethnic-disparities-in-SEER/main/Shiny.csv"

FOOTNOTE_LINES <- c(
  "Model fitted",
  "• Flexible parametric relative survival model fitted on the log cumulative excess hazard scale (stratified by sex) with a period window of 01/01/2021 to 31/12/2022",
  "• Baseline cumulative excess hazard modelled using restricted cubic splines with 3 degrees of freedom.",
  "• Main effects included: race, stage, age, race × stage, stage × age, race × age.",
  "• Age modelled non-linearly using restricted cubic splines with 3 degrees of freedom.",
  "• Age was modelled continuously for the central 95% of the age distribution, while patients younger than 2.5th and older than 97.5th age percentile were assigned the same expected survival as those at these respective cut-off ages.",
  "• Time-dependent effects were included for race, stage, age, race × stage and race × age, each with 2 degrees of freedom.",
  "• Expected survival was incorporated using SEER population life tables stratified by race, year, age and sex.",
  "",
  "Data sources",
  "• Incidence data: Surveillance, Epidemiology, and End Results (SEER) Program (seer.cancer.gov) SEER*Stat Database: Incidence - SEER Research Limited-Field Data, 21 Registries (excl IL), Nov 2024 Sub (2000–2022) - Linked To County Attributes - Time Dependent (1990-2023) Income/Rurality, 1969-2023 Counties, National Cancer Institute, DCCPS, Surveillance Research Program, released April 2025, based on the November 2024 submission",
  "• Expected survival (life tables): Surveillance, Epidemiology, and End Results (SEER) Program (seer.cancer.gov) SEER*Stat Database: Expected Survival - U.S. 1970-2021 by individual year (White, Black, Other (AI/API), Ages 0-99, All races for Other Unspec 1991+ and Unknown), National Cancer Institute, DCCPS, Surveillance Research Program, Cancer Statistics Branch.",
  "",
  "Suggested citation",
  "• To be added (paper in preparation)."
)

# Female-only sites (as per your description)
FEMALE_ONLY_SITES <- c("Breast", "Ovary", "Uterine Corpus")

# Publication/colour-blind-friendly palette (Okabe-Ito)
OKABE_ITO <- c(
  "Black" = "#000000",
  "White" = "#0072B2",
  "Male"  = "#009E73",
  "Female"= "#D55E00"
)

# ---- Load data once at start ----
shiny_data <- readr::read_csv(DATA_URL, show_col_types = FALSE) %>%
  mutate(
    tt = as.numeric(tt),
    PE = as.numeric(PE),
    lci = as.numeric(lci),
    uci = as.numeric(uci)
  )

# Clean/standardise site labels if needed (optional)
# shiny_data$site <- trimws(shiny_data$site)

# choices
SITE_CHOICES <- sort(unique(shiny_data$site))
SEX_CHOICES  <- c("Male", "Female")
RACE_CHOICES <- c("White", "Black")
AGE_CHOICES  <- c("All", "<55", "55-64", "65-74", "75+")

# Defaults you requested
DEFAULTS <- list(
  site = "Lung and bronchus",
  estimate_category = "Standardised relative survival",
  compare_by = "Race/Ethnicity",
  race_multi = c("White", "Black"),
  sex_single = "Male",
  agegroup = "55-64",
  stage = "Regional",
  mediation_effect = "Total causal effect",
  ad_effect = "Total avoidable deaths",
  show_ci = TRUE,
  img_titles = TRUE,
  img_legend = TRUE,
  img_footnote = TRUE,
  img_size = "Small (600×550)"
)

# ---- Helper mappings ----
stage_to_letter <- function(stage) {
  switch(stage,
    "Localised" = "l",
    "Regional"  = "r",
    "Distant"   = "d",
    "r"
  )
}
race_to_letter <- function(race) {
  ifelse(tolower(race) == "white", "w", "b")
}

# Resolve which estimate code(s) should be plotted and which extra codes needed
resolve_estimates <- function(category, stage, mediation_effect, ad_effect) {
  if (category == "Standardised relative survival") {
    # actual codes for survival are: lw, rw, dw, lb, rb, db
    # We will later filter by stage + race
    return(list(plot_codes = c("survival"), extra_codes = character(0)))
  }
  if (category == "Mediation analysis") {
    code <- switch(mediation_effect,
      "Total causal effect" = "tce",
      "Total direct effect" = "tde",
      "Total indirect effect (via stage)" = "tie",
      "tce"
    )
    extra <- character(0)
    # If TCE selected, also pull TIE to compute proportion mediated
    if (code == "tce") extra <- "tie"
    return(list(plot_codes = code, extra_codes = extra))
  }
  # Avoidable deaths
  code <- switch(ad_effect,
    "Total avoidable deaths" = "AD",
    "Avoidable deaths due to other relative survival differences" = "ADa",
    "Avoidable deaths due to stage differences" = "ADb",
    "AD"
  )
  return(list(plot_codes = code, extra_codes = character(0)))
}

# ---- Query string save/restore ----
# Convert input state -> query list
state_to_query <- function(input) {
  list(
    site = input$site,
    est = input$estimate_category,
    cmp = input$compare_by,
    sex = paste(input$sex %||% character(0), collapse = ","),
    race = input$race %||% "",
    racem = paste(input$race_multi %||% character(0), collapse = ","),
    sex1 = input$sex_single %||% "",
    age = input$agegroup,
    stage = input$stage %||% "",
    med = input$mediation_effect %||% "",
    ad  = input$ad_effect %||% "",
    ci  = if (isTRUE(input$show_ci)) "1" else "0"
  )
}

# Apply query list -> update inputs
apply_query <- function(session, q) {
  if (!is.null(q$site)) updateSelectInput(session, "site", selected = q$site)
  if (!is.null(q$est))  updateSelectInput(session, "estimate_category", selected = q$est)
  if (!is.null(q$cmp))  updateRadioButtons(session, "compare_by", selected = q$cmp)

  if (!is.null(q$sex) && nzchar(q$sex)) {
    updateCheckboxGroupInput(session, "sex", selected = strsplit(q$sex, ",", fixed = TRUE)[[1]])
  }
  if (!is.null(q$race) && nzchar(q$race)) updateSelectInput(session, "race", selected = q$race)

  if (!is.null(q$racem) && nzchar(q$racem)) {
    updateCheckboxGroupInput(session, "race_multi", selected = strsplit(q$racem, ",", fixed = TRUE)[[1]])
  }
  if (!is.null(q$sex1) && nzchar(q$sex1)) updateSelectInput(session, "sex_single", selected = q$sex1)

  if (!is.null(q$age)) updateSelectInput(session, "agegroup", selected = q$age)
  if (!is.null(q$stage) && nzchar(q$stage)) updateSelectInput(session, "stage", selected = q$stage)
  if (!is.null(q$med) && nzchar(q$med)) updateSelectInput(session, "mediation_effect", selected = q$med)
  if (!is.null(q$ad)  && nzchar(q$ad))  updateSelectInput(session, "ad_effect", selected = q$ad)
  if (!is.null(q$ci)) updateCheckboxInput(session, "show_ci", value = identical(q$ci, "1"))
}

# ---- UI ----
ui <- fluidPage(
  useShinyjs(),

  tags$head(tags$style(HTML("
    .top-row { display:flex; gap:12px; }
    .top-row .shiny-input-container { flex:1; }
    .icon-btns { display:flex; justify-content:flex-end; gap:10px; margin-top:6px; }
    .hint { display:flex; justify-content:flex-end; align-items:center; gap:8px; margin-top:4px; color:#444; font-size: 0.95em;}
    .footnote-box { margin-top: 18px; padding-top: 10px; border-top: 1px solid #ddd; font-size: 0.95em; }
    .footnote-box ul { margin-bottom: 0; }
  "))),

  titlePanel("Cancer outcomes explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      div(class = "top-row",
          selectInput("site", "Cancer site", choices = SITE_CHOICES, selected = DEFAULTS$site),
          selectInput("estimate_category", "Estimate", choices = c(
            "Standardised relative survival",
            "Mediation analysis",
            "Avoidable deaths"
          ), selected = DEFAULTS$estimate_category)
      ),

      radioButtons(
        "compare_by",
        label = "Compare by:",
        choices = c("Sex", "Race/Ethnicity"),
        selected = DEFAULTS$compare_by,
        inline = TRUE
      ),

      # Compare-by: Sex
      conditionalPanel(
        condition = "input.compare_by == 'Sex'",
        checkboxGroupInput("sex", "Sex", choices = SEX_CHOICES, selected = DEFAULTS$sex_single),
        selectInput("race", "Race/Ethnicity", choices = RACE_CHOICES, selected = "White")
      ),

      # Compare-by: Race
      conditionalPanel(
        condition = "input.compare_by == 'Race/Ethnicity'",
        checkboxGroupInput("race_multi", "Race/Ethnicity", choices = RACE_CHOICES, selected = DEFAULTS$race_multi),
        selectInput("sex_single", "Sex", choices = SEX_CHOICES, selected = DEFAULTS$sex_single)
      ),

      selectInput("agegroup", "Age", choices = AGE_CHOICES, selected = DEFAULTS$agegroup),

      conditionalPanel(
        condition = "input.estimate_category == 'Standardised relative survival'",
        selectInput("stage", "Stage at diagnosis", choices = c("Localised","Regional","Distant"), selected = DEFAULTS$stage)
      ),
      conditionalPanel(
        condition = "input.estimate_category == 'Mediation analysis'",
        selectInput("mediation_effect", "Causal effects", choices = c(
          "Total causal effect",
          "Total direct effect",
          "Total indirect effect (via stage)"
        ), selected = DEFAULTS$mediation_effect)
      ),
      conditionalPanel(
        condition = "input.estimate_category == 'Avoidable deaths'",
        selectInput("ad_effect", "Avoidable deaths", choices = c(
          "Total avoidable deaths",
          "Avoidable deaths due to other relative survival differences",
          "Avoidable deaths due to stage differences"
        ), selected = DEFAULTS$ad_effect)
      ),

      tags$hr(),
      checkboxInput("show_ci", "Show confidence intervals", value = TRUE)
    ),

    mainPanel(
      width = 8,

      div(class="icon-btns",
          actionButton("btn_download_plot", label = NULL, icon = icon("image"), title = "Download graph (PNG)"),
          actionButton("btn_download_data", label = NULL, icon = icon("download"), title = "Download data (CSV)"),
          actionButton("btn_copy_link",   label = NULL, icon = icon("link"), title = "Copy link")
      ),

      uiOutput("title_block"),

      div(class="hint",
          tags$span(icon("mouse-pointer")),
          tags$span("Tap/hover on points for more details")
      ),

      tabsetPanel(
        tabPanel("Graph", plotlyOutput("plot_main", height = "600px")),
        tabPanel("Data table", DTOutput("tbl"))
      ),

      div(class="footnote-box",
          tags$h4("Notes"),
          uiOutput("notes_dynamic"),
          tags$h4("Footnote"),
          tags$ul(
            lapply(FOOTNOTE_LINES, function(x) {
              if (nzchar(x)) tags$li(x) else tags$li(tags$br())
            })
          )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # Restore state from URL query params (if present)
  observeEvent(TRUE, {
    q <- parseQueryString(session$clientData$url_search)
    if (length(q) > 0) apply_query(session, q)
  }, once = TRUE)

  # Enforce: female-only sites cannot compare by Sex; also enforce female sex
  observeEvent(input$site, {
    if (input$site %in% FEMALE_ONLY_SITES) {
      updateRadioButtons(session, "compare_by", selected = "Race/Ethnicity")
      updateSelectInput(session, "sex_single", selected = "Female")
      updateCheckboxGroupInput(session, "sex", selected = "Female")
      disable("compare_by")   # locks the compare_by control
    } else {
      enable("compare_by")
    }
  }, ignoreInit = TRUE)

  # If site becomes female-only while compare_by=Sex, force it
  observe({
    if (input$site %in% FEMALE_ONLY_SITES && input$compare_by == "Sex") {
      updateRadioButtons(session, "compare_by", selected = "Race/Ethnicity")
    }
  })

  # Build a “current selection” description for subtitle
  selection_text <- reactive({
    parts <- c()

    parts <- c(parts, paste0("By ", input$compare_by))

    if (input$compare_by == "Race/Ethnicity") {
      parts <- c(parts, paste(input$race_multi, collapse = " vs "))
      parts <- c(parts, input$sex_single)
    } else {
      parts <- c(parts, paste(input$sex, collapse = " & "))
      parts <- c(parts, input$race)
    }

    parts <- c(parts, paste0("Ages ", input$agegroup))

    if (input$estimate_category == "Standardised relative survival") {
      parts <- c(parts, paste0(input$stage, " stage"))
    }
    paste(parts, collapse = ", ")
  })

  # Title block (title + subtitle)
  output$title_block <- renderUI({
    # main title is: site, estimate by time since diagnosis
    est_title <- input$estimate_category
    tags$div(
      tags$h2(input$site),
      tags$h4(paste0(est_title, " estimate by time since diagnosis")),
      tags$div(style="color:#444; font-size: 1.05em;",
               tags$b("Subtitle: "),
               selection_text()
      )
    )
  })

  # Dynamic notes under plot/table
  output$notes_dynamic <- renderUI({
    notes <- list()
    if (isTRUE(input$show_ci)) {
      notes <- c(notes, list(tags$div("Shaded areas are confidence intervals for the point estimates.")))
    }
    notes <- c(notes, list(
      tags$div("N* represents the number of Black people diagnosed in 2022 for this age group."),
      tags$div("Proportion mediated by stage is calculated as the total indirect effect (via stage) / total causal effect × 100.")
    ))
    do.call(tagList, notes)
  })

  # Determine groups (lines) to plot based on compare_by
  groups_df <- reactive({
    if (input$compare_by == "Race/Ethnicity") {
      tibble(group = input$race_multi, group_type = "Race/Ethnicity")
    } else {
      tibble(group = input$sex, group_type = "Sex")
    }
  })

  # Core filtered dataset for plotting (includes "extra codes" if needed for tooltip calculations)
  data_for_calc <- reactive({
    req(input$site, input$agegroup)

    base <- shiny_data %>%
      filter(site == input$site, agegroup == input$agegroup)

    # Apply sex restriction depending on compare mode
    if (input$compare_by == "Race/Ethnicity") {
      base <- base %>% filter(sex == input$sex_single)
    } else {
      base <- base %>% filter(sex %in% input$sex)
    }

    res <- resolve_estimates(
      category = input$estimate_category,
      stage = input$stage %||% "Regional",
      mediation_effect = input$mediation_effect %||% DEFAULTS$mediation_effect,
      ad_effect = input$ad_effect %||% DEFAULTS$ad_effect
    )

    # For survival, we later select estimate codes by stage+race
    # For mediation/AD, we just filter on estimate codes (plus extra codes if needed).
    if (input$estimate_category == "Standardised relative survival") {
      base
    } else {
      keep_codes <- unique(c(res$plot_codes, res$extra_codes))
      base %>% filter(estimate %in% keep_codes)
    }
  })

  # Data prepared specifically for plotting lines (only the plotted codes)
  data_for_plot <- reactive({
    df <- data_for_calc()

    if (input$estimate_category == "Standardised relative survival") {
      st <- stage_to_letter(input$stage)
      # line identities depend on compare_by:
      if (input$compare_by == "Race/Ethnicity") {
        # compare races; sex fixed
        races <- input$race_multi
        codes <- paste0(st, race_to_letter(races))  # e.g., rw, rb
        df <- df %>%
          filter(estimate %in% codes) %>%
          mutate(group = ifelse(grepl("w$", estimate), "White", "Black"))
      } else {
        # compare sexes; race fixed
        race_fixed <- input$race
        code <- paste0(st, race_to_letter(race_fixed)) # e.g., rw
        df <- df %>%
          filter(estimate == code) %>%
          mutate(group = sex)
      }
      return(df)
    }

    # Mediation / Avoidable deaths
    res <- resolve_estimates(
      category = input$estimate_category,
      stage = input$stage %||% "Regional",
      mediation_effect = input$mediation_effect %||% DEFAULTS$mediation_effect,
      ad_effect = input$ad_effect %||% DEFAULTS$ad_effect
    )

    df <- df %>% filter(estimate %in% res$plot_codes)

    if (input$compare_by == "Race/Ethnicity") {
      df %>% mutate(group = race)
    } else {
      df %>% mutate(group = sex)
    }
  })

  # Compute tooltip dataset for integer time points only
  points_for_tooltip <- reactive({
    dfp <- data_for_plot()
    if (nrow(dfp) == 0) return(dfp)

    pts <- dfp %>% filter(abs(tt - round(tt)) < 1e-9)  # integer years only

    # Add proportion mediated if needed (Mediation + TCE)
    if (input$estimate_category == "Mediation analysis" &&
        (input$mediation_effect %||% "") == "Total causal effect") {

      extra <- data_for_calc() %>%
        filter(estimate %in% c("tce","tie")) %>%
        mutate(group = if (input$compare_by == "Race/Ethnicity") race else sex) %>%
        filter(abs(tt - round(tt)) < 1e-9) %>%
        select(site, sex, agegroup, tt, group, estimate, PE)

      wide <- extra %>%
        pivot_wider(names_from = estimate, values_from = PE)

      pts <- pts %>%
        left_join(wide %>% select(tt, group, tce, tie), by = c("tt","group")) %>%
        mutate(prop_mediated = ifelse(!is.na(tce) & tce != 0, (tie / tce) * 100, NA_real_))
    }

    pts
  })

  # Y-axis label + limits
  y_axis <- reactive({
    if (input$estimate_category == "Standardised relative survival") {
      list(label = "Standardised all-cause survival", limits = c(-0.05, 1))
    } else if (input$estimate_category == "Mediation analysis") {
      list(label = "Difference in all-cause probabilities of death", limits = c(-0.05, 1))
    } else {
      list(label = "Avoidable deaths", limits = NULL) # range-based
    }
  })

  # Build ggplot base
  base_gg <- reactive({
    df <- data_for_plot()
    req(nrow(df) > 0)

    # Choose palette depending on compare mode
    pal <- if (input$compare_by == "Race/Ethnicity") OKABE_ITO[c("White","Black")] else OKABE_ITO[c("Male","Female")]

    g <- ggplot(df, aes(x = tt, y = PE, color = group, fill = group)) +
      geom_line(linewidth = 1) +
      scale_x_continuous("Time since diagnosis (years)", limits = c(0,5), breaks = 0:5) +
      scale_color_manual(values = pal, drop = FALSE) +
      scale_fill_manual(values = pal, drop = FALSE) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), # horizontal grid lines only
        legend.position = "bottom",
        legend.title = element_blank()
      )

    if (isTRUE(input$show_ci)) {
      g <- g + geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.18, colour = NA)
    }

    # y label/limits
    ya <- y_axis()
    g <- g + ylab(ya$label)
    if (!is.null(ya$limits)) {
      g <- g + coord_cartesian(ylim = ya$limits)
    }

    g
  })

  # Build plotly: main line plot + integer-year hover points + (optional) N* inset subplot for Avoidable deaths
  output$plot_main <- renderPlotly({
    df <- data_for_plot()
    req(nrow(df) > 0)

    pts <- points_for_tooltip()

    # Tooltip HTML
    # right-align numbers using a small table
    make_tip <- function(group, tt, PE, lci, uci, show_ci, extra_lines = character(0)) {
      rows <- c(
        sprintf("<tr><td><b>%s</b></td><td style='text-align:right;'><b>%s</b></td></tr>", "Group", group),
        sprintf("<tr><td><b><u>%s</u></b></td><td style='text-align:right;'><b><u>%.0f</u></b></td></tr>", "Time since diagnosis (years)", tt),
        sprintf("<tr><td>%s</td><td style='text-align:right;'>%.3f</td></tr>", "Estimate", PE)
      )
      if (show_ci) {
        rows <- c(rows,
                  sprintf("<tr><td>%s</td><td style='text-align:right;'>%.3f</td></tr>", "Lower 95% CI", lci),
                  sprintf("<tr><td>%s</td><td style='text-align:right;'>%.3f</td></tr>", "Upper 95% CI", uci)
        )
      }
      if (length(extra_lines) > 0) {
        rows <- c(rows, extra_lines)
      }
      paste0("<table>", paste(rows, collapse = ""), "</table>")
    }

    show_ci <- isTRUE(input$show_ci)

    # add extra fields for mediation tce (prop mediated) and avoidable deaths (N*)
    extra_html <- rep("", nrow(pts))
    if (nrow(pts) > 0) {
      extra_list <- vector("list", nrow(pts))
      for (i in seq_len(nrow(pts))) {
        extras <- character(0)
        if (input$estimate_category == "Mediation analysis" &&
            (input$mediation_effect %||% "") == "Total causal effect" &&
            !is.na(pts$prop_mediated[i])) {
          extras <- c(extras,
                      sprintf("<tr><td>%s</td><td style='text-align:right;'>%.1f%%</td></tr>",
                              "Proportion mediated by stage", pts$prop_mediated[i]))
        }
        if (input$estimate_category == "Avoidable deaths") {
          nb <- suppressWarnings(as.numeric(pts$n_black[i]))
          if (!is.na(nb)) {
            extras <- c(extras,
                        sprintf("<tr><td>%s</td><td style='text-align:right;'>%s</td></tr>",
                                "N*", format(round(nb, 0), big.mark = ",")))
          }
        }
        extra_list[[i]] <- extras
      }
      tips <- mapply(
        FUN = make_tip,
        group = pts$group,
        tt = pts$tt,
        PE = pts$PE,
        lci = pts$lci,
        uci = pts$uci,
        MoreArgs = list(show_ci = show_ci),
        extra_lines = extra_list,
        SIMPLIFY = TRUE
      )
      pts$tip <- tips
    } else {
      pts$tip <- character(0)
    }

    g <- base_gg()

    # Overlay integer-year points (so hover only there)
    g <- g + geom_point(
      data = df %>% filter(abs(tt - round(tt)) < 1e-9),
      aes(x = tt, y = PE, color = group),
      size = 2.6
    )

    pl <- ggplotly(g, tooltip = "text") %>% layout(hovermode = "closest")
    # Replace tooltip with our custom tip by adding a marker trace (clean control)
    if (nrow(pts) > 0) {
      pl <- pl %>%
        add_markers(
          data = pts,
          x = ~tt, y = ~PE,
          color = ~group, colors = unname(if (input$compare_by == "Race/Ethnicity") OKABE_ITO[c("White","Black")] else OKABE_ITO[c("Male","Female")]),
          marker = list(size = 9, opacity = 0.0), # invisible markers for tooltip targeting
          text = ~tip,
          hoverinfo = "text",
          showlegend = FALSE
        )
    }

    # Avoidable deaths: add small N* bar chart as a subplot
    if (input$estimate_category == "Avoidable deaths") {
      # Use one N* per group (same repeated) – take first non-missing
      nb <- df %>%
        group_by(group) %>%
        summarise(n_black = suppressWarnings(as.numeric(first(n_black))), .groups = "drop") %>%
        mutate(n_black = ifelse(is.na(n_black), 0, n_black))

      bar <- plot_ly(
        nb, x = ~group, y = ~n_black, type = "bar",
        hovertemplate = paste0("<b>%{x}</b><br>N*: %{y:,}<extra></extra>")
      ) %>%
        layout(
          title = list(text = "N*", font = list(size = 12)),
          margin = list(l=10,r=10,t=30,b=10),
          xaxis = list(title = "", tickfont = list(size=10)),
          yaxis = list(title = "", tickfont = list(size=10))
        )

      pl <- subplot(pl, bar, widths = c(0.78, 0.22), margin = 0.02) %>%
        layout(showlegend = TRUE)
    }

    pl
  })

  # ---- Data table ----
  output$tbl <- renderDT({
    df <- data_for_plot()
    req(nrow(df) > 0)

    out <- df %>%
      select(tt, group, PE, lci, uci) %>%
      arrange(tt, group)

    # Wide table: tt + per-group columns
    if (isTRUE(input$show_ci)) {
      wide <- out %>%
        pivot_wider(
          names_from = group,
          values_from = c(PE, lci, uci),
          names_glue = "{group}__{.value}"
        ) %>%
        arrange(tt)
    } else {
      wide <- out %>%
        select(tt, group, PE) %>%
        pivot_wider(
          names_from = group,
          values_from = PE,
          names_glue = "{group}__PE"
        ) %>%
        arrange(tt)
    }

    # Build a grouped header (2-level)
    groups <- unique(out$group)
    base_cols <- c("tt")
    # per group: PE (+ lci/uci)
    per_group_cols <- if (isTRUE(input$show_ci)) c("PE","lci","uci") else c("PE")

    header_top <- tags$tr(
      tags$th(rowspan = 2, "Time since diagnosis (years)"),
      lapply(groups, function(g) tags$th(colspan = length(per_group_cols), g))
    )
    header_bottom <- tags$tr(
      lapply(groups, function(g) {
        lapply(per_group_cols, function(v) {
          lab <- if (v == "PE") "Estimate" else if (v == "lci") "Lower 95% CI" else "Upper 95% CI"
          tags$th(lab)
        })
      }) |> unlist(recursive = FALSE)
    )

    container <- tags$table(
      class = 'display',
      thead = tags$thead(header_top, header_bottom)
    )

    datatable(
      wide,
      container = container,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(0, "asc"))
      )
    )
  })

  # ---- Download: data (CSV) ----
  output$download_data_csv <- downloadHandler(
    filename = function() {
      safe_site <- gsub("[^A-Za-z0-9]+", "_", input$site)
      paste0("data_", safe_site, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- data_for_plot()
      req(nrow(df) > 0)

      meta1 <- paste0("Cancer site: ", input$site)
      meta2 <- paste0("Estimate: ", input$estimate_category, " by time since diagnosis")
      meta3 <- paste0("By: ", selection_text())

      # table body (same shape as shown in DT)
      out <- df %>%
        select(tt, group, PE, lci, uci) %>%
        arrange(tt, group)

      if (isTRUE(input$show_ci)) {
        wide <- out %>%
          pivot_wider(
            names_from = group,
            values_from = c(PE, lci, uci),
            names_glue = "{group}_{.value}"
          ) %>%
          arrange(tt)
      } else {
        wide <- out %>%
          select(tt, group, PE) %>%
          pivot_wider(
            names_from = group,
            values_from = PE,
            names_glue = "{group}_PE"
          ) %>%
          arrange(tt)
      }

      # write CSV with metadata header + footnote at bottom
      con <- file(file, open = "wt", encoding = "UTF-8")
      on.exit(close(con), add = TRUE)

      writeLines(meta1, con)
      writeLines(meta2, con)
      writeLines(meta3, con)
      writeLines("", con)

      suppressWarnings(write.csv(wide, con, row.names = FALSE))

      writeLines("", con)
      writeLines("Footnote:", con)
      writeLines(FOOTNOTE_LINES, con)
    }
  )

  # ---- Download: plot (PNG) ----
  output$download_plot_png <- downloadHandler(
    filename = function() {
      safe_site <- gsub("[^A-Za-z0-9]+", "_", input$site)
      paste0("plot_", safe_site, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      df <- data_for_plot()
      req(nrow(df) > 0)

      # size
      size_map <- list(
        "Small (600×550)"  = c(600, 550),
        "Medium (900×550)" = c(900, 550),
        "Large (1200×550)" = c(1200, 550)
      )
      wh <- size_map[[ input$img_size %||% "Small (600×550)" ]]
      wpx <- wh[1]; hpx <- wh[2]

      # base plot
      pal <- if (input$compare_by == "Race/Ethnicity") OKABE_ITO[c("White","Black")] else OKABE_ITO[c("Male","Female")]

      g <- ggplot(df, aes(x = tt, y = PE, color = group, fill = group)) +
        geom_line(linewidth = 1) +
        scale_x_continuous("Time since diagnosis (years)", limits = c(0,5), breaks = 0:5) +
        scale_color_manual(values = pal, drop = FALSE) +
        scale_fill_manual(values = pal, drop = FALSE) +
        theme_minimal(base_size = 13) +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = if (isTRUE(input$img_legend)) "bottom" else "none",
          legend.title = element_blank()
        )

      if (isTRUE(input$show_ci)) {
        g <- g + geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.18, colour = NA)
      }

      ya <- y_axis()
      g <- g + ylab(ya$label)
      if (!is.null(ya$limits)) g <- g + coord_cartesian(ylim = ya$limits)

      # titles
      if (isTRUE(input$img_titles)) {
        g <- g + labs(
          title = input$site,
          subtitle = paste0(input$estimate_category, " estimate by time since diagnosis\n", selection_text())
        )
      } else {
        g <- g + labs(title = NULL, subtitle = NULL)
      }

      # footnote
      if (isTRUE(input$img_footnote)) {
        g <- g + labs(caption = paste(
          c(
            if (isTRUE(input$show_ci)) "Shaded areas are confidence intervals for the point estimates." else NULL,
            "N* represents the number of Black people diagnosed in 2022 for this age group.",
            "Proportion mediated by stage = (TIE / TCE) × 100.",
            "",
            FOOTNOTE_LINES
          ),
          collapse = "\n"
        )) +
          theme(plot.caption = element_text(size = 8, hjust = 0))
      }

      # save
      dpi <- 150
      ggsave(file, g, width = wpx/dpi, height = hpx/dpi, dpi = dpi, bg = "white")
    }
  )

  # ---- Download plot modal with preview ----
  observeEvent(input$btn_download_plot, {
    showModal(modalDialog(
      title = "Download graph (PNG)",
      size = "l",
      fluidRow(
        column(
          4,
          h4("Image options"),
          checkboxInput("img_titles", "Show titles", value = TRUE),
          checkboxInput("img_legend", "Show legend", value = TRUE),
          checkboxInput("img_footnote", "Show footnote", value = TRUE),
          br(),
          selectInput("img_size", "Image size", choices = c("Small (600×550)", "Medium (900×550)", "Large (1200×550)"),
                      selected = "Small (600×550)"),
          br(),
          downloadButton("download_plot_png", "Download PNG")
        ),
        column(
          8,
          h4("Preview (layout guide)"),
          plotOutput("preview_layout", height = "360px")
        )
      ),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$preview_layout <- renderPlot({
    par(mar = c(1,1,1,1))
    plot.new()
    # Title area
    rect(0.05, 0.82, 0.95, 0.95)
    text(0.5, 0.885, "Title", cex = 1.1)
    # Main plot area
    rect(0.05, 0.18, 0.72, 0.80)
    text(0.385, 0.49, "Graph", cex = 1.2)
    # Legend area
    rect(0.75, 0.45, 0.95, 0.80)
    text(0.85, 0.625, "Legend", cex = 1.1)
    # Footnote area
    rect(0.05, 0.05, 0.95, 0.15)
    text(0.5, 0.10, "Footnote", cex = 1.1)
  })

  # ---- Download data button triggers actual download ----
  observeEvent(input$btn_download_data, {
    shinyjs::runjs("document.getElementById('download_data_csv').click();")
  })

  # Hidden download button for programmatic click
  outputOptions(output, "download_data_csv", suspendWhenHidden = FALSE)

  # Create a hidden downloadButton in UI dynamically
  # (DT/main UI doesn't include it explicitly)
  insertUI(
    selector = "body",
    where = "beforeEnd",
    ui = tags$div(style="display:none;", downloadButton("download_data_csv", "download"))
  )

  # ---- Copy link ----
  observeEvent(input$btn_copy_link, {
    q <- state_to_query(input)
    # Build query string
    qs <- paste0(
      "?",
      paste0(names(q), "=", vapply(q, URLencode, FUN.VALUE = character(1), reserved = TRUE), collapse = "&")
    )
    # Update browser URL
    updateQueryString(qs, mode = "replace", session = session)

    showModal(modalDialog(
      title = "Link copied",
      tags$p("Your page URL has been updated with the current selections."),
      tags$p("Copy the URL from your browser address bar."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

}

shinyApp(ui, server)

