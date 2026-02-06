library(shiny)
library(shinyjs)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(DT)
library(plotly)
library(patchwork)

# Replace with a GitHub raw URL when the data is hosted in your repo.
data_source <- "Shiny.csv"

load_data <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      site = as.character(site),
      sex = as.character(sex),
      agegroup = as.character(agegroup),
      estimate = as.character(estimate)
    )
}

standardised_map <- function(stage, race) {
  stage_code <- c(
    "Localised" = "l",
    "Regional" = "r",
    "Distant" = "d"
  )[[stage]]
  race_code <- c(
    "White" = "w",
    "Black" = "b"
  )[[race]]
  paste0(stage_code, race_code)
}

estimate_label <- function(category, option) {
  if (category == "Standardised relative survival") {
    return("Standardised relative survival estimate")
  }
  if (category == "Mediation analysis") {
    return(
      c(
        "Total causal effect" = "Total causal effect",
        "Direct effect" = "Total direct effect",
        "Indirect effect (via stage)" = "Total indirect effect (via stage)"
      )[[option]]
    )
  }
  c(
    "Total avoidable deaths" = "Total avoidable deaths",
    "Avoidable deaths due to other relative survival differences" =
      "Avoidable deaths due to other relative survival differences",
    "Avoidable deaths due to stage differences" =
      "Avoidable deaths due to stage differences"
  )[[option]]
}

estimate_code <- function(category, option) {
  if (category == "Mediation analysis") {
    return(
      c(
        "Total causal effect" = "tce",
        "Direct effect" = "tde",
        "Indirect effect (via stage)" = "tie"
      )[[option]]
    )
  }
  if (category == "Avoidable deaths") {
    return(
      c(
        "Total avoidable deaths" = "AD",
        "Avoidable deaths due to other relative survival differences" = "ADa",
        "Avoidable deaths due to stage differences" = "ADb"
      )[[option]]
    )
  }
  NULL
}

comparison_label <- function(compare_by, selections) {
  if (compare_by == "race") {
    return(paste(selections, collapse = ", "))
  }
  paste(selections, collapse = ", ")
}

footnote_text <- paste(
  "Model fitted",
  "• Flexible parametric relative survival model fitted on the log cumulative excess hazard scale (stratified by sex) with a period window of 01/01/2021 to 31/12/2022",
  "• Baseline cumulative excess hazard modelled using restricted cubic splines with 3 degrees of freedom.",
  "• Main effects included: race, stage, age, race × stage, stage × age, race × age.",
  "• Age modelled non-linearly using restricted cubic splines with 3 degrees of freedom.",
  "• Age was modelled continuously for the central 95% of the age distribution, while patients younger than 2.5th and older than 97.5th age percentile were assigned the same expected survival as those at these respective cut-off ages.",
  "• Time-dependent effects were included for race, stage, age, race × stage and race × age, each with 2 degrees of freedom.",
  "• Expected survival was incorporated using SEER population life tables stratified by race, year, age and sex.",
  "Data sources",
  "• Incidence data: Surveillance, Epidemiology, and End Results (SEER) Program (seer.cancer.gov) SEER*Stat Database: Incidence - SEER Research Limited-Field Data, 21 Registries (excl IL), Nov 2024 Sub (2000–2022) - Linked To County Attributes - Time Dependent (1990-2023) Income/Rurality, 1969-2023 Counties, National Cancer Institute, DCCPS, Surveillance Research Program, released April 2025, based on the November 2024 submission",
  "• Expected survival (life tables): Surveillance, Epidemiology, and End Results (SEER) Program (seer.cancer.gov) SEER*Stat Database: Expected Survival - U.S. 1970-2021 by individual year (White, Black, Other (AI/API), Ages 0-99, All races for Other Unspec 1991+ and Unknown), National Cancer Institute, DCCPS, Surveillance Research Program, Cancer Statistics Branch.",
  "Suggested citation",
  "• Citation to be updated once the manuscript title is finalized.",
  sep = "\n"
)

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(
      HTML(
        "
        .top-selects { display: flex; gap: 16px; flex-wrap: wrap; }
        .plot-actions { display: flex; gap: 10px; justify-content: flex-end; }
        .plot-actions .btn { padding: 6px 8px; }
        .footnote { font-size: 0.85em; color: #444; white-space: pre-line; }
        .legend-note { font-size: 0.85em; color: #555; margin-top: 8px; }
        "
      )
    )
  ),
  titlePanel("Ethnic disparities in SEER"),
  div(
    class = "top-selects",
    selectInput(
      "site",
      "Cancer site",
      choices = c(
        "Breast",
        "Colon",
        "Uterine Corpus",
        "Kidney and renal pelvis",
        "Lung and bronchus",
        "Oral cavity and pharynx",
        "Ovary",
        "Pancreas",
        "Rectum"
      ),
      selected = "Lung and bronchus"
    ),
    selectInput(
      "estimate_category",
      "Estimate category",
      choices = c(
        "Standardised relative survival",
        "Mediation analysis",
        "Avoidable deaths"
      ),
      selected = "Standardised relative survival"
    )
  ),
  tabsetPanel(
    id = "compare_by",
    type = "pills",
    selected = "race",
    tabPanel("Compare by: Sex", value = "sex"),
    tabPanel("Compare by: Race/Ethnicity", value = "race")
  ),
  sidebarLayout(
    sidebarPanel(
      uiOutput("compare_controls"),
      uiOutput("conditional_controls"),
      checkboxGroupInput(
        "more_options",
        "More options",
        choices = c("Show confidence intervals" = "show_ci"),
        selected = "show_ci"
      )
    ),
    mainPanel(
      div(
        class = "plot-actions",
        actionButton("download_plot", label = NULL, icon = icon("download"), title = "Download graph"),
        downloadButton("download_data", label = NULL, icon = icon("file-csv"), title = "Download data"),
        actionButton("copy_link", label = NULL, icon = icon("link"), title = "Copy link")
      ),
      uiOutput("plot_title"),
      tags$div(
        style = "margin-bottom: 8px;",
        icon("mouse-pointer"),
        tags$span(" Tap/hover on points for more details")
      ),
      tabsetPanel(
        id = "plot_tab",
        tabPanel(
          "Graph",
          plotlyOutput("estimate_plot", height = "600px"),
          uiOutput("legend_note")
        ),
        tabPanel("Data table", uiOutput("table_panels"))
      ),
      tags$div(class = "footnote", footnote_text)
    )
  )
)

server <- function(input, output, session) {
  data <- reactive({
    load_data(data_source)
  })

  female_only_sites <- c("Breast", "Ovary", "Uterine Corpus")

  observeEvent(input$site, {
    if (input$site %in% female_only_sites && input$compare_by == "sex") {
      updateTabsetPanel(session, "compare_by", selected = "race")
    }
  })

  observeEvent(input$compare_by, {
    if (input$site %in% female_only_sites && input$compare_by == "sex") {
      updateTabsetPanel(session, "compare_by", selected = "race")
    }
  })

  output$compare_controls <- renderUI({
    if (input$compare_by == "sex") {
      tagList(
        checkboxGroupInput(
          "sexes",
          "Sex",
          choices = c("Male", "Female"),
          selected = c("Male", "Female")
        ),
        selectInput(
          "race_single",
          "Race/Ethnicity",
          choices = c("White", "Black"),
          selected = "White"
        ),
        selectInput(
          "agegroup",
          "Age",
          choices = c("All", "<55", "55-64", "65-74", "75+"),
          selected = "55-64"
        )
      )
    } else {
      tagList(
        checkboxGroupInput(
          "races",
          "Race/Ethnicity",
          choices = c("White", "Black"),
          selected = c("White", "Black")
        ),
        selectInput(
          "sex_single",
          "Sex",
          choices = c("Male", "Female"),
          selected = "Male"
        ),
        selectInput(
          "agegroup",
          "Age",
          choices = c("All", "<55", "55-64", "65-74", "75+"),
          selected = "55-64"
        )
      )
    }
  })

  output$conditional_controls <- renderUI({
    if (input$estimate_category == "Standardised relative survival") {
      selectInput(
        "stage",
        "Stage at diagnosis",
        choices = c("Localised", "Regional", "Distant"),
        selected = "Regional"
      )
    } else if (input$estimate_category == "Mediation analysis") {
      selectInput(
        "mediation_effect",
        "Causal effects",
        choices = c(
          "Total causal effect",
          "Direct effect",
          "Indirect effect (via stage)"
        ),
        selected = "Total causal effect"
      )
    } else {
      selectInput(
        "avoidable_effect",
        "Avoidable deaths",
        choices = c(
          "Total avoidable deaths",
          "Avoidable deaths due to other relative survival differences",
          "Avoidable deaths due to stage differences"
        ),
        selected = "Total avoidable deaths"
      )
    }
  })

  filtered_data <- reactive({
    df <- data() %>% filter(site == input$site)

    if (input$compare_by == "sex") {
      df <- df %>%
        filter(
          sex %in% input$sexes,
          agegroup == input$agegroup
        )
      if (input$estimate_category == "Standardised relative survival") {
        codes <- map_chr(input$sexes, ~ standardised_map(input$stage, input$race_single))
        df <- df %>%
          filter(estimate %in% codes) %>%
          mutate(compare = sex)
      } else {
        code <- estimate_code(
          input$estimate_category,
          if (input$estimate_category == "Mediation analysis") {
            input$mediation_effect
          } else {
            input$avoidable_effect
          }
        )
        df <- df %>%
          filter(estimate == code) %>%
          mutate(compare = sex)
      }
    } else {
      df <- df %>%
        filter(
          sex == input$sex_single,
          agegroup == input$agegroup
        )
      if (input$estimate_category == "Standardised relative survival") {
        codes <- map_chr(input$races, ~ standardised_map(input$stage, .x))
        df <- df %>%
          filter(estimate %in% codes) %>%
          mutate(compare = if_else(estimate %in% c("lw", "rw", "dw"), "White", "Black"))
      } else {
        code <- estimate_code(
          input$estimate_category,
          if (input$estimate_category == "Mediation analysis") {
            input$mediation_effect
          } else {
            input$avoidable_effect
          }
        )
        df <- df %>%
          filter(estimate == code) %>%
          mutate(compare = "Black vs White")
      }
    }

    df %>% mutate(tt = as.numeric(tt))
  })

  output$plot_title <- renderUI({
    estimate_text <- estimate_label(
      input$estimate_category,
      if (input$estimate_category == "Mediation analysis") {
        input$mediation_effect
      } else if (input$estimate_category == "Avoidable deaths") {
        input$avoidable_effect
      } else {
        NULL
      }
    )

    by_items <- c(
      if (input$compare_by == "race") "Race" else "Sex",
      if (input$compare_by == "race") input$sex_single else input$race_single,
      paste("Ages", input$agegroup),
      if (input$estimate_category == "Standardised relative survival") {
        paste(input$stage, "stage")
      } else {
        NULL
      }
    )

    tags$div(
      tags$h3(input$site),
      tags$h4(paste(estimate_text, "by time since diagnosis")),
      tags$h5(paste("By", paste(by_items, collapse = ", ")))
    )
  })

  plot_data <- reactive({
    df <- filtered_data()
    tie_df <- NULL

    if (input$estimate_category == "Mediation analysis" &&
      input$mediation_effect == "Total causal effect") {
      tie_df <- data() %>%
        filter(
          site == input$site,
          agegroup == input$agegroup,
          estimate == "tie",
          sex %in% if (input$compare_by == "race") input$sex_single else input$sexes
        ) %>%
        mutate(
          tt = as.numeric(tt),
          compare = if (input$compare_by == "race") "Black vs White" else sex
        )
    }

    df %>%
      group_by(compare) %>%
      mutate(
        tt_int = tt %in% 0:5,
        proportion_mediated = if (!is.null(tie_df)) {
          tie_value <- tie_df$PE[match(paste(compare, tt), paste(tie_df$compare, tie_df$tt))]
          if_else(is.na(tie_value), NA_real_, tie_value / PE * 100)
        } else {
          NA_real_
        },
        tooltip = if_else(
          tt_int,
          paste0(
            "<b>", compare, "</b><br>",
            "<b><u>Time since diagnosis</u></b>: ", sprintf("%.0f", tt), " years<br>",
            "Estimate: ", sprintf("%.3f", PE),
            if ("show_ci" %in% input$more_options) {
              paste0("<br>95% CI: ", sprintf("%.3f", lci), "–", sprintf("%.3f", uci))
            } else {
              ""
            },
            if (!is.null(tie_df)) {
              paste0("<br>Proportion mediated: ", sprintf("%.1f", proportion_mediated), "%")
            } else {
              ""
            },
            if (input$estimate_category == "Avoidable deaths") {
              paste0("<br>N*: ", n_black)
            } else {
              ""
            }
          ),
          ""
        )
      ) %>%
      ungroup()
  })

  build_plot <- reactive({
    df <- plot_data()
    y_label <- switch(
      input$estimate_category,
      "Avoidable deaths" = "Avoidable deaths",
      "Mediation analysis" = "Difference in all-cause probabilities of death",
      "Standardised relative survival" = "Standardised all-cause survival"
    )

    base_plot <- ggplot(df, aes(x = tt, y = PE, color = compare)) +
      geom_line(linewidth = 1) +
      geom_point(
        data = df %>% filter(tt_int),
        aes(text = tooltip),
        size = 2
      ) +
      scale_x_continuous(breaks = 0:5, limits = c(0, 5)) +
      labs(x = "Time since diagnosis (years)", y = y_label, color = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )

    if ("show_ci" %in% input$more_options) {
      base_plot <- base_plot +
        geom_ribbon(
          aes(ymin = lci, ymax = uci, fill = compare),
          alpha = 0.2,
          color = NA
        ) +
        scale_fill_discrete(name = NULL)
    }

    if (input$estimate_category %in% c("Mediation analysis", "Standardised relative survival")) {
      base_plot <- base_plot +
        scale_y_continuous(limits = c(-0.05, 1), breaks = seq(-0.05, 1, by = 0.2))
    }

    if (input$estimate_category == "Avoidable deaths") {
      range_vals <- range(df$PE, na.rm = TRUE)
      base_plot <- base_plot +
        scale_y_continuous(limits = range_vals)
    }

    if (input$estimate_category == "Avoidable deaths") {
      n_black_df <- df %>%
        group_by(compare) %>%
        summarise(n_black = first(n_black), .groups = "drop")

      bar_plot <- ggplot(n_black_df, aes(x = compare, y = n_black, fill = compare)) +
        geom_col(width = 0.6) +
        coord_flip() +
        theme_minimal(base_size = 9) +
        theme(
          axis.title = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          axis.text.y = element_text(size = 8)
        )

      base_plot <- base_plot + inset_element(bar_plot, left = 0.65, bottom = 0.05, right = 0.98, top = 0.35)
    }

    base_plot
  })

  output$estimate_plot <- renderPlotly({
    plot <- build_plot()
    ggplotly(plot, tooltip = "text") %>%
      layout(legend = list(orientation = "h", x = 0.3, y = -0.2))
  })

  output$legend_note <- renderUI({
    note_lines <- c(
      if ("show_ci" %in% input$more_options) {
        "Shaded areas are confidence intervals for the point estimates."
      } else {
        NULL
      },
      "N* represents the number of Black people diagnosed in 2022 for this age group.",
      "Proportion mediated by stage is calculated as the total indirect effect (via stage) / total causal effect * 100."
    )

    tags$p(class = "legend-note", paste(note_lines, collapse = " "))
  })

  output$table_panels <- renderUI({
    df <- filtered_data()
    groups <- unique(df$compare)
    if (length(groups) == 1) {
      DTOutput("table_single")
    } else {
      fluidRow(
        lapply(seq_along(groups), function(i) {
          column(6, DTOutput(paste0("table_", i)))
        })
      )
    }
  })

  observe({
    df <- filtered_data()
    groups <- unique(df$compare)

    if (length(groups) == 1) {
      output$table_single <- renderDT({
        table_df <- df %>%
          select(tt, PE, lci, uci) %>%
          rename(
            `Time since diagnosis (years)` = tt,
            `Point estimate` = PE,
            `Lower 95% CI` = lci,
            `Upper 95% CI` = uci
          )
        if (!"show_ci" %in% input$more_options) {
          table_df <- table_df %>% select(-`Lower 95% CI`, -`Upper 95% CI`)
        }

        table_df %>%
          datatable(options = list(pageLength = 10))
      })
    } else {
      walk(seq_along(groups), function(i) {
        group <- groups[[i]]
        output[[paste0("table_", i)]] <- renderDT({
          table_df <- df %>%
            filter(compare == group) %>%
            select(tt, PE, lci, uci) %>%
            rename(
              `Time since diagnosis (years)` = tt,
              `Point estimate` = PE,
              `Lower 95% CI` = lci,
              `Upper 95% CI` = uci
            )
          if (!"show_ci" %in% input$more_options) {
            table_df <- table_df %>% select(-`Lower 95% CI`, -`Upper 95% CI`)
          }

          table_df %>%
            datatable(
              caption = htmltools::tags$caption(
                style = "caption-side: top; text-align: left;",
                group
              ),
              options = list(pageLength = 10)
            )
        })
      })
    }
  })

  output$download_data <- downloadHandler(
    filename = function() {
      paste0("seer-", gsub(" ", "-", tolower(input$site)), ".csv")
    },
    content = function(file) {
      df <- filtered_data() %>% arrange(compare, tt)
      header_lines <- c(
        paste("Cancer site:", input$site),
        paste("Estimate:", estimate_label(
          input$estimate_category,
          if (input$estimate_category == "Mediation analysis") {
            input$mediation_effect
          } else if (input$estimate_category == "Avoidable deaths") {
            input$avoidable_effect
          } else {
            NULL
          }
        )),
        paste("By:", if (input$compare_by == "race") "Race" else "Sex")
      )

      table_data <- df %>%
        select(compare, tt, PE, lci, uci) %>%
        arrange(compare, tt)

      writeLines(header_lines, file)
      write.table(
        table_data,
        file,
        sep = ",",
        row.names = FALSE,
        col.names = TRUE,
        append = TRUE
      )
      writeLines("", file, sep = "\n", useBytes = TRUE)
      writeLines(footnote_text, file, sep = "\n", useBytes = TRUE)
    }
  )

  observeEvent(input$copy_link, {
    shinyjs::runjs(
      "navigator.clipboard.writeText(window.location.href).then(() => { console.log('Copied'); });"
    )
    showNotification("Link copied to clipboard.", type = "message")
  })

  observeEvent(input$download_plot, {
    showModal(
      modalDialog(
        title = "Download plot options",
        checkboxGroupInput(
          "image_options",
          "Image options",
          choices = c("Show titles" = "show_titles", "Show legend" = "show_legend", "Show footnote" = "show_footnote"),
          selected = c("show_titles", "show_legend", "show_footnote")
        ),
        selectInput(
          "image_size",
          "Image size",
          choices = c("Small (600x550)" = "small", "Medium (900x550)" = "medium", "Large (1200x550)" = "large"),
          selected = "small"
        ),
        plotOutput("download_preview", height = "300px"),
        footer = tagList(
          modalButton("Close"),
          downloadButton("confirm_download", "Download PNG")
        ),
        size = "l"
      )
    )
  })

  output$download_preview <- renderPlot({
    ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.85, ymax = 1, fill = "#d9d9d9") +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.15, ymax = 0.85, fill = "#f2f2f2") +
      annotate("rect", xmin = 0.75, xmax = 1, ymin = 0.15, ymax = 0.6, fill = "#e6e6e6") +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 0.15, fill = "#d9d9d9") +
      annotate("text", x = 0.5, y = 0.93, label = "Title") +
      annotate("text", x = 0.5, y = 0.5, label = "Graph") +
      annotate("text", x = 0.875, y = 0.38, label = "Legend") +
      annotate("text", x = 0.5, y = 0.06, label = "Footnote") +
      theme_void()
  })

  output$confirm_download <- downloadHandler(
    filename = function() {
      paste0("seer-plot-", Sys.Date(), ".png")
    },
    content = function(file) {
      size <- switch(
        input$image_size,
        small = c(600, 550),
        medium = c(900, 550),
        large = c(1200, 550)
      )

      plot <- build_plot()
      if (!"show_legend" %in% input$image_options) {
        plot <- plot + theme(legend.position = "none")
      }
      if ("show_titles" %in% input$image_options) {
        plot <- plot + labs(
          title = input$site,
          subtitle = paste(estimate_label(
            input$estimate_category,
            if (input$estimate_category == "Mediation analysis") {
              input$mediation_effect
            } else if (input$estimate_category == "Avoidable deaths") {
              input$avoidable_effect
            } else {
              NULL
            }
          ), "by time since diagnosis")
        )
      }
      if ("show_footnote" %in% input$image_options) {
        plot <- plot + labs(caption = footnote_text)
      }
      ggsave(file, plot = plot, width = size[1] / 100, height = size[2] / 100, dpi = 100)
    }
  )
}

shinyApp(ui, server)
