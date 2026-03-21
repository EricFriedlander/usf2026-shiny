# tab2_module.R
# Shiny module for Tab 2: Out of Africa Demographic Explorer
#
# Students tune demographic parameters to match simulated allele frequency
# spectra (AFS) to pre-computed "real" data, learning about human history.
#
# Exports:
#   tab2UI(id)     — UI function
#   tab2Server(id) — Server function
#
# Depends on: ooa_simulation.R (must be sourced before this module)

library(shiny)
library(ggplot2)
library(dplyr)

# True hidden parameters (the "answer")
TRUE_T_SPLIT        <- 2400
TRUE_BOTTLENECK     <- 500
TRUE_GROWTH_RATE    <- 2.0

# Display colors
COLOR_AFRICA    <- "#E67E22"   # warm orange
COLOR_NONAFRICA <- "#2980B9"   # blue
COLOR_ANCESTRAL <- "#8E44AD"   # purple
COLOR_GOOD      <- "#27AE60"   # green
COLOR_OK        <- "#F39C12"   # yellow
COLOR_BAD       <- "#E74C3C"   # red


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

#' Tab 2 UI
#'
#' @param id  Shiny module namespace id.
tab2UI <- function(id) {
  ns <- NS(id)

  tagList(
    # --- Parameter sliders --------------------------------------------------
    fluidRow(
      column(
        width = 4,
        sliderInput(
          inputId = ns("t_split"),
          label   = "Time of split out of Africa (generations ago)",
          min = 500, max = 5000, value = 1500, step = 100
        )
      ),
      column(
        width = 4,
        sliderInput(
          inputId = ns("bottleneck_size"),
          label   = "Bottleneck population size",
          min = 100, max = 2000, value = 1000, step = 50
        )
      ),
      column(
        width = 4,
        sliderInput(
          inputId = ns("growth_rate"),
          label   = "Post-split growth rate (% per generation)",
          min = 0, max = 5, value = 1.0, step = 0.1
        )
      )
    ),

    # --- Action buttons -----------------------------------------------------
    fluidRow(
      column(
        width = 12,
        tags$div(
          style = "margin-bottom: 10px;",
          actionButton(
            inputId = ns("simulate"),
            label   = "Simulate!",
            class   = "btn btn-primary btn-lg"
          ),
          actionButton(
            inputId = ns("show_answer"),
            label   = "Show Answer",
            class   = "btn btn-warning"
          )
        )
      )
    ),

    # --- Fit score display --------------------------------------------------
    fluidRow(
      column(
        width = 12,
        uiOutput(ns("fit_display"))
      )
    ),

    # --- Main visualization area --------------------------------------------
    fluidRow(
      # Left: Demography traversal
      column(
        width = 6,
        tags$h5("Population Size History", style = "font-weight: bold; margin-top: 10px;"),
        plotOutput(ns("demography_plot"), height = "350px")
      ),
      # Right: AFS comparison
      column(
        width = 6,
        tags$h5("Allele Frequency Spectrum Comparison", style = "font-weight: bold; margin-top: 10px;"),
        plotOutput(ns("afs_africa_plot"),    height = "160px"),
        plotOutput(ns("afs_nonafrica_plot"), height = "160px")
      )
    ),

    # --- Answer reveal panel ------------------------------------------------
    fluidRow(
      column(
        width = 12,
        uiOutput(ns("answer_panel"))
      )
    )
  )
}


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

#' Tab 2 Server
#'
#' @param id  Shiny module namespace id (must match tab2UI).
tab2Server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Pre-compute "real" AFS at startup with fixed seed
    real_afs <- isolate(
      simulate_afs(TRUE_T_SPLIT, TRUE_BOTTLENECK, TRUE_GROWTH_RATE,
                   n_snps = 200, n_sample = 50, seed = 42)
    )

    # Reactive values
    rv <- reactiveValues(
      sim_result   = NULL,
      show_answer  = FALSE,
      africa_score = NA,
      na_score     = NA
    )

    # -----------------------------------------------------------------------
    # Simulate on button click
    # -----------------------------------------------------------------------
    observeEvent(input$simulate, {
      rv$sim_result <- simulate_afs(
        t_split        = input$t_split,
        bottleneck_size = input$bottleneck_size,
        growth_rate    = input$growth_rate,
        n_snps         = 200,
        n_sample       = 50
      )
      rv$africa_score <- afs_fit_score(
        rv$sim_result$africa_counts, real_afs$africa_counts, 50
      )
      rv$na_score <- afs_fit_score(
        rv$sim_result$nonafrica_counts, real_afs$nonafrica_counts, 50
      )
    })

    # -----------------------------------------------------------------------
    # Show answer
    # -----------------------------------------------------------------------
    observeEvent(input$show_answer, {
      rv$show_answer <- !rv$show_answer
    })

    # -----------------------------------------------------------------------
    # Fit score display
    # -----------------------------------------------------------------------
    output$fit_display <- renderUI({
      req(!is.na(rv$africa_score))

      overall <- round((rv$africa_score + rv$na_score) / 2, 3)
      score_color <- function(s) {
        if (s >= 0.75) COLOR_GOOD
        else if (s >= 0.50) COLOR_OK
        else COLOR_BAD
      }

      tags$div(
        style = paste0(
          "padding: 8px 16px; border-radius: 4px; font-weight: bold; ",
          "background-color: ", score_color(overall), "22; ",
          "border: 2px solid ", score_color(overall), "; ",
          "margin-bottom: 8px;"
        ),
        tags$span(
          paste0(
            "Africa fit: ", round(rv$africa_score * 100), "%  |  ",
            "Non-Africa fit: ", round(rv$na_score * 100), "%  |  ",
            "Overall: ", round(overall * 100), "%"
          )
        )
      )
    })

    # -----------------------------------------------------------------------
    # Answer panel
    # -----------------------------------------------------------------------
    output$answer_panel <- renderUI({
      req(rv$show_answer)
      tags$div(
        style = "margin-top: 12px; padding: 12px; background: #f5f3f7; border-left: 4px solid #533860; border-radius: 4px;",
        tags$strong("True parameters:"),
        tags$ul(
          tags$li(paste("Time of split:", TRUE_T_SPLIT, "generations ago (~",
                        round(TRUE_T_SPLIT * 25 / 1000), "thousand years)")),
          tags$li(paste("Bottleneck size:", TRUE_BOTTLENECK, "individuals")),
          tags$li(paste("Growth rate:", TRUE_GROWTH_RATE, "% per generation"))
        )
      )
    })

    # -----------------------------------------------------------------------
    # Demography plot (updates reactively with sliders)
    # -----------------------------------------------------------------------
    output$demography_plot <- renderPlot({
      t_split   <- input$t_split
      bn_size   <- input$bottleneck_size
      gr        <- input$growth_rate
      t_bottle  <- t_split + 800
      max_time  <- t_split + 1200

      # Build Ne history
      ne_hist <- compute_ne_history(t_split, bn_size, gr, max_time = max_time)

      # Scale dots
      ne_max <- 10000
      max_dots <- 50
      ne_hist$dots_africa    <- ne_to_dots(ne_hist$ne_africa,    ne_max, max_dots)
      ne_hist$dots_nonafrica <- ne_to_dots(ne_hist$ne_nonafrica, ne_max, max_dots)

      # Expand: for each time point and each pop, create one row per dot
      expand_dots <- function(df, col, pop_name, color) {
        do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
          n <- df[[col]][i]
          if (n < 1) return(NULL)
          data.frame(
            gen_ago  = df$gen_ago[i],
            dot_idx  = seq_len(n),
            pop      = pop_name,
            color    = color,
            stringsAsFactors = FALSE
          )
        }))
      }

      # Before split: only one ancestral pop
      df_africa_post   <- ne_hist[ne_hist$gen_ago <= t_split, ]
      df_africa_pre    <- ne_hist[ne_hist$gen_ago >  t_split & ne_hist$gen_ago <= t_bottle, ]
      df_ancestral     <- ne_hist[ne_hist$gen_ago >  t_bottle, ]

      pts_africa_post  <- expand_dots(df_africa_post, "dots_africa", "Africa", COLOR_AFRICA)
      pts_nonafrica    <- expand_dots(ne_hist[ne_hist$gen_ago <= t_split, ],
                                      "dots_nonafrica", "Non-Africa", COLOR_NONAFRICA)
      # Offset non-Africa dots upward so they don't overlap Africa
      if (!is.null(pts_nonafrica) && nrow(pts_nonafrica) > 0)
        pts_nonafrica$dot_idx <- pts_nonafrica$dot_idx + max_dots + 2

      pts_bottleneck   <- expand_dots(df_africa_pre,  "dots_nonafrica", "Bottleneck", COLOR_NONAFRICA)
      if (!is.null(pts_bottleneck) && nrow(pts_bottleneck) > 0)
        pts_bottleneck$dot_idx <- pts_bottleneck$dot_idx + max_dots + 2

      pts_ancestral    <- expand_dots(df_ancestral,   "dots_africa", "Ancestral", COLOR_ANCESTRAL)

      all_pts <- do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0,
                                        list(pts_africa_post, pts_nonafrica,
                                             pts_bottleneck, pts_ancestral)))

      if (is.null(all_pts) || nrow(all_pts) == 0) {
        ggplot() + theme_void() + labs(title = "Adjust sliders to see demography")
      } else {
        ggplot(all_pts, aes(x = -gen_ago, y = dot_idx, color = I(color))) +
          geom_point(size = 1.8, alpha = 0.7) +
          geom_vline(xintercept = -t_split,  linetype = "dashed",
                     color = "grey40", linewidth = 0.6) +
          geom_vline(xintercept = -t_bottle, linetype = "dashed",
                     color = "grey60", linewidth = 0.4) +
          annotate("text", x = -t_split,  y = max_dots * 2.1,
                   label = "Split", hjust = 0.5, size = 3, color = "grey30") +
          annotate("text", x = -t_bottle, y = max_dots * 2.1,
                   label = "Bottleneck\nstart", hjust = 0.5, size = 2.5, color = "grey40") +
          scale_x_continuous(
            name   = "Generations ago",
            breaks = seq(-max_time, 0, by = 500),
            labels = function(x) abs(x)
          ) +
          labs(
            title = paste0("Demography (split: ", t_split, " gen, ",
                           "bottleneck: ", bn_size, ", growth: ", gr, "%)"),
            y     = NULL
          ) +
          theme_minimal(base_size = 11) +
          theme(
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor   = element_blank(),
            plot.title = element_text(size = 10, face = "bold")
          )
      }
    })

    # -----------------------------------------------------------------------
    # AFS plots
    # -----------------------------------------------------------------------
    make_afs_plot <- function(sim_counts, real_counts, pop_name, color, n_sample = 50) {
      bins <- 1:(n_sample - 1)   # exclude 0 (monomorphic) and n_sample (fixed)

      sim_hist  <- tabulate(sim_counts  + 1L, nbins = n_sample + 1)[-c(1, n_sample + 1)]
      real_hist <- tabulate(real_counts + 1L, nbins = n_sample + 1)[-c(1, n_sample + 1)]
      n_snps_sim  <- length(sim_counts)
      n_snps_real <- length(real_counts)

      df <- data.frame(
        count = rep(bins, 2),
        freq  = c(sim_hist  / n_snps_sim,
                  real_hist / n_snps_real),
        type  = rep(c("Simulated", "Real"), each = length(bins))
      )

      ggplot(df, aes(x = count, y = freq, fill = type, alpha = type)) +
        geom_col(position = "identity") +
        scale_fill_manual(values = c("Real" = color, "Simulated" = "#888888")) +
        scale_alpha_manual(values = c("Real" = 0.85, "Simulated" = 0.55)) +
        labs(
          title = paste(pop_name, "AFS"),
          x     = "Derived allele count",
          y     = "Proportion of SNPs",
          fill  = NULL, alpha = NULL
        ) +
        theme_minimal(base_size = 10) +
        theme(
          legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"),
          plot.title      = element_text(size = 10, face = "bold")
        )
    }

    output$afs_africa_plot <- renderPlot({
      req(rv$sim_result)
      make_afs_plot(
        rv$sim_result$africa_counts,
        real_afs$africa_counts,
        "Africa",
        COLOR_AFRICA
      )
    })

    output$afs_nonafrica_plot <- renderPlot({
      req(rv$sim_result)
      make_afs_plot(
        rv$sim_result$nonafrica_counts,
        real_afs$nonafrica_counts,
        "Non-Africa",
        COLOR_NONAFRICA
      )
    })

  })
}
