# tab3_module.R
# Shiny module for Tab 3: Estimating the Mutation Rate.
#
# Students adjust a log-scale mutation-rate slider until a theoretical Beta
# curve matches a pre-computed steady-state allele-frequency histogram.
#
# Exports:
#   tab3UI(id)     — UI function
#   tab3Server(id) — Server function
#
# Depends on: mutation_simulation.R (must be sourced before this module)

library(shiny)
library(ggplot2)

# ---------------------------------------------------------------------------
# Parameters  (N is internal; students never see it)
# ---------------------------------------------------------------------------
TAB3_N     <- 1e6     # effective population size — sets the scale of θ = 4Nμ
TAB3_MU    <- 1e-8    # true mutation rate (hidden from students)
TAB3_NPOPS <- 1000L   # populations sampled for the histogram

# Colours
COLOR_CURVE <- "#E67E22"
COLOR_GOOD  <- "#27AE60"
COLOR_OK    <- "#F39C12"
COLOR_BAD   <- "#E74C3C"

# Helper: format μ as HTML scientific notation, e.g. "3.16 × 10<sup>-9</sup>"
fmt_mu_html <- function(mu) {
  if (!is.finite(mu) || mu <= 0) return("0")
  e <- floor(log10(mu))
  m <- mu / 10^e
  sprintf("%.2g &times; 10<sup>%d</sup>", m, e)
}


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

#' Tab 3 UI
#'
#' @param id  Shiny module namespace id.
tab3UI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      # --- Left: controls ---------------------------------------------------
      column(
        width = 4,

        sliderInput(
          inputId = ns("mu_log10"),
          label   = HTML("Mutation rate &mu; (log<sub>10</sub> scale)"),
          min     = -10,
          max     = -7,
          value   = -9,
          step    = 0.1
        ),

        # Live readout of the actual μ value in scientific notation
        uiOutput(ns("mu_display")),

        uiOutput(ns("fit_display")),

        tags$div(
          style = "margin-top:14px;",
          actionButton(ns("show_answer"), "Show Answer", class = "btn btn-warning")
        ),

        uiOutput(ns("answer_panel"))
      ),

      # --- Right: plot -------------------------------------------------------
      column(
        width = 8,
        plotOutput(ns("main_plot"), height = "420px")
      )
    )
  )
}


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

#' Tab 3 Server
#'
#' @param id  Shiny module namespace id (must match tab3UI).
tab3Server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Pre-compute steady-state histogram data once at startup (fixed seed)
    steady_freqs <- isolate(
      simulate_steady_state(TAB3_N, TAB3_MU, TAB3_NPOPS, seed = 42)
    )

    rv <- reactiveValues(show_answer = FALSE)

    observeEvent(input$show_answer, {
      rv$show_answer <- !rv$show_answer
    })

    # Convenience reactive: the actual μ from the log-scale slider
    mu_r <- reactive({ 10^input$mu_log10 })

    # -----------------------------------------------------------------------
    # μ readout
    # -----------------------------------------------------------------------
    output$mu_display <- renderUI({
      tags$div(
        style = paste0("font-size:1.3em; font-weight:bold; color:#533860;",
                       "text-align:center; margin: 4px 0 14px;"),
        HTML(paste0("μ = ", fmt_mu_html(mu_r())))
      )
    })

    # -----------------------------------------------------------------------
    # Fit score badge
    # -----------------------------------------------------------------------
    output$fit_display <- renderUI({
      score <- mutation_fit_score(steady_freqs, TAB3_N, mu_r())
      pct   <- round(score * 100)
      color <- if (score >= 0.90) COLOR_GOOD
               else if (score >= 0.70) COLOR_OK
               else COLOR_BAD
      tags$div(
        style = paste0(
          "padding:8px 14px; border-radius:4px; font-weight:bold; ",
          "background-color:", color, "22; border:2px solid ", color, ";"
        ),
        paste0("Fit score: ", pct, "%")
      )
    })

    # -----------------------------------------------------------------------
    # Answer reveal
    # -----------------------------------------------------------------------
    output$answer_panel <- renderUI({
      req(rv$show_answer)
      tags$div(
        style = paste0("margin-top:14px; padding:12px; background:#f5f3f7; ",
                       "border-left:4px solid #533860; border-radius:4px;"),
        tags$strong("True mutation rate:"),
        tags$p(
          HTML(paste0("μ = ", fmt_mu_html(TAB3_MU))),
          style = "font-size:1.1em; margin:6px 0 0;"
        )
      )
    })

    # -----------------------------------------------------------------------
    # Main plot
    # -----------------------------------------------------------------------
    output$main_plot <- renderPlot({
      mu    <- mu_r()
      alpha <- 4 * TAB3_N * mu

      # Theoretical density — grid avoids exact boundaries (prevents ±∞ when α < 1)
      p_grid     <- seq(0.005, 0.995, length.out = 400)
      dens_vals  <- dbeta(p_grid, alpha, alpha)
      dens_vals[!is.finite(dens_vals)] <- 0
      density_df <- data.frame(p = p_grid, density = dens_vals)

      freq_df  <- data.frame(freq = steady_freqs)
      hist_obj <- hist(steady_freqs, breaks = 30, plot = FALSE)
      y_ceil   <- max(hist_obj$density, na.rm = TRUE) * 1.5

      ggplot() +
        geom_histogram(
          data      = freq_df,
          aes(x = freq, y = after_stat(density)),
          bins      = 100,
          fill      = "#533860",
          alpha     = 0.65,
          color     = "white",
          linewidth = 0.3
        ) +
        geom_line(
          data      = density_df,
          aes(x = p, y = density),
          color     = COLOR_CURVE,
          linewidth = 1.8
        ) +
        coord_cartesian(ylim = c(0, y_ceil)) +
        scale_x_continuous(
          limits = c(0, 1),
          labels = function(x) paste0(round(x * 100), "%"),
          name   = "Luigi allele frequency"
        ) +
        labs(
          y     = "Density",
          title = sprintf("Steady-state allele frequencies   μ = %.2e", mu)
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title       = element_text(hjust = 0.5, size = 12,
                                          face = "bold", color = "#533860"),
          panel.grid.minor = element_blank()
        )
    })

  })
}
