# tab1_module.R
# Shiny module for Tab 1: Interactive Wright-Fisher Genetic Drift Simulator.
#
# Exports:
#   tab1UI(id)        — UI function
#   tab1Server(id)    — Server function
#
# Depends on: wf_simulation.R (must be sourced before this module)
# Depends on: shinyjs (for enable/disable buttons)

library(shiny)
library(ggplot2)
library(shinyjs)

# ---------------------------------------------------------------------------
# COI brand colors
# ---------------------------------------------------------------------------
COI_RED   <- "#C0392B"   # Mario / allele 0
COI_GREEN <- "#27AE60"   # Luigi / allele 1
COI_GREY  <- "#CCCCCC"   # Unfilled slot


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

#' Tab 1 UI
#'
#' @param id  Shiny module namespace id.
tab1UI <- function(id) {
  ns <- NS(id)

  tagList(
    useShinyjs(),

    # --- Controls row -------------------------------------------------------
    fluidRow(
      column(
        width = 3,
        numericInput(
          inputId  = ns("N"),
          label    = "Population size (N)",
          value    = 20,
          min      = 5,
          max      = 100,
          step     = 1
        )
      ),
      column(
        width = 9,
        tags$div(
          style = "margin-top: 25px;",
          actionButton(
            inputId = ns("roll_one"),
            label   = "\U1F3B2 Roll One Individual",
            class   = "btn btn-primary"
          ),
          actionButton(
            inputId = ns("complete_gen"),
            label   = "Complete Generation",
            class   = "btn btn-warning"
          ),
          actionButton(
            inputId = ns("new_gen"),
            label   = "New Generation",
            class   = "btn btn-success"
          ),
          actionButton(
            inputId = ns("reset"),
            label   = "Reset",
            class   = "btn btn-danger"
          )
        )
      )
    ),

    # --- Status text --------------------------------------------------------
    fluidRow(
      column(
        width = 12,
        tags$div(
          style = "margin-top: 8px; margin-bottom: 4px; font-size: 1.1em;",
          textOutput(ns("status_text"))
        )
      )
    ),

    # --- Plot ---------------------------------------------------------------
    fluidRow(
      column(
        width = 12,
        plotOutput(ns("wf_plot"), height = "400px")
      )
    )
  )
}


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

#' Tab 1 Server
#'
#' @param id  Shiny module namespace id (must match tab1UI).
tab1Server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # -----------------------------------------------------------------------
    # Reactive state
    # -----------------------------------------------------------------------
    rv <- reactiveValues(
      N            = 20,
      current_gen  = NULL,   # complete previous generation (0/1 vector)
      building_gen = NULL,   # partial next generation being assembled
      gen_number   = 1
    )

    # Helper: initialize / reset everything
    do_reset <- function(N = rv$N) {
      rv$N            <- N
      rv$current_gen  <- init_generation(N, p0 = 0.5)
      rv$building_gen <- integer(0)
      rv$gen_number   <- 1
    }

    isolate(do_reset(20))

    # -----------------------------------------------------------------------
    # Reactive expressions
    # -----------------------------------------------------------------------

    n_filled <- reactive({
      length(rv$building_gen)
    })

    gen_complete <- reactive({
      n_filled() >= rv$N
    })

    # -----------------------------------------------------------------------
    # Button enable / disable
    # -----------------------------------------------------------------------
    observe({
      complete <- gen_complete()

      if (complete) {
        shinyjs::disable("roll_one")
        shinyjs::disable("complete_gen")
        shinyjs::enable("new_gen")
      } else {
        shinyjs::enable("roll_one")
        shinyjs::enable("complete_gen")
        shinyjs::disable("new_gen")
      }
    })

    # -----------------------------------------------------------------------
    # Events
    # -----------------------------------------------------------------------

    observeEvent(input$N, {
      req(input$N)
      new_N <- as.integer(input$N)
      if (!is.na(new_N) && new_N >= 5 && new_N <= 100) {
        do_reset(new_N)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$roll_one, {
      req(!gen_complete())
      new_ind <- sample_one_individual(rv$current_gen)
      rv$building_gen <- c(rv$building_gen, new_ind)
    })

    observeEvent(input$complete_gen, {
      req(!gen_complete())
      remaining <- rv$N - n_filled()
      if (remaining > 0) {
        new_inds <- sample(rv$current_gen, remaining, replace = TRUE)
        rv$building_gen <- c(rv$building_gen, new_inds)
      }
    })

    observeEvent(input$new_gen, {
      req(gen_complete())
      rv$current_gen  <- rv$building_gen
      rv$building_gen <- integer(0)
      rv$gen_number   <- rv$gen_number + 1
    })

    observeEvent(input$reset, {
      N_val <- as.integer(input$N)
      if (is.na(N_val) || N_val < 5 || N_val > 100) N_val <- 20
      do_reset(N_val)
    })

    # -----------------------------------------------------------------------
    # Status text
    # -----------------------------------------------------------------------
    output$status_text <- renderText({
      filled <- n_filled()
      N      <- rv$N
      gen    <- rv$gen_number

      if (filled == 0) {
        prop_luigi <- mean(rv$current_gen)
        pct        <- round(prop_luigi * 100, 1)
        sprintf(
          "Generation %d  |  Ready to sample: 0/%d  |  Previous Luigi frequency: %.1f%%",
          gen, N, pct
        )
      } else {
        prop_luigi <- mean(rv$building_gen)
        pct        <- round(prop_luigi * 100, 1)
        sprintf(
          "Generation %d  |  Sampled so far: %d/%d  |  Luigi frequency: %.1f%%",
          gen, filled, N, pct
        )
      }
    })

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    output$wf_plot <- renderPlot({
      N      <- rv$N
      filled <- n_filled()
      gen    <- rv$gen_number

      # Show current completed gen when building_gen is empty,
      # otherwise show partial building_gen with grey slots
      if (filled == 0) {
        display <- rv$current_gen
      } else {
        display <- rep(NA_integer_, N)
        display[seq_len(filled)] <- rv$building_gen
      }

      if (filled == 0) {
        prop_luigi <- mean(rv$current_gen)
        pct        <- round(prop_luigi * 100, 1)
        title_str  <- sprintf(
          "Generation %d \u2014 0/%d sampled (Previous: %.1f%% Luigi)",
          gen, N, pct
        )
      } else {
        prop_luigi <- mean(rv$building_gen)
        pct        <- round(prop_luigi * 100, 1)
        title_str  <- sprintf(
          "Generation %d \u2014 %d/%d sampled (%.1f%% Luigi)",
          gen, filled, N, pct
        )
      }

      cols_per_row <- 10
      df <- data.frame(
        idx    = seq_len(N),
        allele = display,
        col    = ((seq_len(N) - 1) %% cols_per_row) + 1,
        row    = -floor((seq_len(N) - 1) / cols_per_row)
      )

      df$point_fill    <- ifelse(is.na(df$allele), "white",
                           ifelse(df$allele == 1, COI_GREEN, COI_RED))
      df$outline_color <- ifelse(is.na(df$allele), COI_GREY,
                           ifelse(df$allele == 1, "#1E8449", "#922B21"))
      df$alpha_val     <- ifelse(is.na(df$allele), 0.5, 1.0)

      ggplot(df, aes(x = col, y = row)) +
        geom_point(
          aes(fill  = I(point_fill),
              color = I(outline_color),
              alpha = I(alpha_val)),
          shape  = 21,
          size   = 7,
          stroke = 1.5
        ) +
        scale_x_continuous(limits = c(0.5, cols_per_row + 0.5)) +
        labs(title = title_str, x = NULL, y = NULL) +
        theme_void(base_size = 14) +
        theme(
          plot.title       = element_text(hjust = 0.5, size = 14, face = "bold",
                                          margin = margin(b = 12)),
          plot.margin      = margin(10, 20, 10, 20),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA)
        )
    })

  })
}
