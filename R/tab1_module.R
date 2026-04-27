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

    tags$head(
      tags$style(HTML("
        .die-wrap  { text-align:center; padding:8px 0; min-height:110px; }
        .die-face  {
          font-size:52px; font-weight:bold; font-family:monospace;
          display:inline-flex; align-items:center; justify-content:center;
          width:84px; height:84px;
          border:3px solid #533860; border-radius:12px;
          background:white;
          box-shadow:4px 4px 10px rgba(0,0,0,0.25);
        }
        .die-label { font-size:1em; font-weight:bold; margin-top:6px; }
        .coin-wrap  { text-align:center; padding:8px 0; min-height:60px; }
        .coin-face  { font-size:40px; font-weight:bold; display:inline-block; }
        .coin-label { font-size:0.95em; margin-top:4px; color:#555; }
      "))
    ),

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
          actionButton(ns("roll_one"),     "\U1F3B2 Roll One Individual", class = "btn btn-primary"),
          actionButton(ns("complete_gen"), "Complete Generation",         class = "btn btn-warning"),
          actionButton(ns("new_gen"),      "New Generation",              class = "btn btn-success"),
          actionButton(ns("reset"),        "Reset",                       class = "btn btn-danger")
        )
      )
    ),

    # --- Die display --------------------------------------------------------
    fluidRow(
      column(width = 12, htmlOutput(ns("die_display")))
    ),

    # --- Mutation section ---------------------------------------------------
    fluidRow(
      column(
        width = 12,
        checkboxInput(ns("mutation_on"), "Enable Mutation", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s']", ns("mutation_on")),
          fluidRow(
            column(
              width = 4,
              sliderInput(ns("mutation_prob"), "Mutation probability",
                          min = 0, max = 0.5, value = 0.01, step = 0.001)
            ),
            column(
              width = 8,
              tags$div(
                style = "margin-top: 25px;",
                actionButton(ns("sim_mutation"), "Simulate Single Mutation",
                             class = "btn btn-info")
              )
            )
          ),
          htmlOutput(ns("coin_display"))
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

    # --- Main visualization: scrollable tree + frequency plot ---------------
    fluidRow(
      column(
        width = 8,
        tags$h5("Wright-Fisher Tree",
                style = "font-weight:bold; color:#533860; margin-top:10px;"),
        tags$div(
          id    = ns("wf-plot-container"),
          style = "height:520px; overflow-y:auto; border:1px solid #dee2e6; border-radius:4px;",
          uiOutput(ns("wf_plot_ui"))
        )
      ),
      column(
        width = 4,
        tags$h5("Luigi Allele Frequency",
                style = "font-weight:bold; color:#533860; margin-top:10px;"),
        plotOutput(ns("freq_plot"), height = "520px")
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
      N                  = 20,
      current_gen        = NULL,
      building_gen       = NULL,
      gen_number         = 1,
      all_gens           = NULL,    # list: all_gens[[k]] = allele vector for gen k
      all_parent_indices = NULL,    # list: all_parent_indices[[k]] = parent idx for gen k→k+1
      die_index          = NULL,
      die_allele         = NULL,
      parent_indices     = integer(0),
      pending_mut        = FALSE,
      coin_result        = NULL
    )

    do_reset <- function(N = rv$N) {
      gen1                  <- init_generation(N, p0 = 0.5)
      rv$N                  <- N
      rv$current_gen        <- gen1
      rv$all_gens           <- list(gen1)
      rv$all_parent_indices <- list()
      rv$building_gen       <- integer(0)
      rv$gen_number         <- 1
      rv$die_index          <- NULL
      rv$die_allele         <- NULL
      rv$parent_indices     <- integer(0)
      rv$pending_mut        <- FALSE
      rv$coin_result        <- NULL
    }

    isolate(do_reset(20))

    # -----------------------------------------------------------------------
    # Reactive expressions
    # -----------------------------------------------------------------------

    n_filled <- reactive({ length(rv$building_gen) })

    gen_complete <- reactive({ n_filled() >= rv$N })

    # -----------------------------------------------------------------------
    # Button enable / disable
    # -----------------------------------------------------------------------
    observe({
      complete <- gen_complete()
      pending  <- rv$pending_mut
      mut_on   <- isTRUE(input$mutation_on)

      if (complete || pending) shinyjs::disable("roll_one")     else shinyjs::enable("roll_one")
      if (complete || pending) shinyjs::disable("complete_gen") else shinyjs::enable("complete_gen")
      if (complete && !pending) shinyjs::enable("new_gen")      else shinyjs::disable("new_gen")
      if (mut_on && pending)    shinyjs::enable("sim_mutation") else shinyjs::disable("sim_mutation")
    })

    # -----------------------------------------------------------------------
    # Helper: scroll tree container to bottom after plot updates
    # -----------------------------------------------------------------------
    scroll_bottom <- function() {
      shinyjs::runjs(sprintf("
        setTimeout(function(){
          var c = document.getElementById('%s');
          if (c) c.scrollTop = c.scrollHeight;
        }, 500);
      ", session$ns("wf-plot-container")))
    }

    # -----------------------------------------------------------------------
    # Events
    # -----------------------------------------------------------------------

    observeEvent(input$N, {
      req(input$N)
      new_N <- as.integer(input$N)
      if (!is.na(new_N) && new_N >= 5 && new_N <= 100) do_reset(new_N)
    }, ignoreInit = TRUE)

    observeEvent(input$roll_one, {
      req(!gen_complete(), !rv$pending_mut)

      idx               <- sample(length(rv$current_gen), 1L)
      new_ind           <- rv$current_gen[idx]
      rv$building_gen   <- c(rv$building_gen, new_ind)
      rv$parent_indices <- c(rv$parent_indices, idx)
      rv$die_index      <- idx
      rv$die_allele     <- new_ind
      rv$coin_result    <- NULL

      js <- sprintf("
        (function(){
          var el = document.getElementById('%s');
          if (!el) return;
          var max = %d;
          var iv = setInterval(function(){
            el.textContent = Math.floor(Math.random() * max) + 1;
          }, 80);
          setTimeout(function(){ clearInterval(iv); el.textContent = %d; }, 1200);
        })();
      ", session$ns("die-face-el"), rv$N, idx)
      shinyjs::runjs(js)

      if (isTRUE(input$mutation_on)) rv$pending_mut <- TRUE

      scroll_bottom()
    })

    observeEvent(input$sim_mutation, {
      req(rv$pending_mut, isTRUE(input$mutation_on))

      mutated        <- runif(1) < input$mutation_prob
      result         <- if (mutated) "MUTATION" else "NO MUTATION"
      rv$coin_result <- result

      if (mutated && length(rv$building_gen) > 0) {
        n                  <- length(rv$building_gen)
        rv$building_gen[n] <- 1L - rv$building_gen[n]
        rv$die_allele      <- rv$building_gen[n]
      }

      rv$pending_mut <- FALSE

      js <- sprintf("
        (function(){
          var sides = ['MUTATION','NO MUTATION'];
          var el = document.getElementById('%s');
          if (!el) return;
          var i = 0;
          var iv = setInterval(function(){ el.textContent = sides[i++ %% 2]; }, 120);
          setTimeout(function(){ clearInterval(iv); el.textContent = '%s'; }, 1400);
        })();
      ", session$ns("coin-face-el"), result)
      shinyjs::runjs(js)
    })

    observeEvent(input$complete_gen, {
      req(!gen_complete(), !rv$pending_mut)
      remaining <- rv$N - n_filled()
      if (remaining > 0) {
        new_idxs <- sample(length(rv$current_gen), remaining, replace = TRUE)
        new_inds <- rv$current_gen[new_idxs]
        if (isTRUE(input$mutation_on)) {
          mask           <- runif(remaining) < input$mutation_prob
          new_inds[mask] <- 1L - new_inds[mask]
        }
        rv$building_gen   <- c(rv$building_gen, new_inds)
        rv$parent_indices <- c(rv$parent_indices, new_idxs)
      }
      scroll_bottom()
    })

    observeEvent(input$new_gen, {
      req(gen_complete())
      rv$all_parent_indices <- c(rv$all_parent_indices, list(rv$parent_indices))
      rv$all_gens           <- c(rv$all_gens, list(rv$building_gen))
      rv$current_gen        <- rv$building_gen
      rv$building_gen       <- integer(0)
      rv$parent_indices     <- integer(0)
      rv$die_index          <- NULL
      rv$die_allele         <- NULL
      rv$coin_result        <- NULL
      rv$gen_number         <- rv$gen_number + 1
      scroll_bottom()
    })

    observeEvent(input$reset, {
      N_val <- as.integer(input$N)
      if (is.na(N_val) || N_val < 5 || N_val > 100) N_val <- 20
      do_reset(N_val)
    })

    # -----------------------------------------------------------------------
    # Die display
    # -----------------------------------------------------------------------
    output$die_display <- renderUI({
      req(!is.null(rv$die_index))
      allele_label <- if (rv$die_allele == 1) "Luigi (green)" else "Mario (red)"
      allele_color <- if (rv$die_allele == 1) COI_GREEN else COI_RED
      tags$div(
        class = "die-wrap",
        tags$span(id = session$ns("die-face-el"), class = "die-face", rv$die_index),
        tags$div(
          class = "die-label",
          sprintf("Individual #%d  →  ", rv$die_index),
          tags$span(allele_label, style = paste0("color:", allele_color, ";"))
        )
      )
    })

    # -----------------------------------------------------------------------
    # Coin display
    # -----------------------------------------------------------------------
    output$coin_display <- renderUI({
      req(!is.null(rv$coin_result))
      color <- if (rv$coin_result == "MUTATION") "#E74C3C" else COI_GREEN
      tags$div(
        class = "coin-wrap",
        tags$span(
          id    = session$ns("coin-face-el"),
          class = "coin-face",
          rv$coin_result,
          style = paste0("color:", color, ";")
        ),
        tags$div(class = "coin-label", "Coin flip result")
      )
    })

    # -----------------------------------------------------------------------
    # Status text
    # -----------------------------------------------------------------------
    output$status_text <- renderText({
      filled <- n_filled()
      N      <- rv$N
      gen    <- rv$gen_number
      if (filled == 0) {
        pct <- round(mean(rv$current_gen) * 100, 1)
        sprintf("Generation %d  |  Ready to sample: 0/%d  |  Luigi frequency: %.1f%%",
                gen, N, pct)
      } else {
        pct <- round(mean(rv$building_gen) * 100, 1)
        sprintf("Generation %d  |  Sampled so far: %d/%d  |  Luigi frequency: %.1f%%",
                gen, filled, N, pct)
      }
    })

    # -----------------------------------------------------------------------
    # Dynamic plot height
    # -----------------------------------------------------------------------
    output$wf_plot_ui <- renderUI({
      n_gens_total <- length(rv$all_gens) + 1L
      n_rows       <- ceiling(rv$N / 10L)
      block_height <- n_rows + 2L
      px_per_unit  <- 60L
      height_px    <- max(350L, (n_gens_total * block_height + 1L) * px_per_unit)
      plotOutput(session$ns("wf_plot"), height = paste0(height_px, "px"))
    })

    # -----------------------------------------------------------------------
    # Wright-Fisher tree plot (all generations)
    # -----------------------------------------------------------------------
    output$wf_plot <- renderPlot({
      N      <- rv$N
      filled <- n_filled()
      gen    <- rv$gen_number

      cols_per_row  <- 10L
      n_rows        <- ceiling(N / cols_per_row)
      block_height  <- n_rows + 2L
      n_gens_stored <- length(rv$all_gens)   # equals gen_number

      col_idx <- ((seq_len(N) - 1L) %% cols_per_row) + 1L
      row_idx <- -floor((seq_len(N) - 1L) / cols_per_row)

      # --- Dot data for every completed generation -------------------------
      dot_dfs <- lapply(seq_len(n_gens_stored), function(k) {
        data.frame(
          allele = rv$all_gens[[k]],
          col    = col_idx,
          row    = -(k - 1L) * block_height + row_idx
        )
      })

      # --- Dot data for the building generation ----------------------------
      offspring_alleles <- rep(NA_integer_, N)
      if (filled > 0L) offspring_alleles[seq_len(filled)] <- rv$building_gen
      dot_dfs <- c(dot_dfs, list(data.frame(
        allele = offspring_alleles,
        col    = col_idx,
        row    = -n_gens_stored * block_height + row_idx
      )))

      df <- do.call(rbind, dot_dfs)
      df$point_fill    <- ifelse(is.na(df$allele), "white",
                           ifelse(df$allele == 1, COI_GREEN, COI_RED))
      df$outline_color <- ifelse(is.na(df$allele), COI_GREY,
                           ifelse(df$allele == 1, "#1E8449", "#922B21"))
      df$alpha_val     <- ifelse(is.na(df$allele), 0.5, 1.0)

      # --- Arrow data for all completed transitions ------------------------
      arrow_dfs <- lapply(seq_along(rv$all_parent_indices), function(k) {
        pidx  <- rv$all_parent_indices[[k]]
        if (length(pidx) == 0L) return(NULL)
        j_seq <- seq_len(length(pidx))
        data.frame(
          x    = ((j_seq - 1L) %% cols_per_row) + 1L,
          y    = -k * block_height        - floor((j_seq - 1L) / cols_per_row),
          xend = ((pidx  - 1L) %% cols_per_row) + 1L,
          yend = -(k - 1L) * block_height - floor((pidx  - 1L) / cols_per_row)
        )
      })

      # --- Arrow data for the current in-progress transition ---------------
      pidx <- rv$parent_indices
      if (filled > 0L && length(pidx) == filled) {
        j_seq <- seq_len(filled)
        arrow_dfs <- c(arrow_dfs, list(data.frame(
          x    = ((j_seq - 1L) %% cols_per_row) + 1L,
          y    = -n_gens_stored * block_height        - floor((j_seq - 1L) / cols_per_row),
          xend = ((pidx  - 1L) %% cols_per_row) + 1L,
          yend = -(n_gens_stored - 1L) * block_height - floor((pidx  - 1L) / cols_per_row)
        )))
      }

      arrow_df_all <- do.call(rbind, Filter(Negate(is.null), arrow_dfs))
      arrow_layer  <- if (!is.null(arrow_df_all) && nrow(arrow_df_all) > 0) {
        geom_segment(
          data        = arrow_df_all,
          mapping     = aes(x = x, y = y, xend = xend, yend = yend),
          arrow       = arrow(length = unit(0.12, "cm"), type = "closed"),
          color       = "grey50",
          alpha       = 0.55,
          linewidth   = 0.35,
          inherit.aes = FALSE
        )
      } else NULL

      # --- Generation labels -----------------------------------------------
      gen_labels_df <- data.frame(
        x     = 0.5,
        y     = -(seq_len(n_gens_stored + 1L) - 1L) * block_height + 0.65,
        label = paste0("Gen ", seq_len(n_gens_stored + 1L))
      )

      # --- Title -----------------------------------------------------------
      if (filled == 0L) {
        pct       <- round(mean(rv$current_gen) * 100, 1)
        title_str <- sprintf("Generation %d (%.1f%% Luigi)  —  start rolling Generation %d",
                             gen, pct, gen + 1L)
      } else {
        pct       <- round(mean(rv$building_gen) * 100, 1)
        title_str <- sprintf("Sampling Generation %d: %d/%d drawn (%.1f%% Luigi)",
                             gen + 1L, filled, N, pct)
      }

      ggplot(df, aes(x = col, y = row)) +
        arrow_layer +
        geom_point(
          aes(fill  = I(point_fill),
              color = I(outline_color),
              alpha = I(alpha_val)),
          shape = 21, size = 6, stroke = 1.5
        ) +
        geom_text(
          data        = gen_labels_df,
          mapping     = aes(x = x, y = y, label = label),
          hjust       = 0,
          size        = 3,
          color       = "#533860",
          fontface    = "bold",
          inherit.aes = FALSE
        ) +
        scale_x_continuous(limits = c(0.5, cols_per_row + 0.5)) +
        labs(title = title_str, x = NULL, y = NULL) +
        theme_void(base_size = 14) +
        theme(
          plot.title       = element_text(hjust = 0.5, size = 13, face = "bold",
                                          margin = margin(b = 12)),
          plot.margin      = margin(10, 20, 10, 20),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA)
        )
    })

    # -----------------------------------------------------------------------
    # Allele frequency history plot
    # -----------------------------------------------------------------------
    output$freq_plot <- renderPlot({
      n_stored <- length(rv$all_gens)
      freq_df  <- data.frame(
        gen  = seq_len(n_stored),
        freq = sapply(rv$all_gens, mean)
      )

      filled  <- n_filled()
      in_prog <- if (filled > 0L) {
        data.frame(gen = n_stored + 1L, freq = mean(rv$building_gen))
      } else NULL

      p <- ggplot(freq_df, aes(x = gen, y = freq)) +
        geom_hline(yintercept = 0.5, linetype = "dashed",
                   color = "grey60", linewidth = 0.4) +
        geom_line(color = "#533860", linewidth = 1) +
        geom_point(color = "#533860", fill = "white",
                   shape = 21, size = 3, stroke = 1.5) +
        scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, 0.25),
          labels = function(x) paste0(round(x * 100), "%")
        ) +
        scale_x_continuous(
          breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1)
        ) +
        labs(title = "Luigi frequency over time", x = "Generation", y = NULL) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title       = element_text(hjust = 0.5, size = 12,
                                          face = "bold", color = "#533860"),
          panel.grid.minor = element_blank()
        )

      if (!is.null(in_prog)) {
        connector <- rbind(freq_df[nrow(freq_df), , drop = FALSE], in_prog)
        p <- p +
          geom_line(data = connector, linetype = "dashed",
                    color = "grey50", linewidth = 0.8) +
          geom_point(data = in_prog, color = "#533860", fill = "white",
                     shape = 21, size = 3, stroke = 1.5)
      }

      p
    })

  })
}
