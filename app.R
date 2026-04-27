# app.R
# USF 2026 Population Genetics Shiny App
#
# Tab 1: Interactive Wright-Fisher genetic drift simulator
# Tab 2: Out of Africa demographic history explorer
# Tab 3: Estimating the mutation rate

library(shiny)
library(ggplot2)
library(dplyr)
library(shinyjs)

# Source helper modules
source("R/wf_simulation.R")
source("R/ooa_simulation.R")
source("R/mutation_simulation.R")
source("R/tab1_module.R")
source("R/tab2_module.R")
source("R/tab3_module.R")

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- navbarPage(
  title = "Population Genetics Explorer",

  # Custom styling
  header = tags$head(
    tags$style(HTML("
      .navbar { background-color: #533860 !important; }
      .navbar-brand, .navbar-nav > li > a { color: #FFE42A !important; font-weight: bold; }
      .navbar-nav > li.active > a { background-color: #412d5e !important; color: #FFE42A !important; }
      body { font-family: Georgia, serif; }
      h4, h5 { font-family: 'Segoe UI', sans-serif; color: #533860; }
      .btn-primary { background-color: #533860; border-color: #412d5e; }
      .btn-primary:hover { background-color: #412d5e; }
    "))
  ),

  # Tab 1: Wright-Fisher Simulator
  tabPanel(
    title = "\U1F3B2 Wright-Fisher Model",
    fluidPage(
      fluidRow(
        column(
          width = 12,
          tags$h4("Interactive Wright-Fisher Genetic Drift Simulator"),
          tags$p(
            "Build a new generation one individual at a time. ",
            "Each individual randomly inherits an allele from the previous generation.",
            style = "color: #555; font-size: 0.95em;"
          )
        )
      ),
      tab1UI("tab1")
    )
  ),

  # Tab 2: Out of Africa Explorer
  tabPanel(
    title = "\U1F30D Out of Africa",
    fluidPage(
      fluidRow(
        column(
          width = 12,
          tags$h4("Out of Africa: Estimate When Humans Left Africa"),
          tags$p(
            "Tune the demographic parameters below to match the simulated Allele Frequency Spectrum ",
            "to the 'real' data. When your fit score is high, you've found the answer!",
            style = "color: #555; font-size: 0.95em;"
          )
        )
      ),
      tab2UI("tab2")
    )
  ),

  # Tab 3: Estimating Mutation Rate
  tabPanel(
    title = "\U1F9EC Mutation Rate",
    fluidPage(
      fluidRow(
        column(
          width = 12,
          tags$h4("Estimating the Mutation Rate"),
          tags$p(
            "The histogram shows allele frequencies from 1,000 independent populations ",
            "at mutation-drift steady state. Move the μ slider until the orange ",
            "theoretical curve matches the histogram shape. ",
            "When the fit score turns green, you’ve found the answer!",
            style = "color: #555; font-size: 0.95em;"
          )
        )
      ),
      tab3UI("tab3")
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {
  tab1Server("tab1")
  tab2Server("tab2")
  tab3Server("tab3")
}

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------

shinyApp(ui = ui, server = server)
