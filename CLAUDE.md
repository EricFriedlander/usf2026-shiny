# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the app

```r
# From R or RStudio
shiny::runApp()

# Or from the shell
Rscript -e "shiny::runApp()"
```

Required packages: `shiny`, `ggplot2`, `dplyr`, `shinyjs`.

There are no tests or a lint step defined for this project.

## Architecture

This is an educational R Shiny app for a USF 2026 population genetics course. It uses the standard Shiny module pattern: `app.R` is the entry point that wires together two self-contained tabs.

**File layout:**

- `app.R`: Defines the `navbarPage` UI shell and delegates server logic to the two modules. Sources all `R/` files before use.
- `R/wf_simulation.R`: Pure R functions (`init_generation`, `sample_one_individual`, `sample_generation`) for the Wright-Fisher model. No Shiny dependencies.
- `R/tab1_module.R`: Shiny module (`tab1UI` / `tab1Server`) for the interactive Wright-Fisher simulator. Depends on `wf_simulation.R` and `shinyjs`.
- `R/ooa_simulation.R`: Pure R functions (`compute_ne_history`, `ne_to_dots`, `simulate_afs`, `afs_fit_score`) for the Out of Africa demographic model and Allele Frequency Spectrum (AFS). No Shiny dependencies.
- `R/tab2_module.R`: Shiny module (`tab2UI` / `tab2Server`) for the Out of Africa explorer. Depends on `ooa_simulation.R`. The "true" hidden parameters (`TRUE_T_SPLIT=2400`, `TRUE_BOTTLENECK=500`, `TRUE_GROWTH_RATE=2.0`) are constants at the top of this file; the real AFS is pre-computed once at server startup with a fixed seed.

**Key conventions:**

- Simulation logic lives in the pure-function files (`wf_simulation.R`, `ooa_simulation.R`); Shiny reactivity lives only in the module files.
- Brand colors: navbar uses USF purple (`#533860`) and gold (`#FFE42A`). Tab 1 uses COI red (`#C0392B`) for allele 0 (Mario) and COI green (`#27AE60`) for allele 1 (Luigi). Tab 2 uses orange for Africa, blue for non-Africa, purple for ancestral.
- All inputs are namespaced via `ns <- NS(id)` following standard Shiny module conventions.
