library(shiny)
library(tidyverse)
library(plotly)
library(gtools)

# ---------------------------------------------------------------------------
# Load pools once at startup
# ---------------------------------------------------------------------------
data_dir <- tryCatch(
  file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data"),
  error = function(e) file.path(getwd(), "data")
)
mrb_pool <- readRDS(file.path(data_dir, "pool_mrb_summary.rds")) %>%
  filter(polarization <= 4)
srb_pool <- readRDS(file.path(data_dir, "pool_srb_districts.rds"))

MAX_POINTS <- 2500

# Taagepera threshold reference
taag_ref <- data.frame(DM = 1:30) %>%
  mutate(threshold_pct = paste0(round(100 * 0.75 / (DM + 1), 1), "%"),
         y_label = 4.6)

# ---------------------------------------------------------------------------
# Filter definitions (Explore tab)
# ---------------------------------------------------------------------------
pol_choices <- c("All"="all","0–1"="0_1","1–2"="1_2","2–3"="2_3","3–4"="3_4")
pol_limits  <- list(all=c(0,Inf),`0_1`=c(0,1),`1_2`=c(1,2),
                    `2_3`=c(2,3),`3_4`=c(3,4))
all_parties <- sort(unique(c(mrb_pool$n_parties, srb_pool$n_parties)))

# ---------------------------------------------------------------------------
# Seat allocation functions (used in Scenario Designer)
# ---------------------------------------------------------------------------
dhondt <- function(votes, n_seats) {
  if (n_seats == 0L) return(rep(0L, length(votes)))
  votes[is.na(votes) | votes < 0] <- 0
  if (sum(votes) == 0) return(rep(0L, length(votes)))
  seats <- rep(0L, length(votes))
  for (k in seq_len(n_seats)) {
    q <- votes / (seats + 1)
    w <- which.max(q)
    seats[w] <- seats[w] + 1L
  }
  seats
}

sainte_lague <- function(votes, n_seats) {
  if (n_seats == 0L) return(rep(0L, length(votes)))
  votes[is.na(votes) | votes < 0] <- 0
  if (sum(votes) == 0) return(rep(0L, length(votes)))
  seats <- rep(0L, length(votes))
  for (k in seq_len(n_seats)) {
    q <- votes / (2L * seats + 1L)
    w <- which.max(q)
    seats[w] <- seats[w] + 1L
  }
  seats
}

# Modified Sainte-Laguë: first divisor = 1.4, then 3, 5, 7, ...
mod_sainte_lague <- function(votes, n_seats) {
  if (n_seats == 0L) return(rep(0L, length(votes)))
  votes[is.na(votes) | votes < 0] <- 0
  if (sum(votes) == 0) return(rep(0L, length(votes)))
  seats <- rep(0L, length(votes))
  for (k in seq_len(n_seats)) {
    div <- ifelse(seats == 0L, 1.4, 2.0 * seats + 1.0)
    q   <- votes / div
    w   <- which.max(q)
    seats[w] <- seats[w] + 1L
  }
  seats
}

hare_lr <- function(votes, n_seats) {
  n_p <- length(votes)
  if (n_seats == 0L) return(rep(0L, n_p))
  votes[is.na(votes) | votes < 0] <- 0
  tot <- sum(votes)
  if (tot == 0) return(rep(0L, n_p))
  quota <- tot / n_seats
  seats <- floor(votes / quota)
  extra <- as.integer(n_seats - sum(seats))
  if (extra > 0L) {
    rem <- votes - seats * quota
    ord <- order(rem, decreasing = TRUE)
    seats[ord[seq_len(extra)]] <- seats[ord[seq_len(extra)]] + 1L
  }
  as.integer(seats)
}

# Droop quota: tot / (n_seats + 1) — works with vote proportions
droop_lr <- function(votes, n_seats) {
  n_p <- length(votes)
  if (n_seats == 0L) return(rep(0L, n_p))
  votes[is.na(votes) | votes < 0] <- 0
  tot <- sum(votes)
  if (tot == 0) return(rep(0L, n_p))
  quota <- tot / (n_seats + 1L)
  seats <- floor(votes / quota)
  extra <- as.integer(n_seats - sum(seats))
  if (extra > 0L) {
    rem <- votes - seats * quota
    ord <- order(rem, decreasing = TRUE)
    seats[ord[seq_len(extra)]] <- seats[ord[seq_len(extra)]] + 1L
  }
  as.integer(seats)
}

get_alloc_fn <- function(formula) {
  switch(formula,
    dhondt           = dhondt,
    sainte_lague     = sainte_lague,
    mod_sainte_lague = mod_sainte_lague,
    hare             = hare_lr,
    droop            = droop_lr
  )
}

formula_label <- function(formula) {
  switch(formula,
    dhondt           = "D'Hondt",
    sainte_lague     = "Sainte-Laguë",
    mod_sainte_lague = "SL-mod (1.4)",
    hare             = "Hare (LR)",
    droop            = "Droop (LR)"
  )
}

# ---------------------------------------------------------------------------
# Phase 1: Generate base simulations (vote shares + party ideologies only)
# These are the "raw" simulations that Electoral Rules are applied to.
# Only needs to re-run when party system parameters change.
# ---------------------------------------------------------------------------
gen_base_sims <- function(n_sims, np_mode, n_parties_fixed, pol_mode) {
  SC      <- 10.0
  AL      <- 1.5
  N       <- 100L
  beta_sh <- switch(pol_mode, random=2.5, low=6.0, medium=2.5, high=0.8)

  lapply(seq_len(n_sims), function(i) {
    set.seed(i * 1003L)
    n_p   <- if (np_mode == "random")
               sample(2:6, 1, prob = c(0.10, 0.20, 0.30, 0.25, 0.15))
             else n_parties_fixed
    ideol <- rbeta(n_p, beta_sh, beta_sh) * SC
    vs    <- rdirichlet(N, rep(AL, n_p))
    list(n_p=n_p, ideol=ideol, vs=vs)
  })
}

# ---------------------------------------------------------------------------
# Phase 2: Apply electoral rules to base simulations
# Parameters: DM, MAL, malapp_type, threshold, formula
# The malapp assignment uses a separate fixed seed per simulation so it is
# consistent across Refreshes (same districts remain over/under-represented).
# ---------------------------------------------------------------------------
apply_electoral_rules <- function(base_sims, dm_fixed, mal_target, malapp_type,
                                   th_frac, formula) {
  N        <- 100L
  base_pop <- 1000.0
  fn       <- get_alloc_fn(formula)

  res <- vector("list", length(base_sims))

  for (i in seq_along(base_sims)) {
    sim   <- base_sims[[i]]
    n_p   <- sim$n_p
    ideol <- sim$ideol
    vs    <- sim$vs

    IS_d <- as.numeric(vs %*% ideol)

    # Population assignment — seed fixed per sim for consistent Refresh
    if (mal_target < 0.001) {
      pop <- rep(base_pop, N)
    } else {
      set.seed(i * 997L)
      if (malapp_type == "random") {
        idx_over <- sample(N, N %/% 2L)
      } else {
        favor_right <- runif(1) > 0.5
        sorted_idx  <- order(IS_d, decreasing = favor_right)
        idx_over    <- sorted_idx[seq_len(N %/% 2L)]
      }
      pop <- rep((1 + 2 * mal_target) * base_pop, N)
      pop[idx_over] <- (1 - 2 * mal_target) * base_pop
    }

    MAL_act <- 0.5 * sum(abs(rep(1.0/N, N) - pop/sum(pop)))
    IP      <- weighted.mean(IS_d, w = pop)
    IS      <- mean(IS_d)
    MRB     <- IS - IP

    # SRB: seat allocation without threshold
    seats_b <- do.call(rbind, lapply(seq_len(N), function(d) fn(vs[d,], dm_fixed)))
    nat_b   <- colSums(seats_b)
    IR_base <- if (sum(nat_b) > 0) weighted.mean(ideol, w=nat_b) else NA_real_
    SRB     <- if (!is.na(IR_base)) IR_base - IS else NA_real_

    # ThRB: seat allocation with threshold applied per district
    seats_t <- do.call(rbind, lapply(seq_len(N), function(d) {
      vs_f <- ifelse(vs[d,] < th_frac, 0, vs[d,])
      fn(vs_f, dm_fixed)
    }))
    nat_t  <- colSums(seats_t)
    IR_th  <- if (sum(nat_t) > 0) weighted.mean(ideol, w=nat_t) else NA_real_
    ThRB   <- if (!is.na(IR_th) && !is.na(IR_base)) IR_th - IR_base else NA_real_
    TRB    <- if (!is.na(IR_th)) IR_th - IP else NA_real_

    a_MRB  <- abs(MRB)
    a_SRB  <- if (!is.na(SRB))  abs(SRB)  else 0
    a_ThRB <- if (!is.na(ThRB)) abs(ThRB) else 0
    a_tot  <- a_MRB + a_SRB + a_ThRB

    nat_vs  <- colMeans(vs); nat_vs <- nat_vs / sum(nat_vs)
    mean_i  <- sum(nat_vs * ideol)

    res[[i]] <- data.frame(
      n_parties    = n_p,
      MAL          = MAL_act,
      polarization = sum(nat_vs * abs(ideol - mean_i)),
      enep         = 1 / sum(nat_vs^2),
      MRB          = MRB,
      SRB          = SRB,
      ThRB         = ThRB,
      TRB          = TRB,
      BS_MRB       = if (a_tot > 1e-9) a_MRB  / a_tot else NA_real_,
      BS_SRB       = if (a_tot > 1e-9) a_SRB  / a_tot else NA_real_,
      BS_ThRB      = if (a_tot > 1e-9) a_ThRB / a_tot else NA_real_
    )
  }
  do.call(rbind, res)
}

# ---------------------------------------------------------------------------
# ggplot themes
# ---------------------------------------------------------------------------
theme_rep <- theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face="bold", hjust=0.5, size=13),
        plot.subtitle = element_text(hjust=0.5, color="grey45", size=10.5),
        axis.title    = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

theme_rep_dark <- theme_rep +
  theme(panel.background = element_rect(fill="#1e1e2e"),
        plot.background  = element_rect(fill="#1e1e2e"),
        text             = element_text(color="#cdd6f4"),
        axis.text        = element_text(color="#cdd6f4"),
        panel.grid.major = element_line(color="#313244"),
        panel.border     = element_rect(color="#45475a"))

# ---------------------------------------------------------------------------
# CSS (global)
# ---------------------------------------------------------------------------
css_light <- "
  body,.well{background:#fff;color:#212529;}
  .navbar-default{background:#2c3e50!important;border:none!important;}
  .navbar-default .navbar-brand,.navbar-default .navbar-nav>li>a{
    color:#ecf0f1!important;}
  .navbar-default .navbar-nav>.active>a{
    background:#1a252f!important;color:#fff!important;}
  .sidebar-title{font-weight:700;font-size:11px;color:#6c757d;
    text-transform:uppercase;letter-spacing:.05em;margin-bottom:2px;}
  .sc-label{font-weight:600;font-size:12px;color:#495057;margin-bottom:1px;}
  .sc-hint{font-size:11px;color:#868e96;font-style:italic;
    margin-top:1px;margin-bottom:8px;line-height:1.35;}
  .stats-wrap{background:#f8f9fa;border:1px solid #dee2e6;
    border-radius:6px;padding:12px 16px;margin-top:12px;}
  .stats-wrap table{width:100%;font-size:12.5px;border-collapse:collapse;}
  .stats-wrap th{color:#495057;font-weight:600;padding:4px 8px;
    border-bottom:2px solid #dee2e6;text-align:left;}
  .stats-wrap td{padding:3px 8px;}
  .stats-wrap tr:nth-child(even){background:#ffffff;}
  .reset-btn{color:#495057!important;background:#e9ecef!important;
    border:1px solid #ced4da!important;font-size:12px!important;}
  .dl-btn{color:#fff!important;background:#0d6efd!important;
    border:none!important;font-size:12px!important;width:100%;}
  .run-btn{color:#fff!important;background:#198754!important;
    border:none!important;font-size:13px!important;font-weight:600!important;
    width:100%!important;margin-bottom:4px!important;}
  .refresh-btn{color:#fff!important;background:#0d6efd!important;
    border:none!important;font-size:13px!important;font-weight:600!important;
    width:100%!important;margin-bottom:4px!important;}
  .rnd-btn{color:#495057!important;background:#fff3cd!important;
    border:1px solid #ffc107!important;font-size:12px!important;width:100%!important;}
  .dyn-title{font-size:12px;color:#6c757d;font-style:italic;
    margin-bottom:4px;padding-left:2px;}
  .sc-params-box{background:#f8f9fa;border:1px solid #dee2e6;border-radius:8px;
    padding:16px;margin-bottom:12px;}
  .grp-sim{background:#edf7ed;border:1px solid #c3dfc3;border-radius:6px;
    padding:10px 12px;margin-bottom:10px;}
  .grp-rules{background:#fdfbe8;border:1px solid #e8d96b;border-radius:6px;
    padding:10px 12px;margin-bottom:10px;}
  .grp-title{font-size:10px;font-weight:700;text-transform:uppercase;
    letter-spacing:0.06em;margin-bottom:6px;}
  .grp-title-sim{color:#2d6a2d;}
  .grp-title-rules{color:#856404;}
  .history-wrap{background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;
    padding:10px 14px;margin-top:12px;}
  .history-wrap summary{cursor:pointer;font-weight:600;font-size:12px;
    color:#495057;user-select:none;outline:none;}
  .history-wrap summary:hover{color:#212529;}
  .history-tbl{width:100%;border-collapse:collapse;font-size:10px;
    font-family:monospace;margin-top:6px;}
  .history-tbl th{padding:2px 6px;border-bottom:2px solid #dee2e6;
    white-space:nowrap;color:#495057;font-weight:700;text-align:left;}
  .history-tbl td{padding:1px 6px;white-space:nowrap;}
  .history-tbl tr:nth-child(even){background:#ffffff;}
  .history-tbl tr:hover{background:#e8f4f8;}
  .nsims-row{display:flex;align-items:center;gap:8px;margin-top:10px;}
  .nsims-row .form-group{margin-bottom:0;}
  .nsims-row input[type=number]{width:68px!important;height:26px!important;
    font-size:11px!important;padding:2px 4px!important;}
  .nsims-label{font-size:11px;color:#868e96;white-space:nowrap;}
"

css_dark <- "
  body{background:#1e1e2e!important;color:#cdd6f4!important;}
  .well,.sc-params-box{background:#181825!important;border-color:#313244!important;
    color:#cdd6f4!important;}
  .selectize-input,.selectize-dropdown{background:#181825!important;
    color:#cdd6f4!important;border-color:#45475a!important;}
  label,.sidebar-title,.sc-label,.dyn-title{color:#a6adc8!important;}
  .sc-hint{color:#6c7086!important;}
  .stats-wrap{background:#181825!important;border-color:#45475a!important;}
  .stats-wrap th,.stats-wrap td{color:#cdd6f4!important;}
  .stats-wrap tr:nth-child(even){background:#1e1e2e!important;}
  .grp-sim{background:#1a2e1a!important;border-color:#2d4a2d!important;}
  .grp-rules{background:#252210!important;border-color:#4a4020!important;}
  .grp-title-sim{color:#a8d5a8!important;}
  .grp-title-rules{color:#d4c060!important;}
  .history-wrap{background:#181825!important;border-color:#45475a!important;}
  .history-wrap summary{color:#a6adc8!important;}
  .history-tbl th{color:#a6adc8!important;border-color:#45475a!important;}
  .history-tbl td{color:#cdd6f4!important;}
  .history-tbl tr:nth-child(even){background:#1e1e2e!important;}
  .history-tbl tr:hover{background:#24273a!important;}
  .nsims-label{color:#6c7086!important;}
  h2,h3{color:#cba6f7!important;}
  .subtitle{color:#a6adc8!important;}
"

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- navbarPage(
  title = "Unpacking Representation Bias",
  id    = "main_tabs",
  header = tagList(
    tags$head(tags$style(HTML(css_light))),
    uiOutput("dark_css")
  ),

  # ==========================================================================
  # TAB 1: Explore Distributions
  # ==========================================================================
  tabPanel("Explore Distributions",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        checkboxInput("dark_mode", "🌙  Dark mode", value = FALSE),
        hr(style="margin:6px 0;"),
        div(class="sidebar-title", "Bias measure"),
        selectInput("metric", label=NULL,
          choices = c(
            "MRB — Random malapportionment"        = "mrb_random",
            "MRB — Concentrated malapportionment"  = "mrb_concentrated",
            "SRB — Seat Representation Bias"       = "srb",
            "ThRB — Threshold Representation Bias" = "thrb"
          )
        ),
        conditionalPanel(
          condition = "input.metric == 'mrb_concentrated'",
          div(class="sidebar-title", "Direction of malapportionment"),
          selectInput("direction", label=NULL,
            choices = c("All"="all","Leftward (red)"="left","Rightward (blue)"="right"),
            selected = "all"
          )
        ),
        hr(style="margin:8px 0;"),
        div(class="sidebar-title", "Polarization (Dalton index)"),
        selectInput("pol_range", label=NULL, choices=pol_choices, selected="all"),
        div(class="sidebar-title", "Number of parties in the system"),
        checkboxGroupInput("n_parties_sel", label=NULL,
          choices  = setNames(all_parties, paste(all_parties, "parties")),
          selected = all_parties, inline = TRUE
        ),
        hr(style="margin:8px 0;"),
        actionButton("reset_all", "↺  Reset all filters",
                     class="reset-btn", width="100%",
                     style="margin-bottom:6px;"),
        downloadButton("dl_data", "⬇  Download filtered data (.csv)",
                       class="dl-btn"),
        hr(style="margin:8px 0;"),
        uiOutput("n_obs_text")
      ),
      mainPanel(
        width = 9,
        uiOutput("dyn_title"),
        plotlyOutput("main_plot", height="460px"),
        plotOutput("density_strip", height="70px"),
        uiOutput("stats_box")
      )
    )
  ),

  # ==========================================================================
  # TAB 2: Scenario Designer
  # ==========================================================================
  tabPanel("Scenario Designer",
    fluidRow(
      # ---- Left: parameter panel -------------------------------------------
      column(4,
        div(class="sc-params-box",
          h4(style="margin-top:0;font-weight:700;", "Design your electoral system"),
          p(class="sc-hint",
            "Set parameters below and click Run or Refresh. Run generates new
             vote-share simulations; Refresh re-applies electoral rules to the
             same simulations (faster — useful when only changing DM, MAL,
             formula, or threshold)."),

          # ── Simulation parameters (green box) ─────────────────────────────
          div(class="grp-sim",
            div(class="grp-title grp-title-sim",
                "▶  Simulation parameters — Run required to change"),

            div(class="sc-label", "Number of parties"),
            fluidRow(
              column(6,
                radioButtons("np_mode_sc", label=NULL,
                  choices = c("Random"="random","Fixed"="fixed"),
                  selected = "random", inline = TRUE)
              ),
              column(6,
                conditionalPanel(
                  condition = "input.np_mode_sc == 'fixed'",
                  selectInput("np_fixed_sc", label=NULL,
                    choices = setNames(2:6, paste(2:6,"parties")), selected=4)
                )
              )
            ),
            div(class="sc-hint", "Random: drawn from multiparty distribution (2–6)."),

            div(class="sc-label", "Party system polarization"),
            radioButtons("pol_mode_sc", label=NULL,
              choices  = c("Random"="random","Low"="low",
                           "Medium"="medium","High"="high"),
              selected = "random", inline = TRUE),
            div(class="sc-hint",
                "Controls how spread out parties are on the ideological scale (0–10).
                 Low: near centre. High: at the extremes."),

          ),

          # ── Electoral rules (yellow box) ───────────────────────────────────
          div(class="grp-rules",
            div(class="grp-title grp-title-rules",
                "↺  Electoral rules — Refresh OK (same simulations)"),

            div(class="sc-label", "District magnitude (DM, seats per district)"),
            sliderInput("dm_sc", label=NULL, min=1, max=30, value=5, step=1),
            div(class="sc-hint",
                "All 100 districts have the same DM. This isolates the effect of
                 district magnitude from other sources of variation."),

            div(class="sc-label", "Malapportionment (Samuels & Snyder index, MAL)"),
            sliderInput("mal_sc", label=NULL, min=0, max=0.45, value=0.10, step=0.01),
            div(class="sc-hint",
                "MAL = 0: perfect apportionment. MAL = 0.45: extreme disproportion."),
            radioButtons("malapp_type_sc", label=NULL,
              choices  = c("Random"="random","Concentrated"="concentrated"),
              selected = "random"),
            div(class="sc-hint",
                "Random: over-representation without ideological pattern.
                 Concentrated: one ideological bloc systematically over-represented.
                 Same districts stay over-represented when refreshing."),

            div(class="sc-label", "Electoral formula"),
            selectInput("formula_sc", label=NULL,
              choices = c(
                "D'Hondt"                      = "dhondt",
                "Sainte-Laguë"                 = "sainte_lague",
                "Sainte-Laguë modified (1.4)"  = "mod_sainte_lague",
                "Hare (largest remainder)"     = "hare",
                "Droop (largest remainder)"    = "droop"
              ),
              selected = "dhondt"
            ),
            div(class="sc-hint",
                "D'Hondt favours large parties. Sainte-Laguë is more proportional.
                 Modified SL (first divisor 1.4) is used in Nordic countries.
                 Hare and Droop are quota-based largest-remainder methods."),

            div(class="sc-label", "Electoral threshold (%)"),
            sliderInput("th_sc", label=NULL, min=0, max=20, value=3, step=0.5),
            div(class="sc-hint",
                "Applied uniformly to all districts. Parties below threshold
                 are excluded from seat allocation.")
          ),

          hr(style="margin:8px 0;"),

          actionButton("run_sc",      "▶  Run simulations",
                       class="run-btn",     style="margin-bottom:4px;"),
          actionButton("refresh_sc",  "↺  Refresh (same simulations)",
                       class="refresh-btn", style="margin-bottom:4px;"),
          actionButton("randomize_sc","🎲  Randomize all parameters",
                       class="rnd-btn"),
          div(class="nsims-row",
            span(class="nsims-label", "Modify the number of simulations:"),
            numericInput("n_sims_sc", label=NULL, value=300,
                         min=1, max=2000, step=1, width="75px")
          ),

          hr(style="margin:8px 0;"),
          downloadButton("dl_sc", "⬇  Download results (.csv)", class="dl-btn")
        )
      ),

      # ---- Right: results panel --------------------------------------------
      column(8,
        uiOutput("sc_status"),
        plotlyOutput("sc_plot",  height="400px"),
        plotOutput("sc_strip",   height="80px"),
        uiOutput("sc_stats"),
        uiOutput("sc_history")
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  output$dark_css <- renderUI({
    if (isTRUE(input$dark_mode))
      tags$head(tags$style(HTML(css_dark)))
  })

  # ==========================================================================
  # EXPLORE TAB
  # ==========================================================================

  filtered <- reactive({
    pol_lim <- pol_limits[[input$pol_range]]
    sel_p   <- as.integer(input$n_parties_sel)
    if (length(sel_p) == 0) sel_p <- all_parties

    if (input$metric %in% c("mrb_random","mrb_concentrated")) {
      scen <- switch(input$metric, mrb_random="random", mrb_concentrated="concentrated")
      d <- filter(mrb_pool, scenario == scen,
                  n_parties %in% sel_p,
                  polarization >= pol_lim[1], polarization <= pol_lim[2])
      if (input$metric == "mrb_concentrated") {
        dir <- if (!is.null(input$direction)) input$direction else "all"
        if (dir == "left")  d <- filter(d, MRB < 0)
        if (dir == "right") d <- filter(d, MRB > 0)
      }
      d
    } else {
      srb_pool %>% filter(n_parties %in% sel_p,
                          polarization >= pol_lim[1], polarization <= pol_lim[2])
    }
  })

  plot_data <- reactive({
    d <- filtered()
    if (nrow(d) > MAX_POINTS) slice_sample(d, n = MAX_POINTS) else d
  })

  observeEvent(input$reset_all, {
    updateSelectInput(session,        "metric",       selected = "mrb_random")
    updateSelectInput(session,        "direction",    selected = "all")
    updateSelectInput(session,        "pol_range",    selected = "all")
    updateCheckboxGroupInput(session, "n_parties_sel",selected = all_parties)
    updateCheckboxInput(session,      "dark_mode",    value    = FALSE)
  })

  output$dl_data <- downloadHandler(
    filename = function() paste0("repbias_", input$metric, "_",
                                 format(Sys.time(),"%Y%m%d_%H%M"), ".csv"),
    content  = function(f) write.csv(filtered(), f, row.names = FALSE)
  )

  output$n_obs_text <- renderUI({
    n_full <- nrow(filtered()); n_plot <- nrow(plot_data())
    col <- if (n_full == 0) "red" else "#6c757d"
    msg <- if (n_full == 0) "No observations match the selected filters."
           else if (n_full > MAX_POINTS)
             sprintf("Showing %s of %s obs. (random sample).",
                     format(n_plot,big.mark=","), format(n_full,big.mark=","))
           else sprintf("Showing all %s observations.", format(n_full,big.mark=","))
    div(tags$span(style=paste0("color:",col,";font-size:12px;"), msg))
  })

  output$dyn_title <- renderUI({
    pol_lbl <- names(pol_choices)[pol_choices == input$pol_range]
    p_lbl   <- if (length(input$n_parties_sel) == length(all_parties)) "All parties"
               else paste(sort(input$n_parties_sel), collapse=", ")
    lbl <- switch(input$metric,
      mrb_random="MRB (Random)", mrb_concentrated="MRB (Concentrated)",
      srb="SRB", thrb="ThRB")
    div(class="dyn-title",
        paste0(lbl,"  ·  Polarization: ",pol_lbl,"  ·  Parties: ",p_lbl))
  })

  get_theme <- reactive({
    if (isTRUE(input$dark_mode)) theme_rep_dark else theme_rep
  })

  tooltip_mrb <- function(d) paste0(
    "<b>Simulation #",d$Simulacion,"</b><br>",
    "Parties: ",d$n_parties,"  |  Polarization: ",round(d$polarization,2),"<br>",
    "ENEP: ",round(d$enep,2),"<br>",
    "MAL: ",round(d$MAL,3),"<br>",
    "MRB: ",round(d$MRB,3),"  |  SRB: ",round(d$SRB,3),
    "  |  TRB: ",round(d$TRB,3)
  )

  output$main_plot <- renderPlotly({
    d <- plot_data(); req(nrow(d) > 0)
    m <- input$metric; th <- get_theme()

    p <- if (m == "mrb_random") {
      ggplot(d, aes(x=MAL, y=MRB, text=tooltip_mrb(d))) +
        geom_point(alpha=0.35, size=1, color="deepskyblue3") +
        geom_hline(yintercept=0, linetype="dashed", color="grey50") +
        geom_smooth(method="loess", se=FALSE, color="grey30",
                    linewidth=0.9, span=0.6) +
        scale_x_continuous(limits=c(0,0.5)) +
        labs(title="Random Malapportionment and MRB",
             x="Malapportionment index (MAL)",
             y="Malapportionment Representation Bias (MRB)")

    } else if (m == "mrb_concentrated") {
      d <- d %>% mutate(dir_col = ifelse(MRB < 0, "Leftward", "Rightward"))
      ggplot(d, aes(x=MAL, y=MRB, color=dir_col, text=tooltip_mrb(d))) +
        geom_point(alpha=0.35, size=1) +
        geom_hline(yintercept=0, linetype="dashed", color="grey50") +
        geom_smooth(aes(group=1), method="loess", se=FALSE,
                    color="grey30", linewidth=0.9, span=0.6) +
        scale_color_manual(
          values=c("Leftward"="#c0392b","Rightward"="#2980b9"),
          name="Bias direction") +
        scale_x_continuous(limits=c(0,0.5)) +
        labs(title="Concentrated Malapportionment and MRB",
             subtitle="Red = leftward bias  ·  Blue = rightward bias",
             x="Malapportionment index (MAL)",
             y="Malapportionment Representation Bias (MRB)")

    } else if (m == "srb") {
      taag_ann <- taag_ref %>% filter(DM %in% c(1,2,3,5,7,10,15,20,30))
      ggplot(d, aes(x=DM_fixed, y=SRB_distrito_base,
                    text=paste0("<b>District (DM=",DM_fixed,")</b><br>",
                                "SRB: ",round(SRB_distrito_base,3),"<br>",
                                "Parties: ",n_parties,"<br>",
                                "Polarization: ",round(polarization,2)))) +
        geom_jitter(alpha=0.2, size=0.8, color="deepskyblue3", width=0.3) +
        geom_hline(yintercept=0, linetype="dashed", color="grey50") +
        geom_smooth(method="loess", se=FALSE, color="grey30",
                    linewidth=0.9, span=0.5) +
        geom_text(data=taag_ann, aes(x=DM, y=y_label, label=threshold_pct),
                  inherit.aes=FALSE, size=2.8, color="grey55", hjust=0.5) +
        annotate("text", x=1, y=5.05, label="Taagepera threshold:",
                 size=2.8, color="grey55", hjust=0) +
        scale_x_continuous(breaks=c(1,5,10,15,20,25,30)) +
        scale_y_continuous(limits=c(-5,5.2)) +
        labs(title="District Magnitude and Seat Representation Bias (SRB)",
             x="District Magnitude (seats)",
             y="Seat Representation Bias (SRB)")

    } else {
      ggplot(d, aes(x=TH_USED, y=Threshold_Bias,
                    text=paste0("<b>District</b><br>",
                                "Threshold: ",scales::percent(TH_USED,accuracy=0.1),"<br>",
                                "ThRB: ",round(Threshold_Bias,3),"<br>",
                                "Parties: ",n_parties,"<br>",
                                "Polarization: ",round(polarization,2)))) +
        geom_point(alpha=0.25, size=0.8, color="deepskyblue3") +
        geom_hline(yintercept=0, linetype="dashed", color="grey50") +
        geom_smooth(method="loess", se=FALSE, color="grey30",
                    linewidth=0.9, span=0.5) +
        scale_x_continuous(labels=scales::percent, name="Legal Threshold Applied") +
        labs(title="Electoral Threshold and Threshold Representation Bias (ThRB)",
             y="Threshold Representation Bias (ThRB)")
    }

    p <- p + th
    ggplotly(p, tooltip="text") %>%
      layout(hoverlabel=list(bgcolor="white", bordercolor="grey80",
                             font=list(size=12,family="Helvetica")),
             margin=list(t=50)) %>%
      config(displayModeBar=TRUE, displaylogo=FALSE,
             modeBarButtonsToRemove=c("lasso2d","select2d"))
  })

  output$density_strip <- renderPlot(bg="transparent", {
    d <- plot_data(); req(nrow(d) > 0)
    y_var <- switch(input$metric,
      mrb_random="MRB", mrb_concentrated="MRB",
      srb="SRB_distrito_base", thrb="Threshold_Bias")
    y_lbl <- switch(input$metric,
      mrb_random="MRB", mrb_concentrated="MRB", srb="SRB", thrb="ThRB")
    bg_col   <- if (isTRUE(input$dark_mode)) "#1e1e2e" else "white"
    fill_col <- if (isTRUE(input$dark_mode)) "#89b4fa" else "deepskyblue3"
    txt_col  <- if (isTRUE(input$dark_mode)) "#cdd6f4" else "grey30"
    ggplot(d, aes(x=.data[[y_var]])) +
      geom_density(fill=fill_col, color=NA, alpha=0.5) +
      geom_vline(xintercept=mean(d[[y_var]],na.rm=TRUE),
                 linetype="dashed", color=txt_col, linewidth=0.6) +
      scale_x_continuous(name=paste("Distribution of",y_lbl)) +
      theme_void() +
      theme(axis.title.x = element_text(size=10, color=txt_col, hjust=0.5),
            plot.background = element_rect(fill=bg_col, color=NA))
  })

  output$stats_box <- renderUI({
    d <- filtered(); req(nrow(d) > 0)
    y_var <- switch(input$metric,
      mrb_random="MRB", mrb_concentrated="MRB",
      srb="SRB_distrito_base", thrb="Threshold_Bias")
    lbl <- switch(input$metric,
      mrb_random="MRB", mrb_concentrated="MRB", srb="SRB", thrb="ThRB")
    v  <- d[[y_var]]
    mu <- mean(v,na.rm=TRUE); s <- sd(v,na.rm=TRUE)
    fmt <- function(x) formatC(round(x,4), format="f", digits=4)
    rows <- list(
      c("Mean", fmt(mu)),
      c("Mean + 1 SD", fmt(mu+s)),   c("Mean − 1 SD", fmt(mu-s)),
      c("Mean + 2 SD", fmt(mu+2*s)), c("Mean − 2 SD", fmt(mu-2*s)),
      c("Maximum", fmt(max(v,na.rm=TRUE))),
      c("Minimum", fmt(min(v,na.rm=TRUE)))
    )
    tbl_rows <- lapply(rows, function(r)
      tags$tr(tags$td(r[1]),
              tags$td(style="text-align:right;font-family:monospace;", r[2])))
    div(class="stats-wrap",
      div(style="font-weight:700;font-size:12px;color:#6c757d;
                 text-transform:uppercase;letter-spacing:.05em;margin-bottom:6px;",
          paste("Summary statistics —", lbl, "(0–10 scale)")),
      tags$table(
        tags$thead(tags$tr(tags$th("Statistic"),
                           tags$th(style="text-align:right;","Value"))),
        tags$tbody(tbl_rows)
      )
    )
  })

  # ==========================================================================
  # SCENARIO DESIGNER TAB
  # ==========================================================================

  # Stored state: last-used parameters (snapshotted on Run/Refresh)
  sc_params <- reactiveValues(
    dm       = 5L,
    mal      = 0.10,
    malapp   = "random",
    th       = 3.0,
    formula  = "dhondt",
    np_mode  = "random",
    np_fixed = 4L,
    pol      = "random",
    n_sims   = 300L,
    action   = "none",   # "run" or "refresh"
    trigger  = 0L        # increments on every Run/Refresh (for history)
  )

  # Base simulations (vote shares + party ideologies, independent of rules)
  base_sims <- reactiveVal(NULL)

  # Current scenario results
  sc_results_v <- reactiveVal(NULL)

  # History of all runs/refreshes
  scenario_history <- reactiveVal(list())

  # ---- Run button -----------------------------------------------------------
  observeEvent(input$run_sc, {
    sc_params$dm       <- as.integer(input$dm_sc)
    sc_params$mal      <- input$mal_sc
    sc_params$malapp   <- input$malapp_type_sc
    sc_params$th       <- input$th_sc
    sc_params$formula  <- input$formula_sc
    sc_params$np_mode  <- input$np_mode_sc
    sc_params$np_fixed <- as.integer(input$np_fixed_sc)
    sc_params$pol      <- input$pol_mode_sc
    sc_params$n_sims   <- as.integer(input$n_sims_sc)
    sc_params$action   <- "run"

    withProgress(message = sprintf("Running %d simulations…", sc_params$n_sims), value = 0, {
      setProgress(0.10, detail = "Generating vote shares…")
      sims <- gen_base_sims(sc_params$n_sims, sc_params$np_mode,
                            if (sc_params$np_mode=="fixed") sc_params$np_fixed else NULL,
                            sc_params$pol)
      base_sims(sims)

      setProgress(0.55, detail = "Applying electoral rules…")
      results <- apply_electoral_rules(
        sims, sc_params$dm, sc_params$mal, sc_params$malapp,
        sc_params$th / 100, sc_params$formula
      )
      setProgress(1)
      sc_results_v(results)
    })
    sc_params$trigger <- sc_params$trigger + 1L
  })

  # ---- Refresh button -------------------------------------------------------
  observeEvent(input$refresh_sc, {
    if (is.null(base_sims())) {
      showNotification("Run simulations first before refreshing.",
                       type = "warning", duration = 3)
      return()
    }
    sc_params$dm      <- as.integer(input$dm_sc)
    sc_params$mal     <- input$mal_sc
    sc_params$malapp  <- input$malapp_type_sc
    sc_params$th      <- input$th_sc
    sc_params$formula <- input$formula_sc
    sc_params$action  <- "refresh"

    withProgress(message = "Refreshing results…", value = 0.3, {
      results <- apply_electoral_rules(
        base_sims(), sc_params$dm, sc_params$mal, sc_params$malapp,
        sc_params$th / 100, sc_params$formula
      )
      setProgress(1)
      sc_results_v(results)
    })
    sc_params$trigger <- sc_params$trigger + 1L
  })

  # ---- Randomize button -----------------------------------------------------
  observeEvent(input$randomize_sc, {
    dm_val      <- sample(c(1L,2L,3L,5L,7L,10L,15L,20L,30L), 1L)
    mal_val     <- round(runif(1, 0, 0.40), 2)
    th_val      <- round(runif(1, 0, 15), 1)
    np_mode     <- sample(c("random","fixed"), 1)
    np_fix      <- sample(2:6, 1)
    pol_val     <- sample(c("random","low","medium","high"), 1)
    malapp      <- sample(c("random","concentrated"), 1)
    formula_val <- sample(c("dhondt","sainte_lague","mod_sainte_lague","hare","droop"), 1)

    updateSliderInput (session, "dm_sc",         value    = dm_val)
    updateSliderInput (session, "mal_sc",         value    = mal_val)
    updateSliderInput (session, "th_sc",          value    = th_val)
    updateRadioButtons(session, "np_mode_sc",     selected = np_mode)
    updateSelectInput (session, "np_fixed_sc",    selected = np_fix)
    updateRadioButtons(session, "pol_mode_sc",    selected = pol_val)
    updateRadioButtons(session, "malapp_type_sc", selected = malapp)
    updateSelectInput (session, "formula_sc",     selected = formula_val)

    sc_params$dm       <- dm_val
    sc_params$mal      <- mal_val
    sc_params$malapp   <- malapp
    sc_params$th       <- th_val
    sc_params$formula  <- formula_val
    sc_params$np_mode  <- np_mode
    sc_params$np_fixed <- as.integer(np_fix)
    sc_params$pol      <- pol_val
    sc_params$n_sims   <- as.integer(input$n_sims_sc)
    sc_params$action   <- "run"

    withProgress(message = sprintf("Running %d simulations…", sc_params$n_sims), value = 0, {
      setProgress(0.10, detail = "Generating vote shares…")
      sims <- gen_base_sims(sc_params$n_sims, sc_params$np_mode,
                            if (sc_params$np_mode=="fixed") sc_params$np_fixed else NULL,
                            sc_params$pol)
      base_sims(sims)
      setProgress(0.55, detail = "Applying electoral rules…")
      results <- apply_electoral_rules(
        sims, sc_params$dm, sc_params$mal, sc_params$malapp,
        sc_params$th / 100, sc_params$formula
      )
      setProgress(1)
      sc_results_v(results)
    })
    sc_params$trigger <- sc_params$trigger + 1L
  })

  # ---- Update history on every trigger ------------------------------------
  observeEvent(sc_params$trigger, {
    d <- sc_results_v()
    req(!is.null(d) && nrow(d) > 0)

    m_mrb  <- mean(abs(d$MRB),  na.rm = TRUE)
    m_srb  <- mean(abs(d$SRB),  na.rm = TRUE)
    m_thrb <- mean(abs(d$ThRB), na.rm = TRUE)
    m_tot  <- m_mrb + m_srb + m_thrb

    np_lbl   <- if (sc_params$np_mode == "random") "rnd"
                else as.character(sc_params$np_fixed)
    type_lbl <- if (sc_params$malapp == "random") "rnd" else "conc"
    act_lbl  <- if (sc_params$action == "run") "▶" else "↺"
    pct      <- function(x) if (m_tot > 1e-9) sprintf("%.1f%%", 100*x/m_tot) else "—"

    entry <- data.frame(
      `#`      = length(scenario_history()) + 1L,
      Act      = act_lbl,
      Formula  = formula_label(sc_params$formula),
      DM       = sc_params$dm,
      MAL      = sprintf("%.2f", sc_params$mal),
      Type     = type_lbl,
      `TH%`    = sprintf("%.1f", sc_params$th),
      Parties  = np_lbl,
      Pol      = substr(sc_params$pol, 1, 3),
      `|MRB|`  = sprintf("%.3f", m_mrb),
      `|SRB|`  = sprintf("%.3f", m_srb),
      `|ThRB|` = sprintf("%.3f", m_thrb),
      `|TRB|`  = sprintf("%.3f", mean(abs(d$TRB), na.rm=TRUE)),
      `%MRB`   = pct(m_mrb),
      `%SRB`   = pct(m_srb),
      `%ThRB`  = pct(m_thrb),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    hist <- c(list(entry), scenario_history())
    if (length(hist) > 100L) hist <- hist[seq_len(100L)]
    scenario_history(hist)
  }, ignoreInit = TRUE)

  # ---- Status message (shown only before first run) -----------------------
  output$sc_status <- renderUI({
    if (is.null(sc_results_v()))
      div(style="text-align:center;padding:40px;color:#868e96;",
          h4("Set your parameters and click ▶ Run 300 simulations"))
  })

  # ---- Main bar chart -------------------------------------------------------
  output$sc_plot <- renderPlotly({
    req(!is.null(sc_results_v()))
    d <- sc_results_v(); req(nrow(d) > 0)
    th <- get_theme()

    malapp_lbl <- if (sc_params$malapp == "random") "Random" else "Concentrated"
    pol_lbl    <- switch(sc_params$pol,
      random="Random", low="Low", medium="Medium", high="High")
    np_lbl <- if (sc_params$np_mode == "random")
      "Random (2–6)" else paste(sc_params$np_fixed,"parties")
    act_lbl <- if (sc_params$action == "run") "▶ New sims" else "↺ Refreshed"

    m_mrb  <- mean(abs(d$MRB),  na.rm = TRUE)
    m_srb  <- mean(abs(d$SRB),  na.rm = TRUE)
    m_thrb <- mean(abs(d$ThRB), na.rm = TRUE)
    m_tot  <- m_mrb + m_srb + m_thrb

    inst_levels <- c("Malapportionment\n(MRB)",
                     "Seat allocation\n(SRB)",
                     "Threshold\n(ThRB)")

    # Reshape to long format — one row per simulation × institution
    plot_long <- bind_rows(
      data.frame(institution = inst_levels[1], abs_bias = abs(d$MRB),
                 n_parties = d$n_parties, polarization = d$polarization,
                 enep = d$enep, MAL = d$MAL),
      data.frame(institution = inst_levels[2], abs_bias = abs(d$SRB),
                 n_parties = d$n_parties, polarization = d$polarization,
                 enep = d$enep, MAL = d$MAL),
      data.frame(institution = inst_levels[3], abs_bias = abs(d$ThRB),
                 n_parties = d$n_parties, polarization = d$polarization,
                 enep = d$enep, MAL = d$MAL)
    ) %>% mutate(institution = factor(institution, levels = inst_levels))

    means_df <- data.frame(
      institution = factor(inst_levels, levels = inst_levels),
      mean_val    = c(m_mrb, m_srb, m_thrb),
      share_lbl   = if (m_tot > 1e-9)
                      scales::percent(c(m_mrb, m_srb, m_thrb) / m_tot, accuracy = 0.1)
                    else rep("—", 3)
    ) %>%
      left_join(
        plot_long %>%
          group_by(institution) %>%
          summarise(top_val = max(abs_bias, na.rm = TRUE), .groups = "drop"),
        by = "institution"
      )

    fill_vals <- c(
      "Malapportionment\n(MRB)" = "#c0392b",
      "Seat allocation\n(SRB)"  = "#2980b9",
      "Threshold\n(ThRB)"       = "#e67e22"
    )

    p <- ggplot(plot_long,
      aes(x = institution, y = abs_bias, fill = institution,
          text = paste0(
            gsub("\n", " ", as.character(institution)), "<br>",
            "|Bias|: ", round(abs_bias, 3), "<br>",
            "Parties: ", n_parties,
            "  |  Pol: ", round(polarization, 2), "<br>",
            "ENEP: ", round(enep, 2),
            "  |  MAL: ", round(MAL, 3)
          ))) +
      geom_jitter(width = 0.18, alpha = 0.22, size = 0.9, show.legend = FALSE) +
      geom_boxplot(width = 0.35, alpha = 0.55, outlier.shape = NA, color = "grey40",
                   linewidth = 0.6, show.legend = FALSE) +
      geom_point(data = means_df, aes(x = institution, y = mean_val),
                 inherit.aes = FALSE, shape = 18, size = 4.5, color = "grey15") +
      geom_text(data = means_df,
                aes(x = institution, y = top_val, label = share_lbl),
                inherit.aes = FALSE,
                vjust = -0.5, size = 3.0, fontface = "bold", color = "grey15") +
      scale_fill_manual(values = fill_vals) +
      scale_y_continuous(
        name   = "Absolute bias (0–10 scale)",
        expand = expansion(mult = c(0.02, 0.22))) +
      scale_x_discrete(name = NULL) +
      labs(
        title    = "Bias Decomposition — Your Electoral Scenario",
        subtitle = paste0(
          sprintf("%s  ·  %s  ·  DM=%d  ·  MAL=%.2f (%s)  ·  TH=%.1f%%",
            act_lbl, formula_label(sc_params$formula),
            sc_params$dm, sc_params$mal, malapp_lbl, sc_params$th),
          "\n",
          sprintf("Parties: %s  ·  Polarization: %s", np_lbl, pol_lbl))
      ) +
      th +
      theme(
        axis.text     = element_text(size = 9),
        axis.title.y  = element_text(size = 10),
        plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "grey45",
                                     lineheight = 1.3)
      )

    ggplotly(p, tooltip = "text") %>%
      layout(
        showlegend = FALSE,
        hoverlabel = list(bgcolor = "white", bordercolor = "grey80",
                          font = list(size = 12, family = "Helvetica")),
        margin     = list(t = 75, b = 40)
      ) %>%
      config(displayModeBar = FALSE)
  })

  # ---- TRB density strip ---------------------------------------------------
  output$sc_strip <- renderPlot(bg="transparent", {
    req(!is.null(sc_results_v()))
    d <- sc_results_v(); req(nrow(d) > 0)
    bg_col   <- if (isTRUE(input$dark_mode)) "#1e1e2e" else "white"
    fill_col <- if (isTRUE(input$dark_mode)) "#89b4fa" else "#7f8c8d"
    txt_col  <- if (isTRUE(input$dark_mode)) "#cdd6f4" else "grey30"
    ggplot(d, aes(x=abs(TRB))) +
      geom_density(fill=fill_col, color=NA, alpha=0.55) +
      geom_vline(xintercept=mean(abs(d$TRB),na.rm=TRUE),
                 linetype="dashed", color=txt_col, linewidth=0.7) +
      scale_x_continuous(name="Distribution of |TRB| (total bias)") +
      theme_void() +
      theme(axis.title.x = element_text(size=10, color=txt_col, hjust=0.5),
            plot.background = element_rect(fill=bg_col, color=NA))
  })

  # ---- Stats table ---------------------------------------------------------
  output$sc_stats <- renderUI({
    req(!is.null(sc_results_v()))
    d <- sc_results_v(); req(nrow(d) > 0)
    fmt4    <- function(x) formatC(round(x,4), format="f", digits=4)
    fmt_pct <- function(x) scales::percent(x, accuracy=0.1)

    m_mrb  <- mean(abs(d$MRB),  na.rm=TRUE)
    m_srb  <- mean(abs(d$SRB),  na.rm=TRUE)
    m_thrb <- mean(abs(d$ThRB), na.rm=TRUE)
    m_tot  <- m_mrb + m_srb + m_thrb
    safe_pct <- function(x) if (m_tot > 1e-9) fmt_pct(x/m_tot) else "—"

    rows <- list(
      c("── Achieved (simulated) values ──────────", ""),
      c("Mean MAL (Samuels & Snyder)",  fmt4(mean(d$MAL,         na.rm=TRUE))),
      c("Mean polarization (Dalton)",   fmt4(mean(d$polarization, na.rm=TRUE))),
      c("Mean ENEP (Laakso-Taagepera)", fmt4(mean(d$enep,         na.rm=TRUE))),
      c("── Bias components (mean |bias|, 0–10) ──", ""),
      c("Malapportionment |MRB|", fmt4(m_mrb)),
      c("Seat allocation   |SRB|", fmt4(m_srb)),
      c("Threshold        |ThRB|", fmt4(m_thrb)),
      c("Total bias       |TRB|",  fmt4(mean(abs(d$TRB), na.rm=TRUE))),
      c("── Bias shares (% of mean |bias|) ───────", ""),
      c("Share — Malapportionment (MRB)", safe_pct(m_mrb)),
      c("Share — Seat allocation  (SRB)", safe_pct(m_srb)),
      c("Share — Threshold       (ThRB)", safe_pct(m_thrb))
    )
    tbl_rows <- lapply(rows, function(r) {
      if (r[2] == "")
        tags$tr(tags$td(colspan="2",
                        style="font-weight:700;font-size:11px;color:#6c757d;
                               padding-top:8px;", r[1]))
      else
        tags$tr(tags$td(r[1]),
                tags$td(style="text-align:right;font-family:monospace;", r[2]))
    })
    div(class="stats-wrap",
      div(style="font-weight:700;font-size:12px;color:#6c757d;
                 text-transform:uppercase;letter-spacing:.05em;margin-bottom:6px;",
          paste0("Results — ", nrow(d), " simulations  ·  ",
                 formula_label(sc_params$formula))),
      tags$table(
        tags$thead(tags$tr(tags$th("Statistic"),
                           tags$th(style="text-align:right;","Value"))),
        tags$tbody(tbl_rows)
      )
    )
  })

  # ---- Scenario history (collapsible table) --------------------------------
  output$sc_history <- renderUI({
    hist <- scenario_history()
    if (length(hist) == 0L) return(NULL)

    tbl_data  <- do.call(rbind, hist)
    col_names <- names(tbl_data)

    hdr_cells <- lapply(col_names, function(nm) tags$th(nm))

    body_rows <- lapply(seq_len(nrow(tbl_data)), function(r) {
      cells <- lapply(col_names, function(nm)
        tags$td(as.character(tbl_data[r, nm])))
      tags$tr(cells)
    })

    n_runs <- nrow(tbl_data)
    div(class="history-wrap",
      tags$details(
        tags$summary(
          sprintf("Scenario history — %d run%s  (click to expand/collapse)",
                  n_runs, if (n_runs == 1L) "" else "s")
        ),
        div(style="overflow-x:auto; margin-top:4px;",
          tags$table(
            class = "history-tbl",
            tags$thead(tags$tr(hdr_cells)),
            tags$tbody(body_rows)
          )
        )
      )
    )
  })

  # ---- Download scenario results -------------------------------------------
  output$dl_sc <- downloadHandler(
    filename = function() paste0("repbias_scenario_",
                                 format(Sys.time(),"%Y%m%d_%H%M"), ".csv"),
    content  = function(f) {
      d <- sc_results_v()
      req(!is.null(d))
      write.csv(d, f, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
