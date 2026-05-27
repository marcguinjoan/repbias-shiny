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
  /* Tutorial tab */
  .tut-vs-table{border-collapse:collapse;margin-bottom:4px;}
  .tut-vs-table th{font-size:11px;padding:2px 5px;color:#6c757d;font-weight:600;text-align:center;}
  .tut-vs-table td{padding:1px 2px;vertical-align:middle;}
  .tut-table{border-collapse:collapse;font-size:12.5px;width:100%;}
  .tut-table th{background:#f1f3f5;padding:5px 10px;border:1px solid #dee2e6;
    font-weight:600;white-space:nowrap;}
  .tut-table td{padding:4px 10px;border:1px solid #dee2e6;}
  .tut-table tr:nth-child(even){background:#f8f9fa;}
  .metric-card{display:inline-block;min-width:108px;padding:6px 12px;margin:3px 3px 3px 0;
    border-radius:6px;text-align:center;border:1px solid #dee2e6;vertical-align:top;
    background:#f8f9fa;}
  .metric-card .mc-label{font-size:9px;font-weight:700;text-transform:uppercase;
    letter-spacing:.05em;color:#6c757d;line-height:1.4;}
  .metric-card .mc-val{font-size:18px;font-weight:700;font-family:monospace;
    line-height:1.3;color:#495057;}
  .mc-mrb{background:#fde8e8!important;border-color:#f5c6c6!important;}
  .mc-mrb .mc-val{color:#c0392b!important;}
  .mc-srb{background:#e8f0fd!important;border-color:#c6d4f5!important;}
  .mc-srb .mc-val{color:#2980b9!important;}
  .mc-thrb{background:#fdf3e8!important;border-color:#f5d9c6!important;}
  .mc-thrb .mc-val{color:#e67e22!important;}
  .mc-trb{background:#edf7ed!important;border-color:#c3dfc3!important;}
  .mc-trb .mc-val{color:#198754!important;}
  .tut-step{background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;
    padding:12px 16px;margin-bottom:12px;}
  .tut-note{background:#fff9e6;border-left:3px solid #ffc107;padding:8px 12px;
    font-size:12.5px;margin:8px 0;border-radius:0 4px 4px 0;}
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
  .tut-table th{background:#2a2a3d!important;color:#cdd6f4!important;
    border-color:#45475a!important;}
  .tut-table td{border-color:#45475a!important;color:#cdd6f4!important;}
  .tut-table tr:nth-child(even){background:#1e1e2e!important;}
  .tut-step{background:#181825!important;border-color:#45475a!important;}
  .metric-card{background:#181825!important;border-color:#45475a!important;}
  .metric-card .mc-val{color:#cdd6f4!important;}
  .mc-mrb .mc-val{color:#e8a0a0!important;}
  .mc-srb .mc-val{color:#89b4fa!important;}
  .mc-thrb .mc-val{color:#fab387!important;}
  .mc-trb .mc-val{color:#a6e3a1!important;}
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
  ),

  # ==========================================================================
  # TAB 3: Tutorial
  # ==========================================================================
  tabPanel("Tutorial",
    div(style="max-width:980px;margin:0 auto;padding:20px 24px 50px;",

      h2(style="font-weight:700;margin-bottom:4px;",
         "How Representation Bias Measures Work"),
      p(style="color:#6c757d;margin-bottom:22px;",
        "A step-by-step walkthrough of MRB, SRB, ThRB, TRB, and Bias Share using a
         worked example. Adjust the controls on the left to explore how each
         institutional factor shapes the gap between voters and their representatives."),

      fluidRow(
        column(4,
          div(style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:8px;padding:16px;",
            h5(style="font-weight:700;margin-top:0;margin-bottom:10px;",
               "Adjust the Example"),
            p(style="font-size:11.5px;color:#6c757d;margin-bottom:8px;",
              tags$b("Party ideologies (fixed):"),
              " Left = 2 · Centre = 5 · Right = 8  (0–10 scale)"),
            div(class="sc-label", "Vote shares (%)"),
            p(style="font-size:10.5px;color:#868e96;margin:-2px 0 4px;",
              "Rows are auto-normalized to sum to 100."),
            tags$table(class="tut-vs-table",
              tags$thead(tags$tr(
                tags$th(""), tags$th("Left"), tags$th("Centre"), tags$th("Right")
              )),
              tags$tbody(
                tags$tr(
                  tags$td(style="font-size:11px;font-weight:600;padding-right:6px;","Smallia"),
                  tags$td(numericInput("tut_al",NULL,20,0,100,1,width="58px")),
                  tags$td(numericInput("tut_ac",NULL,30,0,100,1,width="58px")),
                  tags$td(numericInput("tut_ar",NULL,50,0,100,1,width="58px"))
                ),
                tags$tr(
                  tags$td(style="font-size:11px;font-weight:600;padding-right:6px;","Mediana"),
                  tags$td(numericInput("tut_bl",NULL,35,0,100,1,width="58px")),
                  tags$td(numericInput("tut_bc",NULL,40,0,100,1,width="58px")),
                  tags$td(numericInput("tut_br",NULL,25,0,100,1,width="58px"))
                ),
                tags$tr(
                  tags$td(style="font-size:11px;font-weight:600;padding-right:6px;","Largua"),
                  tags$td(numericInput("tut_cl",NULL,50,0,100,1,width="58px")),
                  tags$td(numericInput("tut_cc",NULL,35,0,100,1,width="58px")),
                  tags$td(numericInput("tut_cr",NULL,15,0,100,1,width="58px"))
                )
              )
            ),
            br(),
            div(class="sc-label", "Seats per district"),
            p(style="font-size:10.5px;color:#868e96;margin:-2px 0 4px;",
              "Proportional allocation would give 1 / 3 / 6. Equal seats (default) ",
              "over-represent Smallia (small, right-leaning)."),
            fluidRow(
              column(4, numericInput("tut_seats_a","Smallia",3,1,20,1)),
              column(4, numericInput("tut_seats_b","Mediana",3,1,20,1)),
              column(4, numericInput("tut_seats_c","Largua", 3,1,20,1))
            ),
            div(class="sc-label", "Electoral formula"),
            selectInput("tut_formula",NULL,
              choices=c("D'Hondt"="dhondt","Sainte-Laguë"="sainte_lague",
                        "Hare (LR)"="hare","Droop (LR)"="droop"),
              selected="dhondt"),
            div(class="sc-label", "Electoral threshold (%)"),
            sliderInput("tut_th",NULL,min=0,max=30,value=0,step=1),
            div(class="sc-hint",
              "Try 22% to exclude Left from Smallia and see how ThRB shifts IR.")
          )
        ),
        column(8, uiOutput("tut_out"))
      ),

      hr(style="margin:32px 0;"),

      h3(style="font-weight:700;border-bottom:2px solid #dee2e6;padding-bottom:6px;
                margin-bottom:14px;",
         "Interpreting the Bias Share"),
      p("Bias Share measures how much each institutional factor contributes to the",
        tags$em("total distortion"),
        "— the sum of all absolute biases. It is", tags$strong("not"),
        "the share of the net bias (TRB), but of |MRB| + |SRB| + |ThRB|."),
      p("This distinction matters when biases point in",
        tags$strong("opposite directions"),
        ": even if MRB and SRB partially cancel in TRB, both still receive a large
         Bias Share — correctly reflecting that each institution is actively
         distorting representation, just in opposite directions."),

      fluidRow(
        column(6,
          div(class="tut-step",
            div(style="font-weight:700;color:#198754;margin-bottom:8px;",
                "Scenario A — Additive biases (all rightward)"),
            tags$table(class="tut-table",
              tags$thead(tags$tr(
                tags$th("Component"),tags$th("Value"),
                tags$th("Direction"),tags$th("|Bias|"),tags$th("Share")
              )),
              tags$tbody(
                tags$tr(tags$td("MRB"),tags$td("+0.48"),
                  tags$td(style="color:#c0392b;font-weight:600;","rightward ▶"),
                  tags$td("0.48"),tags$td("64.0%")),
                tags$tr(tags$td("SRB"),tags$td("+0.15"),
                  tags$td(style="color:#c0392b;font-weight:600;","rightward ▶"),
                  tags$td("0.15"),tags$td("20.0%")),
                tags$tr(tags$td("ThRB"),tags$td("+0.12"),
                  tags$td(style="color:#c0392b;font-weight:600;","rightward ▶"),
                  tags$td("0.12"),tags$td("16.0%")),
                tags$tr(style="font-weight:700;background:#f1f3f5;",
                  tags$td("TRB"),tags$td("+0.75"),
                  tags$td("rightward ▶"),tags$td("0.75"),tags$td("—"))
              )
            ),
            div(class="tut-note",style="border-left-color:#198754;background:#edf7ed;",
              "All institutions amplify the same rightward bias. |TRB| equals the
               sum of all components. Bias Share identifies malapportionment
               as the dominant factor (64%).")
          )
        ),
        column(6,
          div(class="tut-step",
            div(style="font-weight:700;color:#e67e22;margin-bottom:8px;",
                "Scenario B — Opposing biases (MRB leftward, SRB/ThRB rightward)"),
            tags$table(class="tut-table",
              tags$thead(tags$tr(
                tags$th("Component"),tags$th("Value"),
                tags$th("Direction"),tags$th("|Bias|"),tags$th("Share")
              )),
              tags$tbody(
                tags$tr(tags$td("MRB"),tags$td("−0.35"),
                  tags$td(style="color:#2980b9;font-weight:600;","◀ leftward"),
                  tags$td("0.35"),tags$td("53.8%")),
                tags$tr(tags$td("SRB"),tags$td("+0.20"),
                  tags$td(style="color:#c0392b;font-weight:600;","rightward ▶"),
                  tags$td("0.20"),tags$td("30.8%")),
                tags$tr(tags$td("ThRB"),tags$td("+0.10"),
                  tags$td(style="color:#c0392b;font-weight:600;","rightward ▶"),
                  tags$td("0.10"),tags$td("15.4%")),
                tags$tr(style="font-weight:700;background:#f1f3f5;",
                  tags$td("TRB"),tags$td("−0.05"),
                  tags$td("◀ leftward"),tags$td("0.05"),tags$td("—"))
              )
            ),
            div(class="tut-note",
              "Malapportionment is leftward; seat allocation and threshold are
               rightward — they largely cancel. |TRB| = 0.05 yet total distortion
               is 0.65. Bias Share still sums to 100%, revealing each institution's
               separate contribution despite near-zero net bias.")
          )
        )
      ),
      p(style="font-size:13px;color:#495057;",
        tags$strong("Key takeaway: "),
        "Scenario B looks almost neutral on TRB alone. Bias Share reveals that
         malapportionment and the electoral formula are generating substantial
         opposing distortions that happen to offset each other. Institutions
         that cancel are not 'working well' — they create opposite pressures
         that can decouple when the political context changes."),

      hr(style="margin:32px 0;"),

      h3(style="font-weight:700;border-bottom:2px solid #dee2e6;padding-bottom:6px;
                margin-bottom:14px;",
         "Methodological Notes"),
      fluidRow(
        column(6,
          div(class="tut-step",
            h5(style="font-weight:700;margin-top:0;","Simulation design"),
            tags$ul(style="font-size:13px;line-height:1.75;padding-left:18px;",
              tags$li(tags$strong("Country structure: "),
                "100 districts, each with the same district magnitude (DM). This
                 isolates the pure mechanical effect of DM from other sources of
                 variation."),
              tags$li(tags$strong("Party system: "),
                "2–6 parties per simulation (10/20/30/25/15% probabilities). Party
                 ideologies follow a ", tags$em("Beta(2.5, 2.5) × 10"),
                " distribution — moderate positions are most likely; extreme positions
                 have low probability."),
              tags$li(tags$strong("Vote shares: "),
                "Each district's vote shares are drawn independently from a ",
                tags$em("Dirichlet(1.5)"), " distribution, implying no spatial
                 correlation between districts."),
              tags$li(tags$strong("Seat allocation: "),
                "D'Hondt by default. The electoral threshold is applied uniformly
                 to all districts at the same level before seat allocation.")
            )
          )
        ),
        column(6,
          div(class="tut-step",
            h5(style="font-weight:700;margin-top:0;","Malapportionment scenarios"),
            tags$ul(style="font-size:13px;line-height:1.75;padding-left:18px;",
              tags$li(tags$strong("Random: "),
                "Over-represented districts are selected randomly with no relationship
                 to their ideological composition. The expected MRB is zero."),
              tags$li(tags$strong("Concentrated: "),
                "Either the left-leaning or the right-leaning bloc is systematically
                 over-represented — a logically extreme case that maximises MRB for
                 a given MAL level."),
              tags$li(
                "In practice, malapportionment rarely benefits one ideological bloc
                 consistently. Beramendi, Boix, Guinjoan & Rogers (2021) show that
                 institutional biases are more complex and context-dependent than
                 either scenario implies. The two scenarios therefore bound the
                 realistic distribution of MRB.")
            )
          )
        )
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

  # ==========================================================================
  # TUTORIAL TAB
  # ==========================================================================

  tut_calc <- reactive({
    vs_raw <- matrix(c(
      max(0, input$tut_al), max(0, input$tut_ac), max(0, input$tut_ar),
      max(0, input$tut_bl), max(0, input$tut_bc), max(0, input$tut_br),
      max(0, input$tut_cl), max(0, input$tut_cc), max(0, input$tut_cr)
    ), nrow=3, ncol=3, byrow=TRUE)
    rs <- rowSums(vs_raw)
    req(all(rs > 0))
    vs <- vs_raw / rs

    pop   <- c(100000, 300000, 600000)
    ideol <- c(2.0, 5.0, 8.0)
    seats <- c(max(1L, as.integer(input$tut_seats_a)),
               max(1L, as.integer(input$tut_seats_b)),
               max(1L, as.integer(input$tut_seats_c)))
    th    <- input$tut_th / 100
    fn    <- get_alloc_fn(input$tut_formula)

    IS_d   <- as.numeric(vs %*% ideol)
    IP     <- weighted.mean(IS_d, w=pop)
    IS_nat <- weighted.mean(IS_d, w=seats)
    MRB    <- IS_nat - IP

    seats_b <- do.call(rbind, lapply(1:3, function(d) fn(vs[d,], seats[d])))
    nat_b   <- colSums(seats_b)
    IR_base <- if (sum(nat_b) > 0) weighted.mean(ideol, w=nat_b) else NA_real_
    SRB     <- if (!is.na(IR_base)) IR_base - IS_nat else NA_real_

    seats_t <- do.call(rbind, lapply(1:3, function(d) {
      vf <- ifelse(vs[d,] < th, 0, vs[d,])
      if (sum(vf) == 0) vf <- vs[d,]
      fn(vf, seats[d])
    }))
    nat_t  <- colSums(seats_t)
    IR_th  <- if (sum(nat_t) > 0) weighted.mean(ideol, w=nat_t) else NA_real_
    ThRB   <- if (!is.na(IR_th) && !is.na(IR_base)) IR_th - IR_base else NA_real_
    TRB    <- if (!is.na(IR_th)) IR_th - IP else NA_real_

    aM  <- abs(MRB)
    aS  <- if (!is.na(SRB))  abs(SRB)  else 0
    aT  <- if (!is.na(ThRB)) abs(ThRB) else 0
    tot <- aM + aS + aT

    list(vs=vs, seats=seats, IS_d=IS_d, IP=IP, IS_nat=IS_nat, MRB=MRB,
         seats_b=seats_b, nat_b=nat_b, IR_base=IR_base, SRB=SRB,
         seats_t=seats_t, nat_t=nat_t, IR_th=IR_th, ThRB=ThRB, TRB=TRB,
         th=th, formula=input$tut_formula,
         BS_MRB  = if (tot>1e-9) aM/tot  else NA_real_,
         BS_SRB  = if (tot>1e-9) aS/tot  else NA_real_,
         BS_ThRB = if (tot>1e-9) aT/tot  else NA_real_,
         aM=aM, aS=aS, aT=aT, tot=tot)
  })

  output$tut_out <- renderUI({
    tc <- tut_calc()

    fmt  <- function(x) if (is.na(x)) "—" else sprintf("%.3f", x)
    sgn  <- function(x) if (is.na(x)) "—" else if (x >= 0) sprintf("+%.3f",x) else sprintf("%.3f",x)
    pct  <- function(x) if (is.na(x)) "—" else sprintf("%.1f%%", 100*x)
    dlbl <- function(x) {
      if (is.na(x) || abs(x) < 0.005) "≈ 0"
      else if (x > 0) "rightward ▶" else "◀ leftward"
    }
    dcol <- function(x) {
      if (is.na(x) || abs(x) < 0.005) "#868e96"
      else if (x > 0) "#c0392b" else "#2980b9"
    }

    pop_v <- c(100000, 300000, 600000)
    nm    <- c("Smallia", "Mediana", "Largua")

    # ---- District breakdown table ----
    th_col_lbl <- if (tc$th > 0)
      paste0("After TH=", round(100*tc$th), "% (L|C|R)")
    else "After threshold"

    dist_rows <- lapply(1:3, function(d) {
      vs_pct <- paste0(round(100*tc$vs[d,1]),"% / ",
                       round(100*tc$vs[d,2]),"% / ",
                       round(100*tc$vs[d,3]),"%")
      sb <- paste0(tc$seats_b[d,1],"|",tc$seats_b[d,2],"|",tc$seats_b[d,3])
      st <- if (tc$th > 0)
              paste0(tc$seats_t[d,1],"|",tc$seats_t[d,2],"|",tc$seats_t[d,3])
            else "—"
      tags$tr(
        tags$td(style="font-weight:600;", nm[d]),
        tags$td(format(pop_v[d],big.mark=",")),
        tags$td(tc$seats[d]),
        tags$td(vs_pct),
        tags$td(sprintf("%.2f", tc$IS_d[d])),
        tags$td(style="font-family:monospace;", sb),
        tags$td(style="font-family:monospace;", st)
      )
    })
    nat_row <- tags$tr(style="font-weight:700;background:#f1f3f5;",
      tags$td("National"), tags$td("1,000,000"),
      tags$td(sum(tc$seats)), tags$td("—"),
      tags$td(sprintf("%.2f", tc$IS_nat)),
      tags$td(style="font-family:monospace;",
              paste0(tc$nat_b[1],"|",tc$nat_b[2],"|",tc$nat_b[3])),
      tags$td(style="font-family:monospace;",
              if (tc$th > 0)
                paste0(tc$nat_t[1],"|",tc$nat_t[2],"|",tc$nat_t[3])
              else "—")
    )

    dist_tbl <- div(style="overflow-x:auto;margin-bottom:12px;",
      tags$table(class="tut-table", style="font-size:12px;",
        tags$thead(tags$tr(
          tags$th("District"), tags$th("Population"), tags$th("Seats"),
          tags$th("Votes L/C/R"), tags$th("IS"),
          tags$th("Seats D'Hondt (L|C|R)"), tags$th(th_col_lbl)
        )),
        tags$tbody(c(dist_rows, list(nat_row)))
      )
    )

    # ---- MAL index note ----
    mal_idx <- 0.5 * sum(abs(tc$seats/sum(tc$seats) - pop_v/sum(pop_v)))
    mal_note <- div(style="font-size:11px;color:#6c757d;margin-bottom:10px;",
      sprintf("Malapportionment index (MAL) = %.3f  |  0 = perfectly proportional · 0.5 = extreme", mal_idx)
    )

    # ---- Metric cards ----
    card <- function(lbl, val, cls, use_sgn=TRUE) {
      v_txt <- if (use_sgn) sgn(val) else fmt(val)
      div(class=paste("metric-card", cls),
        div(class="mc-label", lbl),
        div(class="mc-val", v_txt)
      )
    }

    metrics <- div(style="display:flex;flex-wrap:wrap;gap:0;margin:10px 0 8px;",
      card("IP  (voters)",      tc$IP,      "mc-ip",   FALSE),
      card("IS  (districts)",   tc$IS_nat,  "mc-is",   FALSE),
      div(style="width:100%;height:0;flex-basis:100%;"),
      card("MRB = IS − IP",     tc$MRB,     "mc-mrb"),
      card("IR  (legislators)", tc$IR_base, "mc-ir",   FALSE),
      card("SRB = IR − IS",     tc$SRB,     "mc-srb"),
      card("ThRB",              tc$ThRB,    "mc-thrb"),
      div(style="width:100%;height:0;flex-basis:100%;"),
      card("TRB = IR₁ − IP", tc$TRB,  "mc-trb")
    )

    # ---- Bias Share bar ----
    bs_ui <- if (tc$tot > 1e-9) {
      bsm <- max(0, round(100 * tc$BS_MRB))
      bss <- max(0, round(100 * tc$BS_SRB))
      bst <- max(0, 100 - bsm - bss)
      mk_seg <- function(w, bg, lbl) {
        if (w <= 0) return(NULL)
        div(style=paste0("width:",w,"%;background:",bg,";display:flex;align-items:center;
                          justify-content:center;"),
            if (w > 6) span(style="color:#fff;font-size:10px;font-weight:700;", lbl) else NULL)
      }
      div(style="margin-bottom:14px;",
        div(style="font-size:11px;font-weight:700;color:#6c757d;text-transform:uppercase;
                   letter-spacing:.04em;margin-bottom:4px;",
            "Bias Share  —  % of total |distortion|"),
        div(style="display:flex;height:26px;border-radius:4px;overflow:hidden;border:1px solid #dee2e6;",
            mk_seg(bsm, "#c0392b", paste0(bsm,"%")),
            mk_seg(bss, "#2980b9", paste0(bss,"%")),
            mk_seg(bst, "#e67e22", paste0(bst,"%"))
        ),
        div(style="display:flex;font-size:11px;font-weight:600;margin-top:3px;",
          div(style="flex:1;color:#c0392b;", paste0("MRB:  ", pct(tc$BS_MRB))),
          div(style="flex:1;color:#2980b9;", paste0("SRB:  ", pct(tc$BS_SRB))),
          div(style="flex:1;color:#e67e22;", paste0("ThRB: ", pct(tc$BS_ThRB)))
        )
      )
    } else {
      p(style="color:#868e96;font-style:italic;font-size:12px;",
        "All biases are zero — Bias Share is undefined.")
    }

    # ---- Interpretation ----
    ir_lbl <- if (!is.na(tc$IR_th) && tc$th > 0) fmt(tc$IR_th) else fmt(tc$IR_base)
    interp <- div(style="background:#f8f9fa;border:1px solid #dee2e6;border-radius:6px;
                          padding:10px 14px;font-size:12.5px;margin-top:2px;",
      tags$strong("Reading the numbers:"),
      tags$ul(style="margin:4px 0 0;",
        tags$li(
          tags$b(paste0("MRB = ", sgn(tc$MRB))),
          tags$span(style=paste0("color:",dcol(tc$MRB),";font-weight:600;"),
                    paste0(" (", dlbl(tc$MRB), ")")),
          " — IS (district-seat-weighted mean) is ", fmt(tc$IS_nat),
          "; IP (population-weighted mean) is ", fmt(tc$IP), "."
        ),
        tags$li(
          tags$b(paste0("SRB = ", sgn(tc$SRB))),
          tags$span(style=paste0("color:",dcol(tc$SRB),";font-weight:600;"),
                    paste0(" (", dlbl(tc$SRB), ")")),
          " — D'Hondt gives L=", tc$nat_b[1], " / C=", tc$nat_b[2],
          " / R=", tc$nat_b[3], " national seats. IR = ", fmt(tc$IR_base), "."
        ),
        if (tc$th > 0)
          tags$li(
            tags$b(paste0("ThRB = ", sgn(tc$ThRB))),
            tags$span(style=paste0("color:",dcol(tc$ThRB),";font-weight:600;"),
                      paste0(" (", dlbl(tc$ThRB), ")")),
            " — Applying ", round(100*tc$th), "% threshold shifts IR from ",
            fmt(tc$IR_base), " to ", fmt(tc$IR_th), "."
          )
        else
          tags$li("ThRB = 0 — no threshold applied (move the slider to see its effect).")
      )
    )

    tagList(mal_note, dist_tbl, metrics, bs_ui, interp)
  })
}

shinyApp(ui, server)
