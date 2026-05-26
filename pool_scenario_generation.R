# =============================================================================
# Pool C: Scenario Designer pool — 3-institution bias decomposition
#
# Varies per simulation:
#   DM_fixed    — all 100 districts have the same number of seats
#   f_val       — malapportionment severity factor (1–15) → determines MAL
#   TH_applied  — fixed national threshold (0–20%), applied uniformly
#
# Output: 7_Shiny/data/pool_scenario.rds  (~3–8 MB compressed)
# Estimated run time: 30–60 minutes on a multi-core machine
#
# Pool size: 11 DM values × 5,000 sims = 55,000 rows
# =============================================================================

library(pacman)
p_load(tidyverse, gtools, doParallel, foreach)

# ---------------------------------------------------------------------------
# 0. Output directory
# ---------------------------------------------------------------------------
out_dir <- tryCatch({
  file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
}, error = function(e) {
  args      <- commandArgs(trailingOnly = FALSE)
  file_flag <- args[grepl("--file=", args)]
  if (length(file_flag) > 0) {
    file.path(dirname(normalizePath(sub("--file=", "", file_flag))), "data")
  } else {
    file.path(getwd(), "data")
  }
})
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
cat("Output directory:", out_dir, "\n")

# ---------------------------------------------------------------------------
# 1. Constants
# ---------------------------------------------------------------------------
SCALE_FACTOR    <- 10
N_SIM_PER_DM   <- 5000     # sims per DM level → 11 × 5,000 = 55,000 total
N_DISTRICTS     <- 100
DIRICHLET_ALPHA <- 1.5
BETA_SHAPE      <- 2.5
DM_VALUES       <- c(1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30)

# ---------------------------------------------------------------------------
# 2. D'Hondt seat allocation
# ---------------------------------------------------------------------------
distribute_seats_dhondt_multiparty <- function(votes_vector, num_total_seats) {
  num_parties <- length(votes_vector)
  if (num_total_seats == 0) return(rep(0, num_parties))
  votes_vector[is.na(votes_vector) | votes_vector < 0] <- 0
  if (sum(votes_vector) == 0) return(rep(0, num_parties))
  seats_vector <- rep(0, num_parties)
  for (i in 1:num_total_seats) {
    quotients    <- votes_vector / (seats_vector + 1)
    winner_index <- which.max(quotients)
    seats_vector[winner_index] <- seats_vector[winner_index] + 1
  }
  seats_vector
}

# ---------------------------------------------------------------------------
# 3. Single-simulation worker
#    Returns one summary row with all three bias components and their shares.
# ---------------------------------------------------------------------------
run_one_scenario_sim <- function(i, dm_fixed, f_val, th_applied) {
  set.seed(i * 37 + dm_fixed * 13)

  n_p              <- sample(2:6, 1, prob = c(0.1, 0.2, 0.3, 0.25, 0.15))
  party_ideologies <- rbeta(n_p, BETA_SHAPE, BETA_SHAPE) * SCALE_FACTOR

  # Variable district populations (base for MAL calculation)
  pop_base <- sample(100:3000, N_DISTRICTS, replace = TRUE)

  # Random malapportionment: shrink half the districts by f_val, rescale
  sel_idx <- sample(N_DISTRICTS, floor(N_DISTRICTS / 2))
  pop_adj <- ifelse(seq_len(N_DISTRICTS) %in% sel_idx,
                    pop_base / max(f_val, 1), pop_base)
  pop_adj <- pop_adj * (sum(pop_base) / sum(pop_adj))

  # Vote shares (Dirichlet, independent of population)
  vs_mat <- gtools::rdirichlet(N_DISTRICTS, rep(DIRICHLET_ALPHA, n_p))

  # District ideology of voters (IS at district level)
  IS_vec <- as.numeric(vs_mat %*% party_ideologies)

  # ---- MRB --------------------------------------------------------
  # All districts have dm_fixed seats → uniform seat weight → IS = mean(IS_vec)
  IP  <- weighted.mean(IS_vec, w = pop_adj)
  IS  <- mean(IS_vec)
  MRB <- IS - IP

  # MAL index (Samuels & Snyder)
  seats_share <- rep(1 / N_DISTRICTS, N_DISTRICTS)   # equal since dm_fixed
  pop_share   <- pop_adj / sum(pop_adj)
  MAL         <- 0.5 * sum(abs(seats_share - pop_share))

  # ---- SRB (D'Hondt, no threshold) ---------------------------------
  seats_base_list <- lapply(seq_len(N_DISTRICTS), function(d)
    distribute_seats_dhondt_multiparty(vs_mat[d, ], dm_fixed)
  )
  seats_base_mat <- do.call(rbind, seats_base_list)
  nat_seats_base <- colSums(seats_base_mat, na.rm = TRUE)
  total_S_base   <- sum(nat_seats_base)
  IR_base <- if (total_S_base > 0)
    weighted.mean(party_ideologies, w = nat_seats_base) else NA_real_
  SRB <- if (!is.na(IR_base)) IR_base - IS else NA_real_

  # ---- ThRB (fixed national threshold applied at district level) ---
  seats_th_list <- lapply(seq_len(N_DISTRICTS), function(d) {
    vs_filtered <- ifelse(vs_mat[d, ] < th_applied, 0, vs_mat[d, ])
    distribute_seats_dhondt_multiparty(vs_filtered, dm_fixed)
  })
  seats_th_mat <- do.call(rbind, seats_th_list)
  nat_seats_th <- colSums(seats_th_mat, na.rm = TRUE)
  total_S_th   <- sum(nat_seats_th)
  IR_th <- if (total_S_th > 0)
    weighted.mean(party_ideologies, w = nat_seats_th) else NA_real_
  ThRB <- if (!is.na(IR_th) && !is.na(IR_base)) IR_th - IR_base else NA_real_

  TRB <- if (!is.na(IR_th)) IR_th - IP else NA_real_

  # ---- 3-way Bias Shares -------------------------------------------
  abs_MRB  <- abs(MRB)
  abs_SRB  <- if (!is.na(SRB))  abs(SRB)  else 0
  abs_ThRB <- if (!is.na(ThRB)) abs(ThRB) else 0
  abs_total <- abs_MRB + abs_SRB + abs_ThRB

  BS_MRB  <- if (abs_total > 1e-9) abs_MRB  / abs_total else NA_real_
  BS_SRB  <- if (abs_total > 1e-9) abs_SRB  / abs_total else NA_real_
  BS_ThRB <- if (abs_total > 1e-9) abs_ThRB / abs_total else NA_real_

  # ---- Polarization & ENEP -----------------------------------------
  nat_vs <- colMeans(vs_mat)
  nat_vs <- nat_vs / sum(nat_vs)
  mean_i <- sum(nat_vs * party_ideologies)
  pol    <- sum(nat_vs * abs(party_ideologies - mean_i))
  enep   <- 1 / sum(nat_vs^2)

  data.frame(
    DM_fixed     = dm_fixed,
    f            = f_val,
    MAL          = MAL,
    TH_applied   = th_applied,
    n_parties    = n_p,
    polarization = pol,
    enep         = enep,
    MRB          = MRB,
    SRB          = SRB,
    ThRB         = ThRB,
    TRB          = TRB,
    BS_MRB       = BS_MRB,
    BS_SRB       = BS_SRB,
    BS_ThRB      = BS_ThRB
  )
}

# ---------------------------------------------------------------------------
# 4. Parallel setup
# ---------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 3)
cl      <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
cat("Cores:", n_cores, "\n")

parallel::clusterExport(cl, c(
  "N_DISTRICTS", "SCALE_FACTOR", "DIRICHLET_ALPHA", "BETA_SHAPE",
  "distribute_seats_dhondt_multiparty", "run_one_scenario_sim"
))

# ---------------------------------------------------------------------------
# 5. Run all DM levels
# ---------------------------------------------------------------------------
all_results <- list()

for (dm in DM_VALUES) {
  cat("--- DM =", dm, "(N =", N_SIM_PER_DM, ") ---\n")

  set.seed(dm * 777)
  f_vals  <- runif(N_SIM_PER_DM, 1, 15)
  th_vals <- runif(N_SIM_PER_DM, 0, 0.20)

  parallel::clusterExport(cl, c("dm", "f_vals", "th_vals"))

  res <- foreach(
    i         = seq_len(N_SIM_PER_DM),
    .packages = "gtools"
  ) %dopar% {
    run_one_scenario_sim(i, dm, f_vals[i], th_vals[i])
  }

  all_results[[as.character(dm)]] <- do.call(rbind, res)
  cat("  Done. Rows:", nrow(all_results[[as.character(dm)]]), "\n")
}

parallel::stopCluster(cl)

# ---------------------------------------------------------------------------
# 6. Combine and save
# ---------------------------------------------------------------------------
pool_scenario <- do.call(rbind, all_results) %>%
  mutate(DM_fixed = as.integer(DM_fixed))

out_path <- file.path(out_dir, "pool_scenario.rds")
saveRDS(pool_scenario, out_path, compress = "xz")

cat("\nSaved:", out_path, "\n")
cat("Total rows:", nrow(pool_scenario), "\n")
cat("File size:", round(file.size(out_path) / 1e6, 2), "MB\n")
cat("DM values simulated:", sort(unique(pool_scenario$DM_fixed)), "\n")
cat("\nColumn summary:\n")
print(summary(pool_scenario[, c("DM_fixed","MAL","TH_applied","MRB","SRB","ThRB","TRB",
                                 "BS_MRB","BS_SRB","BS_ThRB")]))
