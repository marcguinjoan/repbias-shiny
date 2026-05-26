# =============================================================================
# Pool B: SRB/ThRB district pool — 50k sims × 11 fixed DM values
# Output: 7_Shiny/data/pool_srb_districts.rds  (~50 MB)
# All districts in each simulation have the SAME fixed DM.
# No malapportionment — pure mechanical effect of DM.
# =============================================================================

library(pacman)
p_load(tidyverse, gtools, doParallel, foreach)

set.seed(99)

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
N_SIM_PER_DM   <- 4500    # sims per DM level this batch (4500 × 11 ≈ 50k total)
BATCH_OFFSET   <- 0        # change to 50000 / 150000 for next batches
N_DISTRICTS     <- 100
DISTRICTS_SAVED <- 10
DIRICHLET_ALPHA <- 1.5
BETA_SHAPE      <- 2.5
TH              <- 0.2
DM_VALUES       <- 1:30    # all DM values 1–30

# ---------------------------------------------------------------------------
# 2. D'Hondt function (same as Pool A)
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
    if (length(winner_index) > 1) winner_index <- winner_index[1]
    seats_vector[winner_index] <- seats_vector[winner_index] + 1
  }
  seats_vector
}

# ---------------------------------------------------------------------------
# 3. Worker: one simulation at fixed DM
# ---------------------------------------------------------------------------
run_one_sim_fixed_dm <- function(i, dm_fixed) {
  set.seed(i * 1000 + dm_fixed)   # unique seed per (sim, dm) pair
  n_p             <- sample(2:6, 1, prob = c(0.1, 0.2, 0.3, 0.25, 0.15))
  party_ideologies <- rbeta(n_p, BETA_SHAPE, BETA_SHAPE) * SCALE_FACTOR

  # Vote shares: Dirichlet per district
  vote_shares_matrix <- gtools::rdirichlet(N_DISTRICTS, rep(DIRICHLET_ALPHA, n_p))

  # District-level IS (ideology of voters, pre-conversion)
  IS_vec <- as.numeric(vote_shares_matrix %*% party_ideologies)

  # D'Hondt with fixed DM seats for every district
  seats_list <- lapply(seq_len(N_DISTRICTS), function(d)
    distribute_seats_dhondt_multiparty(vote_shares_matrix[d, ], dm_fixed)
  )
  seats_mat   <- do.call(rbind, seats_list)
  seats_total <- rowSums(seats_mat)

  IR_base_vec <- ifelse(seats_total > 0,
                        as.numeric(seats_mat %*% party_ideologies) / seats_total,
                        NA_real_)
  SRB_vec <- IR_base_vec - IS_vec

  # Threshold bias
  max_share_vec <- apply(vote_shares_matrix, 1, max)
  TH_USED_vec   <- runif(N_DISTRICTS, 0, pmin(max_share_vec, TH))
  below_th      <- sweep(vote_shares_matrix, 1, TH_USED_vec, "<")
  seats_after   <- seats_mat * (!below_th)
  seats_after_t <- rowSums(seats_after)
  IR_th_vec     <- ifelse(seats_after_t > 0,
                          as.numeric(seats_after %*% party_ideologies) / seats_after_t,
                          NA_real_)
  ThBias_vec <- (IR_th_vec - IS_vec) - SRB_vec

  # National polarization and ENEP (equal-weight districts, no pop variation)
  nat_vs  <- colMeans(vote_shares_matrix)[seq_len(n_p)]
  nat_vs  <- nat_vs / sum(nat_vs)
  mean_i  <- sum(nat_vs * party_ideologies)
  pol     <- sum(nat_vs * abs(party_ideologies - mean_i))
  enep    <- 1 / sum(nat_vs^2)

  # Sample DISTRICTS_SAVED random districts
  keep <- sample(N_DISTRICTS, min(DISTRICTS_SAVED, N_DISTRICTS))

  data.frame(
    DM_fixed  = dm_fixed,
    Simulacion = i,
    n_parties  = n_p,
    polarization = pol,
    enep         = enep,
    SRB_distrito_base = SRB_vec[keep],
    TH_USED           = TH_USED_vec[keep],
    Threshold_Bias    = ThBias_vec[keep]
  )
}

# ---------------------------------------------------------------------------
# 4. Parallel setup
# ---------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 3)
cl      <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
cat("Cores:", n_cores, "\n")

# ---------------------------------------------------------------------------
# 5. Run all DM levels
# ---------------------------------------------------------------------------
all_results <- list()

parallel::clusterExport(cl, c(
  "N_DISTRICTS", "DISTRICTS_SAVED", "SCALE_FACTOR",
  "DIRICHLET_ALPHA", "BETA_SHAPE", "TH",
  "distribute_seats_dhondt_multiparty", "run_one_sim_fixed_dm"
))

for (dm in DM_VALUES) {
  cat("--- DM =", dm, "(N =", N_SIM_PER_DM, ", offset =", BATCH_OFFSET, ") ---\n")
  parallel::clusterExport(cl, "dm")

  res <- foreach(
    i         = seq_len(N_SIM_PER_DM),
    .packages = c("gtools")
  ) %dopar% {
    run_one_sim_fixed_dm(i + BATCH_OFFSET, dm)
  }

  all_results[[as.character(dm)]] <- do.call(rbind, res)
  cat("  Done. Rows:", nrow(all_results[[as.character(dm)]]), "\n")
}

parallel::stopCluster(cl)

# ---------------------------------------------------------------------------
# 6. Combine and append
# ---------------------------------------------------------------------------
pool_new <- do.call(rbind, all_results) %>%
  mutate(DM_fixed = as.integer(DM_fixed))

out_path <- file.path(out_dir, "pool_srb_districts.rds")

if (file.exists(out_path)) {
  pool_existing <- readRDS(out_path)
  pool_srb <- bind_rows(pool_existing, pool_new)
  cat("Appended", nrow(pool_new), "rows to existing pool.\n")
} else {
  pool_srb <- pool_new
}

saveRDS(pool_srb, out_path, compress = "xz")
cat("Saved:", out_path, "\n")
cat("Total rows:", nrow(pool_srb), "| File size:", round(file.size(out_path) / 1e6, 2), "MB\n")
cat("DM values:", sort(unique(pool_srb$DM_fixed)), "\n")
print(head(pool_srb))
