# =============================================================================
# Pool A: MRB summary pool — incremental batches (random + concentrated)
# Output: 7_Shiny/data/pool_mrb_summary.rds
# Run with BATCH_OFFSET = 0 (batch 1), 50000 (batch 2), 150000 (batch 3)…
# Each run APPENDS to the existing RDS — no duplicates guaranteed by offset.
# =============================================================================

library(pacman)
p_load(tidyverse, gtools, sf, doParallel, foreach)

set.seed(42)

# ---------------------------------------------------------------------------
# 0. Output directory
# ---------------------------------------------------------------------------
# POOL_MRB_OUT_DIR can be set by a calling script (e.g. master_batch.R) to
# guarantee the correct output path regardless of the active editor document.
out_dir <- if (exists("POOL_MRB_OUT_DIR") && nzchar(POOL_MRB_OUT_DIR)) {
  POOL_MRB_OUT_DIR
} else {
  tryCatch({
    p <- rstudioapi::getSourceEditorContext()$path
    # Reject master_batch.R path — it means we were sourced from there
    if (grepl("master_batch", p, ignore.case = TRUE)) stop("sourced")
    file.path(dirname(p), "data")
  }, error = function(e) {
    args      <- commandArgs(trailingOnly = FALSE)
    file_flag <- args[grepl("--file=", args)]
    if (length(file_flag) > 0) {
      file.path(dirname(normalizePath(sub("--file=", "", file_flag))), "data")
    } else {
      file.path(getwd(), "data")
    }
  })
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
cat("Output directory:", out_dir, "\n")

# ---------------------------------------------------------------------------
# 1. Shared constants (must match Part1_Rsimulations_v2.rmd)
# ---------------------------------------------------------------------------
SCALE_FACTOR   <- 10
N_PER_SCENARIO <- 200000   # sims per scenario this batch (200k × 2 = 400k total)
BATCH_OFFSET   <- 100000   # previous batches used 0 and 50000 → no seed overlap
SEAT_DIVISOR   <- 99
DIRICHLET_ALPHA <- 1.5
BETA_SHAPE     <- 2.5
TH             <- 0.2      # max legal threshold

# ---------------------------------------------------------------------------
# 2. Base grid (100 districts, population 100-3000)
# ---------------------------------------------------------------------------
set.seed(456)
lado_total <- 10000; cellsize <- 1000
grid_base_raw <- st_make_grid(
  st_bbox(c(xmin = -lado_total/2, xmax = lado_total/2,
            ymin = -lado_total/2, ymax = lado_total/2)),
  cellsize = cellsize, crs = 3857
)
grid_base <- grid_base_raw %>%
  st_sf(grid_id = 1:length(.), geometry = .) %>%
  mutate(poblacion = sample(100:3000, n(), replace = TRUE))

geometria_base <- st_geometry(grid_base)

# ---------------------------------------------------------------------------
# 3. Simulation functions (copied from Part1_Rsimulations_v2.rmd)
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

calculate_national_metrics <- function(data, party_ideologies,
                                       escanos_col = "escanos",
                                       pop_col     = "poblacion",
                                       is_col      = "IS_distrito") {
  n_p  <- length(party_ideologies)
  pop  <- data[[pop_col]]
  escs <- data[[escanos_col]]
  IS_d <- data[[is_col]]

  IP  <- weighted.mean(IS_d, w = pop,  na.rm = TRUE)
  IS  <- weighted.mean(IS_d, w = escs, na.rm = TRUE)

  seat_cols <- grep("^escanos_P", names(data), value = TRUE)
  total_seats_by_party <- colSums(data[, seat_cols, drop = FALSE], na.rm = TRUE)
  total_S <- sum(total_seats_by_party)
  IR <- if (total_S > 0) weighted.mean(party_ideologies, w = total_seats_by_party) else NA_real_

  MRB <- IS - IP
  SRB <- IR - IS
  TRB <- IR - IP
  MAL <- 0.5 * sum(abs(escs / sum(escs) - pop / sum(pop)), na.rm = TRUE)

  abs_total <- abs(TRB)
  BIASSHARE_MAL  <- if (!is.na(abs_total) && abs_total > 1e-9) abs(MRB) / (abs(MRB) + abs(SRB)) else NA_real_
  BIASSHARE_SEAT <- if (!is.na(abs_total) && abs_total > 1e-9) abs(SRB) / (abs(MRB) + abs(SRB)) else NA_real_

  list(MAL = MAL, MRB = MRB, SRB = SRB, TRB = TRB,
       BIASSHARE_MAL = BIASSHARE_MAL, BIASSHARE_SEAT = BIASSHARE_SEAT,
       IP = IP, IS = IS, IR = IR)
}

simulate_random_malapp <- function(grid_sim_base, f_value) {
  n_rows  <- nrow(grid_sim_base)
  sel_idx <- sample(n_rows, size = floor(n_rows / 2))
  pop_adj <- ifelse(seq_len(n_rows) %in% sel_idx,
                    grid_sim_base$poblacion_original / max(f_value, 1),
                    grid_sim_base$poblacion_original)
  p_tot <- sum(grid_sim_base$poblacion_original)
  pop_adj <- pop_adj * (p_tot / sum(pop_adj))
  grid_sim_base %>%
    mutate(poblacion = pop_adj,
           escanos   = pmax(1, floor(poblacion_original / SEAT_DIVISOR)))
}

simulate_ideological_malapp <- function(grid_sim_base, f_value,
                                        favor_left_bloc, party_ideologies) {
  n_p <- length(party_ideologies)
  vs_cols <- grep("^vote_share_P", names(grid_sim_base), value = TRUE)[seq_len(n_p)]
  vs_mat  <- as.matrix(st_drop_geometry(grid_sim_base)[, vs_cols])
  cell_ideol <- as.numeric(vs_mat %*% party_ideologies)

  mid <- SCALE_FACTOR / 2
  is_left <- cell_ideol < mid
  pop_adj <- ifelse(
    if (favor_left_bloc) is_left else !is_left,
    grid_sim_base$poblacion_original / max(f_value, 1),
    grid_sim_base$poblacion_original
  )
  p_tot <- sum(grid_sim_base$poblacion_original)
  pop_adj <- pop_adj * (p_tot / sum(pop_adj))
  grid_sim_base %>%
    mutate(poblacion = pop_adj,
           escanos   = pmax(1, floor(poblacion_original / SEAT_DIVISOR)))
}

# ---------------------------------------------------------------------------
# 4. Single-simulation worker (returns one summary row)
# ---------------------------------------------------------------------------
run_one_sim <- function(i, f_val, scenario, grid_base, geometria_base) {
  set.seed(i)
  n_p             <- sample(2:6, 1, prob = c(0.1, 0.2, 0.3, 0.25, 0.15))
  party_ideologies <- rbeta(n_p, BETA_SHAPE, BETA_SHAPE) * SCALE_FACTOR
  vote_shares_matrix <- gtools::rdirichlet(nrow(grid_base), rep(DIRICHLET_ALPHA, n_p))

  grid_sim <- grid_base %>%
    st_drop_geometry() %>%
    mutate(poblacion_original = poblacion)
  for (p in seq_len(n_p))
    grid_sim[[paste0("vote_share_P", p)]] <- vote_shares_matrix[, p]

  # Apply malapportionment
  grid_sim_sf <- st_sf(grid_sim, geometry = geometria_base)
  if (scenario == "random") {
    grid_adj <- simulate_random_malapp(grid_sim_sf, f_val)
  } else {
    favor_left <- runif(1) > 0.5
    grid_adj   <- simulate_ideological_malapp(grid_sim_sf, f_val, favor_left, party_ideologies)
  }

  # D'Hondt seat allocation
  votes_mat <- vote_shares_matrix * grid_adj$poblacion
  seats_list <- mapply(
    distribute_seats_dhondt_multiparty,
    split(as.data.frame(votes_mat), seq(nrow(votes_mat))),
    grid_adj$escanos,
    SIMPLIFY = FALSE
  )
  seats_df <- as.data.frame(do.call(rbind, seats_list))
  colnames(seats_df) <- paste0("escanos_P", seq_len(n_p))

  # District-level IS
  IS_vec <- as.numeric(vote_shares_matrix %*% party_ideologies)

  # Vectorised SRB
  seats_mat   <- as.matrix(seats_df)
  seats_total <- rowSums(seats_mat)
  IR_base_vec <- ifelse(seats_total > 0,
                        as.numeric(seats_mat %*% party_ideologies) / seats_total,
                        NA_real_)

  # Threshold bias
  max_share_vec <- apply(vote_shares_matrix, 1, max)
  TH_USED_vec   <- runif(nrow(grid_adj), 0, pmin(max_share_vec, TH))
  below_th      <- sweep(vote_shares_matrix, 1, TH_USED_vec, "<")
  seats_after   <- seats_mat * (!below_th)
  seats_after_t <- rowSums(seats_after)
  IR_th_vec     <- ifelse(seats_after_t > 0,
                          as.numeric(seats_after %*% party_ideologies) / seats_after_t,
                          NA_real_)
  ThBias_vec    <- (IR_th_vec - IS_vec) - (IR_base_vec - IS_vec)

  # Assemble full district data (temporary, for national metrics)
  sim_data <- st_drop_geometry(grid_adj) %>%
    bind_cols(seats_df) %>%
    mutate(IS_distrito = IS_vec)

  nm <- calculate_national_metrics(
    data             = sim_data,
    party_ideologies = party_ideologies,
    escanos_col      = "escanos",
    pop_col          = "poblacion",
    is_col           = "IS_distrito"
  )

  # Polarization and ENEP (national, population-weighted)
  nat_vs <- colSums(vote_shares_matrix * grid_adj$poblacion) / sum(grid_adj$poblacion)
  nat_vs <- nat_vs[seq_len(n_p)]
  nat_vs <- nat_vs / sum(nat_vs)
  mean_i <- sum(nat_vs * party_ideologies)

  data.frame(
    Simulacion        = i,
    scenario          = scenario,
    f                 = f_val,
    n_parties         = n_p,
    MAL               = nm$MAL,
    MRB               = nm$MRB,
    SRB               = nm$SRB,
    TRB               = nm$TRB,
    BIASSHARE_MAL     = nm$BIASSHARE_MAL,
    BIASSHARE_SEAT    = nm$BIASSHARE_SEAT,
    TH_USED_mean      = mean(TH_USED_vec, na.rm = TRUE),
    Threshold_Bias_mean = mean(ThBias_vec, na.rm = TRUE),
    polarization      = sum(nat_vs * abs(party_ideologies[seq_len(n_p)] - mean_i)),
    enep              = 1 / sum(nat_vs^2)
  )
}

# ---------------------------------------------------------------------------
# 5. Parallel setup
# ---------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 3)
cl      <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
cat("Cores:", n_cores, "\n")

# Export everything once — avoids the "already exporting" warning on 2nd foreach
parallel::clusterExport(cl, c(
  "grid_base", "geometria_base",
  "SCALE_FACTOR", "SEAT_DIVISOR", "DIRICHLET_ALPHA", "BETA_SHAPE", "TH",
  "distribute_seats_dhondt_multiparty", "calculate_national_metrics",
  "simulate_random_malapp", "simulate_ideological_malapp", "run_one_sim"
))

# ---------------------------------------------------------------------------
# 6. Run Pool A — random scenario
# ---------------------------------------------------------------------------
cat("--- Pool A: RANDOM scenario (N =", N_PER_SCENARIO, ", offset =", BATCH_OFFSET, ") ---\n")
set.seed(42 + BATCH_OFFSET)
f_random <- runif(N_PER_SCENARIO, 1, 15)
parallel::clusterExport(cl, "f_random")

results_random <- foreach(
  i         = seq_len(N_PER_SCENARIO),
  .packages = c("sf", "dplyr", "gtools")
) %dopar% {
  run_one_sim(i + BATCH_OFFSET, f_random[i], "random", grid_base, geometria_base)
}
pool_random <- do.call(rbind, results_random)
cat("Random done. Rows:", nrow(pool_random), "\n")

# ---------------------------------------------------------------------------
# 7. Run Pool A — concentrated scenario
# ---------------------------------------------------------------------------
cat("--- Pool A: CONCENTRATED scenario (N =", N_PER_SCENARIO, ", offset =", BATCH_OFFSET, ") ---\n")
set.seed(123 + BATCH_OFFSET)
f_ideol <- runif(N_PER_SCENARIO, 1, 15)
parallel::clusterExport(cl, "f_ideol")

results_ideol <- foreach(
  i         = seq_len(N_PER_SCENARIO),
  .packages = c("sf", "dplyr", "gtools")
) %dopar% {
  run_one_sim(i + N_PER_SCENARIO + BATCH_OFFSET, f_ideol[i], "concentrated", grid_base, geometria_base)
}
pool_ideol <- do.call(rbind, results_ideol)
cat("Concentrated done. Rows:", nrow(pool_ideol), "\n")

parallel::stopCluster(cl)

# ---------------------------------------------------------------------------
# 8. Combine and save
# ---------------------------------------------------------------------------
pool_new <- bind_rows(pool_random, pool_ideol) %>%
  mutate(scenario = factor(scenario, levels = c("random", "concentrated")))

out_path <- file.path(out_dir, "pool_mrb_summary.rds")

# Append to existing pool if it exists
if (file.exists(out_path)) {
  pool_existing <- readRDS(out_path)
  pool_mrb <- bind_rows(pool_existing, pool_new)
  cat("Appended", nrow(pool_new), "rows to existing pool.\n")
} else {
  pool_mrb <- pool_new
}

saveRDS(pool_mrb, out_path, compress = "xz")
cat("Saved:", out_path, "\n")
cat("Total rows:", nrow(pool_mrb), "| File size:", round(file.size(out_path) / 1e6, 2), "MB\n")
print(head(pool_mrb))
