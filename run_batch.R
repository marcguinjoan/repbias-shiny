# =============================================================================
# run_batch.R — Run one batch of Pool A (MRB) + Pool B (SRB/ThRB) in sequence
#
# HOW TO USE:
#   Batch 1 (first run):  BATCH_OFFSET <- 0      in both scripts → run this file
#   Batch 2 (next run):   set BATCH_OFFSET <- 50000 in both scripts → run again
#   Batch 3:              BATCH_OFFSET <- 150000, etc.
#
# Each pool saves independently: if Pool B fails, Pool A is already on disk.
# =============================================================================

cat("\n========================================\n")
cat("BATCH START:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# Detect script directory (works in RStudio and via Rscript)
script_dir <- tryCatch({
  dirname(rstudioapi::getSourceEditorContext()$path)
}, error = function(e) {
  args      <- commandArgs(trailingOnly = FALSE)
  file_flag <- args[grepl("--file=", args)]
  if (length(file_flag) > 0) {
    dirname(normalizePath(sub("--file=", "", file_flag)))
  } else {
    getwd()
  }
})

# ---------------------------------------------------------------------------
# POOL A — MRB (50k total: 25k random + 25k concentrated)
# ---------------------------------------------------------------------------
cat(">>> POOL A: MRB simulation starting...\n\n")
t0 <- proc.time()
source(file.path(script_dir, "pool_mrb_generation.R"), local = FALSE)
elapsed_a <- round((proc.time() - t0)["elapsed"] / 60, 1)
cat("\n>>> POOL A done in", elapsed_a, "min\n\n")

# ---------------------------------------------------------------------------
# POOL B — SRB / ThRB (50k total: 4500 × 11 DM values)
# ---------------------------------------------------------------------------
cat(">>> POOL B: SRB/ThRB simulation starting...\n\n")
t1 <- proc.time()
source(file.path(script_dir, "pool_srb_generation.R"), local = FALSE)
elapsed_b <- round((proc.time() - t1)["elapsed"] / 60, 1)
cat("\n>>> POOL B done in", elapsed_b, "min\n\n")

cat("========================================\n")
cat("BATCH COMPLETE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total time:", round(elapsed_a + elapsed_b, 1), "min\n")
cat("========================================\n")
