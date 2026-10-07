#Set working directory to "coastalspruce" GitHub repo

################################################################################
################### NEEDLE SPECTRAL ANALYSIS ############################
################################################################################

######################### JOIN WC TO NP SPECTRA ##############################
# Define input file for np_spectra
infile <- "./data/branch_experiment/np_spectra.csv"
# Read CSV back into R
np_spectra <- read.csv(infile,
                       row.names = 1,          # keep the tree IDs as row names
                       check.names = FALSE,    # preserve original column names
                       stringsAsFactors = FALSE)
# Input file for np_experiment_data
infile2 <- "./data/branch_experiment/np_experiment_data.csv"
# Read CSV back into R
np_experiment <- read.csv(infile2,
                       check.names = FALSE,    # preserve original column names
                       stringsAsFactors = FALSE)

# Create key field in np_spectra
# Extract IDs from row names
ids <- rownames(np_spectra)

# Pull out tree number (digits after "Tree")
tree_num <- sub(".*Tree(\\d+).*", "\\1", ids)

# Pull out round number (digits after "R")
round_num <- sub(".*_R(\\d+)_.*", "\\1", ids)

# Combine into key field (tree_round)
np_spectra$KeyField <- paste0(tree_num, "_", round_num)

# Create key field in np_experiment
np_experiment$KeyField <- paste0(np_experiment$Tree, "_", np_experiment$Round)

# Quick check
head(np_experiment[, c("Tree","Round","KeyField")])

# Do the join!!
# Preserve np_spectra rownames through merge
np_spectra$._ID <- rownames(np_spectra)

np_spectra_joined <- merge(
  np_spectra,
  np_experiment,
  by = "KeyField",
  all.x = TRUE,      # left join: keep all spectra
  sort = FALSE
)

# Restore original row names and drop helper
rownames(np_spectra_joined) <- np_spectra_joined$._ID
np_spectra_joined$._ID <- NULL

# Quick diagnostics
cat("Rows in np_spectra:        ", nrow(np_spectra), "\n")
cat("Rows in joined (should match np_spectra): ", nrow(np_spectra_joined), "\n")
if (nrow(np_experiment) > 1) {
  first_meta_col <- setdiff(names(np_experiment), "KeyField")[1]
  cat("Unmatched spectra rows:     ",
      sum(is.na(np_spectra_joined[[ first_meta_col ]])), "\n")
}

####################### ADD VEGETATION INDICES TO NP SPECTRA JOINED ####################
#Use np_spectra_joined here

# 1) Identify wavelength columns and order numerically
is_wl_col      <- function(nm) startsWith(nm, "nm_")
wl_from_names  <- function(nms) as.numeric(gsub("_", ".", sub("^nm_", "", nms), fixed = TRUE))

wl_cols <- which(is_wl_col(names(np_spectra_joined)))
if (length(wl_cols) == 0) stop("No 'nm_' wavelength columns found in np_spectra_joined.")

wl_all <- wl_from_names(names(np_spectra_joined)[wl_cols])
o      <- order(wl_all)
wl     <- wl_all[o]                        # numeric wavelengths (integers in your case)
X      <- as.matrix(np_spectra_joined[, wl_cols[o], drop = FALSE])  # spectra matrix

# 2) Ensure required exact bands exist (integer wavelengths)
need <- c(531, 570, 550, 670, 680, 700, 703, 704, 705, 709, 720, 753, 754, 755, 800)
need_names <- paste0("nm_", need)
have_names <- names(np_spectra_joined)[wl_cols[o]]
missing <- setdiff(need_names, have_names)
if (length(missing) > 0) stop("Missing required bands: ", paste(missing, collapse = ", "))

# Helper to fetch columns by exact band name
col_ix <- function(nm) match(paste0("nm_", nm), have_names)

# Extract band columns (vectors over all rows)
R531 <- X[, col_ix(531)]
R570 <- X[, col_ix(570)]
R550 <- X[, col_ix(550)]
R670 <- X[, col_ix(670)]
R680 <- X[, col_ix(680)]
R700 <- X[, col_ix(700)]
R704 <- X[, col_ix(704)]
R709 <- X[, col_ix(709)]
R720 <- X[, col_ix(720)]
R754 <- X[, col_ix(754)]
R800 <- X[, col_ix(800)]
R703 <- X[, col_ix(703)]
R705 <- X[, col_ix(705)]
R753 <- X[, col_ix(753)]
R755 <- X[, col_ix(755)]
eps  <- .Machine$double.eps

# 3) Compute indices (vectorized)
PRI   <- (R531 - R570) / pmax(R531 + R570, eps)
NDVI  <- (R800 - R680) / pmax(R800 + R680, eps)
GNDVI  <- (R800 - R550) / pmax(R800 + R550, eps)
NDRE  <- (R800 - R720) / pmax(R800 + R720, eps)

TCARI <- 3 * ((R700 - R670) - 0.2 * (R700 - R550) * (R700 / pmax(R670, eps)))
OSAVI <- (1 + 0.16) * (R800 - R670) / pmax(R800 + R670 + 0.16, eps)
TCARIOSAVI <- TCARI / pmax(OSAVI, eps)

D1_754 <- (R755 - R753) / 2
D1_704 <- (R705 - R703) / 2
Datt3  <- D1_754 / pmax(D1_704, eps)

# 4) Bind back to np_spectra_joined (clean names, no suffixes)
vi_df <- data.frame(
  PRI        = PRI,
  NDVI       = NDVI,
  NDRE       = NDRE,
  TCARIOSAVI = TCARIOSAVI,
  Datt3      = Datt3,
  GNDVI      = GNDVI,
  check.names = FALSE
)

np_spectra_joined <- cbind(np_spectra_joined, vi_df)

# Define output path
outfile <- "./data/branch_experiment/np_spectra_VIs_manu.csv"

# Write CSV with row names preserved
write.csv(np_spectra_joined,
          file = outfile,
          row.names = TRUE)

################# CALCULATE MEDIAN SPECTRA/INDICES PER NP AND ROUND (AGGREGATE) ############
# Read back in np_spectra_joined (if needed)
# Define input path
infile <- "./data/branch_experiment/np_spectra_VIs_manu.csv"

# Read CSV back into R, keeping row names
np_spectra_joined <- read.csv(infile,
                              row.names = 1,
                              check.names = FALSE,
                              stringsAsFactors = FALSE)

library(dplyr)

# Columns to take first value from (non-spectral metadata)
meta_cols <- c("Tree", "Round", "Date", "Time", "Military_Time",
               "Temp_C", "RH", "Fresh_mass", "Dry_mass", "WC", "Lights", "Notes")

# Columns to take median of (spectral + VIs)
spectral_cols <- c(grep("^nm_", names(np_spectra_joined), value = TRUE),
                   c("PRI", "NDVI", "NDRE", "TCARIOSAVI", "Datt3", "GNDVI"))

# Aggregate
np_spectra_agg <- np_spectra_joined %>%
  group_by(KeyField) %>%
  summarise(
    across(all_of(meta_cols), first),
    across(all_of(spectral_cols), median, na.rm = TRUE),
    .groups = "drop"
  )

# Remove trees 6 and 7 (needles not from canopy red spruce branches)
np_spectra_agg <- np_spectra_agg %>%
  filter(!Tree %in% c(6, 7))

# Export
write.csv(np_spectra_agg,
          file = "./data/branch_experiment/np_spectra_agg.csv",
          row.names = FALSE)

cat("Aggregated dataframe dimensions:", nrow(np_spectra_agg), "rows x", ncol(np_spectra_agg), "cols\n")

######################## PLSR ANALYSIS FOR NP SPECTRA! ##############################
# Read back in np_spectra_agg (if needed)
# Define input path
infile <- "./data/branch_experiment/np_spectra_agg.csv"

# Read CSV back into R, keeping row names
np_spectra_agg <- read.csv(infile,
                              row.names = 1,
                              check.names = FALSE,
                              stringsAsFactors = FALSE)

library(pls)
library(ggplot2)

# ---- SETTINGS ----
SPEC_MIN   <- 398      # lower wavelength bound (nm). 398 / 1000 = drone-equivalent run
SPEC_MAX   <- 1000     # upper wavelength bound (nm). 350 / 2500 = full spectrometer run
NCOMP_MAX  <- 5        # max components (n = 30 samples)
SCALE      <- FALSE    # TRUE = unit-variance scaling of each wavelength; FALSE = mean-centering only
PREPROCESS <- "none"   # "none", "snv" (standard normal variate), or "d1" (first derivative)
RUN_PERMTEST <- FALSE  # TRUE to run the (within-pile) permutation test
N_PERMS    <- 999      # permutations (999 is a good start; 10,000 for publication)
SEED       <- 42

AXIS_LABEL  <- "WC"
TITLE_LABEL <- "Needle Water Content"

OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE)
run_tag <- paste0(SPEC_MIN, "-", SPEC_MAX, "_", PREPROCESS, "_",
                  ifelse(SCALE, "scaled", "centered"))

# ---- DATA PREP ----
spec_cols_all <- grep("^nm_[0-9]+$", names(np_spectra_agg), value = TRUE)
wl_all <- as.numeric(sub("nm_", "", spec_cols_all))
keep   <- wl_all >= SPEC_MIN & wl_all <= SPEC_MAX
if (sum(keep) < 10) stop("Fewer than 10 wavelengths fall inside SPEC_MIN-SPEC_MAX.")
ord       <- order(wl_all[keep])
spec_cols <- spec_cols_all[keep][ord]
wl        <- wl_all[keep][ord]

needed <- c(spec_cols, "WC", "Tree", "Round")
dat <- np_spectra_agg[complete.cases(np_spectra_agg[, needed]), ]

X_raw <- as.matrix(dat[, spec_cols])
y     <- dat$WC
g     <- factor(dat$Tree)
rnd   <- factor(dat$Round, levels = sort(unique(dat$Round)))

# Preprocessing is applied spectrum-by-spectrum (row-wise), so it cannot leak
# information between training and test samples.
preprocess_spectra <- function(X, wl, method) {
  if (method == "snv") {
    X <- t(apply(X, 1, function(r) (r - mean(r)) / sd(r)))
  } else if (method == "d1") {
    dw <- diff(wl)
    X  <- sweep(t(apply(X, 1, diff)), 2, dw, "/")  # n x (p-1) slope per nm
    wl <- wl[-1] - dw / 2                        # midpoint wavelengths
  } else if (method != "none") {
    stop("PREPROCESS must be 'none', 'snv', or 'd1'")
  }
  colnames(X) <- wl
  list(X = X, wl = wl)
}
pp    <- preprocess_spectra(X_raw, wl, PREPROCESS)
X_mat <- pp$X
wl    <- pp$wl

cat("Modelling WC | range:", min(wl), "-", max(wl), "nm |",
    ncol(X_mat), "predictors | preprocess:", PREPROCESS,
    "| scaled:", SCALE, "\n")
cat("n =", length(y), "samples in", nlevels(g), "piles\n")
print(table(Pile = g, Round = rnd))

# ---- HELPER FUNCTIONS ----
make_df <- function(X, y = NULL) {
  d <- data.frame(row = seq_len(nrow(X)))
  if (!is.null(y)) d$y <- y
  d$X <- X
  d
}

fit_pls <- function(X, y, ncomp) {
  plsr(y ~ X, ncomp = ncomp, data = make_df(X, y),
       scale = SCALE, validation = "none")
}

# Predictions for several component numbers at once -> n x length(ncomps) matrix
predict_pls <- function(model, Xnew, ncomps) {
  p <- predict(model, ncomp = ncomps, newdata = make_df(Xnew))
  matrix(p, nrow = nrow(Xnew), ncol = length(ncomps))
}

rmsep_by_comp <- function(pred, obs) sqrt(colMeans((pred - obs)^2))

# Leave-one-pile-out CV. Returns out-of-fold predictions for 1..K components
# plus a baseline (training-fold mean) used for Q2.
lopo_cv <- function(X, y, g, K) {
  g    <- droplevels(g)
  pred <- matrix(NA_real_, length(y), K)
  base <- rep(NA_real_, length(y))
  for (lv in levels(g)) {
    te <- g == lv
    tr <- !te
    m  <- fit_pls(X[tr, , drop = FALSE], y[tr], K)
    pred[te, ] <- predict_pls(m, X[te, , drop = FALSE], 1:K)
    base[te]   <- mean(y[tr])
  }
  list(pred = pred, base = base)
}

# Nested CV: for each held-out pile, choose ncomp by leave-one-pile-out within
# the remaining piles, refit, then predict the held-out pile.
nested_lopo <- function(X, y, g, K) {
  g    <- droplevels(g)
  pred <- rep(NA_real_, length(y))
  base <- rep(NA_real_, length(y))
  k_chosen <- setNames(rep(NA_integer_, nlevels(g)), levels(g))
  for (lv in levels(g)) {
    te <- g == lv
    tr <- !te
    inner  <- lopo_cv(X[tr, , drop = FALSE], y[tr], g[tr], K)
    k_best <- which.min(rmsep_by_comp(inner$pred, y[tr]))
    m <- fit_pls(X[tr, , drop = FALSE], y[tr], k_best)
    pred[te] <- predict_pls(m, X[te, , drop = FALSE], k_best)[, 1]
    base[te] <- mean(y[tr])
    k_chosen[lv] <- k_best
  }
  list(pred = pred, base = base, k_chosen = k_chosen)
}

calc_metrics <- function(obs, pred, base) {
  rmsep <- sqrt(mean((obs - pred)^2))
  c(RMSEP = rmsep,
    R2    = if (sd(pred) > 0 && sd(obs) > 0) cor(obs, pred)^2 else NA_real_,
    Q2    = 1 - sum((obs - pred)^2) / sum((obs - base)^2),
    Bias  = mean(pred - obs),
    RPD   = sd(obs) / rmsep)
}

# VIP scores (works for any ncomp, including 1)
vip_func <- function(model, ncomp) {
  W  <- model$loading.weights[, 1:ncomp, drop = FALSE]
  Q  <- as.vector(model$Yloadings)[1:ncomp]
  T  <- model$scores[, 1:ncomp, drop = FALSE]
  SS <- Q^2 * colSums(T^2)
  W_norm <- sweep(W, 2, sqrt(colSums(W^2)), "/")
  as.vector(sqrt(nrow(W) * (W_norm^2 %*% SS) / sum(SS)))
}

# ---- 1. LEAVE-ONE-PILE-OUT CV (component-number diagnostic) ----
cv_simple  <- lopo_cv(X_mat, y, g, NCOMP_MAX)
rmsep_curve <- rmsep_by_comp(cv_simple$pred, y)
n_opt <- which.min(rmsep_curve)
cat("\nRMSEP by number of components (leave-one-pile-out):\n")
print(round(setNames(rmsep_curve, 1:NCOMP_MAX), 4))
cat("Optimal ncomp (lowest RMSEP on all 30 samples):", n_opt, "\n")

# ---- 2. NESTED CV (headline performance) ----
# NOTE: choosing ncomp and reporting its error on the same CV is slightly
# optimistic, so the headline metrics below come from nested CV.
nested <- nested_lopo(X_mat, y, g, NCOMP_MAX)
m_nested <- calc_metrics(y, nested$pred, nested$base)
cat("\nComponents chosen in each outer fold:\n"); print(nested$k_chosen)
cat("\n--- Nested leave-one-pile-out metrics (headline) ---\n")
print(round(m_nested, 4))

m_simple <- calc_metrics(y, cv_simple$pred[, n_opt], cv_simple$base)
cat("\n--- Non-nested LOPO metrics at ncomp =", n_opt, "(diagnostic) ---\n")
print(round(m_simple, 4))
cat("(If R2 and Q2 differ a lot, the model has offset/slope errors.)\n")

# Per-pile metrics (descriptive only: 4-7 points each)
pile_metrics <- do.call(rbind, lapply(levels(g), function(lv) {
  i <- g == lv
  data.frame(Pile = lv, n = sum(i),
             RMSEP = sqrt(mean((y[i] - nested$pred[i])^2)),
             R2 = if (sd(nested$pred[i]) > 0) cor(y[i], nested$pred[i])^2 else NA,
             Bias = mean(nested$pred[i] - y[i]))
}))
cat("\n--- Per-pile metrics (nested predictions; descriptive) ---\n")
print(pile_metrics, digits = 3, row.names = FALSE)

# ---- 3. FINAL MODEL ON ALL DATA + VIP ----
pls_final  <- fit_pls(X_mat, y, n_opt)
vip_scores <- vip_func(pls_final, n_opt)
stopifnot(length(vip_scores) == length(wl))
vip_df <- data.frame(Wavelength = wl, VIP = vip_scores)

# ---- 4. OPTIONAL PERMUTATION TEST (within-pile, full nested procedure) ----
p_val <- NA
if (RUN_PERMTEST) {
  cat("\nRunning permutation test (", N_PERMS, "permutations)...\n")
  set.seed(SEED)
  permute_within <- function(y, g) {
    yp <- y
    for (lv in levels(g)) {
      i <- which(g == lv)
      yp[i] <- y[i[sample.int(length(i))]]
    }
    yp
  }
  perm_q2 <- numeric(N_PERMS)
  pb <- txtProgressBar(min = 0, max = N_PERMS, style = 3)
  for (i in seq_len(N_PERMS)) {
    yp <- permute_within(y, g)
    nn <- nested_lopo(X_mat, yp, g, NCOMP_MAX)
    perm_q2[i] <- calc_metrics(yp, nn$pred, nn$base)["Q2"]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  p_val <- (1 + sum(perm_q2 >= m_nested["Q2"])) / (N_PERMS + 1)
  cat("\nPermutation p-value (statistic = nested Q2):", p_val, "\n")
}

# ---- 5. PLOTS ----
plot_df <- data.frame(Observed = y, Predicted = nested$pred,
                      Residual = nested$pred - y, Pile = g, Round = rnd)
shape_vals <- rep_len(c(16, 17, 15, 18, 3, 4, 8), nlevels(rnd))
lims <- range(c(plot_df$Observed, plot_df$Predicted))

p_obs <- ggplot(plot_df, aes(Observed, Predicted, color = Pile)) +
  geom_abline(slope = 1, intercept = 0, color = "red",
              linetype = "dashed", linewidth = 0.9) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2") +
  coord_equal(xlim = lims, ylim = lims) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, size = 3.8,
           label = paste0("Q² = ", round(m_nested["Q2"], 3),
                          "\nR² = ", round(m_nested["R2"], 3),
                          "\nRMSEP = ", round(m_nested["RMSEP"], 3),
                          "\nRPD = ", round(m_nested["RPD"], 3))) +
  labs(title = paste0(TITLE_LABEL, " (", SPEC_MIN, "-", SPEC_MAX, " nm)",
                      "\nObserved vs Predicted (nested leave-one-pile-out CV)"),
       x = paste("Observed", AXIS_LABEL), y = paste("Predicted", AXIS_LABEL)) +
  theme_bw(base_size = 14)
print(p_obs)

p_resid <- ggplot(plot_df, aes(Observed, Residual, color = Pile,)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = shape_vals) +
  labs(title = paste("Residuals (predicted - observed) -", TITLE_LABEL),
       x = paste("Observed", AXIS_LABEL), y = "Residual") +
  theme_bw(base_size = 14)
print(p_resid)

p_rmsep <- ggplot(data.frame(Components = 1:NCOMP_MAX, RMSEP = rmsep_curve),
                  aes(Components, RMSEP)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  geom_point(data = data.frame(Components = n_opt, RMSEP = rmsep_curve[n_opt]),
             color = "red", size = 4, shape = 1, stroke = 1.5) +
  labs(title = paste("RMSEP vs Components (leave-one-pile-out) -", TITLE_LABEL),
       x = "Number of components", y = "RMSEP") +
  theme_bw(base_size = 14)
print(p_rmsep)

p_vip <- ggplot(vip_df, aes(Wavelength, VIP)) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 0.9) +
  labs(title = paste("VIP Scores -", TITLE_LABEL,
                     paste0("(", n_opt, " components)")),
       x = ifelse(PREPROCESS == "d1", "Wavelength (nm, midpoints)", "Wavelength (nm)"),
       y = "VIP Score") +
  theme_bw(base_size = 14)
print(p_vip)

if (RUN_PERMTEST) {
  p_perm <- ggplot(data.frame(perm_q2 = perm_q2), aes(perm_q2)) +
    geom_histogram(bins = 30, fill = "lightgrey", color = "white") +
    geom_vline(xintercept = m_nested["Q2"], color = "red", linewidth = 0.9) +
    annotate("text", x = m_nested["Q2"], y = Inf, hjust = 1.1, vjust = 2,
             size = 3.8, color = "red",
             label = paste0("Observed Q² = ", round(m_nested["Q2"], 3),
                            "\np = ", signif(p_val, 3))) +
    labs(title = paste("Within-pile Permutation Test -", TITLE_LABEL),
         x = "Permuted nested Q²", y = "Frequency") +
    theme_bw(base_size = 14)
  print(p_perm)
}

# ---- 6. RESULTS TABLE ----
results_row <- data.frame(
  Range_nm        = paste0(SPEC_MIN, "-", SPEC_MAX),
  Preprocess      = PREPROCESS,
  Scaled          = SCALE,
  n_predictors    = length(wl),
  n_samples       = length(y),
  nComp_optimal   = n_opt,
  nComp_nested_median = median(nested$k_chosen),
  RMSEP_nested    = m_nested["RMSEP"],
  R2_nested       = m_nested["R2"],
  Q2_nested       = m_nested["Q2"],
  Bias_nested     = m_nested["Bias"],
  RPD_nested      = m_nested["RPD"],
  RMSEP_LOPO      = m_simple["RMSEP"],
  Q2_LOPO         = m_simple["Q2"],
  Perm_p          = p_val,
  row.names = NULL
)
print(results_row)

# ---- 7. EXPORT ----
f <- function(stem, ext) file.path(OUT_DIR, paste0("np_PLSR_", stem, "_", run_tag, ".", ext))

write.csv(results_row, f("results", "csv"), row.names = FALSE)
write.csv(pile_metrics, f("pile_metrics", "csv"), row.names = FALSE)
write.csv(data.frame(Pile = g, Round = rnd, Observed = y,
                     Predicted_nested = nested$pred),
          f("obs_vs_pred", "csv"), row.names = FALSE)
write.csv(data.frame(Components = 1:NCOMP_MAX, RMSEP = rmsep_curve),
          f("RMSEP_by_comp", "csv"), row.names = FALSE)
write.csv(vip_df, f("VIP", "csv"), row.names = FALSE)

ggsave(f("obs_vs_pred", "png"), p_obs,   width = 7,  height = 7, dpi = 300)
ggsave(f("residuals", "png"),   p_resid, width = 7,  height = 5, dpi = 300)
ggsave(f("RMSEP_by_comp", "png"), p_rmsep, width = 6, height = 5, dpi = 300)
ggsave(f("VIP", "png"),         p_vip,   width = 10, height = 5, dpi = 300)
if (RUN_PERMTEST) ggsave(f("permtest", "png"), p_perm, width = 7, height = 6, dpi = 300)


####################### PLOTTING INDEX vs WC IN NP (ALL SAMPLES) ###########################
# Read back in np_spectra_joined (if needed)
# Define input path
infile <- "./data/branch_experiment/np_spectra_VIs_manu.csv"

# Read CSV back into R, keeping row names
np_spectra_joined <- read.csv(infile,
                              row.names = 1,
                              check.names = FALSE,
                              stringsAsFactors = FALSE)

# USER SETTINGS 
INDEX        <- "PRI"          # "PRI","NDVI","NDRE","TCARIOSAVI","Datt3","CARI","Boochs"
X_AXIS       <- "WC"           # "WC" or "INDEX"  (the other will be Y)
STAT         <- "MEDIAN"       # "MEDIAN" or "MEAN" for within Tree×Round aggregation
SHOW_FIT     <- FALSE           # draw a single overall linear fit?
POINT_LABELS <- FALSE          # label each point with the round number?
#

stopifnot(exists("np_spectra_joined"))
need <- c("KeyField", "WC", INDEX)
if (!all(need %in% names(np_spectra_joined))) {
  stop("np_spectra_joined must contain columns: ", paste(need, collapse = ", "))
}

# --------- Dynamic plot labels based on INDEX and X_AXIS ----------
X_AXIS <- toupper(X_AXIS)
if (!X_AXIS %in% c("WC","INDEX")) stop("X_AXIS must be 'WC' or 'INDEX'")

if (X_AXIS == "WC") {
  XLAB <- "Water Content"
  YLAB <- INDEX
  PLOT_TITLE <- sprintf("%s vs Water Content (Tree × Round medians)", INDEX)
} else {
  XLAB <- INDEX
  YLAB <- "Water Content"
  PLOT_TITLE <- sprintf("Water Content vs %s (Tree × Round medians)", INDEX)
}

# -------------------- Parse Tree/Round ---------------------
kf <- as.character(np_spectra_joined$KeyField)
parts <- do.call(rbind, strsplit(kf, "_", fixed = TRUE))
TreeNum  <- as.integer(parts[, 1])
RoundNum <- as.integer(parts[, 2])

# Build working frame with generic 'Index' column
df <- data.frame(
  KeyField = kf,
  TreeNum  = TreeNum,
  RoundNum = RoundNum,
  Index    = np_spectra_joined[[INDEX]],
  WC       = np_spectra_joined$WC
)

# Filter valid rows
df <- df[is.finite(df$Index) & is.finite(df$WC) & !is.na(df$TreeNum) & !is.na(df$RoundNum), , drop = FALSE]
if (nrow(df) == 0) stop("No valid rows after filtering.")

# Aggregate to one row per Tree×Round (keeps name 'agg' for your LME code)
agg_fun <- switch(toupper(STAT),
                  "MEDIAN" = function(x) median(x, na.rm = TRUE),
                  "MEAN"   = function(x) mean(x,   na.rm = TRUE),
                  stop("STAT must be 'MEDIAN' or 'MEAN'"))
agg <- aggregate(cbind(Index, WC) ~ TreeNum + RoundNum + KeyField, data = df, FUN = agg_fun)

# Colors for trees
trees <- sort(unique(agg$TreeNum))
cols  <- setNames(rainbow(length(trees)), trees)

# Choose axes
x <- if (X_AXIS == "WC") agg$WC else agg$Index
y <- if (X_AXIS == "WC") agg$Index else agg$WC

# Plot with legend outside right
op <- par(mar = c(5, 4, 4, 10), xpd = NA); on.exit(par(op), add = TRUE)

plot(x, y, pch = 19,
     col = cols[as.character(agg$TreeNum)],
     xlab = XLAB, ylab = YLAB, main = PLOT_TITLE)

if (POINT_LABELS) {
  text(x, y, labels = agg$RoundNum, pos = 3, cex = 0.8)
}

if (SHOW_FIT && nrow(agg) >= 2 && all(is.finite(x)) && all(is.finite(y))) {
  fit <- lm(y ~ x)
  abline(fit, lwd = 2)
  r2 <- summary(fit)$r.squared
  usr <- par("usr")
  text(x = usr[1] + 0.02 * diff(usr[1:2]),
       y = usr[4] - 0.05 * diff(usr[3:4]),
       labels = paste0("R² = ", sprintf("%.2f", r2)),
       adj = c(0, 1))
}

legend("topright",
       inset = c(-0.15, 0),
       legend = trees, title = "Spruce",
       col = cols[as.character(trees)],
       pch = 19, bty = "n", cex = 0.9)


####################### LINEAR MIXED EFFECTS MODEL #############################
#install.packages("lme4")
library(lme4)
#install.packages("car")
library(car)

# Feed results from above OR if necessary read back in:
agg <- read.csv("./data/branch_experiment/np_spectra_agg.csv",
                check.names = FALSE,
                stringsAsFactors = FALSE)

#Index <- "PRI"

# Normality tests: QQ and Shapiro-Wilk
par(mfrow = c(1, 2))  # 1 row, 2 columns

# Index
qqnorm(agg$PRI, main = "QQ Plot of Index")
qqline(agg$PRI, col = "red", lwd = 2)

# WC
qqnorm(agg$WC, main = "QQ Plot of WC")
qqline(agg$WC, col = "red", lwd = 2)

par(mfrow = c(1, 1))  # reset layout

# Shapiro-Wilk tests for normality
shapiro_PRI <- shapiro.test(agg$PRI)
shapiro_WC  <- shapiro.test(agg$WC)

shapiro_PRI
shapiro_WC

# Fit a linear mixed-effects model:
#   Response: Index - PRI, NDVI, TCARI/OSAVI, etc.
#   Fixed effect: WC (common slope across trees)
#   Random effect: random intercept for each TreeNum (tree-specific baseline Index)
agg_model = lmer(Index~WC+(1|Tree),data=agg)
summary(agg_model)
anova(agg_model)
Anova(agg_model) # car version of the anova

# Null model with random intercepts only (no WC effect)
null_model = lmer(Index~1+(1|TreeNum),data=agg)
# Likelihood-ratio test comparing models (refitted with ML):
# Tests whether adding WC significantly improves model fit
anova(agg_model, null_model)

# Coniditional (tree level) predictions
# Includes each tree’s random intercept (same slope, tree-specific intercepts)
agg$pred = predict(agg_model)

plot(agg$WC, agg$pred, pch = 19,
     col = cols[as.character(agg$TreeNum)],
     xlab = XLAB, ylab = YLAB, main = "Predicted CARI vs Water Content by Tree (Random Intercept LMM)")
legend("topright",
       inset = c(-0.15, 0),      # push outside; adjust if needed
       legend = trees,           # just the numbers
       title  = "Spruce",
       col = cols[as.character(trees)],
       pch = 19, bty = "n", cex = 0.9)

agg$pred_overall = predict(agg_model, re.form = NA)

# Compute both marginal and conditional R²
#install.packages("MuMIn")
#library(MuMIn)
r.squaredGLMM(agg_model)

# Population-level predictions:
# Uses fixed effects only (no random intercepts); single overall line
plot(agg$WC, agg$pred_overall, pch = 19,
     col = cols[as.character(agg$TreeNum)],
     xlab = XLAB, ylab = YLAB, main = "Predicted Overall PRI vs Water Content - Fixed Intercept")
legend("topright",
       inset = c(-0.15, 0),      # push outside; adjust if needed
       legend = trees,           # just the numbers
       title  = "Spruce",
       col = cols[as.character(trees)],
       pch = 19, bty = "n", cex = 0.9)

# Compute both marginal and conditional R²
install.packages("MuMIn")
library(MuMIn)
r.squaredGLMM(agg_model)

# Write to .csv
outfile <- "G:/Branch_Experiment/r_outputs/np_agg_medians.csv"
write.csv(agg, file = outfile, row.names = FALSE)


################# PLOT SPECTRAL PROFILES FROM ONE TREE OVER TIME #######################
setwd("G:/Branch_Experiment/sprucedry")
np_spectra <- read.csv("G:/Branch_Experiment/r_outputs/np_spectra.csv", row.names = 1)

# USER SETTINGS
TREE       <- "Tree1"                    # e.g., "Tree1", "Tree2", ...
ROUNDS     <- c("R1","R2","R3","R5")     # c(1,3,5) or c("R1","R3","R5"); NULL/empty = all rounds
WL_RANGE   <- c(400, 2500)               # wavelength window in nm; set to NULL for full range
PLOT_TITLE <- "Needle Pile Spectra Across Rounds - Tree 1, Median"  # "" to auto-generate a title
STAT       <- "MEDIAN"                     # "MEDIAN", "MEAN" "MIN", or "MAX"

stopifnot(exists("np_spectra"))

## ----------- Select rows for the chosen tree & parse rounds -----------
all_rows <- rownames(np_spectra)
pat <- paste0("^", TREE, "_NeedlePile_R\\d+_")
ix_tree <- grepl(pat, all_rows, ignore.case = FALSE)
if (!any(ix_tree)) stop("No rows matched pattern: ", pat)

tree_rows <- all_rows[ix_tree]

# Extract the round number after "_R"
get_round <- function(s) as.integer(sub(".*_R(\\d+)_.*", "\\1", s))
round_id <- vapply(tree_rows, get_round, integer(1), USE.NAMES = FALSE)

## ----------- Apply optional ROUNDS filter -----------
normalize_rounds <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  if (is.character(x)) as.integer(gsub("^[Rr]", "", x)) else as.integer(x)
}
rounds_requested <- normalize_rounds(ROUNDS)

available_rounds <- sort(unique(round_id))
if (!is.null(rounds_requested)) {
  missing <- setdiff(rounds_requested, available_rounds)
  if (length(missing) > 0) {
    warning("Requested rounds not found for ", TREE, ": R", paste(missing, collapse = ", R"))
  }
  keep <- round_id %in% rounds_requested
  if (!any(keep)) stop("No rows for requested rounds. Available for ", TREE, ": R", paste(available_rounds, collapse = ", R"))
  tree_rows <- tree_rows[keep]
  round_id  <- round_id[keep]
}
## ----------- Prepare wavelength axis (numeric) -----------
wl_chr <- sub("^nm_", "", colnames(np_spectra))
wl_num <- as.numeric(gsub("_", ".", wl_chr, fixed = TRUE))
o <- order(wl_num)
wl_num <- wl_num[o]

## ----------- Choose aggregation function -----------
STAT <- toupper(trimws(STAT))
agg_fun <- switch(
  STAT,
  "MEDIAN" = function(x) median(x, na.rm = TRUE),
  "MEAN"   = function(x) mean(x,   na.rm = TRUE),
  "MIN"    = function(x) min(x,    na.rm = TRUE),
  "MAX"    = function(x) max(x,    na.rm = TRUE),
  stop("STAT must be one of: 'MEDIAN', 'MEAN', 'MIN', 'MAX'")
)
## ----------- Compute per-round column-wise aggregates -----------
rounds_to_compute <- sort(unique(round_id))
agg_list <- lapply(rounds_to_compute, function(r) {
  rows_r <- tree_rows[round_id == r]
  # aggregate across all replicates in the round, per wavelength column
  apply(np_spectra[rows_r, , drop = FALSE], 2, agg_fun)
})
agg_mat <- do.call(rbind, agg_list)[, o, drop = FALSE]
rownames(agg_mat) <- paste0("R", rounds_to_compute)
colnames(agg_mat) <- paste0("nm_", wl_num)

## ----------- Apply wavelength window (if requested) -----------
if (!is.null(WL_RANGE)) {
  WL_RANGE <- sort(as.numeric(WL_RANGE))
  idx_wl <- wl_num >= WL_RANGE[1] & wl_num <= WL_RANGE[2]
  if (!any(idx_wl)) stop("No wavelengths within requested range: ", paste(WL_RANGE, collapse = "-"), " nm")
  wl_num  <- wl_num[idx_wl]
  agg_mat <- agg_mat[, idx_wl, drop = FALSE]
}

## ----------- Plot -----------
cols <- grDevices::rainbow(nrow(agg_mat))
ymin <- min(agg_mat, na.rm = TRUE)
ymax <- max(agg_mat, na.rm = TRUE)

auto_title <- paste(TREE, sprintf("round-%s spectra", tolower(STAT)))
main_title <- if (nzchar(PLOT_TITLE)) PLOT_TITLE else auto_title

plot(wl_num, agg_mat[1, ], type = "l", lwd = 2, col = cols[1],
     xlab = "Wavelength (nm)",
     ylab = "Reflectance (fraction)",
     main = main_title,
     ylim = c(ymin, ymax))

if (nrow(agg_mat) > 1) {
  for (i in 2:nrow(agg_mat)) lines(wl_num, agg_mat[i, ], lwd = 2, col = cols[i])
}
legend("topright", inset = c(-0.15, 0),legend = rownames(agg_mat), lwd = 2, col = cols, bty = "n")
