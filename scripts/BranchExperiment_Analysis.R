#setwd("G:/Branch_Experiment/sprucedry")

################################################################################
################### BRANCH/NEEDLE SPECTRAL ANALYSIS ############################
################################################################################

################### CREATE DATAFRAME OF *NEEDLE* REFLECTANCE #####################
## ---------- helper: read one SED reflectance
read_sed_reflect <- function(f) {
  lns <- readLines(f, warn = FALSE)
  i <- grep("^Data:\\s*$", lns)
  if (length(i) != 1) stop("Could not find 'Data:' section in ", f)
  hdr <- strsplit(lns[i + 1], "\t")[[1]]
  dat <- read.table(text = paste(lns[(i + 2):length(lns)], collapse = "\n"),
                    sep = "\t", header = FALSE, col.names = hdr, check.names = FALSE)
  wl   <- as.numeric(dat$Wvl)
  refl <- as.numeric(dat[["Reflect. %"]]) / 100  # 0–1
  data.frame(wl = wl, refl = refl)
}

# list files
sed_files <- list.files(pattern = "^Tree.*\\.sed$", ignore.case = TRUE)
if (length(sed_files) == 0) stop("No Tree*.sed files found in working directory.")

# read all
dfs <- lapply(sed_files, read_sed_reflect)
wls <- lapply(dfs, function(d) d$wl)

## FAST PATH: all wavelength vectors identical?
same_grid <- all(vapply(wls, function(w) identical(w, wls[[1]]), logical(1)))

if (same_grid) {
  # stack directly (no resampling)
  wl_vec <- wls[[1]]
  mat <- do.call(rbind, lapply(dfs, function(d) d$refl))
  rownames(mat) <- tools::file_path_sans_ext(basename(sed_files))
  # clean column names like nm_350, nm_351, ...
  colnames(mat) <- paste0("nm_", format(wl_vec, trim = TRUE))
  np_spectra <- as.data.frame(mat, stringsAsFactors = FALSE)
} else {
  # SAFE PATH: use exact intersection of wavelengths (no interpolation)
  common_wl <- Reduce(intersect, wls)
  if (length(common_wl) < 50) {
    warning("Very few common wavelengths across files; check instrument settings.")
  }
  # Keep common wavelengths in original (sorted) order from the first file
  common_wl <- sort(common_wl)
  mat <- do.call(rbind, lapply(dfs, function(d) {
    idx <- match(common_wl, d$wl)
    d$refl[idx]
  }))
  rownames(mat) <- tools::file_path_sans_ext(basename(sed_files))
  colnames(mat) <- paste0("nm_", gsub("\\.", "_", format(common_wl, trim = TRUE)))
  np_spectra <- as.data.frame(mat, stringsAsFactors = FALSE)
  
  # Optional: tell you how many wavelengths were dropped (purely informational)
  dropped <- length(wls[[1]]) - length(common_wl)
  message("Used exact intersection of wavelengths (no interpolation). Dropped ",
          dropped, " bands from the first file's grid.")
}

## quick sanity checks
dim(np_spectra)        # rows = files, cols = wavelengths kept
np_spectra[1, 1:10]    # first row, first few bands
summary(rowMeans(np_spectra, na.rm = TRUE))

# Define output folder
outdir <- "G:/Branch_Experiment/r_outputs"

# Write CSV (row names = your Tree file IDs)
outfile <- file.path(outdir, "np_spectra.csv")
write.csv(np_spectra, file = outfile, row.names = TRUE)

#################### OUTLIERS IN NP SPECTRA #############################
# Visual envelope, does anything poke out?
#--- helpers ---
is_wl_col <- function(nm) startsWith(nm, "nm_")
wl_num_of <- function(nms) as.numeric(gsub("_",".", sub("^nm_","", nms), fixed=TRUE))

tree <- "Tree2"             # choose
round_pat <- "^R(\\d+)$"    # we’ll compute all rounds present

rows_tree <- grepl(paste0("^", tree, "_NeedlePile_R\\d+_"), rownames(np_spectra))
X <- np_spectra[rows_tree, is_wl_col(colnames(np_spectra)), drop=FALSE]
wl <- wl_num_of(colnames(X))

# per round, draw envelope
get_round <- function(s) as.integer(sub(".*_R(\\d+)_.*", "\\1", s))
rids <- get_round(rownames(X))
rounds <- sort(unique(rids))

op <- par(mfrow=c(ceiling(length(rounds)/2), 2), mar=c(4,4,2,1))
for (r in rounds) {
  M <- as.matrix(X[rids==r, , drop=FALSE])
  med <- apply(M, 2, median, na.rm=TRUE)
  q05 <- apply(M, 2, quantile, probs=0.05, na.rm=TRUE)
  q95 <- apply(M, 2, quantile, probs=0.95, na.rm=TRUE)
  
  plot(wl, med, type="l", lwd=2, xlab="nm", ylab="Reflectance",
       main=paste(tree, "R", r, sep=""))
  polygon(c(wl, rev(wl)), c(q05, rev(q95)),
          border=NA, col=adjustcolor("gray70", alpha.f=0.5))
  lines(wl, med, lwd=2)
  # overlay each spectrum (thin) to spot outliers
  apply(M, 1, function(v) lines(wl, v, lwd=0.7, col=adjustcolor("black", 0.5)))
}
par(op)

# Robust per-band z scores
M <- as.matrix(X)  # all replicates for the chosen tree (all rounds)
band_median <- apply(M, 2, median, na.rm=TRUE)
band_mad    <- apply(M, 2, mad, constant=1.4826, na.rm=TRUE)  # ~sd under normality
robust_z    <- sweep(sweep(M, 2, band_median, "-"), 2, pmax(band_mad, .Machine$double.eps), "/")

# row-level rule: % of bands with |z| > 3
prop_high <- apply(abs(robust_z) > 3, 1, mean, na.rm=TRUE)
flag_rows <- names(prop_high[prop_high > 0.1])   # >10% bands extreme; tune threshold
flag_rows

if (length(flag_rows)) {
  matplot(wl, t(M[flag_rows, , drop=FALSE]), type="l", lty=1, xlab="nm", ylab="Reflectance",
          main="Flagged spectra", col=adjustcolor("red", 0.6))
  lines(wl, band_median, lwd=2)
}

# PCA Analysis
# keep a stable VNIR core (e.g., 450–900 nm) to avoid edge noise
keep <- wl >= 450 & wl <= 900
Msub <- M[, keep, drop=FALSE]
Msub <- Msub[, colSums(is.na(Msub)) < nrow(Msub)]  # drop all-NA bands

pca <- prcomp(Msub, center=TRUE, scale.=FALSE)
scores <- pca$x[, 1:5, drop=FALSE]   # first few PCs capture shape/level
center <- colMeans(scores)
covS   <- cov(scores)
d2     <- mahalanobis(scores, center, covS)        # Hotelling T²-ish
cut    <- quantile(d2, 0.99, na.rm=TRUE)           # 99th percentile
pca_outliers <- rownames(Msub)[d2 > cut]
pca_outliers

# Spike detection with residuals
# Simple Savitzky–Golay style using stats::filter as a quick proxy
k <- 11  # odd window
w <- rep(1/k, k)
smoothM <- t(apply(M, 1, function(v) stats::filter(v, w, sides=2)))
resid_sd <- apply((M - smoothM)^2, 1, function(x) sqrt(mean(x, na.rm=TRUE)))
spiky <- names(resid_sd[resid_sd > quantile(resid_sd, 0.99, na.rm=TRUE)])
spiky

# Test with VIs
wl_idx <- function(target) which.min(abs(wl - target))
i531 <- wl_idx(531); i570 <- wl_idx(570)
i680 <- wl_idx(680); i800 <- wl_idx(800)

PRI   <- (M[, i531] - M[, i570]) / (M[, i531] + M[, i570])
NDVI  <- (M[, i800] - M[, i680]) / (M[, i800] + M[, i680])

df_vi <- data.frame(
  id = rownames(M),
  round = paste0("R", get_round(rownames(M))),
  PRI = PRI, NDVI = NDVI
)
boxplot(PRI ~ round, df_vi, main=paste(tree, "PRI by round"))
boxplot(NDVI ~ round, df_vi, main=paste(tree, "NDVI by round"))

######################### JOIN WC TO NP SPECTRA ##############################
# Define input file for np_spectra
infile <- "G:/Branch_Experiment/r_outputs/np_spectra.csv"
# Read CSV back into R
np_spectra <- read.csv(infile,
                       row.names = 1,          # keep the tree IDs as row names
                       check.names = FALSE,    # preserve original column names
                       stringsAsFactors = FALSE)
# Input file for np_experiment_data
infile2 <- "G:/Branch_Experiment/np_experiment_data.csv"
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
## Activate hsdar if available (optional Boochs)
hsdar_ok <- requireNamespace("hsdar", quietly = TRUE)
if (hsdar_ok) {
  library(hsdar)  # will attach; skip if not installed
} else {
  message("Package 'hsdar' not installed; Boochs will be set to NA.")
}

stopifnot(exists("np_spectra_joined"))

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
need <- c(531, 570, 550, 670, 680, 700, 704, 709, 720, 754, 800)
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
eps  <- .Machine$double.eps

# 3) Compute indices (vectorized)
PRI   <- (R531 - R570) / pmax(R531 + R570, eps)
NDVI  <- (R800 - R680) / pmax(R800 + R680, eps)
NDRE  <- (R800 - R720) / pmax(R800 + R720, eps)

TCARI <- 3 * ((R700 - R670) - 0.2 * (R700 - R550) * (R700 / pmax(R670, eps)))
OSAVI <- (1 + 0.16) * (R800 - R670) / pmax(R800 + R670 + 0.16, eps)
TCARIOSAVI <- TCARI / pmax(OSAVI, eps)

Datt3 <- (R754 - R704) / pmax(R754 + R704, eps)

# Practical CARI variant (baseline in 550–700 region)
CARI <- abs((R700 - R670 - 0.2 * (R700 - R550)) * (R670 / pmax(R700, eps)))

# 4) Boochs via hsdar (optional)
Boochs <- rep(NA_real_, nrow(X))
if (hsdar_ok) {
  # Build speclib once (rows = samples, columns = wavelengths)
  colnames(X) <- paste0("nm_", wl)  # ensure names align with wl vector
  sl <- hsdar::speclib(X, wavelength = wl)
  # Try Boochs, then Boochs2, then Boochs1
  for (idx_name in c("Boochs", "Boochs2", "Boochs1")) {
    res <- try(hsdar::vegindex(sl, index = idx_name), silent = TRUE)
    if (!inherits(res, "try-error")) {
      # robust numeric extraction
      if (isS4(res) && "vi" %in% slotNames(res)) {
        Boochs <- as.numeric(slot(res, "vi"))
      } else if (is.data.frame(res) || is.matrix(res)) {
        Boochs <- as.numeric(res[, 1])
      } else {
        Boochs <- as.numeric(res)
      }
      break
    }
  }
}

# 5) Bind back to np_spectra_joined (clean names, no suffixes)
vi_df <- data.frame(
  PRI = PRI,
  NDVI = NDVI,
  NDRE = NDRE,
  TCARIOSAVI = TCARIOSAVI,
  Datt3 = Datt3,
  CARI = CARI,
  Boochs = Boochs,
  check.names = FALSE
)

np_spectra_joined <- cbind(np_spectra_joined, vi_df)

# Define output path
outfile <- "G:/Branch_Experiment/r_outputs/np_spectra_VIs.csv"

# Write CSV with row names preserved
write.csv(np_spectra_joined,
          file = outfile,
          row.names = TRUE)

####################### INDEX vs WC IN NP (ALL SAMPLES) ###########################
# Read back in np_spectra_joined (if needed)
# Define input path
infile <- "G:/Branch_Experiment/r_outputs/np_spectra_VIs.csv"

# Read CSV back into R, keeping row names
np_spectra_joined <- read.csv(infile,
                              row.names = 1,
                              check.names = FALSE,
                              stringsAsFactors = FALSE)

# USER SETTINGS 
INDEX        <- "CARI"          # "PRI","NDVI","NDRE","TCARIOSAVI","Datt3","CARI","Boochs"
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
#agg <- read.csv("G:/Branch_Experiment/r_outputs/np_agg_medians.csv",
                #check.names = FALSE,
                #stringsAsFactors = FALSE)

# Normality tests: QQ and Shapiro-Wilk
par(mfrow = c(1, 2))  # 1 row, 2 columns

# Index
qqnorm(agg$Index, main = "QQ Plot of Index")
qqline(agg$Index, col = "red", lwd = 2)

# WC
qqnorm(agg$WC, main = "QQ Plot of WC")
qqline(agg$WC, col = "red", lwd = 2)

par(mfrow = c(1, 1))  # reset layout

# Shapiro-Wilk tests for normality
shapiro_Index <- shapiro.test(agg$Index)
shapiro_WC  <- shapiro.test(agg$WC)

shapiro_Index
shapiro_WC

# Fit a linear mixed-effects model:
#   Response: Index - PRI, NDVI, TCARI/OSAVI, etc.
#   Fixed effect: WC (common slope across trees)
#   Random effect: random intercept for each TreeNum (tree-specific baseline Index)
agg_model = lmer(Index~WC+(1|TreeNum),data=agg)
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

######################### R1 INDICES PER TREE #####################################
## --- Inputs ---
csv  <- "G:/Branch_Experiment/r_outputs/np_spectra.csv"
np_spectra <- read.csv(csv, row.names = 1, check.names = FALSE)

## --- Identify wavelength columns and order them numerically ---
is_wl_col <- function(nm) startsWith(nm, "nm_")
wl_from_names <- function(nms) as.numeric(gsub("_", ".", sub("^nm_", "", nms), fixed = TRUE))

wl_cols <- which(is_wl_col(colnames(np_spectra)))
stopifnot(length(wl_cols) > 0)

wl      <- wl_from_names(colnames(np_spectra)[wl_cols])
o       <- order(wl)
wl      <- wl[o]
X       <- as.matrix(np_spectra[, wl_cols[o], drop = FALSE])

## --- Helpers to grab nearest band and compute indices for one spectrum vector ---
nearest_idx <- function(target_nm) which.min(abs(wl - target_nm))
val_at <- function(v, nm) v[nearest_idx(nm)]

compute_indices <- function(v) {
  # bands used
  R531 <- val_at(v, 531); R570 <- val_at(v, 570)
  R550 <- val_at(v, 550); R670 <- val_at(v, 670)
  R680 <- val_at(v, 680); R681 <- val_at(v, 681)
  R700 <- val_at(v, 700); R709 <- val_at(v, 709)
  R720 <- val_at(v, 720); R754 <- val_at(v, 754)
  R800 <- val_at(v, 800)
  
  # indices (with small eps to avoid divide-by-zero)
  eps <- .Machine$double.eps
  
  PRI  <- (R531 - R570) / pmax(R531 + R570, eps)
  NDVI <- (R800 - R680) / pmax(R800 + R680, eps)
  
  # red-edge family
  NDRE <- (R800 - R720) / pmax(R800 + R720, eps)
  MTCI <- (R754 - R709) / pmax(R709 - R681, eps)
  TCARI <- 3 * ((R700 - R670) - 0.2 * (R700 - R550) * (R700 / pmax(R670, eps)))
  OSAVI <- (1 + 0.16) * (R800 - R670) / pmax(R800 + R670 + 0.16, eps)
  TCARI_OSAVI <- TCARI / pmax(OSAVI, eps)
  
  c(PRI = PRI, NDVI = NDVI, NDRE = NDRE, MTCI = MTCI, `TCARI/OSAVI` = TCARI_OSAVI)
}

## --- Compute for every row and append to np_spectra ---
idx_mat <- t(apply(X, 1, compute_indices))
idx_df  <- as.data.frame(idx_mat, stringsAsFactors = FALSE)

# Bind new columns at the end; keep original object name
np_spectra <- cbind(np_spectra, idx_df)

## (Optional) save updated table
# write.csv(np_spectra, "G:/Branch_Experiment/r_outputs/np_spectra_with_VIs.csv", row.names = TRUE)

############### BOX PLOT OF ROUND 1 SPECTRA #####################################
# USER SETTINGS
INDEX <- "NDRE.nm_800"   # choose from names in np_spectra: "PRI", "NDVI", "NDRE", "MTCI", "TCARI/OSAVI"
ROUND <- 1       # which round to subset (e.g. 1, 2, 3...)

## --- Parse tree and round from row names ---
get_tree  <- function(s) sub("^([^_]+)_.*$", "\\1", s)
get_round <- function(s) as.integer(sub(".*_R(\\d+)_.*", "\\1", s))

rn   <- rownames(np_spectra)
tree <- get_tree(rn)
rnd  <- get_round(rn)

## --- Build dataframe with selected index ---
if (!INDEX %in% colnames(np_spectra)) {
  stop("Column ", INDEX, " not found in np_spectra. Available: ", 
       paste(colnames(np_spectra), collapse = ", "))
}

df <- data.frame(Tree = tree, Round = rnd, Value = np_spectra[[INDEX]])

## --- Subset to selected round ---
df_sub <- df[df$Round == ROUND, ]

## --- Boxplot ---
boxplot(Value ~ Tree, data = df_sub,
        main = paste("R", ROUND, INDEX, "values per Tree"),
        xlab = "Tree",
        ylab = INDEX,
        outline = TRUE,   # show outliers
        col = "lightblue", border = "darkblue")

################# SPECTRAL PROFILES FROM ONE TREE OVER TIME #######################
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

###############################################################################
################# LEAF CLIP SPECTRA!!! CREATING THE DF #######################
## ---------- helper: read one SED reflectance (0–1) ----------
read_sed_reflect <- function(f) {
  lns <- readLines(f, warn = FALSE)
  i <- grep("^Data:\\s*$", lns)
  if (length(i) != 1) stop("Could not find 'Data:' section in ", f)
  hdr <- strsplit(lns[i + 1], "\t")[[1]]
  dat <- read.table(text = paste(lns[(i + 2):length(lns)], collapse = "\n"),
                    sep = "\t", header = FALSE, col.names = hdr, check.names = FALSE)
  wl   <- as.numeric(dat$Wvl)
  refl <- as.numeric(dat[["Reflect. %"]]) / 100  # convert % to fraction
  data.frame(wl = wl, refl = refl)
}

## ---------- list S1_Ba / S2_Ba files ----------
sed_files <- list.files(pattern = "^(S1_Ba|S2_Ba).*\\.sed$", ignore.case = TRUE)
if (length(sed_files) == 0) stop("No S1_Ba/S2_Ba .sed files found in the working directory.")

## ---------- read all ----------
dfs <- lapply(sed_files, read_sed_reflect)
wls <- lapply(dfs, function(d) d$wl)

## ---------- FAST PATH: identical wavelength grid? ----------
same_grid <- all(vapply(wls, function(w) identical(w, wls[[1]]), logical(1)))

if (same_grid) {
  # Stack directly (no resampling)
  wl_vec <- wls[[1]]
  mat <- do.call(rbind, lapply(dfs, function(d) d$refl))
  rownames(mat) <- tools::file_path_sans_ext(basename(sed_files))
  # clean column names like nm_350, nm_351, nm_351_5, ...
  colnames(mat) <- paste0("nm_", gsub("\\.", "_", format(wl_vec, trim = TRUE)))
  lc_spectra <- as.data.frame(mat, stringsAsFactors = FALSE)
} else {
  # SAFE PATH: exact intersection of wavelengths (still no interpolation)
  common_wl <- Reduce(intersect, wls)
  if (length(common_wl) < 50) {
    warning("Very few common wavelengths across files; check instrument settings.")
  }
  common_wl <- sort(common_wl)
  mat <- do.call(rbind, lapply(dfs, function(d) {
    idx <- match(common_wl, d$wl)
    d$refl[idx]
  }))
  rownames(mat) <- tools::file_path_sans_ext(basename(sed_files))
  colnames(mat) <- paste0("nm_", gsub("\\.", "_", format(common_wl, trim = TRUE)))
  lc_spectra <- as.data.frame(mat, stringsAsFactors = FALSE)
  
  message("Used exact intersection of wavelengths (no interpolation). Kept ",
          length(common_wl), " bands.")
}

## ---------- quick sanity checks ----------
dim(lc_spectra)        # rows = files, cols = wavelengths
lc_spectra[1, 1:10]    # first row, first few bands

# Define output folder
outdir <- "G:/Branch_Experiment/r_outputs"

# Write CSV (row names = your Tree file IDs)
outfile <- file.path(outdir, "lc_spectra.csv")
write.csv(np_spectra, file = outfile, row.names = TRUE)

#################### CLEAN LEAF CLIP SPECTRA ###############################
lc_spectra <- lc_spectra[!(rownames(lc_spectra) %in% 
                             c("S1_Ba_R1_P5_00005", "S1_Ba_R1_P5_00006")), ]

rownames(lc_spectra)[rownames(lc_spectra) == "S2_Ba_R5_P2_00004"] <- 
  "S2_Ba_R5_P3_00004"

# Define output folder
outdir <- "G:/Branch_Experiment/r_outputs"

# Write CSV (row names = your Tree file IDs)
outfile <- file.path(outdir, "lc_spectra.csv")
write.csv(np_spectra, file = outfile, row.names = TRUE)

# QC: Distance-to-median + NDVI sanity #
# INPUT: lc_spectra (rows = spectra IDs like "S1_Ba_R1_P5_00005", cols = nm_* bands)

# USER SETTINGS ==
WL_CORE      <- c(400, 900)   # distance computed on this wavelength window (nm)
NDVI_MIN_OK  <- 0.20          # NDVI below this is flagged (tune: 0.15–0.30 common)
MAD_K        <- 2.5             # flag distances > median(dist) + K*MAD within group
SAVE_CSV     <- FALSE
OUTFILE      <- "G:/Branch_Experiment/r_outputs/lc_spectra_QC_flags.csv"
#

stopifnot(exists("lc_spectra"))

# Helpers
is_wl_col <- function(nm) startsWith(nm, "nm_")
wl_num_of <- function(nms) as.numeric(gsub("_", ".", sub("^nm_", "", nms), fixed = TRUE))
get_tree  <- function(s) sub("^([^_]+)_.*$", "\\1", s)                  # "S1_Ba" -> tree label block
get_round <- function(s) as.integer(sub(".*_R(\\d+)_.*", "\\1", s))     # R number

# Identify wavelength columns, order them numerically
wl_cols <- which(is_wl_col(colnames(lc_spectra)))
stopifnot(length(wl_cols) > 0)
wl_all  <- wl_num_of(colnames(lc_spectra)[wl_cols])
ord     <- order(wl_all)
wl_all  <- wl_all[ord]
X_all   <- as.matrix(lc_spectra[, wl_cols[ord], drop = FALSE])

# Restrict to WL_CORE window for distance calc (avoids noisy edges)
keep_wl <- wl_all >= WL_CORE[1] & wl_all <= WL_CORE[2]
stopifnot(any(keep_wl))
wl_core <- wl_all[keep_wl]
X_core  <- X_all[, keep_wl, drop = FALSE]

# Nearest-band helpers for indices
nearest_idx <- function(target_nm) which.min(abs(wl_all - target_nm))
val_at <- function(v, nm) v[nearest_idx(nm)]

# Parse grouping fields
ids   <- rownames(lc_spectra)
tree  <- get_tree(ids)   # e.g., "S1_Ba" and "S2_Ba"
round <- get_round(ids)  # numeric 1,2,3,...

# Compute NDVI per row (using nearest 800 & 680 nm in your grid)
R800 <- X_all[, nearest_idx(800)]
R680 <- X_all[, nearest_idx(680)]
eps  <- .Machine$double.eps
NDVI <- (R800 - R680) / pmax(R800 + R680, eps)
ndvi_flag <- NDVI < NDVI_MIN_OK | !is.finite(NDVI)

# Distance-to-median per (tree, round)
grp <- paste(tree, round, sep = " | ")
ug  <- unique(grp)

dist_vec <- numeric(nrow(X_core))
dist_vec[] <- NA_real_

for (g in ug) {
  idx <- which(grp == g)
  if (length(idx) < 2) {
    # With only 1 spectrum, distance=0 by definition
    dist_vec[idx] <- 0
    next
  }
  # robust group median spectrum (per wavelength)
  med_spec <- apply(X_core[idx, , drop = FALSE], 2, median, na.rm = TRUE)
  # Euclidean distance (you can swap for cosine if you prefer shape-only)
  d <- sqrt(rowSums((X_core[idx, , drop = FALSE] - rep(med_spec, each = length(idx)))^2, na.rm = TRUE))
  dist_vec[idx] <- d
}

# Robust cutoff using MAD within each group
dist_flag <- logical(length(dist_vec))
for (g in ug) {
  idx <- which(grp == g)
  d   <- dist_vec[idx]
  if (all(!is.finite(d))) next
  med_d <- median(d, na.rm = TRUE)
  mad_d <- mad(d, constant = 1.4826, na.rm = TRUE)  # ~sd under normality
  thr   <- med_d + MAD_K * mad_d
  dist_flag[idx] <- d > thr
}

# Assemble QC table
qc <- data.frame(
  id          = ids,
  tree        = tree,
  round       = round,
  NDVI        = NDVI,
  ndvi_flag   = ndvi_flag,
  dist_core   = dist_vec,
  dist_flag   = dist_flag,
  suspect     = ndvi_flag | dist_flag,
  stringsAsFactors = FALSE
)

# (Optional) save flags
if (isTRUE(SAVE_CSV)) {
  write.csv(qc, OUTFILE, row.names = FALSE)
  message("Wrote QC flags to: ", OUTFILE)
}

# Quick summary
cat("\nQC summary:\n")
print(table(qc$round, qc$suspect, dnn = c("Round", "Suspect")))
cat("\nTotal suspect:", sum(qc$suspect), "of", nrow(qc), "rows\n")

# Example: drop suspects or keep only clean
lc_spectra_clean <- lc_spectra[!qc$suspect, , drop = FALSE]
lc_spectra_flagged <- merge(lc_spectra, qc[, c("id","suspect","ndvi_flag","dist_flag","NDVI")],
                             by.x = 0, by.y = "id", all.x = TRUE, sort = FALSE)
rownames(lc_spectra_flagged) <- lc_spectra_flagged$Row.names; lc_spectra_flagged$Row.names <- NULL


################# ADD WATER CONTENT DATA TO LEAF SPECTRA #########################
lc_experiment_data <-read.csv("G:/Branch_Experiment/lc_experiment_data.csv")

# JOIN lc_experiment_data → lc_spectra_clean (by prefix)
stopifnot(exists("lc_spectra_clean"), exists("lc_experiment_data"))

# Extract prefix key like "S1_Ba_R1" from IDs such as "S1_Ba_R1_P3_00001"
extract_prefix <- function(id) sub("^((?:[^_]+)_[^_]+_R\\d+).*", "\\1", id)

# Preserve ID as a proper column for the join and to restore row order later
lc_spectra_clean$._ID <- rownames(lc_spectra_clean)
lc_spectra_clean$KeyPrefix <- extract_prefix(rownames(lc_spectra_clean))

# Expect lc_experiment_data to have a column named "Keyfield" with values like "S1_Ba_R1"
if (!"Keyfield" %in% colnames(lc_experiment_data)) {
  stop("lc_experiment_data must contain a column named 'Keyfield'.")
}

# Left-join (many spectra to one experiment row)
lc_spectra_joined <- merge(
  lc_spectra_clean,
  lc_experiment_data,
  by.x = "KeyPrefix",
  by.y = "Keyfield",
  all.x = TRUE,
  sort = FALSE
)

# Restore rownames and drop helper columns (keep KeyPrefix if you find it handy)
rownames(lc_spectra_joined) <- lc_spectra_joined$._ID
lc_spectra_joined$._ID <- NULL
# (Optionally) drop KeyPrefix:
# lc_spectra_joined$KeyPrefix <- NULL

# Optional: quick sanity check
# head(lc_spectra_joined[, c("KeyPrefix", setdiff(names(lc_experiment_data), "Keyfield")), drop = FALSE])

# Define output folder
outdir <- "G:/Branch_Experiment/r_outputs"

# Write CSV (row names = your Tree file IDs)
outfile <- file.path(outdir, "lc_spectra_clean.csv")
write.csv(lc_spectra_joined, file = outfile, row.names = TRUE) # Output is cleaned and joined lc_spectra

############## CALCULATE SPECTRAL INDICIES FOR LC SPECTRA ########################
# APPEND SPECTRAL INDICES TO lc_spectra_joined
stopifnot(exists("lc_spectra_joined"))

# Identify wavelength columns and order by numeric wavelength
is_wl_col <- function(nm) startsWith(nm, "nm_")
wl_from_names <- function(nms) as.numeric(gsub("_", ".", sub("^nm_", "", nms), fixed = TRUE))

wl_cols <- which(is_wl_col(colnames(lc_spectra_joined)))
stopifnot(length(wl_cols) > 0)

wl_all <- wl_from_names(colnames(lc_spectra_joined)[wl_cols])
ord    <- order(wl_all)
wl     <- wl_all[ord]
X      <- as.matrix(lc_spectra_joined[, wl_cols[ord], drop = FALSE])

# Helpers
nearest_idx <- function(target_nm) which.min(abs(wl - target_nm))
val_at      <- function(v, nm) v[nearest_idx(nm)]
eps         <- .Machine$double.eps

# Safely extract a numeric from hsdar::vegindex() outputs of varying classes
.safe_num <- function(x) {
  tryCatch({
    if (isS4(x)) {
      sn <- slotNames(x)
      if ("vi" %in% sn) return(as.numeric(slot(x, "vi")))
      if ("si" %in% sn) return(as.numeric(slot(x, "si")))
      # fall through if unknown S4 layout
    }
    if (is.data.frame(x) || is.matrix(x)) return(as.numeric(x[1, 1]))
    as.numeric(x)[1]
  }, error = function(e) NA_real_)
}

# Try Boochs flavors via hsdar; return NA if not available
safe_boochs <- function(v) {
  if (!requireNamespace("hsdar", quietly = TRUE)) return(NA_real_)
  re <- matrix(v, nrow = 1); colnames(re) <- paste0("nm_", gsub("\\.", "_", as.character(wl)))
  sl <- hsdar::speclib(re, wavelength = wl)
  for (nm_idx in c("Boochs", "Boochs2", "Boochs1")) {
    res <- try(hsdar::vegindex(sl, index = nm_idx), silent = TRUE)
    if (!inherits(res, "try-error")) {
      val <- .safe_num(res)
      if (is.finite(val)) return(val)
    }
  }
  NA_real_
}

# Compute indices for one spectrum vector (v)
compute_indices <- function(v) {
  # Band picks (nearest)
  R531 <- val_at(v, 531); R570 <- val_at(v, 570)
  R550 <- val_at(v, 550); R670 <- val_at(v, 670)
  R680 <- val_at(v, 680); R700 <- val_at(v, 700)
  R704 <- val_at(v, 704); R709 <- val_at(v, 709)
  R720 <- val_at(v, 720); R754 <- val_at(v, 754)
  R800 <- val_at(v, 800); R681 <- val_at(v, 681)
  
  # Core indices
  PRI  <- (R531 - R570) / pmax(R531 + R570, eps)
  NDVI <- (R800 - R680) / pmax(R800 + R680, eps)
  
  # Red-edge family
  NDRE <- (R800 - R720) / pmax(R800 + R720, eps)
  MTCI <- (R754 - R709) / pmax(R709 - R681, eps)  # helpful diagnostic, not requested for join
  
  # TCARI / OSAVI and their ratio
  TCARI <- 3 * ((R700 - R670) - 0.2 * (R700 - R550) * (R700 / pmax(R670, eps)))
  OSAVI <- (1 + 0.16) * (R800 - R670) / pmax(R800 + R670 + 0.16, eps)
  TCARI_OSAVI <- TCARI / pmax(OSAVI, eps)
  
  # Datt3 (common form)
  Datt3 <- (R754 - R704) / pmax(R754 + R704, eps)
  
  # CARI (baseline across ~550–700 nm; practical variant)
  CARI <- abs((R700 - R670 - 0.2 * (R700 - R550)) * (R670 / pmax(R700, eps)))
  
  Boochs <- safe_boochs(v)
  
  c(
    PRI = PRI,
    NDVI = NDVI,
    `TCARI/OSAVI` = TCARI_OSAVI,
    Datt3 = Datt3,
    CARI = CARI,
    Boochs = Boochs,
    NDRE = NDRE
  )
}

# Apply to all rows and bind back
idx_mat <- t(apply(X, 1, compute_indices))
idx_df  <- as.data.frame(idx_mat, stringsAsFactors = FALSE)

# Attach to lc_spectra_joined
lc_spectra_joined <- cbind(lc_spectra_joined, idx_df)

# Heads-up if Boochs was not computed anywhere
if (all(is.na(lc_spectra_joined$Boochs))) {
  message("Boochs not computed (hsdar unavailable or index unsupported for your setup). ",
          "Column filled with NA; other indices are computed.")
}

colnames(lc_spectra_joined) <- sub("\\.nm_\\d+(?:_\\d+)?$", "", colnames(lc_spectra_joined), perl = TRUE)

colnames(lc_spectra_joined)[colnames(lc_spectra_joined) == "TCARI/OSAVI"] <- "TCARIOSAVI"

# Save to .csv after editing names
names(lc_spectra_joined)[names(lc_spectra_joined) == ""] <- "Unnamed"

# Rename the known columns
names(lc_spectra_joined)[names(lc_spectra_joined) == "Military.Time"]   <- "Military_Time"
names(lc_spectra_joined)[names(lc_spectra_joined) == "Fresh.Mass..g."]  <- "Fresh_Mass"
names(lc_spectra_joined)[names(lc_spectra_joined) == "Dry.Mass..g."]    <- "Dry_Mass"
names(lc_spectra_joined)[names(lc_spectra_joined) == "X.R.WC"]          <- "WC"
names(lc_spectra_joined)[names(lc_spectra_joined) == "Tree.."]          <- "Tree"
names(lc_spectra_joined)[names(lc_spectra_joined) == "X.."]             <- "Date"
names(lc_spectra_joined)[names(lc_spectra_joined) == "Branch.."]        <- "Branch"
names(lc_spectra_joined)[names(lc_spectra_joined) == "TCPchamber."]     <- "TCP_chamber"
names(lc_spectra_joined)[names(lc_spectra_joined) == "spectra."]        <- "spectra"
names(lc_spectra_joined)[names(lc_spectra_joined) == "Sampling.Round"]  <- "Sampling_Round"

write.csv(lc_spectra_joined, "G:/Branch_Experiment/r_outputs/lc_spectra_VIs.csv", row.names = TRUE)


############### PLOT INDICES PER ROUND ################################
lc_spectra_joined <- read.csv("G:/Branch_Experiment/r_outputs/lc_spectra_VIs.csv",
                              check.names = FALSE,    # keep original column names
                              stringsAsFactors = FALSE)

## ================= USER SETTINGS ===
INDEX      <- "Boochs"                 # one of: "PRI","TCARIOSAVI","NDVI","Boochs","NDRE","CARI","Datt3"
ROUNDS     <- c(1, 3, 5)            # which rounds to include (numeric), e.g., c(1,3,5) or c(1,5)
PLOT_TITLE <- "Boochs by Sampling Round (leaf-clip measurements)"
#

stopifnot(exists("lc_spectra_joined"))

# Require KeyPrefix to parse tree/round
if (!"KeyPrefix" %in% names(lc_spectra_joined)) {
  stop("lc_spectra_joined must contain a 'KeyPrefix' column like 'S1_Ba_R1'.")
}
if (!INDEX %in% names(lc_spectra_joined)) {
  stop("Index column '", INDEX, "' not found in lc_spectra_joined.")
}

# Parse Tree (S1/S2) and Round (integer) from KeyPrefix
kp <- as.character(lc_spectra_joined$KeyPrefix)
tree_code <- sub("^([^_]+)_.*", "\\1", kp)                        # "S1_Ba_R1" -> "S1"
round_num <- suppressWarnings(as.integer(sub(".*_R(\\d+).*", "\\1", kp)))  # -> 1

# Build plotting data
df <- data.frame(
  Tree  = factor(ifelse(tree_code == "S1", "Spruce 1",
                        ifelse(tree_code == "S2", "Spruce 2", NA_character_)),
                 levels = c("Spruce 1","Spruce 2")),
  Round = factor(paste0("R", round_num), levels = paste0("R", ROUNDS)),
  Value = lc_spectra_joined[[INDEX]]
)

# Filter to requested rounds and finite values
df <- df[!is.na(df$Tree) & !is.na(df$Round) &
           df$Round %in% paste0("R", ROUNDS) &
           is.finite(df$Value), , drop = FALSE]
if (nrow(df) == 0) stop("No data after filtering. Check INDEX and ROUNDS.")

# Shared y-axis across panels for fair comparison
ylim_all <- range(df$Value, finite = TRUE)

# Plot: two panels (Spruce 1, Spruce 2)
op <- par(mfrow = c(1, 2), mar = c(5,4,3,1), oma = c(0,0,3,0))
on.exit(par(op), add = TRUE)

for (spruce in c("Spruce 1","Spruce 2")) {
  subdf <- df[df$Tree == spruce, , drop = FALSE]
  if (nrow(subdf) > 0) {
    boxplot(Value ~ Round, data = subdf,
            main = spruce, xlab = "Round", ylab = INDEX,
            outline = TRUE, ylim = ylim_all)
  } else {
    # keep panel space even if no data
    plot.new(); title(main = spruce)
    mtext("No data for selected rounds", side = 3, line = 0.5, cex = 0.9)
  }
}

mtext(PLOT_TITLE, outer = TRUE, cex = 1.2, line = 1)
################ SCATTERPLOT: WC vs PRI ########################
################ PRI ~ Water Content (Boxplots) ################
## User settings
PLOT_TITLE <- "Leaf Clip Spectra, PRI vs Water Content"
FACET_BY_TREE <- FALSE          # TRUE = one panel per spruce; FALSE = combined
TREES <- c("S1","S2")          # which trees to include (by KeyPrefix code)
ROUNDS <- c(1,3,5)             # which rounds to include (numeric)
POINTS <- FALSE                 # overlay jittered points on boxes

stopifnot(exists("lc_spectra_joined"))
needed <- c("KeyPrefix","PRI","WC")
if (!all(needed %in% names(lc_spectra_joined))) {
  stop("lc_spectra_joined must contain columns: ", paste(needed, collapse=", "))
}

kp <- as.character(lc_spectra_joined$KeyPrefix)
Tree_code  <- sub("^([^_]+)_.*", "\\1", kp)                          # "S1_Ba_R1" -> "S1"
Round_num  <- suppressWarnings(as.integer(sub(".*_R(\\d+).*", "\\1", kp)))

df <- data.frame(
  Tree  = factor(ifelse(Tree_code == "S1","Spruce 1",
                        ifelse(Tree_code == "S2","Spruce 2", NA_character_)),
                 levels=c("Spruce 1","Spruce 2")),
  Round = Round_num,
  WC    = lc_spectra_joined$WC,
  PRI   = lc_spectra_joined$PRI
)

# filter
df <- df[!is.na(df$Tree) &
           df$Round %in% ROUNDS &
           Tree_code %in% TREES &
           is.finite(df$PRI) & is.finite(df$WC), , drop=FALSE]
if (nrow(df) == 0) stop("No rows after filtering; check TREES/ROUNDS.")

# make WC a factor ordered by numeric value so boxes appear sorted by WC
wc_levels <- sort(unique(df$WC))
df$WCf <- factor(df$WC, levels = wc_levels, labels = format(wc_levels))

# layout
if (FACET_BY_TREE) {
  trees <- levels(droplevels(df$Tree))
  op <- par(mfrow=c(1, length(trees)), mar=c(6,4,3,1), oma=c(0,0,3,0))
  on.exit(par(op), add=TRUE)
  ylim_all <- range(df$PRI, finite=TRUE)
  for (tr in trees) {
    subdf <- df[df$Tree == tr, , drop=FALSE]
    if (nrow(subdf)==0) { plot.new(); title(tr); next }
    boxplot(PRI ~ WCf, data=subdf, ylim=ylim_all,
            xlab="Water Content (%)", ylab="PRI", main=tr, las=2, outline=TRUE)
    if (POINTS) stripchart(PRI ~ WCf, data=subdf, vertical=TRUE, method="jitter",
                           pch=19, col=adjustcolor("black", 0.5), add=TRUE)
  }
  mtext(PLOT_TITLE, outer=TRUE, cex=1.2, line=1)
} else {
  op <- par(mar=c(6,4,3,1))
  on.exit(par(op), add=TRUE)
  boxplot(PRI ~ WCf, data=df, xlab="Water Content", ylab="PRI",
          main=PLOT_TITLE, las=2, outline=TRUE)
  if (POINTS) stripchart(PRI ~ WCf, data=df, vertical=TRUE, method="jitter",
                         pch=19, col=adjustcolor("black", 0.5), add=TRUE)
}

################ PRI vs Water Content (all points) ################
# USER SETTINGS
PLOT_TITLE <- "PRI vs Relative Water Content (all replicates)"
COLOR_BY_TREE <- TRUE   # TRUE = color S1 vs S2 differently
POINT_ALPHA   <- 0.5    # transparency for overlapping points

stopifnot(exists("lc_spectra_joined"))
needed <- c("KeyPrefix","PRI","X.R.WC")
if (!all(needed %in% names(lc_spectra_joined))) {
  stop("lc_spectra_joined must contain columns: ", paste(needed, collapse=", "))
}

# parse tree codes (S1/S2 -> Spruce labels)
kp <- as.character(lc_spectra_joined$KeyPrefix)
tree_code <- sub("^([^_]+)_.*", "\\1", kp)
Tree <- factor(ifelse(tree_code=="S1","Spruce 1",
                      ifelse(tree_code=="S2","Spruce 2", NA_character_)),
               levels=c("Spruce 1","Spruce 2"))

df <- data.frame(Tree=Tree,
                 PRI = lc_spectra_joined$PRI,
                 WC  = lc_spectra_joined$X.R.WC)

df <- df[is.finite(df$PRI) & is.finite(df$WC) & !is.na(df$Tree), , drop=FALSE]
if (nrow(df) == 0) stop("No valid rows with PRI and WC.")

# plot
cols <- c("Spruce 1"="darkgreen","Spruce 2"="purple")
plot(df$WC, df$PRI,
     xlab="Relative Water Content (%)",
     ylab="PRI",
     main=PLOT_TITLE,
     pch=19,
     col=if (COLOR_BY_TREE) adjustcolor(cols[df$Tree], alpha.f=POINT_ALPHA)
     else adjustcolor("black", alpha.f=POINT_ALPHA))

if (COLOR_BY_TREE) legend("topright", legend=levels(df$Tree),
                          col=cols, pch=19, bty="n")

################# LEAF CLIP SPECTRA!!! PLOTTING OUT SPECTRAL PROFILE #######################

# Uses lc_spectra_joined and parses rounds from the "Unnamed" ID column
lc_spectra_joined <- read.csv("G:/Branch_Experiment/r_outputs/lc_spectra_VIs.csv",
                              check.names = FALSE,    # keep original column names
                              stringsAsFactors = FALSE)

# User settings
S_CHOICE     <- "S2"                     # "S1" or "S2"
NORMALIZE    <- FALSE                    # TRUE = scale each round-median to [0,1]
MEDIAN_RANGE <- c(400, 1000)             # nm used to compute round medians
VIEW_RANGE   <- c(500, 700)             # nm shown in plot; set NULL for full
ROUNDS       <- c(1, 3, 5)            # which rounds to include
PLOT_TITLE   <- "Leaf Clip Spectra by Round, Spruce 2 - 500 to 700 nm"
#

stopifnot(exists("lc_spectra_joined"))
if (!all(c("Unnamed","KeyPrefix") %in% names(lc_spectra_joined))) {
  stop("Expected columns 'Unnamed' and 'KeyPrefix' in lc_spectra_joined.")
}

# --- Select rows for the chosen spruce using the Unnamed measurement IDs ---
ids <- as.character(lc_spectra_joined$Unnamed)
pat <- paste0("^", S_CHOICE, "_Ba_")        # e.g., "^S1_Ba_"
ix  <- grepl(pat, ids, ignore.case = FALSE)
if (!any(ix)) stop("No rows matched pattern: ", pat)

spruce_ids <- ids[ix]

# Parse round number from Unnamed (e.g., S1_Ba_R3_P2_00001 -> 3)
get_round_num <- function(s) as.integer(sub(".*_R(\\d+)_.*", "\\1", s))
round_id <- vapply(spruce_ids, get_round_num, integer(1), USE.NAMES = FALSE)

# --- Prepare wavelength axis (numeric) from nm_* columns ---
wl_cols <- which(startsWith(names(lc_spectra_joined), "nm_"))
if (length(wl_cols) == 0) stop("No 'nm_' wavelength columns found.")
wl_chr <- sub("^nm_", "", names(lc_spectra_joined)[wl_cols])
wl_num <- suppressWarnings(as.numeric(gsub("_", ".", wl_chr, fixed = TRUE)))
o <- order(wl_num)
wl_num <- wl_num[o]; wl_cols <- wl_cols[o]

# Indices for median computation and plotting window
idx_median <- wl_num >= MEDIAN_RANGE[1] & wl_num <= MEDIAN_RANGE[2]
if (!any(idx_median)) stop("No wavelengths in MEDIAN_RANGE.")

idx_view <- if (is.null(VIEW_RANGE)) rep(TRUE, length(wl_num)) else
  (wl_num >= VIEW_RANGE[1] & wl_num <= VIEW_RANGE[2])
if (!any(idx_view)) stop("No wavelengths in VIEW_RANGE.")

# --- Compute per-round column-wise medians (using only MEDIAN_RANGE bands) ---
rounds_found <- sort(unique(round_id))
rounds_to_use <- intersect(rounds_found, ROUNDS)
if (length(rounds_to_use) == 0) stop("No matching rounds found in ROUNDS setting.")

# Grab the numeric matrix of reflectance for selected rows, ordered by wavelength
X_sel <- as.matrix(lc_spectra_joined[ix, wl_cols, drop = FALSE])

med_list <- lapply(rounds_to_use, function(r) {
  rows_r <- which(round_id == r)
  med_full <- apply(X_sel[rows_r, , drop = FALSE], 2, function(x) median(x, na.rm = TRUE))
  med_full[idx_median]
})

# Build matrix aligned to wavelengths in MEDIAN_RANGE
wl_num_med <- wl_num[idx_median]
med_mat <- do.call(rbind, med_list)
rownames(med_mat) <- paste0("R", rounds_to_use)
colnames(med_mat) <- paste0("nm_", wl_num_med)

# --- Optional normalization (shape-only) ---
if (NORMALIZE) {
  rng <- apply(med_mat, 1, function(v) diff(range(v, na.rm = TRUE)))
  mn  <- apply(med_mat, 1, function(v) min(v, na.rm = TRUE))
  med_mat <- (med_mat - mn[col(med_mat)]) / pmax(rng[col(med_mat)], .Machine$double.eps)
}

# --- Apply VIEW_RANGE for plotting ---
idx_view_med <- idx_view[idx_median]
wl_plot  <- wl_num_med[idx_view_med]
plot_mat <- med_mat[, idx_view_med, drop = FALSE]

# --- Plot ---
cols <- grDevices::rainbow(nrow(plot_mat))
ymin <- min(plot_mat, na.rm = TRUE)
ymax <- max(plot_mat, na.rm = TRUE)

plot(wl_plot, plot_mat[1, ], type = "l", lwd = 2, col = cols[1],
     xlab = "Wavelength (nm)",
     ylab = if (NORMALIZE) "Scaled reflectance (0–1)" else "Reflectance (fraction)",
     main = PLOT_TITLE,
     ylim = c(ymin, ymax))

if (nrow(plot_mat) > 1) {
  for (i in 2:nrow(plot_mat)) lines(wl_plot, plot_mat[i, ], lwd = 2, col = cols[i])
}
legend("topleft", legend = rownames(plot_mat), lwd = 2, col = cols, bty = "n")
