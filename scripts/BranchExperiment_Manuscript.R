#Set working directory to "coastalspruce" GitHub repo

################################################################################
################### BRANCH/NEEDLE SPECTRAL ANALYSIS ############################
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
outfile <- "./data/branch_experiment/np_spectra_VIs.csv"

# Write CSV with row names preserved
write.csv(np_spectra_joined,
          file = outfile,
          row.names = TRUE)

####################### INDEX vs WC IN NP (ALL SAMPLES) ###########################
# Read back in np_spectra_joined (if needed)
# Define input path
infile <- "./data/branch_experiment/np_spectra_VIs.csv"

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
#agg <- read.csv(./data/branch_experiment/np_agg_medians.csv",
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
