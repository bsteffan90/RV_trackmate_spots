############################################################
# P101 Infection Dynamics - FULL R-ONLY PIPELINE (FINAL)
# Processes all TrackMate "spots" CSVs directly in R
# Output: track_features_master.csv + summary_by_treatment.csv
############################################################

#### ───────────────────────────────────────────────────────
#### User settings
#### ───────────────────────────────────────────────────────
project_dir <- "C:/Users/bsteffan/Desktop/P101Project"  # <— change if needed
plate_date  <- "2026.03.01"
donor_cells <- "B02-25"

# Analysis parameters
INTENSITY_COL <- "MEAN_INTENSITY_CH1"  # set to "TOTAL_INTENSITY_CH1" if desired
THRESHOLD_DFF <- 3.0                   # onset threshold in ΔF/F0
K_CONSEC      <- 2                     # sustained frames above/below threshold
SMOOTH_K      <- 5                     # rolling median window (frames)
MIN_TRACK_PTS <- 4                     # drop very short tracks

#### ───────────────────────────────────────────────────────
#### Package bootstrap (auto-install if missing)
#### ───────────────────────────────────────────────────────
req_pkgs <- c("tidyverse","zoo")
inst <- installed.packages()[,"Package"]
to_install <- setdiff(req_pkgs, inst)
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

#### ───────────────────────────────────────────────────────
#### 0) Start clean & set working directory
#### ───────────────────────────────────────────────────────
rm(list = setdiff(ls(), c("project_dir","plate_date","donor_cells",
                          "INTENSITY_COL","THRESHOLD_DFF","K_CONSEC",
                          "SMOOTH_K","MIN_TRACK_PTS")))
gc()
setwd(project_dir)
cat("Working dir:", getwd(), "\n")

#### ───────────────────────────────────────────────────────
#### 1) Locate files and derive series number (1..60)
#### ───────────────────────────────────────────────────────
files <- list.files(pattern = "p101-2026\\.03\\.01_spots_series\\d+\\.csv")
stopifnot(length(files) == 60)

series_no <- as.numeric(stringr::str_extract(files, "(?<=series)\\d+"))
stopifnot(all(!is.na(series_no)))

key <- tibble(file = files, series = series_no) %>%
  arrange(series)

#### ───────────────────────────────────────────────────────
#### 2) Build well/treatment key from series mapping
####    (5 fields per well)
#### ───────────────────────────────────────────────────────
series_map <- tribble(
  ~series_range,           ~well, ~treatment,
  1:5,                     "A1", "AA_INDO",
  6:10,                    "A2", "AA_PGE2_IL13",
  11:15,                   "A3", "AA_PGE2_IL13_INDO",
  16:20,                   "A4", "AA_PGE2_IL13_INDO",
  21:25,                   "B1", "AA",
  26:30,                   "B2", "AA_PGE2_IL13",
  31:35,                   "B3", "AA_IL13_INDO",
  36:40,                   "B4", "AA_PGE2_INDO",
  41:45,                   "C1", "AA_IL13",
  46:50,                   "C2", "AA_PGE2",
  51:55,                   "C3", "AA_PGE2_IL13",
  56:60,                   "C4", "AA_PGE2_INDO"
)

assign_meta <- function(s) {
  row <- series_map %>% dplyr::filter(purrr::map_lgl(series_range, ~ s %in% .x))
  tibble(well = row$well, treatment = row$treatment)
}

meta <- purrr::map_dfr(key$series, assign_meta)

key <- dplyr::bind_cols(key, meta) %>%
  mutate(
    cells      = donor_cells,
    plate_date = plate_date,
    field_idx  = ((series - 1) %% 5) + 1,
    tech_rep   = field_idx,      # 1..5 within well
    bio_rep    = plate_date
  )

stopifnot(nrow(key) == 60)
print(head(key, 10))

#### ───────────────────────────────────────────────────────
#### 3) Per-file processing function
#### ───────────────────────────────────────────────────────

process_file <- function(file) {
  
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  needed <- c("TRACK_ID","POSITION_T","FRAME", INTENSITY_COL)
  if (!all(needed %in% names(df))) return(NULL)
  
  df <- df %>%
    mutate(across(c(POSITION_T, FRAME, !!INTENSITY_COL),
                  ~ suppressWarnings(as.numeric(.)))) %>%
    filter(!is.na(TRACK_ID), !is.na(POSITION_T), !is.na(.data[[INTENSITY_COL]])) %>%
    arrange(TRACK_ID, POSITION_T)
  
  if (nrow(df) == 0) return(NULL)
  
  bg <- df %>%
    group_by(POSITION_T) %>%
    summarise(bg = quantile(.data[[INTENSITY_COL]], 0.10), .groups = "drop")
  
  df <- left_join(df, bg, by = "POSITION_T") %>%
    mutate(F_corr = pmax(.data[[INTENSITY_COL]] - bg, 0))
  
  results <- df %>%
    group_by(TRACK_ID) %>%
    group_modify(~{
      
      d <- .x
      if (nrow(d) < MIN_TRACK_PTS) {
        return(tibble(
          t_on_min     = NA_real_,
          t_peak_min   = NA_real_,
          dFF_peak     = NA_real_,
          auc_corr     = NA_real_,
          duration_min = NA_real_,
          t_clear_min  = NA_real_,
          cleared      = NA_integer_
        )[0, ])
      }
      
      t_min <- d$POSITION_T / 60
      Fc    <- d$F_corr
      F0    <- median(head(Fc, 3), na.rm = TRUE)
      dFF   <- (Fc - F0) / (F0 + 1e-6)
      
      dFF_s <- zoo::rollmedian(
        dFF,
        k = min(SMOOTH_K, max(3, sum(!is.na(dFF)) | 1)),
        fill = NA
      )
      
      above <- dFF_s >= THRESHOLD_DFF
      rle1  <- rle(above)
      
      t_on <- NA_real_; pos <- 1
      for (i in seq_along(rle1$lengths)) {
        if (rle1$values[i] && rle1$lengths[i] >= K_CONSEC) {
          t_on <- t_min[pos]; break
        }
        pos <- pos + rle1$lengths[i]
      }
      
      pk_idx   <- if (all(is.na(dFF_s))) NA_integer_ else which.max(dFF_s)
      t_peak   <- if (is.na(pk_idx)) NA_real_ else t_min[pk_idx]
      dFF_peak <- if (is.na(pk_idx)) NA_real_ else dFF_s[pk_idx]
      
      tibble(
        t_on_min     = t_on,
        t_peak_min   = t_peak,
        dFF_peak     = dFF_peak,
        auc_corr     = 0,
        duration_min = 0,
        t_clear_min  = NA_real_,
        cleared      = 0L
      )
    }) %>%
    ungroup()
  
  results
}

cat(">>> ENTERING SECTION 4 <<<\n")

#### ───────────────────────────────────────────────────────
#### 4) Batch across all 60 files → master table
#### ───────────────────────────────────────────────────────

all_tracks_list <- purrr::map2(
  key$file,
  seq_len(nrow(key)),
  ~{
    message("Processing: ", .x)
    
    feats <- process_file(.x)
    
    if (!is.null(feats) && nrow(feats) > 0) {
      dplyr::bind_cols(
        feats,
        key[.y, c(
          "file", "treatment", "well", "cells",
          "plate_date", "field_idx", "bio_rep", "tech_rep"
        )]
      )
    } else {
      NULL   # safe here
    }
  }
)

all_tracks <- dplyr::bind_rows(all_tracks_list)

#### ───────────────────────────────────────────────────────
#### 5) Quick summaries (sanity check)
#### ───────────────────────────────────────────────────────
summary_by_treatment <- all_tracks %>%
  filter(!is.na(t_on_min)) %>%
  group_by(treatment) %>%
  summarise(
    n_tracks      = n(),
    onset_mean    = mean(t_on_min, na.rm = TRUE),
    peak_mean     = mean(dFF_peak, na.rm = TRUE),
    auc_mean      = mean(auc_corr, na.rm = TRUE),
    duration_mean = mean(duration_min, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_by_treatment)
readr::write_csv(summary_by_treatment, "summary_by_treatment.csv")
cat("Saved: summary_by_treatment.csv\n")

#### ───────────────────────────────────────────────────────
#### 6) (Optional) tiny QC plots – uncomment to save PNGs
#### ───────────────────────────────────────────────────────
# if (nrow(all_tracks)) {
#   p_onset <- ggplot(all_tracks %>% filter(!is.na(t_on_min)),
#                     aes(x = treatment, y = t_on_min, fill = treatment)) +
#     geom_violin(alpha = 0.4) + geom_boxplot(width = 0.15, outlier.shape = NA) +
#     theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = "Onset by treatment", x = "Treatment", y = "Minutes")
#   ggsave("onset_by_treatment.png", p_onset, width = 8, height = 5, dpi = 300)
# }
############################################################
# END OF FILE
############################################################
# ---- helpers ----
safe_rollmedian <- function(x, k = 5) {
  n <- length(x)
  if (n < 3 || all(is.na(x))) return(x)
  # choose an odd window <= length, otherwise rollmedian() gets unhappy
  k_eff <- min(k, n - (1 - n %% 2))
  if (is.na(k_eff) || k_eff < 3) return(x)
  sm <- zoo::rollmedian(x, k = k_eff, fill = NA, align = "center")
  # replace NA edges with original values so we never propagate NAs into logic
  idx_na <- is.na(sm)
  if (any(idx_na)) sm[idx_na] <- x[idx_na]
  sm
}

safe_auc <- function(t_min, y) {
  if (length(t_min) < 2) return(0)
  # trapezoid with 2‑point rolling mean; when length==2 it's a single area
  mean_2 <- zoo::rollmean(y, 2, fill = NA, align = "left")
  sum(diff(t_min) * mean_2, na.rm = TRUE)
}

process_file <- function(file) {
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  needed <- c("TRACK_ID","POSITION_T","FRAME", INTENSITY_COL)
  if (!all(needed %in% names(df))) {
    warning("Missing columns in file: ", file)
    return(NULL)
  }
  
  df <- df %>%
    mutate(across(c(POSITION_T, FRAME, !!INTENSITY_COL),
                  ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(!is.na(TRACK_ID), !is.na(POSITION_T), !is.na(.data[[INTENSITY_COL]])) %>%
    arrange(TRACK_ID, POSITION_T)
  
  if (nrow(df) == 0) return(NULL)
  
  # frame‑wise background (10th percentile)
  bg <- df %>%
    group_by(POSITION_T) %>%
    summarise(bg = quantile(.data[[INTENSITY_COL]], 0.10), .groups = "drop")
  
  df <- left_join(df, bg, by = "POSITION_T") %>%
    mutate(F_corr = pmax(.data[[INTENSITY_COL]] - bg, 0))
  
  # ---- per‑track features ----
  results <- df %>%
    group_by(TRACK_ID) %>%
    group_modify(~{
      d <- .x
      if (nrow(d) < MIN_TRACK_PTS) return(NULL)
      
      t_min <- d$POSITION_T / 60
      Fc    <- d$F_corr
      
      F0  <- median(head(Fc, 3), na.rm = TRUE)
      dFF <- (Fc - F0) / (F0 + 1e-6)
      
      # robust smoothing that never errors on short vectors
      dFF_s <- safe_rollmedian(dFF, SMOOTH_K)
      
      # onset: first ≥ K_CONSEC consecutive frames >= THRESHOLD_DFF
      above <- dFF_s >= THRESHOLD_DFF
      r1    <- rle(above)
      t_on  <- NA_real_; pos <- 1
      for (i in seq_along(r1$lengths)) {
        if (r1$values[i] && r1$lengths[i] >= K_CONSEC) { t_on <- t_min[pos]; break }
        pos <- pos + r1$lengths[i]
      }
      
      # peak on smoothed curve
      pk_idx   <- which.max(dFF_s)
      t_peak   <- t_min[pk_idx]
      dFF_peak <- dFF_s[pk_idx]
      
      # metrics after onset
      if (!is.na(t_on)) {
        post      <- d %>% filter(POSITION_T/60 >= t_on)
        t_post    <- post$POSITION_T/60
        Fc_post   <- post$F_corr
        
        # safe post‑onset smoothing
        dFF_post_raw <- (Fc_post - F0)/(F0+1e-6)
        dFF_post     <- safe_rollmedian(dFF_post_raw, SMOOTH_K)
        
        above_post <- dFF_post >= THRESHOLD_DFF
        above_post[is.na(above_post)] <- FALSE
        
        # clearance: first ≥ K_CONSEC consecutive below threshold
        r2     <- rle(!above_post)
        t_clear <- NA_real_; pos2 <- 1
        for (j in seq_along(r2$lengths)) {
          if (r2$values[j] && r2$lengths[j] >= K_CONSEC) { t_clear <- t_post[pos2]; break }
          pos2 <- pos2 + r2$lengths[j]
        }
        
        duration <- if (any(above_post)) max(t_post[above_post]) - min(t_post[above_post]) else 0
        auc      <- safe_auc(t_post, Fc_post)
      } else {
        t_clear <- NA_real_; duration <- 0; auc <- 0
      }
      
      tibble(
        t_on_min     = t_on,
        t_peak_min   = t_peak,
        dFF_peak     = dFF_peak,
        auc_corr     = auc,
        duration_min = duration,
        t_clear_min  = t_clear,
        cleared      = ifelse(is.na(t_clear), 0, 1)
      )
    }) %>%
    ungroup()
  
  results
}
test_feats <- process_file(key$file[1])
nrow(test_feats); head(test_feats)