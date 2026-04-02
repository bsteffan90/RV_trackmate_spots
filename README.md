# TrackMate Infection Dynamics Pipeline – README

## Overview
This R script processes fluorescence microscopy data from TrackMate "spots" CSV files to extract single-cell calcium dynamics features. It analyzes 60 image series across 12 treatment conditions, computing onset times, peak intensities, and response durations for each tracked cell.

**Outputs:**
- `track_features_master.csv` – Per-track metrics (7 features × ~5,785 tracks)
- `summary_by_treatment.csv` – Treatment-level summaries

---

## Quick Start

1. **Set user parameters** (lines 7–14):
   ```r
   project_dir <- "C:/path/to/your/data"  # Folder containing 60 CSV files
   plate_date  <- "2026.02.27"             # Must match filename pattern
   donor_cells <- "B02-25"                 # Cell line identifier
   ```

2. **Run the script** in RStudio or R console. Dependencies auto-install.

3. **Check outputs** in `project_dir`:
   - `track_features_master.csv` (all tracks)
   - `summary_by_treatment.csv` (aggregated stats)

---

## Analysis Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `INTENSITY_COL` | `MEAN_INTENSITY_CH1` | Fluorescence metric (or `TOTAL_INTENSITY_CH1`) |
| `THRESHOLD_DFF` | 3.0 | ΔF/F₀ threshold for response onset |
| `K_CONSEC` | 2 | Frames above threshold to confirm onset |
| `SMOOTH_K` | 5 | Rolling median window (frames) |
| `MIN_TRACK_PTS` | 4 | Minimum points per track (shorter tracks dropped) |

---

## Pipeline Sections

### Section 0: Setup
- Clears workspace, sets working directory, loads `tidyverse` + `zoo`

### Section 1: File Discovery
- Locates 60 CSV files matching pattern `p100-YYYY.MM.DD_spots_seriesN.csv`
- Extracts series numbers (1–60)

### Section 2: Metadata Mapping
- Maps series 1–5 → Well A1, series 6–10 → Well A2, etc.
- Assigns treatment labels (e.g., "AA_INDO", "AA_PGE2_IL13")
- Creates metadata table with well, treatment, field index, replicates

### Section 3: Per-File Processing (`process_file()`)
For each CSV:
1. Reads TrackMate spots data
2. Validates required columns
3. **Background correction:** 10th percentile intensity per timepoint
4. **Per-track analysis:**
   - Computes ΔF/F₀ (normalized fluorescence change)
   - Smooths with rolling median
   - Detects **onset** (first sustained rise ≥ threshold)
   - Finds **peak** intensity and time
   - Calculates **duration** (time above threshold)
   - Computes **AUC** (area under curve)

### Section 4: Batch Processing
- Loops through all 60 files
- Combines results into master table with metadata

### Section 5: Summary Statistics
- Groups by treatment
- Computes mean onset, peak, AUC, duration per condition
- Saves summary CSV

### Section 6: Optional Plots
- Uncomment to generate violin/box plots of onset by treatment

---

## Output Columns

**`track_features_master.csv`:**
- `t_on_min` – Onset time (minutes)
- `t_peak_min` – Peak time (minutes)
- `dFF_peak` – Peak ΔF/F₀
- `auc_corr` – Area under curve (background-corrected)
- `duration_min` – Response duration (minutes)
- `t_clear_min` – Time when signal clears
- `cleared` – Binary (1 = response detected, 0 = no response)
- `treatment`, `well`, `cells`, `plate_date`, `field_idx`, `bio_rep`, `tech_rep`

---

## Troubleshooting

**"replacement has length zero"** → Check that `|` operators are `||` (logical OR, not bitwise) in `rollmedian()` calls.

**Missing columns** → Verify CSV has `TRACK_ID`, `POSITION_T`, `FRAME`, and intensity column.

**No output** → Check file naming matches pattern exactly (date, series numbers).

---

## Dependencies
- `tidyverse` (dplyr, readr, purrr, stringr)
- `zoo` (rolling statistics)

Auto-installed on first run.