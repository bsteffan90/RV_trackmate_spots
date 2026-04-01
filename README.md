# RV_trackmate_spots
Using spots input file from trackmate to get infection dynamic output
Below is a **clean, reusable README** you can drop directly into your repository (for example as `README.md`).  
It is written as if this pipeline will be used long‑term, shared with lab members, and referenced later for methods or grants.

***

# Live Imaging Infection Dynamics Pipeline (TrackMate → R)

## Overview

This repository contains an **R‑only pipeline** for processing live‑cell imaging data tracked with **FIJI TrackMate**, extracting **single‑cell infection dynamics**, and generating standardized outputs suitable for **downstream statistical and biological analysis**.

The pipeline is designed to be:

*   ✅ **Deterministic**: identical inputs + parameters produce identical outputs
*   ✅ **Scalable**: can be run independently on many experiments
*   ✅ **Reproducible**: outputs include all metadata needed to combine datasets
*   ✅ **Decoupled**: raw data processing is separated from biological analysis

This pipeline should be run **once per experiment / plate**, after which the resulting CSVs can be combined across **technical and biological replicates** using a separate analysis script.

***

## Intended use

This pipeline is appropriate for:

*   Live imaging of infection or activation reporters
*   Single‑cell fluorescence time series tracked with TrackMate
*   Experiments with:
    *   multiple treatments
    *   multiple imaging fields (technical replicates)
    *   multiple plates, timepoints, or donors (biological replicates)

It is **not** intended to run statistical tests or combine experiments directly.  
Those steps belong in a separate **analysis layer**.

***

## Input data requirements

### 1. TrackMate CSV files

*   One CSV file per imaging field
*   File names must encode a **series index**, e.g.:

<!---->

    p101-2026.03.01_spots_series1.csv
    p101-2026.03.01_spots_series2.csv
    …

*   Required columns (from TrackMate):
    *   `TRACK_ID`
    *   `POSITION_T`
    *   `FRAME`
    *   intensity column (e.g. `MEAN_INTENSITY_CH1`)

Each file represents **one field of view**.

***

### 2. Experiment metadata (embedded in script)

The pipeline assigns metadata including:

*   treatment
*   well
*   plate date
*   donor / cell source
*   technical replicate (field index)
*   biological replicate (plate or donor)

> **Note:** For long‑term reuse, metadata assignment can be refactored into an external `metadata.csv`, but the current version uses an explicit mapping defined in the script.

***

## Processing steps (what the pipeline does)

For each TrackMate CSV file, the pipeline performs the following steps:

### 1️⃣ Data cleaning

*   Enforces numeric time and intensity values
*   Drops incomplete or malformed rows
*   Sorts by `TRACK_ID` and time

### 2️⃣ Frame‑wise background subtraction

*   Computes the **10th percentile intensity per frame**
*   Subtracts background from each spot
*   Corrected signal: `F_corr`

### 3️⃣ Per‑track normalization (ΔF/F₀)

*   Baseline `F₀` = median of the first 3 frames per track
*   Signal expressed as:

<!---->

    ΔF/F₀ = (F_corr − F₀) / F₀

### 4️⃣ Robust smoothing

*   Rolling median smoothing (user‑defined window)
*   Safe for short tracks and missing values
*   Preserves original data at boundaries

### 5️⃣ Infection event detection

*   **Infection onset**:
    *   First sustained (≥ K frames) crossing of ΔF/F₀ threshold
*   **Peak response**:
    *   Maximum smoothed ΔF/F₀ value
*   **Clearance**:
    *   First sustained return below threshold
*   **Duration**:
    *   Time spent above threshold
*   **Integrated burden**:
    *   Area under corrected fluorescence curve (AUC)

Each **TRACK\_ID** produces exactly **one feature row**.

***

## Output files

The pipeline produces two CSV files per experiment:

### 1️⃣ `track_features_master.csv` (primary output)

**One row per track (single‑cell level)**.

Key columns:

| Column         | Description                         |
| -------------- | ----------------------------------- |
| `t_on_min`     | Time of infection onset (minutes)   |
| `t_peak_min`   | Time of maximal response            |
| `dFF_peak`     | Peak ΔF/F₀ (infection burden)       |
| `auc_corr`     | Integrated corrected signal         |
| `duration_min` | Time above infection threshold      |
| `t_clear_min`  | Time of clearance (if applicable)   |
| `cleared`      | Binary indicator of resolution      |
| `treatment`    | Experimental condition              |
| `well`         | Plate well                          |
| `cells`        | Donor / cell source                 |
| `plate_date`   | Experiment date                     |
| `field_idx`    | Imaging field (technical replicate) |
| `tech_rep`     | Technical replicate index           |
| `bio_rep`      | Biological replicate identifier     |

This file is the **interface between processing and biological analysis**.

***

### 2️⃣ `summary_by_treatment.csv` (quick sanity summary)

Contains coarse, per‑treatment summaries such as:

*   number of infected tracks
*   mean time to onset
*   mean peak response
*   mean duration

This file is useful for:

*   quick sanity checks
*   confirming that processing behaved as expected

It is **not** intended to replace full analysis.

***

## How this pipeline fits into a larger analysis workflow

### ✅ Recommended workflow

1.  **Run this pipeline independently for each experiment**
    *   Produces one `track_features_master.csv` per plate / date

2.  **Store outputs in experiment‑specific folders**
        experiment_1/outputs/track_features_master.csv
        experiment_2/outputs/track_features_master.csv

3.  **Use a separate analysis script** to:
    *   load multiple `track_features_master.csv` files
    *   combine biological and technical replicates
    *   perform statistics and visualization

This separation:

*   prevents pseudoreplication
*   avoids unnecessary re‑processing
*   allows new experiments to be added seamlessly

***

## Analysis philosophy (important)

*   **Cells are not independent biological replicates**
*   Technical replicates (fields) should be collapsed at the **well level**
*   Statistical tests should be performed on **biological replicates**
*   Infection should be analyzed using **multiple orthogonal metrics**:
    *   infection incidence
    *   onset kinetics
    *   peak burden
    *   duration and clearance

***

## User‑defined parameters

The following parameters control infection detection and should be documented for every run:

```r
INTENSITY_COL  # fluorescence metric
THRESHOLD_DFF  # infection threshold (ΔF/F₀)
K_CONSEC       # sustained frames for onset/clearance
SMOOTH_K       # rolling median window
MIN_TRACK_PTS  # minimum track length
```

If these change, the pipeline should be rerun.

***

## Reproducibility notes

*   The pipeline is deterministic
*   All metadata required for downstream combination is embedded in outputs
*   The pipeline should **not** be modified between experiments that will be compared
*   Parameter values should be recorded alongside outputs

***

## Summary

This pipeline converts **raw TrackMate tracking data** into **single‑cell infection phenotypes** in a reproducible, scalable way.  
It is intentionally designed to **separate data processing from biological analysis**, enabling robust combination of technical and biological replicates across experiments.

***

