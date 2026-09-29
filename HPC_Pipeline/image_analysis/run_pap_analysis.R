# =============================================================================
# run_pap_analysis.R -- cropped six-well PAP plate images -> PAP feature tables.
#
# Run from this directory (local workstation, not the cluster; Fiji needs a
# desktop session on macOS or Windows):
#   PAP_CROP_DIR=/path/to/pap_cropped Rscript run_pap_analysis.R
#
# Reads  : $PAP_CROP_DIR/{parent,monoclonal,selected,passage}/*_crop.png
#          pap_roiset.zip                        fixed six-well ROI template
# Writes : $PAP_WORK_DIR/<set>/pap_feature_results.csv + qc/   PAPArea
#          ../../SourceData/<set>_pap_feature_results.csv  <- deposited as SourceData
#                                         (only when PAP_WRITE_SOURCEDATA=true)
#
# The cropped images were made from the original plate photographs with
#   AutoCrop(photoDir = <photos>, outputDir = <cropped>, plate = "six")
# Set PAP_PHOTO_DIR (same four sub-folders, *.jpg) to repeat that step first;
# the crops are then written to $PAP_WORK_DIR/<set>/cropped/.
#
# By default SourceData/ is not modified: each regenerated table is compared
# with the deposited one and reported as identical or different.
#
# Software used for the deposited tables:
#   diskImageR fork   pak::pak("tony27786/diskImageR")
#   Fiji              ImageJ 1.54p
#   R                 4.6.0
# =============================================================================

library(diskImageR)

get_script_dir <- function() {
  args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(args)) {
    return(dirname(normalizePath(sub("^--file=", "", args[1]), mustWork = TRUE)))
  }
  if (file.exists("run_pap_analysis.R")) return(normalizePath("."))
  stop("Run run_pap_analysis.R with Rscript, or from its own directory.")
}
SCRIPT_DIR <- get_script_dir()
SOURCE_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "SourceData"), mustWork = TRUE)

# =========================
# Settings
# =========================
CROP_DIR  <- Sys.getenv("PAP_CROP_DIR", file.path(SCRIPT_DIR, "pap_cropped"))
PHOTO_DIR <- Sys.getenv("PAP_PHOTO_DIR", "")   # optional: start from the photographs
WORK_DIR  <- Sys.getenv("PAP_WORK_DIR", file.path(SCRIPT_DIR, "pap_work"))
FIJI_APP  <- Sys.getenv("FIJI_APP", "/Applications/Fiji/Fiji.app")
ROI_ZIP   <- file.path(SCRIPT_DIR, "pap_roiset.zip")
WRITE_SOURCEDATA <- tolower(Sys.getenv("PAP_WRITE_SOURCEDATA", "false")) %in% c("1", "true", "yes")

# =========================
# Image sets
# =========================
# One image per six-well plate, named <sample>_crop.png; 000blank is a
# medium-only plate and stays in the table (the figure scripts drop it).
# selected and passage share the same 000blank photograph.
#   parent      <isolate>H / <isolate>L   two plates per isolate (PAPDual layout)
#                 L: TL..BR = 0, 0.5, 1, 2, 4, 8 ug/mL FLC
#                 H: TL..BR = 0, 8, 16, 32, 64, 128 ug/mL FLC
#   selected    <isolate>01   one plate per isolate (PAPSingle layout)
#   monoclonal  <isolate>02     TL..BR = 0, 8, 12, 16, 24, 32 ug/mL FLC
#   passage     <isolate>03
# The well -> FLC mapping is applied in Plot_Codes/, not here.
pap_sets <- data.frame(
  set      = c("parent", "monoclonal", "selected", "passage"),
  n_images = c(41L, 21L, 21L, 21L),
  stringsAsFactors = FALSE
)

wells <- c("TL", "TM", "TR", "BL", "BM", "BR")
feature_cols <- c("image", "well", "roi_area", "raw_mean", "raw_median", "raw_intden",
                  "mean_gray", "median_gray", "sd_gray", "min_gray", "max_gray",
                  "intden_gray", "cv_gray", "dark_frac_40", "dark_frac_60", "dark_frac_80")

# =========================
# Run one image set
# =========================
run_pap_set <- function(set, n_images) {
  message("\n=== ", set, " ===")
  proj_dir <- file.path(WORK_DIR, set)

  # 1. Optional: crop each photograph to the plate. Start from an empty
  #    folder, because PAPArea measures every image in it.
  if (nzchar(PHOTO_DIR)) {
    crop_dir <- file.path(WORK_DIR, set, "cropped")
    unlink(crop_dir, recursive = TRUE)
    AutoCrop(photoDir = file.path(PHOTO_DIR, set), outputDir = crop_dir,
             imageJLoc = FIJI_APP, plate = "six")
  } else {
    crop_dir <- file.path(CROP_DIR, set)
  }

  crops <- list.files(crop_dir)
  if (length(crops) != n_images || !all(grepl("_crop\\.png$", crops))) {
    stop(set, ": expected ", n_images, " *_crop.png images (and nothing else) in ", crop_dir,
         ", found ", length(crops), " files")
  }
  if (!"000blank_crop.png" %in% crops) stop(set, ": 000blank_crop.png is missing")

  # 2. Measure the six wells with the fixed ROI template (pap.ijm).
  #    PAPSingle()/PAPDual() call the same PAPArea() step with run_macro = TRUE;
  #    their normalization options only change the R objects they return,
  #    not pap_feature_results.csv.
  PAPArea(inputDir = crop_dir, projectDir = proj_dir, roiZip = ROI_ZIP,
          imageJLoc = FIJI_APP, overwrite = TRUE)
  out_csv <- file.path(proj_dir, "pap_feature_results.csv")
  if (!file.exists(out_csv)) stop(set, ": PAPArea did not write ", out_csv)

  d <- read.csv(out_csv, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(feature_cols, names(d))
  if (length(missing)) stop(set, ": output is missing ", paste(missing, collapse = ", "))
  wells_per_image <- tapply(d$well, d$image, function(w) setequal(w, wells) && length(w) == 6)
  if (length(wells_per_image) != n_images || !all(wells_per_image)) {
    stop(set, ": every image must have exactly the wells ", paste(wells, collapse = "/"))
  }

  # 3. Compare with (or replace) the deposited table.
  deposited <- file.path(SOURCE_DIR, paste0(set, "_pap_feature_results.csv"))
  status <- if (!file.exists(deposited)) {
    "no deposited file"
  } else if (unname(tools::md5sum(out_csv)) == unname(tools::md5sum(deposited))) {
    "identical"
  } else {
    "DIFFERENT"
  }
  if (WRITE_SOURCEDATA) {
    file.copy(out_csv, deposited, overwrite = TRUE)
    status <- paste(status, "-> written to SourceData")
  }
  data.frame(set = set, images = n_images, rows = nrow(d), vs_SourceData = status)
}

summary_tbl <- do.call(rbind, Map(run_pap_set, pap_sets$set, pap_sets$n_images))
rownames(summary_tbl) <- NULL
message("\n=== Summary ===")
print(summary_tbl)
