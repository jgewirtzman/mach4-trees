# =============================================================================
# 01_import.R
# Subset the LI-7810 .data file to the measurement dates, then import with
# goFlux::import.LI7810(). Saves the imported data frame to RData.
#
# The raw file is ~163 MB / ~950k 1-Hz rows spanning 2026-05-04..05-27.
# We keep only the field-campaign dates (built from the native file so the
# LI-COR header is preserved for import.LI7810).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))

# --- Build a date-subset of the native .data file ----------------------------

if (!file.exists(licor_subset)) {
  message("Subsetting LI-7810 file to: ", paste(measurement_dates, collapse = ", "))
  # Preserve the 7 header lines (Model..Timezone + DATAH + DATAU), then append
  # DATA rows whose DATE (column 7, tab-separated) is a measurement date.
  date_alt <- paste(measurement_dates, collapse = "|")
  awk_prog <- sprintf(
    'NR<=7 {print; next} $1=="DATA" && $7 ~ /^(%s)$/ {print}', date_alt)
  # Use shQuote to handle the path with spaces / the awk program safely.
  cmd <- sprintf("awk -F'\\t' %s %s > %s",
                 shQuote(awk_prog), shQuote(licor_raw), shQuote(licor_subset))
  status <- system(cmd)
  if (status != 0) stop("awk subsetting failed (status ", status, ")")
  message("Wrote subset: ", licor_subset)
} else {
  message("Using existing subset: ", licor_subset)
}

n_lines <- as.integer(strsplit(trimws(system(
  sprintf("wc -l < %s", shQuote(licor_subset)), intern = TRUE)), "\\s+")[[1]][1])
message("Subset rows (incl. header): ", n_lines)

# --- Import with goFlux::import.LI7810() --------------------------------------

message("Importing with import.LI7810() ...")
licor <- import.LI7810(
  inputfile   = licor_subset,
  date.format = "ymd",
  timezone    = tz_data,
  prec        = li7810_prec
)

message("Imported ", nrow(licor), " rows, ", ncol(licor), " columns")
message("Columns: ", paste(names(licor), collapse = ", "))
message("Time range: ", format(min(licor$POSIX.time)), " to ",
        format(max(licor$POSIX.time)), " (", tz_data, ")")

# Quick gas sanity (non-NA portion)
for (g in c("CO2dry_ppm", "CH4dry_ppb", "H2O_ppm")) {
  v <- licor[[g]]
  message(sprintf("  %-11s n=%d  range %.1f .. %.1f  (%d NA)",
                  g, sum(!is.na(v)), min(v, na.rm = TRUE),
                  max(v, na.rm = TRUE), sum(is.na(v))))
}

save(licor, file = file.path(rdata_dir, "licor_imported.RData"))
message("Saved: ", file.path(rdata_dir, "licor_imported.RData"))
message("Proceed to 02_build_auxfile.R")
