# =============================================================================
# 02_build_auxfile.R
# Parse the field-notes CSV into a goFlux auxfile (one row per measurement).
#
# Each measured tree row carries paired, semicolon-delimited lists:
#   fluxStartNotes_2026_05 : HH:MM start time of each height's measurement
#   heights_2026_05        : stem height (cm) of each measurement
# plus a single measurement date (date_2026_05). time[i] pairs with height[i].
#
# Output auxfile columns required by goFlux: UniqueID, start.time, Vtot, Area,
# Tcham, Pcham (+ obs.length and metadata we carry through for the final merge).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))

# --- Read field notes --------------------------------------------------------

fn <- readr::read_csv(field_csv, show_col_types = FALSE,
                      na = c("", "NA", "MISSING"))
message("Field notes: ", nrow(fn), " rows")

# Keep only rows that were actually measured (have start times AND heights).
measured <- fn %>%
  filter(!is.na(fluxStartNotes_2026_05), !is.na(heights_2026_05),
         !is.na(date_2026_05)) %>%
  mutate(meas_date = as.Date(gsub("_", "-", date_2026_05)))

message("Measured trees: ", nrow(measured),
        "  (", paste(measured$TreeID, collapse = ", "), ")")

# --- Explode paired start-time / height lists into long format ---------------

long <- measured %>%
  mutate(
    .time = strsplit(fluxStartNotes_2026_05, "\\s*;\\s*"),
    .hgt  = strsplit(as.character(heights_2026_05), "\\s*;\\s*")
  ) %>%
  rowwise() %>%
  mutate(.n = {
    if (length(.time) != length(.hgt))
      warning("Tree ", TreeID, ": ", length(.time), " times but ",
              length(.hgt), " heights")
    min(length(.time), length(.hgt))
  }) %>%
  ungroup() %>%
  select(Transect, TreeID, species, status_2026_03, Description,
         meas_date, .time, .hgt, .n) %>%
  # one row per (tree, measurement index)
  mutate(.idx = map(.n, seq_len)) %>%
  unnest(c(.idx)) %>%
  rowwise() %>%
  mutate(
    time_str  = .time[[.idx]],
    height_cm = suppressWarnings(as.numeric(.hgt[[.idx]]))
  ) %>%
  ungroup() %>%
  select(-.time, -.hgt, -.n)

# --- Build POSIX start time (date + HH:MM, UTC) ------------------------------

long <- long %>%
  mutate(
    start.time = as.POSIXct(
      paste(format(meas_date), time_str),
      format = "%Y-%m-%d %H:%M", tz = tz_data),
    UniqueID = paste0(TreeID, "_", height_cm, "cm")
  ) %>%
  arrange(TreeID, start.time)

# --- Auxfile = start times + the chamber geometry / T / P that goFlux needs ---
# Only START times were recorded. The display window length (obs_length) and
# shoulder are processing parameters applied in 03 via obs.win(); they are NOT
# stored here. True start/end are set interactively in 03_manual_id.R.

aux <- long %>%
  transmute(
    UniqueID,
    start.time,
    Vtot  = Vtot_L,
    Area  = chamber_area_cm2,
    Tcham = Tcham_default_C,
    Pcham = Pcham_default_kPa
  )

# --- Field metadata carried through to the final merge -----------------------

field_meta <- long %>%
  transmute(
    UniqueID, Transect, TreeID, species,
    status_2026_03, Description,
    meas_date, time_str, height_cm,
    start.time
  )

# --- Report & save -----------------------------------------------------------

message("\n=== Auxfile (", nrow(aux), " measurements) ===")
print(as.data.frame(field_meta %>%
        select(UniqueID, meas_date, time_str, height_cm)),
      row.names = FALSE)

save(aux, field_meta, file = file.path(rdata_dir, "auxfile.RData"))
write.csv(aux, file.path(results_dir, "auxfile.csv"), row.names = FALSE)
write.csv(field_meta, file.path(results_dir, "field_metadata.csv"), row.names = FALSE)
message("\nSaved auxfile + field metadata to RData/ and results/")
message("Proceed to 03_manual_id.R (interactive).")
