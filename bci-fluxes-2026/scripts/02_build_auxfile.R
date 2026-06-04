# =============================================================================
# 02_build_auxfile.R  (BCI)
# Parse the BCI field-notes CSV into a goFlux auxfile (one row per measurement).
#
# BCI field notes have ONE row per tree with TWO start-time columns:
#   Start_time_2m   : HH:MM start of the 2 m measurement   -> height_cm = 200
#   Start_time_Base : HH:MM start of the Base measurement  -> height_cm = 0  ("Base")
# We melt these into long format -> 2 measurements per tree.
#
# Output auxfile columns required by goFlux: UniqueID, start.time, Vtot, Area,
# Tcham, Pcham. field_meta carries species/height through to the final merge.
# Only START times are recorded; true start/end are set interactively in 03.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/bci-fluxes-2026/scripts/00_setup.R"))

# --- Read field notes --------------------------------------------------------

fn <- readr::read_csv(field_csv, show_col_types = FALSE,
                      na = c("", "NA", "MISSING"))
message("Field notes: ", nrow(fn), " tree rows")

# --- Melt the two start-time columns into long format ------------------------

long <- fn %>%
  tidyr::pivot_longer(
    cols      = c(Start_time_2m, Start_time_Base),
    names_to  = "position",
    values_to = "time_str"
  ) %>%
  filter(!is.na(time_str)) %>%
  mutate(
    meas_date    = lubridate::mdy(Date),   # handles 2- or 4-digit year (5/12/26)
    height_label = dplyr::recode(position,
                                 Start_time_2m   = "2m",
                                 Start_time_Base = "Base"),
    height_cm    = dplyr::recode(position,
                                 Start_time_2m   = 200,
                                 Start_time_Base = 0),
    TreeID       = paste0("BCI", Tree),
    species      = Species,
    # date + HH:MM in the instrument clock; field_offset_hours handles a Panama
    # vs New_York mismatch if needed (default 0).
    start.time   = as.POSIXct(paste(format(meas_date), time_str),
                              format = "%Y-%m-%d %H:%M", tz = tz_data) +
                   field_offset_hours * 3600,
    UniqueID     = paste0(TreeID, "_", height_label)
  ) %>%
  arrange(TreeID, start.time)

# --- Auxfile = start times + chamber geometry / T / P that goFlux needs -------

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
    UniqueID, TreeID, species,
    meas_date, time_str, height_label, height_cm,
    start.time
  )

# --- Report & save -----------------------------------------------------------

message("\n=== Auxfile (", nrow(aux), " measurements) ===")
print(as.data.frame(field_meta %>%
        select(UniqueID, species, meas_date, time_str, height_label)),
      row.names = FALSE)

save(aux, field_meta, file = file.path(rdata_dir, "auxfile.RData"))
write.csv(aux, file.path(results_dir, "auxfile.csv"), row.names = FALSE)
write.csv(field_meta, file.path(results_dir, "field_metadata.csv"), row.names = FALSE)
message("\nSaved auxfile + field metadata to RData/ and results/")
message("Proceed to 03_manual_id.R (interactive).")
