# =============================================================================
# 02_build_auxfile.R  (BASIC trees)
# Parse the transcribed basic-tree datasheet master (all_basic_sheets.csv) into
# a goFlux auxfile. Each tree row carries up to 3 heights (H1/H2/H3, cm) with a
# matching flux start time (flux30Start/flux80Start/flux160Start). We pair them
# POSITION-WISE (so a skipped height keeps the others aligned) and emit one
# measurement per (tree, height) that has a flux start time.
#
# Auxfile columns goFlux needs: UniqueID, start.time, Vtot, Area, Tcham, Pcham.
# Only START times are recorded; true start/end are set interactively in 03.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-basic-2026/scripts/00_setup.R"))

m <- read.csv(field_csv, stringsAsFactors = FALSE, encoding = "UTF-8")
m <- m[!grepl("no data", m$Note, ignore.case = TRUE), ]
message("Datasheet rows with data: ", nrow(m))

# --- Normalize the per-sheet date to a Date ----------------------------------
norm_date <- function(s) {
  s <- trimws(as.character(s))
  out <- suppressWarnings(lubridate::ymd(stringr::str_sub(s, 1, 10)))
  md  <- stringr::str_match(s, "(\\d{1,2})/(\\d{1,2})")
  need <- is.na(out) & !is.na(md[, 1])
  out[need] <- lubridate::make_date(2026, as.integer(md[need, 2]), as.integer(md[need, 3]))
  out
}
m$meas_date <- norm_date(m$date)

# --- Explode the 3 height positions, paired with their flux start times -------
extract_time <- function(v) stringr::str_extract(trimws(v), "\\d{1,2}:\\d{2}")
extract_h    <- function(v) suppressWarnings(as.numeric(stringr::str_extract(v, "\\d+\\.?\\d*")))
hcols <- c("H1", "H2", "H3"); fcols <- c("flux30Start", "flux80Start", "flux160Start")
defh  <- c(40, 80, 160)

recs <- lapply(1:3, function(i) {
  t <- extract_time(m[[fcols[i]]])
  h <- extract_h(m[[hcols[i]]]); h[is.na(h)] <- defh[i]
  keep <- !is.na(t)
  data.frame(TreeID = m$TreeID[keep], species = m$species[keep],
             meas_date = m$meas_date[keep], height_cm = h[keep],
             time_str = t[keep], stringsAsFactors = FALSE)
})
long <- do.call(rbind, recs)

long <- long %>%
  mutate(start.time = as.POSIXct(paste(format(meas_date), time_str),
                                 format = "%Y-%m-%d %H:%M", tz = tz_data),
         UniqueID   = paste0(TreeID, "_", height_cm, "cm")) %>%
  filter(!is.na(start.time)) %>%
  arrange(TreeID, start.time)

# Flag (and de-duplicate) any repeated UniqueID
dups <- long$UniqueID[duplicated(long$UniqueID)]
if (length(dups)) {
  warning("Duplicate UniqueIDs (kept first): ", paste(unique(dups), collapse = ", "))
  long <- long[!duplicated(long$UniqueID), ]
}

# --- Auxfile + field metadata ------------------------------------------------
aux <- long %>%
  transmute(UniqueID, start.time, Vtot = Vtot_L, Area = chamber_area_cm2,
            Tcham = Tcham_default_C, Pcham = Pcham_default_kPa)
field_meta <- long %>%
  transmute(UniqueID, TreeID, species, meas_date, time_str, height_cm, start.time)

message("\n=== Auxfile: ", nrow(aux), " measurements across ",
        length(unique(long$TreeID)), " trees ===")
print(as.data.frame(field_meta %>% count(meas_date, name = "measurements")), row.names = FALSE)

save(aux, field_meta, file = file.path(rdata_dir, "auxfile.RData"))
write.csv(aux, file.path(results_dir, "auxfile.csv"), row.names = FALSE)
write.csv(field_meta, file.path(results_dir, "field_metadata.csv"), row.names = FALSE)
message("Saved auxfile + field metadata. Proceed to 03_manual_id.R (interactive).")
