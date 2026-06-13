# =============================================================================
# 10b_add_basic_isotopes.R   (incremental, non-destructive)
# The May Isotope Sampling Sheet (Trees-Actual) lists basic-tree isotope vials
# that were not originally flagged in 10_. These trees DO have processed basic
# flux traces; here we flag the new (tree, height) points as isotope samples and
# extract their t0/t1 concentrations (basic protocol: t1 = 5 min, no t2), WITHOUT
# recomputing or overwriting any row already reviewed via 12_ (iso_review kept).
# New rows get iso_review = "" so 12_isotope_click_adjust.R (DATASET="basic")
# picks them up for clicking. Run with Rscript.
#
# New points (t1 only, per the sheet):
#   NE59 H2/80 ; N88 H1/40,H2/80 ; NW16 H1/40,H2/80 ; S10 H1/40,H2/80 ; SW22 H1/40,H2/80
# =============================================================================

suppressMessages({ library(dplyr) })
base <- "/Users/jongewirtzman/My Drive/Research/mach4-trees"
ok  <- function(v) is.finite(v) & v > 0 & v < 1e6
amb <- function(mid, uid, gas) {
  m <- mid[mid$UniqueID == uid & mid$Etime > -90 & mid$Etime < -15 & ok(mid[[gas]]), ]
  if (nrow(m) < 3) return(NA_real_); mean(m[[gas]]) }
at_etime <- function(mid, uid, et_s, gas) {
  m <- mid[mid$UniqueID == uid & ok(mid[[gas]]), ]; if (!nrow(m)) return(NA_real_)
  i <- which.min(abs(m$Etime - et_s)); if (abs(m$Etime[i] - et_s) > 20) NA_real_ else m[[gas]][i] }
hbin <- function(h) cut(h, c(0,55,120,300,600,1e4), labels = c("40","80","160","5m","top"))

load(file.path(base, "transect-basic-2026/RData/manID.RData")); mb <- manID
csv <- file.path(base, "transect-basic-2026/results/transect_basic_2026_fluxes_isotopes.csv")
fb  <- read.csv(csv, stringsAsFactors = FALSE, encoding = "UTF-8")
fb$bin <- as.character(hbin(fb$height_cm))

new <- tibble::tribble(
  ~TreeID, ~h_named,
  "NE59",80,
  "N88",40, "N88",80,
  "NW16",40, "NW16",80,
  "S10",40,  "S10",80,
  "SW22",40, "SW22",80) %>%
  mutate(bin = as.character(hbin(h_named)))

added <- character(0)
for (k in seq_len(nrow(new))) {
  sel <- which(toupper(fb$TreeID) == toupper(new$TreeID[k]) & fb$bin == new$bin[k])
  if (!length(sel)) { message("  no basic flux row for ", new$TreeID[k], " bin ", new$bin[k]); next }
  i <- sel[1]
  if (isTRUE(fb$isotope[i]) && !is.na(fb$iso_review[i]) && fb$iso_review[i] != "") {
    message("  already reviewed, skipping: ", fb$UniqueID[i]); next }
  uid <- fb$UniqueID[i]
  fb$isotope[i]        <- TRUE
  fb$iso_t1_min[i]     <- 5;  fb$iso_t2_min[i] <- NA
  fb$iso_t0_min[i]     <- NA
  fb$iso_t0_CO2_ppm[i] <- amb(mb, uid, "CO2dry_ppm");        fb$iso_t0_CH4_ppb[i] <- amb(mb, uid, "CH4dry_ppb")
  fb$iso_t1_CO2_ppm[i] <- at_etime(mb, uid, 5*60, "CO2dry_ppm"); fb$iso_t1_CH4_ppb[i] <- at_etime(mb, uid, 5*60, "CH4dry_ppb")
  fb$iso_t2_CO2_ppm[i] <- NA; fb$iso_t2_CH4_ppb[i] <- NA
  fb$iso_review[i]     <- ""                       # -> 12_ will offer it for clicking
  fb$iso_flag[i] <- if (is.finite(fb$iso_t1_CO2_ppm[i]) && is.finite(fb$iso_t0_CO2_ppm[i]) &&
                        fb$iso_t1_CO2_ppm[i] <= fb$iso_t0_CO2_ppm[i])
                      "t1 not above t0 (short window / low flux?)" else ""
  added <- c(added, uid)
}
fb$bin <- NULL
write.csv(fb, csv, row.names = FALSE, fileEncoding = "UTF-8")

cat("Added", length(added), "new basic isotope measurements:\n  ", paste(added, collapse = ", "), "\n")
cat("Total basic isotope measurements now:", sum(fb$isotope),
    " (reviewed:", sum(fb$isotope & fb$iso_review != ""), "/ to click:", sum(fb$isotope & fb$iso_review == ""), ")\n")
print(fb %>% filter(UniqueID %in% added) %>%
      transmute(UniqueID, CO2_t0 = round(iso_t0_CO2_ppm), CO2_t1 = round(iso_t1_CO2_ppm),
                CH4_t0 = round(iso_t0_CH4_ppb), CH4_t1 = round(iso_t1_CH4_ppb), iso_flag) %>%
      as.data.frame(), row.names = FALSE)
