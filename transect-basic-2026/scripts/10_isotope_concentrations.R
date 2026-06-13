# =============================================================================
# 10_isotope_concentrations.R
# For every flux with an isotope sample, record the CO2/CH4 concentration at
# t0 (pre-closure ambient), t1, and t2. Updates BOTH datasets:
#   BASIC    (transect-basic-2026): t1 = 5 min into the flux; t2 = 5 + stated;
#            samples only at the named heights.
#   CLIMBING (transect-trees-2026): t1/t2 are the exact clock times in
#            data/isotope_times.csv (all climbed trees, all heights).
# t0 is always the ambient just before the chamber closed.
# Run with Rscript. Writes *_fluxes_isotopes.csv in each project's results/.
# =============================================================================

suppressMessages({ library(dplyr); library(tidyr); library(stringr); library(lubridate) })
base <- "/Users/jongewirtzman/My Drive/Research/mach4-trees"

# --- concentration helpers (drop LI-7810 diagnostic/garbage rows: 0<gas<1e6) -
ok <- function(v) is.finite(v) & v > 0 & v < 1e6
amb <- function(mid, uid, gas) {                       # pre-closure ambient
  m <- mid[mid$UniqueID == uid & mid$Etime > -90 & mid$Etime < -15 & ok(mid[[gas]]), ]
  if (nrow(m) < 3) return(NA_real_); mean(m[[gas]])
}
at_etime <- function(mid, uid, et_s, gas) {            # nearest Etime (basic)
  if (is.na(et_s)) return(NA_real_)
  m <- mid[mid$UniqueID == uid & ok(mid[[gas]]), ]; if (!nrow(m)) return(NA_real_)
  i <- which.min(abs(m$Etime - et_s))
  if (abs(m$Etime[i] - et_s) > 20) NA_real_ else m[[gas]][i]
}
at_time <- function(licor, t0, gas) {                  # nearest POSIX time (climbing)
  if (is.na(t0)) return(NA_real_)
  m <- licor[ok(licor[[gas]]), ]
  i <- which.min(abs(as.numeric(m$POSIX.time) - as.numeric(t0)))
  if (abs(as.numeric(m$POSIX.time[i]) - as.numeric(t0)) > 20) NA_real_ else m[[gas]][i]
}
hbin <- function(h) cut(h, c(0,55,120,300,600,1e4), labels = c("40","80","160","5m","top"))

# =============================================================================
# BASIC trees
# =============================================================================
load(file.path(base, "transect-basic-2026/RData/manID.RData"))            # manID
mb <- manID
fb <- read.csv(file.path(base, "transect-basic-2026/results/transect_basic_2026_fluxes_with_mdf.csv"),
               stringsAsFactors = FALSE, encoding = "UTF-8")
# Basic isotope samples (named heights, t1 = 5 min, t2 = 5 + stated). The last
# block (NE59 H2 .. SW22) was reconciled from the May Isotope Sampling Sheet
# (Trees-Actual); see 10b_add_basic_isotopes.R, which adds them incrementally so
# manual reviews from 12_ are preserved. NOTE: a from-scratch run of THIS script
# regenerates defaults for ALL rows and overwrites manual (clicked) reviews —
# re-review in 12_ afterward, or use 10b_ to add new samples without clobbering.
biso <- tibble::tribble(
  ~TreeID, ~h_named, ~t1_min, ~t2_min,
  "S14b",40,5,NA, "S14b",80,5,NA, "NW46a",40,5,NA, "NW46a",80,5,NA,
  "E20",40,5,NA,  "E20",80,5,NA,  "W23a",40,5,NA,  "W23a",80,5,NA,
  "W10",40,5,NA,  "W10",80,5,NA,  "NE14",40,5,7,   "NE59",40,5,7,
  "SW61",40,5,NA, "SW61",80,5,NA, "SW61",160,5,NA,
  "NE59",80,5,NA, "N88",40,5,NA,  "N88",80,5,NA,   "NW16",40,5,NA,
  "NW16",80,5,NA, "S10",40,5,NA,  "S10",80,5,NA,   "SW22",40,5,NA, "SW22",80,5,NA) %>%
  mutate(bin = as.character(hbin(h_named)))
fb$bin <- as.character(hbin(fb$height_cm))
fb <- fb %>% left_join(biso, by = c("TreeID","bin"))
fb$isotope <- !is.na(fb$t1_min)
fb[c("iso_t0_CO2_ppm","iso_t0_CH4_ppb","iso_t1_CO2_ppm","iso_t1_CH4_ppb",
     "iso_t2_CO2_ppm","iso_t2_CH4_ppb")] <- NA_real_
for (i in which(fb$isotope)) {
  uid <- fb$UniqueID[i]
  fb$iso_t0_CO2_ppm[i] <- amb(mb, uid, "CO2dry_ppm"); fb$iso_t0_CH4_ppb[i] <- amb(mb, uid, "CH4dry_ppb")
  fb$iso_t1_CO2_ppm[i] <- at_etime(mb, uid, fb$t1_min[i]*60, "CO2dry_ppm")
  fb$iso_t1_CH4_ppb[i] <- at_etime(mb, uid, fb$t1_min[i]*60, "CH4dry_ppb")
  fb$iso_t2_CO2_ppm[i] <- at_etime(mb, uid, fb$t2_min[i]*60, "CO2dry_ppm")
  fb$iso_t2_CH4_ppb[i] <- at_etime(mb, uid, fb$t2_min[i]*60, "CH4dry_ppb")
}
fb$iso_flag <- ifelse(fb$isotope & is.finite(fb$iso_t1_CO2_ppm) & is.finite(fb$iso_t0_CO2_ppm) &
                        fb$iso_t1_CO2_ppm <= fb$iso_t0_CO2_ppm,
                      "t1 not above t0 (short window / low flux?)", NA_character_)
fb <- fb %>% rename(iso_t1_min = t1_min, iso_t2_min = t2_min) %>% select(-h_named, -bin)
write.csv(fb, file.path(base, "transect-basic-2026/results/transect_basic_2026_fluxes_isotopes.csv"), row.names = FALSE)
cat("BASIC: isotope measurements =", sum(fb$isotope), "/", nrow(fb), "\n")

# =============================================================================
# CLIMBING trees
# =============================================================================
load(file.path(base, "transect-trees-2026/RData/manID.RData"));  mc <- manID
load(file.path(base, "transect-trees-2026/RData/licor_imported.RData"));  Lc <- licor
fc <- read.csv(file.path(base, "transect-trees-2026/results/transect_trees_2026_fluxes.csv"),
               stringsAsFactors = FALSE, encoding = "UTF-8")
it <- read.csv(file.path(base, "transect-trees-2026/data/isotope_times.csv"),
               stringsAsFactors = FALSE, check.names = TRUE)
names(it)[1] <- "tree"
itl <- it %>%
  pivot_longer(-c(tree,date), names_to = c(".value","bin"),
               names_pattern = "(t1|t2)\\.(.*)") %>%
  mutate(bin = recode(bin, X5m = "5m"),
         tree = recode(tree, SW46a = "NW46a"),
         d = lubridate::mdy(date),
         t1_clk = suppressWarnings(as.POSIXct(paste(d, str_extract(t1, "\\d{1,2}:\\d{2}")), format="%Y-%m-%d %H:%M", tz="UTC")),
         t2_clk = suppressWarnings(as.POSIXct(paste(d, str_extract(t2, "\\d{1,2}:\\d{2}")), format="%Y-%m-%d %H:%M", tz="UTC")))
fc$tree <- sub("^T-", "", fc$TreeID); fc$bin <- as.character(hbin(fc$height_cm))
fc <- fc %>% left_join(itl %>% select(tree, bin, t1_clk, t2_clk), by = c("tree","bin"))
fc$isotope <- !is.na(fc$t1_clk) | !is.na(fc$t2_clk)
fc[c("iso_t0_CO2_ppm","iso_t0_CH4_ppb","iso_t1_CO2_ppm","iso_t1_CH4_ppb",
     "iso_t2_CO2_ppm","iso_t2_CH4_ppb","iso_t1_min","iso_t2_min")] <- NA_real_
fc$iso_flag <- NA_character_
# Climbing protocol: t1 = 3 min, t2 = 5 min into the flux. The clock times in
# isotope_times.csv are the field cross-check; when a clock falls OUTSIDE the
# flux window (corrupt/mis-recorded entry) we fall back to the protocol elapsed
# time measured from this flux's own start, and flag it.
PROT_T1_MIN <- 3; PROT_T2_MIN <- 5
for (i in which(fc$isotope)) {
  uid <- fc$UniqueID[i]
  st <- unique(mc$start.time_corr[mc$UniqueID == uid])[1]
  en <- unique(mc$end.time_corr[mc$UniqueID == uid])[1]
  fc$iso_t0_CO2_ppm[i] <- amb(mc, uid, "CO2dry_ppm"); fc$iso_t0_CH4_ppb[i] <- amb(mc, uid, "CH4dry_ppb")
  # Valid sampling clock: at/after chamber closure (start) and within a generous
  # 12-min ceiling. The clicked window (end.time_corr) is only the linear slope
  # portion; the chamber stays closed and logging continues past it, so a 5-6 min
  # isotope sample legitimately lands after `en`. Only clocks BEFORE closure (or
  # absurdly late) are the genuine mis-records.
  inwin <- function(t) !is.na(t) && !is.na(st) &&
    as.numeric(t) >= as.numeric(st) - 30 && as.numeric(t) <= as.numeric(st) + 12*60
  flags <- character(0)
  # t1
  if (inwin(fc$t1_clk[i])) { tt1 <- fc$t1_clk[i]
  } else { tt1 <- if (!is.na(st)) st + PROT_T1_MIN*60 else NA
    flags <- c(flags, "t1 clock outside flux window -> used 3min protocol") }
  # t2
  if (inwin(fc$t2_clk[i])) { tt2 <- fc$t2_clk[i]
  } else { tt2 <- if (!is.na(st)) st + PROT_T2_MIN*60 else NA
    flags <- c(flags, if (is.na(fc$t2_clk[i])) "t2 clock missing -> used 5min protocol"
                      else "t2 clock outside flux window -> used 5min protocol") }
  fc$iso_t1_CO2_ppm[i] <- at_time(Lc, tt1, "CO2dry_ppm"); fc$iso_t1_CH4_ppb[i] <- at_time(Lc, tt1, "CH4dry_ppb")
  fc$iso_t2_CO2_ppm[i] <- at_time(Lc, tt2, "CO2dry_ppm"); fc$iso_t2_CH4_ppb[i] <- at_time(Lc, tt2, "CH4dry_ppb")
  fc$iso_t1_min[i] <- if (!is.na(st) && !is.na(tt1)) round(as.numeric(difftime(tt1, st, units="mins")),1) else NA
  fc$iso_t2_min[i] <- if (!is.na(st) && !is.na(tt2)) round(as.numeric(difftime(tt2, st, units="mins")),1) else NA
  if (!is.na(fc$iso_t2_min[i]) && !is.na(fc$iso_t1_min[i]) && fc$iso_t2_min[i] < fc$iso_t1_min[i])
    flags <- c(flags, "t2<t1")
  fc$iso_flag[i] <- paste(flags, collapse = "; ")
}
fc <- fc %>% mutate(iso_t1_clock = format(t1_clk, "%H:%M"), iso_t2_clock = format(t2_clk, "%H:%M")) %>% select(-tree, -bin, -t1_clk, -t2_clk)
write.csv(fc, file.path(base, "transect-trees-2026/results/transect_trees_2026_fluxes_isotopes.csv"), row.names = FALSE)
cat("CLIMBING: isotope measurements =", sum(fc$isotope), "/", nrow(fc), "\n")

# --- quick sanity: concentrations should rise t0 -> t1 -> t2 -----------------
cat("\n=== BASIC isotope samples (CO2 ppm: t0->t1->t2) ===\n")
print(fb %>% filter(isotope) %>% transmute(UniqueID, t1=iso_t1_min, t2=iso_t2_min,
      CO2_t0=round(iso_t0_CO2_ppm), CO2_t1=round(iso_t1_CO2_ppm), CO2_t2=round(iso_t2_CO2_ppm),
      CH4_t0=round(iso_t0_CH4_ppb), CH4_t1=round(iso_t1_CH4_ppb)) %>% as.data.frame(), row.names=FALSE)
cat("\n=== CLIMBING isotope samples (sample of 12) ===\n")
print(fc %>% filter(isotope) %>% transmute(UniqueID, t1c=iso_t1_clock, t2c=iso_t2_clock, t1=iso_t1_min, t2=iso_t2_min,
      CO2_t0=round(iso_t0_CO2_ppm), CO2_t1=round(iso_t1_CO2_ppm), CO2_t2=round(iso_t2_CO2_ppm), flag=iso_flag) %>%
      head(12) %>% as.data.frame(), row.names=FALSE)
