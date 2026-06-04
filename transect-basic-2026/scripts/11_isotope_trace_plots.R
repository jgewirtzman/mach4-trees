# =============================================================================
# 11_isotope_trace_plots.R
# One multi-page PDF per dataset: for every flux that has an isotope sample,
# plot the CO2 (top) and CH4 (bottom) trace and mark t0/t1/t2 ON the trace line
# (time + concentration). t0 = pre-closure ambient; t1/t2 = the sampling times
# recorded in 10_isotope_concentrations.R. The clicked flux window (linear
# slope portion) is shaded for context.
#   -> transect-basic-2026/plots/isotope_traces_basic.pdf
#   -> transect-trees-2026/plots/isotope_traces_climbing.pdf
# Run with Rscript.
# =============================================================================

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
base <- "/Users/jongewirtzman/My Drive/Research/mach4-trees"
ok <- function(v) is.finite(v) & v > 0 & v < 1e6
T0_MIN <- -50/60   # representative position of the ambient (pre-closure) marker

# Some upstream CSVs serialized non-ASCII species under a non-UTF8 locale, so a
# "ú" became the literal text "<U+00FA>". Convert any such escape back.
unescape_u <- function(x) {
  vapply(x, function(s) {
    if (is.na(s) || !grepl("<U\\+[0-9A-Fa-f]{4,6}>", s)) return(s)
    m <- gregexpr("<U\\+([0-9A-Fa-f]{4,6})>", s)[[1]]
    for (cp in regmatches(s, gregexpr("<U\\+([0-9A-Fa-f]{4,6})>", s))[[1]]) {
      hex <- sub("<U\\+([0-9A-Fa-f]+)>", "\\1", cp)
      s <- sub(cp, intToUtf8(strtoi(hex, 16L)), s, fixed = TRUE)
    }
    enc2utf8(s)
  }, character(1), USE.NAMES = FALSE)
}

# Build the multi-page PDF for one dataset -----------------------------------
make_pdf <- function(manID, iso, out_pdf, label) {
  iso <- iso[iso$isotope %in% c(TRUE, "TRUE"), ]
  iso <- iso[order(iso$UniqueID), ]
  message(label, ": ", nrow(iso), " isotope measurements -> ", basename(out_pdf))
  grDevices::cairo_pdf(out_pdf, width = 7.5, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off())

  for (k in seq_len(nrow(iso))) {
    r   <- iso[k, ]
    uid <- r$UniqueID
    tr  <- manID[manID$UniqueID == uid, ]
    if (!nrow(tr)) next
    win_end_min <- suppressWarnings(unique(tr$obs.length_corr)[1] / 60)

    # long-format trace (CO2 + CH4) ------------------------------------------
    trL <- bind_rows(
      data.frame(emin = tr$Etime/60, value = tr$CO2dry_ppm, gas = "CO2 (ppm)")[ok(tr$CO2dry_ppm), ],
      data.frame(emin = tr$Etime/60, value = tr$CH4dry_ppb, gas = "CH4 (ppb)")[ok(tr$CH4dry_ppb), ])

    # markers: one row per (gas, stage) with a time + concentration ----------
    # use the manually-adjusted t0 time when present, else the default position
    t0_min <- if (!is.null(r$iso_t0_min) && is.finite(r$iso_t0_min)) r$iso_t0_min else T0_MIN
    mk <- bind_rows(
      data.frame(gas = "CO2 (ppm)", stage = c("t0","t1","t2"),
                 emin  = c(t0_min, r$iso_t1_min, r$iso_t2_min),
                 value = c(r$iso_t0_CO2_ppm, r$iso_t1_CO2_ppm, r$iso_t2_CO2_ppm)),
      data.frame(gas = "CH4 (ppb)", stage = c("t0","t1","t2"),
                 emin  = c(t0_min, r$iso_t1_min, r$iso_t2_min),
                 value = c(r$iso_t0_CH4_ppb, r$iso_t1_CH4_ppb, r$iso_t2_CH4_ppb)))
    mk <- mk[is.finite(mk$emin) & is.finite(mk$value), ]
    unit <- ifelse(mk$gas == "CO2 (ppm)", "ppm", "ppb")
    mk$lab <- sprintf("%s: %s %s", mk$stage,
                      formatC(round(mk$value), format = "d", big.mark = ","), unit)

    cols <- c(t0 = "#1f78b4", t1 = "#e31a1c", t2 = "#ff7f00")
    ttl  <- sprintf("%s   |   %s   %s cm", uid, unescape_u(r$species), r$height_cm)
    sub  <- if (!is.na(r$iso_flag) && r$iso_flag != "")
              paste0("flag: ", r$iso_flag) else
              sprintf("t1 = %.1f min,  t2 = %s min into flux",
                      r$iso_t1_min, ifelse(is.na(r$iso_t2_min),"—",sprintf("%.1f",r$iso_t2_min)))

    p <- ggplot(trL, aes(emin, value)) +
      annotate("rect", xmin = 0, xmax = ifelse(is.finite(win_end_min), win_end_min, 0),
               ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.45) +
      geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
      geom_line(colour = "grey35", linewidth = 0.5) +
      geom_vline(data = mk[mk$stage != "t0", ], aes(xintercept = emin, colour = stage),
                 linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
      geom_point(data = mk, aes(colour = stage), size = 3) +
      ggrepel::geom_text_repel(data = mk, aes(label = lab, colour = stage),
                               size = 3, show.legend = FALSE, seg.color = NA,
                               direction = "y", nudge_y = 0) +
      facet_wrap(~ gas, ncol = 1, scales = "free_y") +
      scale_colour_manual(values = cols, name = NULL) +
      labs(x = "time relative to chamber closure (min)", y = NULL,
           title = ttl, subtitle = sub) +
      theme_classic(base_size = 11) +
      theme(strip.background = element_blank(),
            strip.text = element_text(face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 8, colour = "grey30"),
            legend.position = "top")
    print(p)
  }
}

# BASIC ----------------------------------------------------------------------
load(file.path(base, "transect-basic-2026/RData/manID.RData")); mB <- manID
fb <- read.csv(file.path(base, "transect-basic-2026/results/transect_basic_2026_fluxes_isotopes.csv"),
               stringsAsFactors = FALSE, encoding = "UTF-8")
dir.create(file.path(base, "transect-basic-2026/plots"), showWarnings = FALSE)
make_pdf(mB, fb, file.path(base, "transect-basic-2026/plots/isotope_traces_basic.pdf"), "BASIC")

# CLIMBING -------------------------------------------------------------------
load(file.path(base, "transect-trees-2026/RData/manID.RData")); mC <- manID
fc <- read.csv(file.path(base, "transect-trees-2026/results/transect_trees_2026_fluxes_isotopes.csv"),
               stringsAsFactors = FALSE, encoding = "UTF-8")
dir.create(file.path(base, "transect-trees-2026/plots"), showWarnings = FALSE)
make_pdf(mC, fc, file.path(base, "transect-trees-2026/plots/isotope_traces_climbing.pdf"), "CLIMBING")

message("Done.")
