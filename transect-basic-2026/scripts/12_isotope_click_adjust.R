# =============================================================================
# 12_isotope_click_adjust.R   *** RUN IN AN INTERACTIVE R SESSION ***
# Manually correct isotope sampling times by clicking on the flux trace, the
# same way click.peak2 works. For each isotope measurement the stacked CO2 (top)
# + CH4 (bottom) trace pops up with the CURRENT t0/t1/t2 markers drawn. You then
# click the CORRECT sampling points, LEFT -> RIGHT, in time order:
#       t0 (ambient, before closure),  t1,  [t2 if present]
# It is fully click-driven (no keyboard): after your Nth click each point snaps
# to the nearest trace sample, the markers flash green ("SAVED"), the row is
# written, and it AUTO-ADVANCES to the next trace. N is shown in the title.
# The *_fluxes_isotopes.csv files are rewritten after every trace.
#
# Per-trace controls (all in the popup window):
#   - click N points to set the new times -> saves & advances automatically
#   - press Esc (no clicks) -> keep current values, advance
#
# Resume: rows already reviewed (iso_review != "") are skipped. To re-do a few,
# set REDO_IDS <- c("T-W10_40cm", ...); to redo everything, REDO_ALL <- TRUE.
#
# Choose the dataset below, run the whole file, then work through the pop-ups.
# =============================================================================

DATASET  <- "climbing"     # "climbing" or "basic"
REDO_ALL <- FALSE          # TRUE = review every isotope trace again
REDO_IDS <- c("T-E20_160cm")   # specific UniqueIDs to force-review (reset to character(0) when done)

if (!interactive()) stop("Run this in an interactive R session (RStudio / R.app).")
suppressMessages({ library(dplyr) })
base <- "/Users/jongewirtzman/My Drive/Research/mach4-trees"
ok <- function(v) is.finite(v) & v > 0 & v < 1e6
unescape_u <- function(s) {
  if (length(s) != 1 || is.na(s) || !grepl("<U\\+[0-9A-Fa-f]{4,6}>", s)) return(s)
  for (cp in regmatches(s, gregexpr("<U\\+[0-9A-Fa-f]{4,6}>", s))[[1]]) {
    hex <- sub("<U\\+([0-9A-Fa-f]+)>", "\\1", cp)
    s <- sub(cp, intToUtf8(strtoi(hex, 16L)), s, fixed = TRUE)
  }
  enc2utf8(s)
}

cfg <- list(
  climbing = list(
    manID = file.path(base, "transect-trees-2026/RData/manID.RData"),
    csv   = file.path(base, "transect-trees-2026/results/transect_trees_2026_fluxes_isotopes.csv")),
  basic = list(
    manID = file.path(base, "transect-basic-2026/RData/manID.RData"),
    csv   = file.path(base, "transect-basic-2026/results/transect_basic_2026_fluxes_isotopes.csv"))
)[[DATASET]]

load(cfg$manID)                       # manID
fx <- read.csv(cfg$csv, stringsAsFactors = FALSE, encoding = "UTF-8")
# ensure adjustable columns exist
for (cc in c("iso_t0_min","iso_review")) if (!cc %in% names(fx))
  fx[[cc]] <- if (cc == "iso_review") "" else NA_real_
fx$iso_review[is.na(fx$iso_review)] <- ""

# --- dedicated interactive window (NOT the RStudio plot pane) ----------------
# locator() is unreliable in the RStudio Plots pane; open a native window so
# clicks register. Reused across traces; reopened if the user closes it.
ISO_DEV <- NULL
ensure_dev <- function() {
  if (is.null(ISO_DEV) || !(ISO_DEV %in% dev.list())) {
    dev.new(noRStudioGD = TRUE, width = 8, height = 7.2)
    ISO_DEV <<- dev.cur()
  } else dev.set(ISO_DEV)
}

# --- nearest valid trace sample to a time (minutes) -> time + both gases ------
snap <- function(tr, emin) {
  v  <- tr[ok(tr$CO2dry_ppm) & ok(tr$CH4dry_ppb), ]
  if (!nrow(v)) return(NULL)
  i  <- which.min(abs(v$Etime/60 - emin))
  list(min = round(v$Etime[i]/60, 2), CO2 = v$CO2dry_ppm[i], CH4 = v$CH4dry_ppb[i])
}

# --- draw the stacked trace with a set of markers ----------------------------
draw_trace <- function(tr, row, marks, col_pts, title_extra = "") {
  win_end <- suppressWarnings(unique(tr$obs.length_corr)[1] / 60)
  v  <- tr[ok(tr$CO2dry_ppm) & ok(tr$CH4dry_ppb), ]
  xr <- range(v$Etime/60)
  op <- par(mfrow = c(2, 1), mar = c(2.4, 4.2, 2.6, 1), oma = c(2.5, 0, 2.2, 0))
  on.exit(par(op), add = TRUE)
  panel <- function(yv, ylab, mk_y, marker_vals) {
    plot(v$Etime/60, v[[yv]], type = "l", col = "grey30", lwd = 1.3,
         xlim = xr, xlab = "", ylab = ylab)
    if (is.finite(win_end)) rect(0, par("usr")[3], win_end, par("usr")[4],
                                 col = adjustcolor("grey80", 0.35), border = NA)
    lines(v$Etime/60, v[[yv]], col = "grey30", lwd = 1.3)
    abline(v = 0, col = "grey60")
    for (j in seq_len(nrow(marks))) {
      if (!is.finite(marks$min[j]) || !is.finite(marker_vals[j])) next
      abline(v = marks$min[j], lty = 2, col = col_pts[j])
      points(marks$min[j], marker_vals[j], pch = 19, col = col_pts[j], cex = 1.5)
      text(marks$min[j], marker_vals[j], labels = marks$stage[j],
           pos = 3, col = col_pts[j], font = 2, cex = 0.9)
    }
  }
  panel("CO2dry_ppm", "CO2 (ppm)", , marks$CO2)
  panel("CH4dry_ppb", "CH4 (ppb)", , marks$CH4)
  mtext("time relative to chamber closure (min)", side = 1, outer = TRUE, line = 1)
  mtext(sprintf("%s   |   %s   %s cm   %s", row$UniqueID,
                unescape_u(row$species), row$height_cm, title_extra),
        side = 3, outer = TRUE, line = 0.6, font = 2, cex = 1.05)
}

# --- which isotope rows to review --------------------------------------------
idx <- which(fx$isotope %in% c(TRUE, "TRUE"))
todo <- if (REDO_ALL) idx else idx[fx$iso_review[idx] == "" | fx$UniqueID[idx] %in% REDO_IDS]
message(DATASET, ": ", length(todo), " isotope traces to review (",
        length(idx) - length(todo), " already reviewed).")

stage_cols <- c(t0 = "#1f78b4", t1 = "#e31a1c", t2 = "#ff7f00")

for (ii in todo) {
  row <- fx[ii, ]
  tr  <- manID[manID$UniqueID == row$UniqueID, ]
  if (!nrow(tr)) { message("  no trace for ", row$UniqueID, " - skipped"); next }

  has_t2 <- is.finite(row$iso_t2_min) || is.finite(row$iso_t2_CO2_ppm)
  stages <- if (has_t2) c("t0","t1","t2") else c("t0","t1")
  N <- length(stages)

  # current markers (t0 shown near pre-closure if no stored time)
  cur <- data.frame(
    stage = c("t0","t1","t2"),
    min   = c(ifelse(is.finite(row$iso_t0_min), row$iso_t0_min, -50/60), row$iso_t1_min, row$iso_t2_min),
    CO2   = c(row$iso_t0_CO2_ppm, row$iso_t1_CO2_ppm, row$iso_t2_CO2_ppm),
    CH4   = c(row$iso_t0_CH4_ppb, row$iso_t1_CH4_ppb, row$iso_t2_CH4_ppb),
    stringsAsFactors = FALSE)
  cur <- cur[cur$stage %in% stages, ]

  # Fully click-driven: click N points -> snap, flash result, save, auto-advance.
  # No keyboard confirmation (the popup steals focus, so readline would hang).
  ensure_dev()                                           # draw in the popup window
  draw_trace(tr, row, cur, stage_cols[cur$stage],
             sprintf("[click %d pts L->R: %s   (Esc = keep current)]", N, paste(stages, collapse=", ")))
  loc <- try(locator(n = N, type = "p", pch = 3, col = "forestgreen", lwd = 2), silent = TRUE)
  nclick <- if (inherits(loc, "try-error") || is.null(loc)) 0 else length(loc$x)

  if (nclick == 0) {                                     # Esc -> keep current, advance
    fx$iso_review[ii] <- "kept"
    message(sprintf("  %-16s kept current", row$UniqueID))
  } else if (nclick != N) {                              # partial -> leave for a redo pass
    message(sprintf("  %-16s got %d clicks, expected %d -> kept current, revisit via REDO_IDS",
                    row$UniqueID, nclick, N))
  } else {                                               # N clicks -> snap + save + advance
    o   <- order(loc$x)                                  # chronological
    new <- cur
    for (k in seq_len(N)) {
      sp <- snap(tr, loc$x[o][k])
      if (is.null(sp)) next
      new$min[k] <- sp$min; new$CO2[k] <- sp$CO2; new$CH4[k] <- sp$CH4
    }
    draw_trace(tr, row, new, rep("forestgreen", N), "SAVED")    # flash the result
    Sys.sleep(0.45)
    g <- function(st, col) new[new$stage == st, col][1]
    fx$iso_t0_min[ii]     <- g("t0","min")
    fx$iso_t0_CO2_ppm[ii] <- g("t0","CO2"); fx$iso_t0_CH4_ppb[ii] <- g("t0","CH4")
    fx$iso_t1_min[ii]     <- g("t1","min")
    fx$iso_t1_CO2_ppm[ii] <- g("t1","CO2"); fx$iso_t1_CH4_ppb[ii] <- g("t1","CH4")
    if (has_t2) {
      fx$iso_t2_min[ii]     <- g("t2","min")
      fx$iso_t2_CO2_ppm[ii] <- g("t2","CO2"); fx$iso_t2_CH4_ppb[ii] <- g("t2","CH4")
    }
    fx$iso_review[ii] <- "manual"
    fx$iso_flag[ii]   <- trimws(paste("manually adjusted",
                                      ifelse(is.na(row$iso_flag), "", row$iso_flag)))
    message(sprintf("  %-16s saved: t0=%.2f t1=%.2f%s min", row$UniqueID,
                    g("t0","min"), g("t1","min"),
                    if (has_t2) sprintf(" t2=%.2f", g("t2","min")) else ""))
  }
  write.csv(fx, cfg$csv, row.names = FALSE, fileEncoding = "UTF-8")   # checkpoint
}

if (!is.null(ISO_DEV) && ISO_DEV %in% dev.list()) dev.off(ISO_DEV)   # close popup
message("\nDone. Reviewed: ", sum(fx$iso_review != "" & fx$isotope %in% c(TRUE,"TRUE")),
        " / ", length(idx), "  (saved -> ", basename(cfg$csv), ")")
message("Re-run 11_isotope_trace_plots.R to regenerate the PDF with the adjusted times.")
