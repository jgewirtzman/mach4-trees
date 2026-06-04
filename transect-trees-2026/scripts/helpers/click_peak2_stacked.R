# =============================================================================
# helpers/click_peak2_stacked.R
# click.peak2.stacked(): a drop-in variant of goFlux::click.peak2() that shows
# TWO stacked panels -- CO2 (top) and CH4 (bottom) -- sharing one time axis, so
# you can judge each measurement window using both gas traces at once.
#
# You click the start and end on the TOP (CO2) panel; the same window is applied
# to both gases (goFlux is run per-gas afterwards using the shared `flag`).
# Uses split.screen() so the click coordinates map to the CO2 panel correctly.
#
# Output is identical in structure to click.peak2() (adds flag, Etime,
# start.time_corr, end.time_corr, obs.length_corr), so 04_flux_calculation.R
# works unchanged.
#
# Based on goFlux::click.peak2() (Qepanna/goFlux); plotting changed to stacked.
# =============================================================================

suppressMessages({ library(goFlux); library(dplyr) })

click.peak2.stacked <- function(ow.list,
                                 gastype   = "CO2dry_ppm",   # clicked (top panel)
                                 gastype2  = "CH4dry_ppb",   # shown (bottom panel)
                                 seq       = NULL,
                                 ref.secs  = NULL,           # typical-length guide (dashed), s
                                 sleep     = 3,
                                 plot.lim  = c(300, 5000),    # CO2 y clip
                                 plot.lim2 = c(1500, 100000), # CH4 y clip
                                 warn.length = 60,
                                 save.plots  = NULL,
                                 width = 7, height = 6,
                                 abline = TRUE, abline_corr = TRUE) {

  if (is.null(seq)) seq <- 1:length(ow.list)

  # y-limits like click.peak2: data range + 5% pad, clipped to plot.lim
  ylim_fun <- function(v, lim) {
    ymax <- max(v, na.rm = TRUE); ymin <- min(v, na.rm = TRUE)
    yd <- ymax - ymin; lo <- ymin - yd * 0.05; hi <- ymax + yd * 0.05
    c(ifelse(lo < lim[1], lim[1], lo), ifelse(hi > lim[2], lim[2], hi))
  }
  sleeploop <- function(x) { p <- proc.time(); Sys.sleep(x); proc.time() - p }

  # --- stacked drawing used for both the live view and the saved validation ---
  draw_stacked <- function(x1, x2, y1, y2, yl1, yl2, uid, top.title,
                           v_blue = NULL, v_red = NULL, v_guide = NULL,
                           time.axis = FALSE, col1 = "grey30", col2 = "grey30") {
    vlines <- function() {
      if (!is.null(v_blue))  abline(v = v_blue,  col = "blue")
      if (!is.null(v_guide)) abline(v = v_guide, col = "grey55", lty = 2)
      if (!is.null(v_red))   abline(v = v_red,   col = "red")
    }
    scr <- split.screen(c(2, 1))
    # TOP: CO2
    screen(scr[1]); par(mar = c(2.2, 4.6, 3.0, 1))
    plot(y1 ~ x1, col = col1, xlab = "", ylab = gastype, xaxt = "n", ylim = yl1)
    vlines()
    if (time.axis) axis.POSIXct(1, at = pretty(x1, 12), format = "%H:%M:%S") else axis(1)
    title(main = top.title, font.main = 2)
    # BOTTOM: CH4
    screen(scr[2]); par(mar = c(3.2, 4.6, 1.6, 1))
    plot(y2 ~ x2, col = col2, xlab = if (time.axis) "Time" else "Etime (s)",
         ylab = gastype2, xaxt = "n", ylim = yl2)
    vlines()
    if (time.axis) axis.POSIXct(1, at = pretty(x2, 12), format = "%H:%M:%S") else axis(1)
    scr
  }

  default.par <- par(no.readonly = TRUE)
  on.exit(par(default.par)); on.exit(Sys.unsetenv("TZ"), add = TRUE)
  ow.list.name <- deparse(substitute(ow.list))
  ow.corr.ls <- vector("list", length(ow.list))

  for (ow in seq) {
    g1 <- Reduce("c", ow.list[[ow]][, gastype])    # CO2
    g2 <- Reduce("c", ow.list[[ow]][, gastype2])   # CH4
    time.meas <- Reduce("c", ow.list[[ow]][, "POSIX.time"])
    uid <- unique(ow.list[[ow]]$UniqueID)
    start.time <- Reduce("c", unique(ow.list[[ow]][, "start.time"]))
    obs.length <- Reduce("c", unique(ow.list[[ow]][, "obs.length"]))
    end.time   <- start.time + obs.length
    yl1 <- ylim_fun(g1, plot.lim); yl2 <- ylim_fun(g2, plot.lim2)
    Sys.setenv(TZ = attr(time.meas, "tzone"))

    rownum <- NULL; identify.error <- NULL
    tryCatch({
      dev.new(noRStudioGD = TRUE, width = width, height = height)
      scr <- draw_stacked(time.meas, time.meas, g1, g2, yl1, yl2, uid,
                          paste0(uid, "   —  click START then END here (CO2)"),
                          v_blue  = if (abline) start.time else NULL,         # recorded start
                          v_guide = if (!is.null(ref.secs)) start.time + ref.secs else NULL,
                          time.axis = TRUE)
      screen(scr[1], new = FALSE)      # reactivate CO2 panel -> identify coords
      rownum <- identify(time.meas, g1, pos = FALSE, n = 2, plot = TRUE,
                         atpen = FALSE, offset = 0.5, tolerance = 0.25)
    }, error = function(e) identify.error <<- e)
    try(close.screen(all.screens = TRUE), silent = TRUE)
    if (!is.null(dev.list())) dev.off()

    if (!is.null(identify.error)) {
      warning(ow.list.name, "[[", ow, "]] ", uid, ": ", identify.error[[1]],
              call. = FALSE); next
    }
    if (length(rownum) < 2) rownum <- c(1, 1)

    start.time_corr <- time.meas[rownum[1]]
    end.time_corr   <- time.meas[rownum[2]]
    flux.flag <- which(dplyr::between(time.meas, start.time_corr, end.time_corr))
    d <- ow.list[[ow]] %>%
      mutate(flag = if_else(row_number() %in% flux.flag, 1, 0),
             Etime = as.numeric(POSIX.time - start.time_corr, units = "secs"),
             start.time_corr = start.time_corr,
             end.time_corr   = end.time_corr,
             obs.length_corr = as.numeric(end.time_corr - start.time_corr, units = "secs"))
    ow.corr.ls[[ow]] <- d

    n_flag <- sum(d$flag == 1)
    if (n_flag < warn.length)
      warning(ow.list.name, "[[", ow, "]] ", uid, ": only ", n_flag,
              " obs (< warn.length=", warn.length, ")", call. = FALSE)
    else
      message(ow.list.name, "[[", ow, "]] ", uid, ": good window (", n_flag, " obs)")

    # --- live validation flash (stacked, red = selected window) --------------
    et <- d$Etime; fl <- d$flag
    v_red <- c(0, max(et[fl == 1]))
    if (!is.null(sleep) && sleep > 0) {
      dev.new(noRStudioGD = TRUE, width = width, height = height)
      draw_stacked(et, et, g1, g2, yl1, yl2, uid,
                   paste0(uid, "   (red = selected window)"),
                   v_red = v_red, col1 = fl + 1, col2 = fl + 1)
      sleeploop(sleep)
      try(close.screen(all.screens = TRUE), silent = TRUE)
      if (!is.null(dev.list())) dev.off()
    }
  }

  # --- save all validation panels to one stacked PDF ---------------------------
  if (!is.null(save.plots)) {
    done <- which(!vapply(ow.corr.ls, is.null, logical(1)))
    if (length(done) > 0) {
      pdf(file = save.plots, width = width, height = height)
      for (ow in done) {
        d <- ow.corr.ls[[ow]]
        g1 <- d[[gastype]]; g2 <- d[[gastype2]]; et <- d$Etime; fl <- d$flag
        draw_stacked(et, et, g1, g2,
                     ylim_fun(g1, plot.lim), ylim_fun(g2, plot.lim2),
                     unique(d$UniqueID),
                     paste0(unique(d$UniqueID), "   (red = selected window)"),
                     v_red = c(0, max(et[fl == 1])), col1 = fl + 1, col2 = fl + 1)
        close.screen(all.screens = TRUE)
      }
      dev.off()
    }
  }

  ow.corr.ls <- ow.corr.ls[!vapply(ow.corr.ls, is.null, logical(1))]
  as.data.frame(dplyr::bind_rows(lapply(ow.corr.ls, as.data.frame)))
}
