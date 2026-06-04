# =============================================================================
# 08_raincloud.R
# Distribution of fluxes by height bin as a raincloud: ridge (half-density) +
# box + jittered points, stacked with LOW heights at the BOTTOM and HIGH heights
# at the TOP. Points with |flux| >= MDF are SOLID, < MDF are OPEN, where MDF is
# the per-measurement Allan-deviation (Christiansen) definition.
# CH4 uses an asinh x-scale so the low-flux end is readable.
# Can be run with Rscript (non-interactive).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))
suppressMessages({ library(ggplot2); library(ggdist); library(scales); library(dplyr) })
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

# --- MDF: per-measurement Allan deviation (Christiansen). 95% conf. -----------
mdf_level <- "chr95"   # switch to "chr99" / "chr90" if desired

dat <- read.csv(file.path(results_dir, "transect_trees_2026_fluxes_with_mdf.csv"),
                stringsAsFactors = FALSE, encoding = "UTF-8")

# --- Height bins, ordered base (bottom) -> canopy (top) ----------------------
hb_levels <- c("Base (40–50 cm)", "80 cm", "160 cm",
               "Mid (400–500 cm)", "Canopy (800–1000 cm)")
dat$height_bin <- factor(case_when(
  dat$height_cm <= 50                          ~ hb_levels[1],
  dat$height_cm == 80                          ~ hb_levels[2],
  dat$height_cm == 160                         ~ hb_levels[3],
  dat$height_cm >= 400 & dat$height_cm <= 500  ~ hb_levels[4],
  dat$height_cm >= 800                         ~ hb_levels[5]
), levels = hb_levels)   # first level -> bottom of y axis

# --- Detection vs the per-measurement Allan MDF ------------------------------
det_lab <- c("≥ MDF", "< MDF")               # solid, open
detect <- function(flux, mdf) factor(ifelse(abs(flux) >= mdf, det_lab[1], det_lab[2]),
                                     levels = det_lab)
dat$CH4_detect <- detect(dat$CH4_best.flux, dat[[paste0("CH4_MDF_", mdf_level)]])
dat$CO2_detect <- detect(dat$CO2_best.flux, dat[[paste0("CO2_MDF_", mdf_level)]])

shape_vals <- c(19, 1);            names(shape_vals) <- det_lab
col_vals   <- c("black", "grey55"); names(col_vals)  <- det_lab

theme_rc <- theme_classic(base_size = 12) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        plot.subtitle = element_text(size = 9, colour = "grey30"))

sub <- paste0("solid = |flux| ≥ MDF, open = < MDF   (MDF = per-measurement ",
              "Allan deviation, ", sub("chr", "", mdf_level), "% conf.)")

# asinh axis that shows the low-flux end while compressing the high base values.
# Default: compact K/M labels. Pass labs = label_pow10 for power-of-ten ticks.
compact_lab <- scales::label_number(scale_cut = scales::cut_short_scale(), accuracy = 1)
asinh_x <- function(brks, labs = compact_lab) scale_x_continuous(
  trans = "asinh", breaks = brks, labels = labs)

# Power-of-ten tick labels (10^x), signed: 0 stays "0", negatives -> -10^x.
label_pow10 <- function(b) parse(text = vapply(b, function(x) {
  if (is.na(x)) return(NA_character_)
  if (x == 0)  return("0")
  paste0(if (x < 0) "-" else "", "10^", round(log10(abs(x))))
}, character(1)))

raincloud <- function(flux_col, detect_col, xlab, brks = waiver(), use_asinh = FALSE,
                      xmin = -1, xlabs = compact_lab) {
  p <- ggplot(dat, aes(x = .data[[flux_col]], y = height_bin)) +
    # ridge: half-density per bin, sitting above each row
    stat_halfeye(adjust = 1, width = 0.7, .width = 0, justification = -0.18,
                 height = 0.7, point_colour = NA, fill = "grey85",
                 normalize = "groups") +
    # box per bin (narrow), nudged just below the ridge baseline
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, colour = "grey30",
                 position = position_nudge(y = -0.05)) +
    # jittered raindrops below the box; solid/open by MDF
    geom_point(aes(shape = .data[[detect_col]], colour = .data[[detect_col]]),
               position = position_jitter(width = 0, height = 0.06, seed = 1),
               size = 2.1, stroke = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    scale_shape_manual(values = shape_vals, drop = FALSE) +
    scale_colour_manual(values = col_vals, drop = FALSE) +
    labs(x = xlab, y = NULL, subtitle = sub) + theme_rc
  if (use_asinh) p <- p + asinh_x(brks, xlabs)
  # extend the flux axis down to xmin (e.g. -1) without dropping any points
  p + coord_cartesian(xlim = c(xmin, NA))
}

# Dot-density "cloud" variant: replaces the kernel ridge with stacked dots
# (ggdist::stat_dots), which reads better than a smooth density at small n.
# Dots are coloured/shaped solid (>= MDF) vs open (< MDF) like the raindrops.
raincloud_dots <- function(flux_col, detect_col, xlab, brks = waiver(),
                           use_asinh = FALSE, xmin = -1, xlabs = compact_lab) {
  p <- ggplot(dat, aes(x = .data[[flux_col]], y = height_bin)) +
    stat_dots(aes(shape = .data[[detect_col]], colour = .data[[detect_col]],
                  fill = .data[[detect_col]]),
              side = "top", justification = -0.05, scale = 0.75,
              dotsize = 1.0, stackratio = 1.05, layout = "hex", overflow = "compress") +
    geom_boxplot(width = 0.10, outlier.shape = NA, alpha = 0.5, colour = "grey30",
                 position = position_nudge(y = -0.04)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    scale_shape_manual(values = shape_vals, drop = FALSE) +
    scale_colour_manual(values = col_vals, drop = FALSE) +
    scale_fill_manual(values = col_vals, drop = FALSE) +
    labs(x = xlab, y = NULL, subtitle = sub) + theme_rc
  if (use_asinh) p <- p + asinh_x(brks, xlabs)
  p + coord_cartesian(xlim = c(xmin, NA))
}

ch4_lab <- expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
co2_lab <- expression(CO[2]~flux~(mu*mol~m^{-2}~s^{-1}))

# --- CH4 mass flux: ug CH4 m-2 h-1 -------------------------------------------
# 1 nmol CH4 = 16.04 ng = 0.01604 ug; nmol m-2 s-1 -> ug m-2 h-1 multiply by
# 0.01604 ug/nmol * 3600 s/h = 57.744.
ch4_molar_mass <- 16.04
nmol_s_to_ug_h <- ch4_molar_mass * 1e-3 * 3600    # = 57.744
dat$CH4_flux_ug_m2_h <- dat$CH4_best.flux * nmol_s_to_ug_h
ch4_mass_lab <- expression(CH[4]~flux~(mu*g~m^{-2}~h^{-1}))

save_both <- function(p, name, w = 8, h = 5.5) {
  ggsave(file.path(plots_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(plots_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 200, device = ragg::agg_png)
}

ch4_brks <- c(-1, 0, 1, 5, 25, 100, 400)

# --- Kernel-density ridge raincloud (axis extended to -1) ---------------------
save_both(raincloud("CH4_best.flux", "CH4_detect", ch4_lab,
                    brks = ch4_brks, use_asinh = TRUE),
          "CH4_raincloud_by_heightbin")
save_both(raincloud("CO2_best.flux", "CO2_detect", co2_lab), "CO2_raincloud_by_heightbin")

# --- Dot-density cloud variant ------------------------------------------------
save_both(raincloud_dots("CH4_best.flux", "CH4_detect", ch4_lab,
                         brks = ch4_brks, use_asinh = TRUE),
          "CH4_dotcloud_by_heightbin")
save_both(raincloud_dots("CO2_best.flux", "CO2_detect", co2_lab),
          "CO2_dotcloud_by_heightbin")

# --- CH4 in mass flux units (ug m-2 h-1) -------------------------------------
# Same data/detection as above, just rescaled; axis extended to the -1 nmol/s
# equivalent (= -57.7 ug/h). Round mass breaks, compact K labels.
ch4_brks_ug <- c(-100, -10, -1, 0, 1, 10, 100, 1e3, 1e4)
xmin_ug     <- -nmol_s_to_ug_h
save_both(raincloud("CH4_flux_ug_m2_h", "CH4_detect", ch4_mass_lab,
                    brks = ch4_brks_ug, use_asinh = TRUE, xmin = xmin_ug,
                    xlabs = label_pow10),
          "CH4_raincloud_by_heightbin_ug_h")
save_both(raincloud_dots("CH4_flux_ug_m2_h", "CH4_detect", ch4_mass_lab,
                         brks = ch4_brks_ug, use_asinh = TRUE, xmin = xmin_ug,
                         xlabs = label_pow10),
          "CH4_dotcloud_by_heightbin_ug_h")

summ <- dat %>% group_by(height_bin) %>%
  summarise(n = n(), CH4_ge_MDF = sum(CH4_detect == det_lab[1]),
            CO2_ge_MDF = sum(CO2_detect == det_lab[1]), .groups = "drop")
message("Detection (>= per-measurement Allan MDF, ", mdf_level, ") by height bin:")
print(as.data.frame(summ), row.names = FALSE)
message("Saved raincloud plots to: ", plots_dir)
