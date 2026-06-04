# =============================================================================
# 07_plots.R
# Flux vs stem height plots for CO2 and CH4. Style follows
# whole_tree_flux/.../08_ch4_height_plot.R: coord_flip so height is the vertical
# axis, dashed zero line, points + loess, faceted by tree, plus an asinh-scaled
# version for the wide-range / near-zero CH4 fluxes.
# Can be run with Rscript (non-interactive).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))
suppressMessages({ library(ggplot2); library(scales) })

# Render accented species names (Macucú, Ingá, ...) correctly even under a C locale
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

dat <- read.csv(file.path(results_dir, "transect_trees_2026_fluxes.csv"),
                stringsAsFactors = FALSE, encoding = "UTF-8") %>%
  mutate(species = enc2utf8(species),
         Height_m   = height_cm / 100,
         Tree_label = enc2utf8(paste0(TreeID, " — ", species)))

# --- Common theme & layers ---------------------------------------------------

plot_theme <- theme_classic(base_size = 13) +
  theme(strip.text = element_text(size = 10),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.9, "lines"))

# Exponential fit per group: flux = a * exp(k * height), fit as log-linear on
# the positive fluxes (robust; exp can't take <=0). Returns a fine height grid
# with fitted values for geom_line. Groups with < 2 positive points are skipped.
exp_fit <- function(df, yvar, by) {
  do.call(rbind, lapply(split(df, df[[by]]), function(g) {
    gp <- g[is.finite(g[[yvar]]) & g[[yvar]] > 0, ]
    if (nrow(gp) < 2 || length(unique(gp$Height_m)) < 2) return(NULL)
    m <- try(lm(log(gp[[yvar]]) ~ gp$Height_m), silent = TRUE)
    if (inherits(m, "try-error") || any(is.na(coef(m)))) return(NULL)
    # extend the curve from ground (0 m) up through the canopy so it spans the
    # full decay (peak at the base, asymptoting toward zero with height)
    hg <- seq(0, max(g$Height_m), length.out = 80)
    data.frame(setNames(list(g[[by]][1]), by),
               Height_m = hg, fit = exp(coef(m)[1] + coef(m)[2] * hg),
               check.names = FALSE)
  }))
}

# Power-law fit per group: flux = a * height^b, fit as log-log linear on the
# positive fluxes. Beats the exponential here (AIC 175 vs 189, adj R2 0.66 vs
# 0.51): the steep super-exponential rise toward the base is a straight line in
# log-log space. The curve DIVERGES at h = 0 (can't draw to the ground), so the
# grid runs from the lowest measured height up; it asymptotes to ~0 with height.
pow_fit <- function(df, yvar, by) {
  do.call(rbind, lapply(split(df, df[[by]]), function(g) {
    gp <- g[is.finite(g[[yvar]]) & g[[yvar]] > 0 & g$Height_m > 0, ]
    if (nrow(gp) < 2 || length(unique(gp$Height_m)) < 2) return(NULL)
    m <- try(lm(log(gp[[yvar]]) ~ log(gp$Height_m)), silent = TRUE)
    if (inherits(m, "try-error") || any(is.na(coef(m)))) return(NULL)
    hg <- exp(seq(log(min(gp$Height_m)), log(max(gp$Height_m)), length.out = 160))
    data.frame(setNames(list(g[[by]][1]), by),
               Height_m = hg, fit = exp(coef(m)[1] + coef(m)[2] * log(hg)),
               check.names = FALSE)
  }))
}

# lines between points (grey) + exponential fit (coloured) + points; no loess.
prof_layers <- function(ylab, fit_df, fit_col) list(
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40"),
  geom_line(colour = "grey60", linewidth = 0.4),                       # connect points
  if (!is.null(fit_df) && nrow(fit_df) > 0)
    geom_line(data = fit_df, aes(Height_m, fit), inherit.aes = FALSE,
              colour = fit_col, linewidth = 0.9),                      # exponential fit
  geom_point(size = 2.6, alpha = 0.9),
  labs(x = "Height (m)", y = ylab),
  plot_theme, coord_flip()
)
ch4_col <- "#2E8B57"; co2_col <- "#D2691E"

co2_lab <- expression(CO[2]~flux~(mu*mol~m^{-2}~s^{-1}))
ch4_lab <- expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))

save_both <- function(p, name, w = 11, h = 8) {
  # cairo_pdf + ragg::agg_png both render UTF-8 (accented species names) reliably
  ggsave(file.path(plots_dir, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(plots_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300, device = ragg::agg_png)
}

# =============================================================================
# 1. Faceted by tree
# =============================================================================

ch4_fit <- exp_fit(dat, "CH4_best.flux", "Tree_label")
ch4_facet <- ggplot(dat, aes(Height_m, CH4_best.flux, group = Tree_label)) +
  prof_layers(ch4_lab, ch4_fit, ch4_col) +
  facet_wrap(~ Tree_label, scales = "free_x")
save_both(ch4_facet, "CH4_flux_by_height_byTree")

ch4_facet_asinh <- ch4_facet + scale_y_continuous(trans = "asinh")
save_both(ch4_facet_asinh, "CH4_flux_by_height_byTree_asinh")

co2_fit <- exp_fit(dat, "CO2_best.flux", "Tree_label")
co2_facet <- ggplot(dat, aes(Height_m, CO2_best.flux, group = Tree_label)) +
  prof_layers(co2_lab, co2_fit, co2_col) +
  facet_wrap(~ Tree_label, scales = "free_x")
save_both(co2_facet, "CO2_flux_by_height_byTree")

# =============================================================================
# 2. All trees combined (one profile, coloured by tree)
# =============================================================================

combo_theme <- theme(legend.position = "right")
ch4_fit_t <- exp_fit(dat, "CH4_best.flux", "TreeID")
co2_fit_t <- exp_fit(dat, "CO2_best.flux", "TreeID")

ch4_all <- ggplot(dat, aes(Height_m, CH4_best.flux, colour = TreeID)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_line(aes(group = TreeID), alpha = 0.35) +                         # connect points
  geom_line(data = ch4_fit_t, aes(Height_m, fit, colour = TreeID),       # exp fit (dashed)
            inherit.aes = FALSE, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 2.6, alpha = 0.9) +
  labs(x = "Height (m)", y = ch4_lab, colour = "Tree") +
  plot_theme + combo_theme + coord_flip()
save_both(ch4_all, "CH4_flux_by_height_allTrees", w = 8, h = 6)
save_both(ch4_all + scale_y_continuous(trans = "asinh"),
          "CH4_flux_by_height_allTrees_asinh", w = 8, h = 6)

co2_all <- ggplot(dat, aes(Height_m, CO2_best.flux, colour = TreeID)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_line(aes(group = TreeID), alpha = 0.35) +
  geom_line(data = co2_fit_t, aes(Height_m, fit, colour = TreeID),
            inherit.aes = FALSE, linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 2.6, alpha = 0.9) +
  labs(x = "Height (m)", y = co2_lab, colour = "Tree") +
  plot_theme + combo_theme + coord_flip()
save_both(co2_all, "CO2_flux_by_height_allTrees", w = 8, h = 6)

# =============================================================================
# 3. Mean profile across trees: mean +/- SE per HEIGHT BIN, with a POWER-LAW
#    fit (flux = a*h^b; better than exponential, AIC 175 vs 189) using ALL
#    individual measurements. The power law captures the super-exponential rise
#    toward the base; it asymptotes to ~0 with height (but diverges at h=0, so
#    the curve is drawn from the lowest measured height up, not to the ground).
# =============================================================================

# Height bins (actual heights vary within each nominal level)
hb_levels <- c("Base (40-50 cm)", "80 cm", "160 cm",
               "Mid (400-500 cm)", "Canopy (800-1000 cm)")
dat$height_bin <- factor(dplyr::case_when(
  dat$height_cm <= 50                          ~ hb_levels[1],
  dat$height_cm == 80                          ~ hb_levels[2],
  dat$height_cm == 160                         ~ hb_levels[3],
  dat$height_cm >= 400 & dat$height_cm <= 500  ~ hb_levels[4],
  dat$height_cm >= 800                         ~ hb_levels[5]
), levels = hb_levels)

# mean +/- SE per bin, positioned at the bin's mean actual height
mean_prof <- dat %>%
  group_by(height_bin) %>%
  summarise(Height_m = mean(height_cm, na.rm = TRUE) / 100,
            CH4_mean = mean(CH4_best.flux, na.rm = TRUE),
            CH4_se   = sd(CH4_best.flux, na.rm = TRUE) / sqrt(sum(!is.na(CH4_best.flux))),
            n = n(), .groups = "drop")

# Power-law fit on ALL individual points (flux = a*h^b; drawn from lowest
# measured height up, asymptoting toward zero with height)
dat$grp <- "all"
ch4_all_fit <- pow_fit(dat, "CH4_best.flux", "grp")

ch4_mean <- ggplot(mean_prof, aes(Height_m, CH4_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  { if (!is.null(ch4_all_fit))
      geom_line(data = ch4_all_fit, aes(Height_m, fit), inherit.aes = FALSE,
                colour = ch4_col, linewidth = 0.9) } +                  # power-law fit (all points)
  geom_errorbar(aes(ymin = CH4_mean - CH4_se, ymax = CH4_mean + CH4_se),
                width = 0.15, colour = "grey50") +
  geom_line(colour = "grey60", linewidth = 0.4) +                       # connect bin means
  geom_point(size = 3, colour = ch4_col) +
  labs(x = "Height (m)", y = ch4_lab,
       title = expression("Mean CH"[4]*" profile (± SE), power-law fit  flux = 3.0·h"^{-2.5})) +
  plot_theme + coord_flip()
save_both(ch4_mean, "CH4_flux_by_height_mean", w = 6.5, h = 6)

message("Saved height-profile plots to: ", plots_dir)
message("  CH4_flux_by_height_byTree[_asinh], CO2_flux_by_height_byTree")
message("  CH4/CO2_flux_by_height_allTrees[_asinh]")
message("  CH4_flux_by_height_mean")
