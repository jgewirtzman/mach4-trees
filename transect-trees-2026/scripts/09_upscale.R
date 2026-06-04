# =============================================================================
# 09_upscale.R
# Upscale per-height CH4 chamber fluxes to a whole-stem emission, integrating
# from the WATER SURFACE (h = 0, the physical lower boundary in this flooded
# forest) to the top of the measured stem, in 20-cm height bins (= chamber
# height). The base bins (below the lowest chamber, ~0-40 cm) are unmeasured and
# dominate the integral, so we compute the whole-stem flux under several base
# treatments and report the spread.
#
# Output is VERTICALLY-INTEGRATED FLUX per unit stem circumference:
#   VIF = sum_bins  F(bin) * 0.20 m      [nmol m^-1 s^-1]   (per m of circumference)
# Multiply by stem circumference (pi * diameter) for a whole-stem rate. No
# diameters were recorded, so a nominal cylinder is used for the whole-stem
# numbers (edit stem_diam_m); taper would make the base dominate even more.
#
# Methods (each = interior model + base-bin treatment):
#   powerlaw_full   : per-tree power law  F=a*h^b   for ALL bins (10 cm -> top)
#   exp_full        : per-tree exponential F=a*e^kh  for ALL bins
#   data_plateau    : measured pts (log-linear interp) interior; base = lowest
#                     measured flux held flat to the waterline      [central]
#   data_meas_only  : measured interior; base bins = 0              [lower bound]
#   powerlaw_capped : power law interior+base, every bin capped at the tree's
#                     max measured flux (physical ceiling)          [tamed upper]
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))
suppressMessages({ library(dplyr); library(ggplot2); library(scales) })
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

bin_w        <- 0.20                         # m, = chamber height
nmol_s_to_ug_h <- 16.04 * 1e-3 * 3600        # 57.744 : nmol m-2 s-1 -> ug m-2 h-1
stem_diam_m  <- 0.20                         # <-- nominal stem diameter (NO data); edit
circ_m       <- pi * stem_diam_m             # stem circumference (m)

dat <- read.csv(file.path(results_dir, "transect_trees_2026_fluxes.csv"),
                stringsAsFactors = FALSE, encoding = "UTF-8") %>%
  mutate(Height_m = height_cm / 100) %>%
  filter(is.finite(CH4_best.flux))

# --- per-tree flux predictor for each method ---------------------------------
# returns a function f(h_vector) -> flux at those heights, for one tree
make_predictors <- function(g) {
  gp  <- g[g$CH4_best.flux > 0, ]                  # positive pts for log fits
  hmin <- min(g$Height_m); hmax <- max(g$Height_m)
  fmin <- g$CH4_best.flux[which.min(g$Height_m)]   # lowest-chamber flux
  fmax <- max(g$CH4_best.flux)
  # log-linear interpolation across measured pts (floor tiny/neg for the log)
  oh <- order(g$Height_m)
  xi <- g$Height_m[oh]; yi <- pmax(g$CH4_best.flux[oh], 1e-4)
  interp <- approxfun(xi, log(yi), rule = 2)       # rule=2 -> flat beyond ends
  pw <- if (nrow(gp) >= 3 && length(unique(gp$Height_m)) >= 3)
          lm(log(CH4_best.flux) ~ log(Height_m), data = gp) else NULL
  ex <- if (nrow(gp) >= 3 && length(unique(gp$Height_m)) >= 3)
          lm(log(CH4_best.flux) ~ Height_m, data = gp) else NULL
  list(
    hmin = hmin, hmax = hmax,
    powerlaw_full   = function(h) if (is.null(pw)) rep(NA, length(h)) else
                        exp(predict(pw, data.frame(Height_m = h))),
    exp_full        = function(h) if (is.null(ex)) rep(NA, length(h)) else
                        exp(predict(ex, data.frame(Height_m = h))),
    data_plateau    = function(h) exp(interp(h)),                  # flat base via rule=2
    data_meas_only  = function(h) ifelse(h < hmin, 0, exp(interp(h))),
    powerlaw_capped = function(h) if (is.null(pw)) rep(NA, length(h)) else
                        pmin(exp(predict(pw, data.frame(Height_m = h))), fmax)
  )
}

methods <- c("powerlaw_full","exp_full","data_plateau","data_meas_only","powerlaw_capped")

# --- integrate each tree x method --------------------------------------------
rows <- list()
for (tid in unique(dat$TreeID)) {
  g  <- dat[dat$TreeID == tid, ]
  pr <- make_predictors(g)
  mids <- seq(bin_w/2, pr$hmax, by = bin_w)        # bin midpoints: 0.10, 0.30, ...
  base_idx <- mids < pr$hmin                       # bins below the lowest chamber
  for (m in methods) {
    f <- pr[[m]](mids)
    if (all(is.na(f))) next
    VIF      <- sum(f * bin_w, na.rm = TRUE)        # nmol m^-1 s^-1
    VIF_base <- sum(f[base_idx] * bin_w, na.rm = TRUE)
    rows[[length(rows)+1]] <- data.frame(
      TreeID = tid, method = m,
      VIF_nmol_m_s = VIF,
      base_share = ifelse(VIF != 0, VIF_base / VIF, NA),
      whole_stem_ug_h = VIF * circ_m * nmol_s_to_ug_h)   # cylinder assumption
  }
}
up <- bind_rows(rows)
up$method <- factor(up$method, levels = methods)

# --- summaries ---------------------------------------------------------------
tot <- up %>% group_by(method) %>%
  summarise(n_trees = n(),
            mean_VIF = mean(VIF_nmol_m_s),
            median_VIF = median(VIF_nmol_m_s),
            mean_base_share = mean(base_share, na.rm = TRUE),
            mean_stem_ug_h = mean(whole_stem_ug_h), .groups = "drop")

write.csv(up,  file.path(results_dir, "ch4_upscaling_by_tree_method.csv"), row.names = FALSE)
write.csv(tot, file.path(results_dir, "ch4_upscaling_method_summary.csv"), row.names = FALSE)

cat("\n=== Whole-stem CH4 by base treatment (per unit circumference) ===\n")
cat(sprintf("(nominal stem diameter %.2f m for the ug/h column; base = below lowest chamber)\n\n",
            stem_diam_m))
print(as.data.frame(tot %>% mutate(across(where(is.numeric), ~round(.,2)))), row.names = FALSE)

cat("\n=== Per-tree vertically-integrated flux (nmol m^-1 s^-1) ===\n")
wide <- up %>% select(TreeID, method, VIF_nmol_m_s) %>%
  tidyr::pivot_wider(names_from = method, values_from = VIF_nmol_m_s)
print(as.data.frame(wide %>% mutate(across(where(is.numeric), ~round(.,1)))), row.names = FALSE)

# --- plot: whole-stem emission by tree, grouped by method --------------------
p <- ggplot(up, aes(TreeID, whole_stem_ug_h, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_y_continuous(trans = "log10", labels = label_number(),
                     breaks = c(1,10,100,1e3,1e4,1e5)) +
  annotation_logticks(sides = "l") +
  labs(x = NULL, y = expression("Whole-stem CH"[4]*" emission ("*mu*g~h^{-1}*", log; cylinder D=0.2 m)"),
       fill = "base treatment",
       title = "CH4 upscaling: sensitivity to the unmeasured base",
       subtitle = "integrated waterline -> top in 20-cm bins; bars within a tree = the five methods") +
  theme_classic(base_size = 12) + theme(legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(plots_dir, "CH4_upscaling_by_method.png"), p, width = 10, height = 6,
       dpi = 200, device = ragg::agg_png)
ggsave(file.path(plots_dir, "CH4_upscaling_by_method.pdf"), p, width = 10, height = 6,
       device = cairo_pdf)
message("\nSaved: results/ch4_upscaling_by_tree_method.csv, ch4_upscaling_method_summary.csv")
message("Saved: plots/CH4_upscaling_by_method.{png,pdf}")
