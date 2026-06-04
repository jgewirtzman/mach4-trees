# =============================================================================
# 07_plots.R  (BCI)
# Paired Base-vs-2m comparison of CO2 and CH4 fluxes (slopegraph): one line per
# tree connecting its Base and 2 m measurement, coloured by species. CH4 points
# are SOLID if |flux| >= MDF and OPEN if < MDF (per-measurement Allan / Christiansen).
# Can be run with Rscript (non-interactive).
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/bci-fluxes-2026/scripts/00_setup.R"))
suppressMessages({ library(ggplot2); library(dplyr); library(scales) })
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)

mdf_level <- "chr95"   # per-measurement Allan deviation (Christiansen), 95%

dat <- read.csv(file.path(results_dir, "bci_2026_fluxes_with_mdf.csv"),
                stringsAsFactors = FALSE, encoding = "UTF-8") %>%
  mutate(species = enc2utf8(species),
         height  = factor(height_label, levels = c("Base", "2m"),
                          labels = c("Base", "2 m")))

det_lab <- c("≥ MDF", "< MDF")              # solid, open
dat$CH4_detect <- factor(ifelse(abs(dat$CH4_best.flux) >=
                                  dat[[paste0("CH4_MDF_", mdf_level)]],
                                det_lab[1], det_lab[2]), levels = det_lab)
shape_vals <- c(19, 1); names(shape_vals) <- det_lab

sp_cols <- c("Simarouba amara" = "#1B9E77", "Heisteria concinna" = "#7570B3")

theme_p <- theme_classic(base_size = 12) +
  theme(legend.position = "top", legend.title = element_blank(),
        plot.title = element_text(face = "bold"))

# CH4 mass flux: 1 nmol CH4 = 0.01604 ug; nmol m-2 s-1 -> ug m-2 h-1 (x3600)
dat$CH4_flux_ug_h <- dat$CH4_best.flux * 16.04e-3 * 3600
mdf_sub <- paste0("solid = |flux| ≥ MDF, open = < MDF (per-measurement Allan, ",
                  sub("chr", "", mdf_level), "%)")

# --- CO2 panel ---------------------------------------------------------------
co2_lab <- expression(CO[2]~flux~(mu*mol~m^{-2}~s^{-1}))
p_co2 <- ggplot(dat, aes(height, CO2_best.flux, group = TreeID, colour = species)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = sp_cols) +
  labs(x = NULL, y = co2_lab, title = "CO2") + theme_p

# --- CH4 panel builder (solid >= MDF, open < MDF) ----------------------------
ch4_panel <- function(yvar, ylab) {
  ggplot(dat, aes(height, .data[[yvar]], group = TreeID)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
    geom_line(aes(colour = species), linewidth = 0.7, alpha = 0.8) +
    geom_point(aes(colour = species, shape = CH4_detect), size = 3, stroke = 0.9) +
    scale_colour_manual(values = sp_cols) +
    scale_shape_manual(values = shape_vals, drop = FALSE) +
    labs(x = NULL, y = ylab, title = "CH4", subtitle = mdf_sub) +
    theme_p + theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
}
ch4_lab    <- expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
ch4_ug_lab <- expression(CH[4]~flux~(mu*g~m^{-2}~h^{-1}))
p_ch4    <- ch4_panel("CH4_best.flux", ch4_lab)
p_ch4_ug <- ch4_panel("CH4_flux_ug_h", ch4_ug_lab)

save_both <- function(p, name, w = 6.5, h = 5.5) {
  ggsave(file.path(plots_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(plots_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 200, device = ragg::agg_png)
}

# height on x (upright slopegraph)
save_both(p_co2,    "BCI_base_vs_2m_CO2")
save_both(p_ch4,    "BCI_base_vs_2m_CH4")        # nmol m-2 s-1
save_both(p_ch4_ug, "BCI_base_vs_2m_CH4_ug")     # ug m-2 s-1

# height on the VERTICAL axis (coord_flip): Base at bottom, 2 m at top
flip <- coord_flip()
save_both(p_co2    + flip, "BCI_base_vs_2m_CO2_byheight")
save_both(p_ch4    + flip, "BCI_base_vs_2m_CH4_byheight")
save_both(p_ch4_ug + flip, "BCI_base_vs_2m_CH4_ug_byheight")

if (have_patchwork) {
  ann <- function(p) p + patchwork::plot_annotation(
    title = "BCI stem fluxes: Base vs 2 m (paired by tree)",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  wp <- function(a, b) patchwork::wrap_plots(a, b, ncol = 2)
  save_both(ann(wp(p_co2, p_ch4)),               "BCI_base_vs_2m",             w = 11, h = 5.5)
  save_both(ann(wp(p_co2, p_ch4_ug)),            "BCI_base_vs_2m_ug",          w = 11, h = 5.5)
  save_both(ann(wp(p_co2 + flip, p_ch4 + flip)), "BCI_base_vs_2m_byheight",    w = 11, h = 5)
  save_both(ann(wp(p_co2 + flip, p_ch4_ug + flip)), "BCI_base_vs_2m_ug_byheight", w = 11, h = 5)
  message("Saved combined: plots/BCI_base_vs_2m{,_ug}{,_byheight}.{png,pdf}")
}

# --- requested standalone: CH4 (ug/h), height on y, fixed -200..+400 axis ----
# Heisteria = orange, Simarouba = yellow. Note: BCI CH4 is ~0 on this scale.
sp_cols2 <- c("Heisteria concinna" = "#E69F00",   # orange
              "Simarouba amara"    = "#F0E442")   # yellow
p_req <- ggplot(dat, aes(CH4_flux_ug_h, height, group = TreeID)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey50") +
  geom_line(aes(colour = species), linewidth = 0.7, alpha = 0.85) +
  geom_point(aes(colour = species, shape = CH4_detect), size = 3, stroke = 0.9) +
  scale_colour_manual(values = sp_cols2) +
  scale_shape_manual(values = shape_vals, drop = FALSE) +
  coord_cartesian(xlim = c(-200, 400)) +
  labs(x = ch4_ug_lab, y = NULL, title = "CH4", subtitle = mdf_sub) +
  theme_p + theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
save_both(p_req, "BCI_CH4_ug_byheight_axis200_400", w = 7, h = 5)
message("Saved: plots/BCI_CH4_ug_byheight_axis200_400.{png,pdf}")

# --- same figure with an INSET zooming to the actual data range --------------
if (have_patchwork) {
  rng <- range(dat$CH4_flux_ug_h); pad <- diff(rng) * 0.15
  inset <- ggplot(dat, aes(CH4_flux_ug_h, height, group = TreeID)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_line(aes(colour = species), linewidth = 0.6, alpha = 0.85) +
    geom_point(aes(colour = species, shape = CH4_detect), size = 2, stroke = 0.8) +
    scale_colour_manual(values = sp_cols2) +
    scale_shape_manual(values = shape_vals, drop = FALSE) +
    coord_cartesian(xlim = c(rng[1] - pad, rng[2] + pad)) +
    labs(x = expression(mu*g~m^{-2}~h^{-1}), y = NULL, title = "zoom: actual range") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, face = "italic"),
          plot.background = element_rect(colour = "grey40", fill = "white"))
  p_main <- p_req +
    annotate("curve", x = 20, y = 1.55, xend = 150, yend = 2.05,
             curvature = -0.3, linewidth = 0.4, colour = "grey45",
             arrow = arrow(length = unit(0.18, "cm"), ends = "first")) +
    annotate("text", x = 150, y = 2.12, label = "values ≈ 0", size = 3,
             hjust = 0.4, colour = "grey35")
  p_inset <- p_main +
    patchwork::inset_element(inset, left = 0.46, bottom = 0.46, right = 0.99, top = 0.93)
  save_both(p_inset, "BCI_CH4_ug_byheight_axis200_400_inset", w = 7.5, h = 5)
  message("Saved: plots/BCI_CH4_ug_byheight_axis200_400_inset.{png,pdf}")
}

# --- paired differences (Base - 2m) ------------------------------------------
wide <- dat %>%
  select(TreeID, species, height_label, CO2_best.flux, CH4_best.flux) %>%
  tidyr::pivot_wider(names_from = height_label,
                     values_from = c(CO2_best.flux, CH4_best.flux))
cat("\nPaired Base - 2m differences:\n")
cat(sprintf("  CO2: mean Base-2m = %+.2f umol m-2 s-1 (Base higher in %d/%d trees)\n",
            mean(wide$CO2_best.flux_Base - wide$`CO2_best.flux_2m`),
            sum(wide$CO2_best.flux_Base > wide$`CO2_best.flux_2m`), nrow(wide)))
cat(sprintf("  CH4: mean Base-2m = %+.3f nmol m-2 s-1\n",
            mean(wide$CH4_best.flux_Base - wide$`CH4_best.flux_2m`)))
message("Saved: plots/BCI_base_vs_2m_CO2.{png,pdf}, BCI_base_vs_2m_CH4.{png,pdf}")
