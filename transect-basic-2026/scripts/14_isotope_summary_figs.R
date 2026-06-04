# =============================================================================
# 14_isotope_summary_figs.R
# Two presentation summary figures for the isotope sampling campaign, built
# from results/sample_inventory_isotopes.csv:
#   FIG 1  sampling coverage  - vials by tree x height, per campaign
#   FIG 2  what we captured   - CH4 (and CO2) concentration vs height & timepoint
# Run with Rscript.  -> results/figs/isotope_summary_{coverage,concentration}.{pdf,png}
# =============================================================================

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
base <- "/Users/jongewirtzman/My Drive/Research/mach4-trees"
inv  <- read.csv(file.path(base, "results/sample_inventory_isotopes.csv"),
                 stringsAsFactors = FALSE, encoding = "UTF-8")
inv$height_m  <- inv$height_cm / 100
inv$campaign  <- ifelse(inv$dataset == "climbing", "Climbed (intensive)", "Ground (basic)")
tp_cols <- c(t0 = "#1f78b4", t1 = "#e31a1c", t2 = "#ff7f00")
figdir  <- file.path(base, "results/figs"); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
save_both <- function(p, name, w, h) {
  ggsave(file.path(figdir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, device = ragg::agg_png)
}

n_vial <- nrow(inv)
n_meas <- inv %>% distinct(dataset, TreeID, height_cm) %>% nrow()
n_tree <- length(unique(inv$TreeID))

# ===== FIG 1: sampling coverage ==============================================
meas <- inv %>%
  group_by(campaign, TreeID, height_cm, height_m) %>%
  summarise(n_vial = n(), n_tp = n_distinct(timepoint), .groups = "drop")
# order trees within each campaign by the max height reached
ord <- meas %>% group_by(campaign, TreeID) %>% summarise(mx = max(height_m), .groups = "drop") %>%
  arrange(campaign, mx)
meas$TreeID <- factor(meas$TreeID, levels = unique(ord$TreeID))

p1 <- ggplot(meas, aes(height_m, TreeID)) +
  geom_line(aes(group = TreeID), colour = "grey80", linewidth = 0.5) +
  geom_point(aes(size = n_tp, colour = campaign)) +
  facet_grid(campaign ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(breaks = c(0.4, 1, 2, 4, 6, 8, 10)) +
  scale_size_continuous(range = c(2.2, 4.6), breaks = c(2, 3)) +
  scale_colour_manual(values = c("Climbed (intensive)" = "#1b7837", "Ground (basic)" = "#762a83"),
                      guide = "none") +
  labs(x = "Sampling height (m)", y = NULL, size = "timepoints\n(t0,t1[,t2])",
       title = "Isotope sampling coverage",
       subtitle = sprintf("%d trees · %d chamber measurements · %d vials (t0 ambient, t1, t2)",
                          n_tree, n_meas, n_vial)) +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        panel.grid.major.x = element_line(colour = "grey92"))
save_both(p1, "isotope_summary_coverage", w = 8, h = 6)

# ===== FIG 2: concentrations captured ========================================
gl <- inv %>%
  select(campaign, TreeID, height_m, timepoint, CO2_ppm, CH4_ppm) %>%
  pivot_longer(c(CH4_ppm, CO2_ppm), names_to = "gas", values_to = "conc") %>%
  mutate(gas = recode(gas, CH4_ppm = "CH[4]~(ppm)", CO2_ppm = "CO[2]~(ppm)"))

p2 <- ggplot(gl, aes(height_m, conc, colour = timepoint)) +
  geom_point(position = position_jitter(width = 0.06, height = 0), size = 1.9, alpha = 0.85) +
  facet_wrap(~ gas, scales = "free_y", labeller = label_parsed) +
  scale_colour_manual(values = tp_cols, name = "timepoint",
                      labels = c(t0 = "t0 (ambient)", t1 = "t1", t2 = "t2")) +
  scale_x_continuous(breaks = c(0.4, 1, 2, 4, 6, 8, 10)) +
  labs(x = "Sampling height (m)", y = "concentration in vial",
       title = "What the isotope vials captured",
       subtitle = "Strong CH4 build-up at the stem base; near-ambient toward the canopy") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10, colour = "grey30"), legend.position = "top")

# CH4 panel reads best on a log axis (spans ambient ~2 to ~200 ppm); CO2 linear.
# Use a single sqrt transform that tames both without hiding the ambient floor.
p2 <- p2 + scale_y_continuous(trans = "sqrt")
save_both(p2, "isotope_summary_concentration", w = 9, h = 5)

cat("Saved 2 figures to", figdir, "\n  isotope_summary_coverage.[pdf/png]\n  isotope_summary_concentration.[pdf/png]\n")
