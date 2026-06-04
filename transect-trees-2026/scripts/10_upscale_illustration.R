# =============================================================================
# 10_upscale_illustration.R
# Illustrate the upscaling problem: plot mean +/- SE CH4 flux per height bin
# (the data), then overlay each base-treatment model across height. All models
# agree where there are data and FAN OUT by orders of magnitude below the lowest
# chamber (the unmeasured 0-30 cm base, shaded), which is the part that dominates
# the whole-stem integral.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/scripts/00_setup.R"))
suppressMessages({ library(dplyr); library(ggplot2); library(scales) })
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

dat <- read.csv(file.path(results_dir, "transect_trees_2026_fluxes.csv"),
                stringsAsFactors = FALSE, encoding = "UTF-8") %>%
  mutate(Height_m = height_cm / 100) %>% filter(is.finite(CH4_best.flux))

# --- mean +/- SE per height bin ----------------------------------------------
hb_levels <- c("Base (40-50 cm)","80 cm","160 cm","Mid (400-500 cm)","Canopy (800-1000 cm)")
dat$height_bin <- factor(case_when(
  dat$height_cm <= 50                         ~ hb_levels[1],
  dat$height_cm == 80                         ~ hb_levels[2],
  dat$height_cm == 160                        ~ hb_levels[3],
  dat$height_cm >= 400 & dat$height_cm <= 500 ~ hb_levels[4],
  dat$height_cm >= 800                        ~ hb_levels[5]), levels = hb_levels)
binm <- dat %>% group_by(height_bin) %>%
  summarise(H = mean(Height_m), F = mean(CH4_best.flux),
            SE = sd(CH4_best.flux)/sqrt(n()), n = n(), .groups = "drop") %>%
  filter(is.finite(H))

# --- pooled fits on positive points ------------------------------------------
pos <- dat %>% filter(CH4_best.flux > 0)
pw <- lm(log(CH4_best.flux) ~ log(Height_m), data = pos)   # power law
ex <- lm(log(CH4_best.flux) ~ Height_m,      data = pos)   # exponential
a <- exp(coef(pw)[1]); b <- coef(pw)[2]
Fpow <- function(h) a*h^b
Fexp <- function(h) exp(coef(ex)[1] + coef(ex)[2]*h)
cap  <- max(binm$F)                                        # physical ceiling = hottest bin mean
# plateau: log-linear interp of bin means, held flat below the lowest bin
ip <- approxfun(binm$H, log(binm$F), rule = 2)
Fplt <- function(h) exp(ip(h))

h_low_edge <- 0.30                                         # lowest chamber lower edge (30 cm)
hg  <- exp(seq(log(0.02), log(max(binm$H)), length.out = 400))
mods <- bind_rows(
  data.frame(h = hg, F = Fpow(hg),               model = "power law (full)"),
  data.frame(h = hg, F = pmin(Fpow(hg), cap),    model = "power law (capped)"),
  data.frame(h = hg, F = Fplt(hg),               model = "plateau"),
  data.frame(h = hg, F = Fexp(hg),               model = "exponential")
)
# measured-only line: interior interpolation drawn ONLY where there are data
mo <- data.frame(h = hg[hg >= min(binm$H)], model = "measured only (base = 0)")
mo$F <- Fplt(mo$h)
mods <- bind_rows(mods, mo)
mod_levels <- c("power law (full)","power law (capped)","plateau",
                "exponential","measured only (base = 0)")
mods$model <- factor(mods$model, levels = mod_levels)
cols <- setNames(c("#D7191C","#FDAE61","#2C7BB6","#7B3294","grey55"), mod_levels)
ltys <- setNames(c("solid","solid","solid","solid","21"), mod_levels)

# predicted flux in the 0-20 cm base bin (h = 0.10 m) for each model
base_pts <- data.frame(model = factor(mod_levels[1:4], levels = mod_levels),
                       h = 0.10,
                       F = c(Fpow(0.10), min(Fpow(0.10), cap), Fplt(0.10), Fexp(0.10)))
base_pts$lab <- format(round(base_pts$F), big.mark = ",")

p <- ggplot() +
  annotate("rect", ymin = 0, ymax = h_low_edge, xmin = 1e-3, xmax = Inf,
           fill = "grey85", alpha = 0.5) +
  annotate("text", x = 1.2e-2, y = 0.20,
           label = "unmeasured base", size = 3, hjust = 0, fontface = "italic", colour = "grey35") +
  annotate("text", x = 300, y = 6.2,
           label = "diamonds = predicted flux\nin the 0-20 cm base bin",
           size = 3.1, hjust = 0, fontface = "italic", colour = "grey25") +
  geom_path(data = mods, aes(F, h, colour = model, linetype = model), linewidth = 0.9) +
  geom_errorbarh(data = binm, aes(y = H, xmin = pmax(F-SE,1e-3), xmax = F+SE),
                 height = 0.12, colour = "black") +
  geom_point(data = binm, aes(F, H), size = 3, colour = "black") +
  # explicit base-bin predictions: the >100x fan
  geom_point(data = base_pts, aes(F, h, colour = model), size = 3.4, shape = 18) +
  geom_text(data = base_pts, aes(F, h, colour = model, label = lab),
            vjust = -1.1, size = 3.1, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = cols) + scale_linetype_manual(values = ltys) +
  scale_x_log10(breaks = c(1e-2,1e-1,1,10,100,1e3,1e4,1e5),
                labels = label_number(scale_cut = cut_short_scale())) +
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim = c(1e-2, 1e5), ylim = c(0, max(binm$H)+0.2)) +
  labs(x = expression("CH"[4]~flux~(nmol~m^{-2}~s^{-1})*", log"),
       y = "Height above water (m)", colour = NULL, linetype = NULL,
       title = "Why the whole-stem flux hinges on the unmeasured base",
       subtitle = "Points = mean ± SE per height bin. Models agree over the data, diverge >100x below the lowest chamber.") +
  theme_classic(base_size = 12) + theme(legend.position = "right")

ggsave(file.path(plots_dir, "CH4_upscale_models_illustration.png"), p,
       width = 9, height = 6.5, dpi = 200, device = ragg::agg_png)
ggsave(file.path(plots_dir, "CH4_upscale_models_illustration.pdf"), p,
       width = 9, height = 6.5, device = cairo_pdf)

# --- LINEAR flux-axis version ------------------------------------------------
# Power law (full) base prediction = 1,001 dwarfs every measurement (all < 85);
# its curve runs off the right edge (clipped). Capped/plateau base = 68, exp = 7.
x_hi <- 1080
p_lin <- ggplot() +
  annotate("rect", ymin = 0, ymax = h_low_edge, xmin = -Inf, xmax = Inf,
           fill = "grey85", alpha = 0.5) +
  annotate("text", x = 540, y = 6.2,
           label = "diamonds = predicted flux\nin the 0-20 cm base bin",
           size = 3.1, hjust = 0, fontface = "italic", colour = "grey25") +
  annotate("text", x = x_hi, y = 1.5, label = "power law (full)\nruns off scale →",
           size = 3, hjust = 1, colour = cols[["power law (full)"]]) +
  geom_path(data = mods, aes(F, h, colour = model, linetype = model), linewidth = 0.9) +
  geom_errorbarh(data = binm, aes(y = H, xmin = F-SE, xmax = F+SE),
                 height = 0.12, colour = "black") +
  geom_point(data = binm, aes(F, H), size = 3, colour = "black") +
  geom_point(data = base_pts, aes(F, h, colour = model), size = 3.4, shape = 18) +
  geom_text(data = base_pts, aes(F, h, colour = model, label = lab),
            vjust = -1.1, size = 3.1, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = cols) + scale_linetype_manual(values = ltys) +
  scale_x_continuous(breaks = seq(0, 1000, 250), labels = label_number()) +
  coord_cartesian(xlim = c(-30, x_hi), ylim = c(0, max(binm$H)+0.2)) +
  labs(x = expression("CH"[4]~flux~(nmol~m^{-2}~s^{-1})),
       y = "Height above water (m)", colour = NULL, linetype = NULL,
       title = "Why the whole-stem flux hinges on the unmeasured base",
       subtitle = "Linear flux axis: the full power-law base prediction (1,001) dwarfs every measurement (all < 85).") +
  theme_classic(base_size = 12) + theme(legend.position = "right")

ggsave(file.path(plots_dir, "CH4_upscale_models_illustration_linear.png"), p_lin,
       width = 9, height = 6.5, dpi = 200, device = ragg::agg_png)
ggsave(file.path(plots_dir, "CH4_upscale_models_illustration_linear.pdf"), p_lin,
       width = 9, height = 6.5, device = cairo_pdf)

cat(sprintf("Pooled power law: flux = %.2f * h^%.2f ;  cap = %.1f ;  exp slope k = %.2f\n",
            a, b, cap, coef(ex)[2]))
cat("Model flux at the 0-20 cm base bin (h = 0.10 m), nmol m-2 s-1:\n")
cat(sprintf("  power law(full)=%.0f   capped=%.0f   plateau=%.0f   exp=%.0f\n",
            Fpow(0.10), min(Fpow(0.10),cap), Fplt(0.10), Fexp(0.10)))
message("Saved: plots/CH4_upscale_models_illustration{,_linear}.{png,pdf}")
