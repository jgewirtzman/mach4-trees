# =============================================================================
# 09_plots_dbh_remeasure.R  (BASIC trees)
# DBH: interpolate datasheet diameters (D30/D80/D160) to 130 cm. Adds DBH (+
# species) to the mixed model, a DBH-vs-flux plot, a climbed-vs-remeasured
# scatter (basic vs climbing fluxes for shared trees), and a linear-scale
# flux-by-height with multiple candidate models. Run with Rscript.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-basic-2026/scripts/08_plots.R"))  # reuse d, species cleaning, helpers
suppressMessages({ library(ggplot2); library(dplyr); library(scales) })

# --- DBH per tree: interpolate measured diameters to 130 cm ------------------
master <- read.csv(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/data/datasheets/parsed/all_basic_sheets.csv"),
  stringsAsFactors = FALSE, encoding = "UTF-8")
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
dbh_of <- function(H, D) { ok <- is.finite(H) & is.finite(D)
  if (sum(ok) == 0) return(NA_real_); if (sum(ok) == 1) return(D[ok])
  approx(H[ok], D[ok], xout = 130, rule = 2)$y }
master$DBH <- mapply(function(h1,h2,h3,d1,d2,d3)
  dbh_of(c(num(h1),num(h2),num(h3)), c(num(d1),num(d2),num(d3))),
  master$H1,master$H2,master$H3, master$D30,master$D80,master$D160)
dbh_lut <- master %>% filter(!grepl("no data", Note)) %>%
  group_by(TreeID) %>% summarise(DBH = dplyr::first(DBH), .groups = "drop")
d <- d %>% left_join(dbh_lut, by = "TreeID")
cat("trees with DBH:", sum(!is.na(unique(d[c("TreeID","DBH")])$DBH)),
    " | DBH range:", round(range(d$DBH, na.rm=TRUE),1), "cm\n")

# ===== 6) DBH vs FLUX ========================================================
p6 <- ggplot(d %>% filter(!is.na(DBH)), aes(DBH, CH4_best.flux)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey70") +
  geom_point(aes(colour = factor(height_cm)), alpha = 0.8, size = 1.9) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.6) +
  scale_y_continuous(trans = "asinh", breaks = c(0,1,5,25,100,400)) +
  labs(x = "DBH (cm, interpolated to 130 cm)", y = ch4_lab, colour = "height (cm)",
       title = "CH4 vs DBH", subtitle = "do larger-diameter stems emit more?") +
  theme_b + theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
save_both(p6, "06_CH4_vs_DBH", w = 7, h = 5.5)

# ===== 7) MIXED MODEL with DBH + species (random) ============================
if (have_lme4) {
  suppressMessages(library(lme4))
  dm <- d %>% filter(!is.na(DBH)) %>%
    mutate(aCH4 = asinh(CH4_best.flux), hz = scale(height_cm)[,1],
           dbhz = scale(DBH)[,1], todz = scale(tod)[,1])
  m <- try(lmer(aCH4 ~ hz + dbhz + todz + transect + (1|species/TreeID), data = dm), silent = TRUE)
  if (inherits(m, "try-error") || isSingular(m))
    m <- lmer(aCH4 ~ hz + dbhz + todz + transect + (1|species) + (1|TreeID), data = dm)
  fe <- summary(m)$coefficients
  cat("\n=== mixed model: asinh(CH4) ~ height + DBH + ToD + transect + (1|species/tree) ===\n")
  print(round(fe, 3))
  cf <- data.frame(term = rownames(fe), est = fe[,1], se = fe[,2]) %>% filter(term != "(Intercept)") %>%
    mutate(term = recode(term, hz="height (per SD)", dbhz="DBH (per SD)", todz="time of day (per SD)",
                         `transectN–S`="N–S vs E–W", `transectNE–SW`="NE–SW vs E–W", `transectNW–SE`="NW–SE vs E–W"),
           lo = est-1.96*se, hi = est+1.96*se)
  p7 <- ggplot(cf, aes(est, reorder(term, est))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = "grey50") +
    geom_point(size = 2.6, colour = "#1f4e79") +
    labs(x = "effect on asinh(CH4 flux)  (±95% CI)", y = NULL,
         title = "Adjusted effects on CH4 (+ DBH, species random effect)",
         subtitle = "height + DBH + time-of-day + transect; (1|species/tree)") +
    theme_b + theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
  save_both(p7, "07_CH4_adjusted_effects_DBH", w = 8, h = 4.8)
}

# ===== 8) REMEASURE: climbed (1st) vs basic (remeasure) ======================
cl <- read.csv(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-trees-2026/results/transect_trees_2026_fluxes.csv"),
  stringsAsFactors = FALSE, encoding = "UTF-8") %>%
  mutate(TreeID = sub("^T-", "", TreeID))
hbin <- function(h) cut(h, c(0,55,120,9999), labels = c("base (~40)","80 cm","160 cm"))
clb <- cl %>% filter(height_cm <= 200) %>% transmute(TreeID, bin = hbin(height_cm),
          CH4_climb = CH4_best.flux, CO2_climb = CO2_best.flux)
bab <- d   %>% filter(height_cm <= 200) %>% transmute(TreeID, bin = hbin(height_cm),
          CH4_basic = CH4_best.flux, CO2_basic = CO2_best.flux)
rem <- inner_join(clb, bab, by = c("TreeID","bin"))
cat("\nremeasure-matched points:", nrow(rem), "across", length(unique(rem$TreeID)), "trees\n")
lim <- range(c(rem$CH4_climb, rem$CH4_basic), na.rm = TRUE)
p8 <- ggplot(rem, aes(CH4_climb, CH4_basic, colour = TreeID, shape = bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2.8) +
  scale_x_continuous(trans = "asinh", breaks = c(0,1,5,25,100,400)) +
  scale_y_continuous(trans = "asinh", breaks = c(0,1,5,25,100,400)) +
  labs(x = expression("CH"[4]*" climbed (1st)  nmol m"^-2*" s"^-1),
       y = expression("CH"[4]*" basic (remeasure)"), colour = "tree", shape = "height",
       title = "Remeasured trees: climbed vs basic CH4",
       subtitle = "dashed = 1:1; 6 shared trees at matched heights") +
  theme_b + theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
save_both(p8, "08_remeasure_climb_vs_basic", w = 7.5, h = 6)

# ===== 9) FLUX BY HEIGHT, LINEAR scale, RAW-scale (NLS) models ===============
# Fit power/exp/linear all on the RAW response so (a) they track the arithmetic
# means on a linear axis and (b) their AICs are directly comparable.
g  <- d %>% filter(is.finite(CH4_best.flux))
gp <- g %>% filter(CH4_best.flux > 0)
mp <- g %>% group_by(Height_m) %>%
  summarise(mu = mean(CH4_best.flux), se = sd(CH4_best.flux)/sqrt(n()), .groups = "drop")
hg <- seq(min(g$Height_m), max(g$Height_m), length.out = 160)
lp <- lm(log(CH4_best.flux) ~ log(Height_m), gp)   # log fits -> NLS start values
le <- lm(log(CH4_best.flux) ~ Height_m, gp)
ctl <- nls.control(maxiter = 1000, warnOnly = TRUE)
m_lin <- lm(CH4_best.flux ~ Height_m, g)
m_pow <- try(nls(CH4_best.flux ~ a*Height_m^b,      g, start = list(a = exp(coef(lp)[1]), b = coef(lp)[2]), control = ctl), silent = TRUE)
m_exp <- try(nls(CH4_best.flux ~ a*exp(k*Height_m), g, start = list(a = exp(coef(le)[1]), k = coef(le)[2]), control = ctl), silent = TRUE)
fl <- list("linear  a+bh" = coef(m_lin)[1] + coef(m_lin)[2]*hg)
if (!inherits(m_pow,"try-error")) fl[["power law  a·h^b (NLS)"]]    <- coef(m_pow)["a"]*hg^coef(m_pow)["b"]
if (!inherits(m_exp,"try-error")) fl[["exponential  a·e^(kh) (NLS)"]] <- coef(m_exp)["a"]*exp(coef(m_exp)["k"]*hg)
fits <- bind_rows(lapply(names(fl), function(nm) data.frame(Height_m = hg, fit = fl[[nm]], model = nm)))
aics <- c(linear = AIC(m_lin))
if (!inherits(m_pow,"try-error")) aics["power"]       <- AIC(m_pow)
if (!inherits(m_exp,"try-error")) aics["exponential"] <- AIC(m_exp)
cat("\nCH4 ~ height RAW-scale model AIC (comparable; lower=better):\n"); print(round(aics,1))
sub9 <- paste0("raw-scale NLS fits (track the means); AIC  ",
               paste(sprintf("%s %.0f", names(aics), aics), collapse = "  "),
               "  -> ", names(which.min(aics)), " best")
cols9 <- c("power law  a·h^b (NLS)"="#2E8B57","exponential  a·e^(kh) (NLS)"="#D2691E","linear  a+bh"="#6A5ACD")
p9 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey70") +
  geom_point(data = g, aes(Height_m, CH4_best.flux), colour = "grey70", alpha = 0.5, size = 1.3) +
  geom_line(data = fits, aes(Height_m, fit, colour = model), linewidth = 1) +
  geom_errorbar(data = mp, aes(Height_m, ymin = mu-se, ymax = mu+se), width = 0.04, colour = "black") +
  geom_point(data = mp, aes(Height_m, mu), size = 3, colour = "black") +
  scale_colour_manual(values = cols9) +
  labs(x = "Height (m)", y = ch4_lab, colour = "model",
       title = "CH4 by height — linear scale, raw-scale (NLS) models", subtitle = sub9) +
  coord_flip(ylim = c(min(g$CH4_best.flux), max(mp$mu+mp$se)*1.1)) + theme_b +
  theme(plot.subtitle = element_text(size = 8, colour = "grey30"), legend.position = "top")
save_both(p9, "09_CH4_by_height_linear_models", w = 7.5, h = 6)

message("\nSaved plots 06..09 to: ", plots_dir)
