# =============================================================================
# water_sonde_profile.R
# Multiparameter sonde depth profile (Tower water column, 5/27/26 ~6pm), from
# the transcribed datasheet data/datasheets/parsed/IMG_3173.csv. Plots each
# parameter vs depth (surface at top). Transcription is uncertain in places:
# values flagged "?" are drawn as open/grey points, and the ORP column (logged
# as a drifting range, e.g. "169-134?") is shown as its range with a midpoint.
# Run with Rscript.  -> results/figs/water_sonde_profile.{pdf,png}
# =============================================================================

suppressMessages({ library(dplyr); library(tidyr); library(stringr); library(ggplot2) })
base <- "/Users/jongewirtzman/My Drive/Research/mach4-trees"
raw  <- read.csv(file.path(base, "data/datasheets/parsed/IMG_3173.csv"),
                 skip = 1, stringsAsFactors = FALSE, na.strings = c("", "NA"),
                 check.names = FALSE, fill = TRUE)
raw$Depth_cm <- suppressWarnings(as.numeric(raw$Depth_cm)) # was character (mis-aligned blank rows)
raw  <- raw[!is.na(raw$Depth_cm), ]                        # drop blank (no-data) rows

# --- parse one column -> long rows with value, range (lo/hi), uncertain -------
# extract UNSIGNED numbers: the "-" in an ORP range ("75-104") is a separator,
# not a minus sign, and every parameter on this sheet is non-negative
nums <- function(x) as.numeric(str_extract_all(x, "\\d+\\.?\\d*")[[1]])
parse_col <- function(col, label) {
  do.call(rbind, lapply(seq_len(nrow(raw)), function(i) {
    s <- as.character(raw[[col]][i]); n <- nums(s); n <- n[is.finite(n)]
    if (!length(n)) return(NULL)
    data.frame(Depth_cm = raw$Depth_cm[i], parameter = label,
               value = mean(n), lo = min(n), hi = max(n),
               uncertain = grepl("\\?", s), stringsAsFactors = FALSE)
  }))
}
# facet order + readable labels. ORP rows logged as "X-X" are meaned by parse_col
# (value = mean of the numbers); uS_cm is converted to practical salinity below.
specs <- tribble(
  ~col,            ~label,
  "pH",            "pH",
  "mV_ORP",        "ORP (mV)",
  "DO_ppm",        "DO (ppm)",
  "uS_cm",         "Salinity (PSU)",
  "FNU_turbidity", "Turbidity (FNU)",
  "degC",          "Temp (C)",
  "DO_pct",        "DO (%)"
)
prof <- bind_rows(Map(parse_col, specs$col, specs$label))
prof$parameter <- factor(prof$parameter, levels = specs$label)

# specific conductance (uS/cm at 25C) -> practical salinity (PSU); near-fresh
# here, so linear seawater ratio (35 PSU = 53087 uS/cm at 25C) is fine
ec2psu <- function(uScm) uScm / 53087 * 35
prof <- prof %>% mutate(
  value = ifelse(parameter == "Salinity (PSU)", ec2psu(value), value),
  lo    = ifelse(parameter == "Salinity (PSU)", ec2psu(lo),    lo),
  hi    = ifelse(parameter == "Salinity (PSU)", ec2psu(hi),    hi))

# force each panel to span at least one whole unit (Temp shown as-is)
lims <- tribble(
  ~parameter,        ~value,
  "pH",              6,    "pH",              7,
  "ORP (mV)",       -200,  "ORP (mV)",        200,
  "DO (ppm)",        0,    "DO (ppm)",        1,
  "Salinity (PSU)",  0,    "Salinity (PSU)",  5,
  "Turbidity (FNU)", 0,    "Turbidity (FNU)", 20,
  "DO (%)",          0,    "DO (%)",          10
) %>% mutate(parameter = factor(parameter, levels = specs$label), Depth_cm = 0)

p <- ggplot(prof, aes(value, Depth_cm)) +
  geom_blank(data = lims, aes(value, Depth_cm), inherit.aes = FALSE) +
  geom_path(aes(group = parameter), colour = "grey55", linewidth = 0.5) +
  geom_point(aes(fill = uncertain), shape = 21, size = 2.6, colour = "grey20") +
  facet_wrap(~ parameter, scales = "free_x", nrow = 2) +
  scale_y_reverse(breaks = seq(0, 275, 50)) +
  scale_fill_manual(values = c("FALSE" = "#1f6fb4", "TRUE" = "white"),
                    labels = c("FALSE" = "transcribed", "TRUE" = "uncertain (?)"),
                    name = NULL) +
  labs(x = NULL, y = "Depth (cm)") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), legend.position = "top")

figdir <- file.path(base, "results/figs"); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(figdir, "water_sonde_profile.pdf"), p, width = 10, height = 6, device = cairo_pdf)
ggsave(file.path(figdir, "water_sonde_profile.png"), p, width = 10, height = 6, dpi = 200, device = ragg::agg_png)
cat("Saved water_sonde_profile.{pdf,png} to", figdir, "\n")
cat("Depths:", paste(range(prof$Depth_cm), collapse = "-"), "cm |",
    length(unique(prof$parameter)), "parameters\n")
