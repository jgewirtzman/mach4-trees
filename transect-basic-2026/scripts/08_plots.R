# =============================================================================
# 08_plots.R  (BASIC trees)
# Exploratory flux plots: flux-by-height, species x height, slope-by-species,
# transect/time-of-day (model-adjusted), plus CH4-vs-CO2. CH4 nmol m-2 s-1,
# CO2 umol m-2 s-1. Run with Rscript.
# =============================================================================

source(file.path(
  "/Users/jongewirtzman/My Drive/Research/mach4-trees/transect-basic-2026/scripts/00_setup.R"))
suppressMessages({ library(ggplot2); library(dplyr); library(scales) })
try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
have_lme4 <- requireNamespace("lme4", quietly = TRUE)

d <- read.csv(file.path(results_dir, "transect_basic_2026_fluxes_with_mdf.csv"),
              stringsAsFactors = FALSE, encoding = "UTF-8")

# --- Consolidate species spelling variants -----------------------------------
canon <- c(
  "Araca Peua"="Araçá Peua","Araça Peua"="Araçá Peua","Araçá Peua"="Araçá Peua",
  "Macuru"="Macucú","Macucú"="Macucú",
  "Mata Folha Miuda"="Mata mata folha miúda","Mata Mata Folha Miuda"="Mata mata folha miúda","Mata Mata folha miuda"="Mata mata folha miúda",
  "Uruçurana"="Urucurana","Urucurana"="Urucurana",
  "Uxirana de tuxinho"="Uxirana","Uxirana de uxinho"="Uxirana",
  "Piranhieira"="Piranheira","Piranheira (struck through)"="Piranheira","Piranheira"="Piranheira")
d$species <- enc2utf8(ifelse(is.na(canon[d$species]) , d$species, canon[d$species]))

# --- Covariates --------------------------------------------------------------
arm <- gsub("[0-9].*", "", d$TreeID)
d$transect <- dplyr::recode(arm, E="E–W", W="E–W", N="N–S", S="N–S",
                            NE="NE–SW", SW="NE–SW", NW="NW–SE", SE="NW–SE")
d$Height_m <- d$height_cm / 100
d$tod <- as.numeric(format(as.POSIXct(d$start.time, tz = tz_data), "%H")) +
         as.numeric(format(as.POSIXct(d$start.time, tz = tz_data), "%M"))/60
ch4_lab <- expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
co2_lab <- expression(CO[2]~flux~(mu*mol~m^{-2}~s^{-1}))
sp_n <- d %>% group_by(species) %>% summarise(trees=n_distinct(TreeID), .groups="drop")
multi <- sp_n$species[sp_n$trees >= 2]

save_both <- function(p,name,w=7.5,h=6){
  ggsave(file.path(plots_dir,paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf)
  ggsave(file.path(plots_dir,paste0(name,".png")),p,width=w,height=h,dpi=200,device=ragg::agg_png)
}
asinh_x <- scale_x_continuous(trans="asinh", breaks=c(-2,0,2,10,50,200,500))
theme_b <- theme_classic(base_size=12)

# --- power-law fit on all positive points (flux = a*h^b) ---------------------
powfit <- function(df,yv){ g<-df[is.finite(df[[yv]])&df[[yv]]>0&df$Height_m>0,]
  if(nrow(g)<3) return(NULL); m<-lm(log(g[[yv]])~log(g$Height_m))
  hg<-exp(seq(log(min(g$Height_m)),log(max(g$Height_m)),length.out=80))
  data.frame(Height_m=hg, fit=exp(coef(m)[1]+coef(m)[2]*log(hg))) }

# ===== 1) FLUX BY HEIGHT (vertical) ==========================================
mk_height <- function(yv,ylab,col,use_asinh=FALSE){
  mp <- d %>% group_by(Height_m) %>%
    summarise(mu=mean(.data[[yv]],na.rm=T), se=sd(.data[[yv]],na.rm=T)/sqrt(n()), .groups="drop")
  fit <- powfit(d,yv)
  p <- ggplot(d, aes(Height_m, .data[[yv]])) +
    geom_hline(yintercept=0,linetype="dashed",linewidth=0.3,colour="grey60") +
    geom_line(aes(group=TreeID),colour="grey80",linewidth=0.3,alpha=0.6) +
    geom_point(colour=col,alpha=0.35,size=1.4) +
    { if(!is.null(fit)) geom_line(data=fit,aes(Height_m,fit),colour=col,linewidth=1) } +
    geom_errorbar(data=mp,aes(Height_m,ymin=mu-se,ymax=mu+se),width=0.04,inherit.aes=F,colour="black") +
    geom_point(data=mp,aes(Height_m,mu),inherit.aes=F,size=3,colour="black") +
    labs(x="Height (m)",y=ylab,title=paste(deparse(substitute(yv)))) + theme_b + coord_flip()
  if(use_asinh) p <- p + scale_y_continuous(trans="asinh",breaks=c(0,1,5,25,100,400))
  p }
save_both(mk_height("CH4_best.flux",ch4_lab,"#2E8B57",TRUE)+labs(title="CH4 flux by height"),"01_CH4_by_height")
save_both(mk_height("CO2_best.flux",co2_lab,"#D2691E")+labs(title="CO2 flux by height"),"01_CO2_by_height")

# ===== 2) SPECIES x HEIGHT (multi-tree species) ==============================
ds <- d %>% filter(species %in% multi)
p2 <- ggplot(ds, aes(Height_m, CH4_best.flux)) +
  geom_hline(yintercept=0,linetype="dashed",linewidth=0.3,colour="grey60") +
  geom_line(aes(group=TreeID),colour="grey70",linewidth=0.3) +
  geom_point(aes(colour=species),size=1.6,show.legend=FALSE) +
  facet_wrap(~species, scales="free_x") +
  scale_y_continuous(trans="asinh",breaks=c(0,1,5,25,100,400)) +
  labs(x="Height (m)",y=ch4_lab,title="CH4 by height, per species (>=2 trees)") +
  theme_b + theme(strip.background=element_blank(),strip.text=element_text(size=8)) + coord_flip()
save_both(p2,"02_CH4_species_x_height",w=10,h=8)

# ===== 3) SLOPE BY SPECIES (log-linear CH4 decline) ==========================
slopes <- lapply(multi, function(s){ g<-d[d$species==s & d$CH4_best.flux>0 & is.finite(d$CH4_best.flux),]
  if(nrow(g)<4 || length(unique(g$Height_m))<2) return(NULL)
  m<-lm(log(CH4_best.flux)~Height_m,data=g); ci<-confint(m)[2,]
  data.frame(species=s,n=nrow(g),slope=coef(m)[2],lo=ci[1],hi=ci[2]) })
sl <- bind_rows(slopes) %>% arrange(slope) %>% mutate(species=factor(species,levels=species))
p3 <- ggplot(sl, aes(slope,species)) +
  geom_vline(xintercept=0,linetype="dashed",colour="grey60") +
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=0.25,colour="grey50") +
  geom_point(aes(size=n),colour="#2E8B57") +
  labs(x="CH4 height-decline slope  d[ln CH4]/d[height_m]  (more negative = steeper)",
       y=NULL,size="n meas",title="CH4 vertical-decline slope by species") + theme_b
save_both(p3,"03_slope_by_species",w=8,h=5.5)

# ===== 4) TRANSECT / TIME-OF-DAY (model-adjusted) ============================
# raw descriptive
p4a <- ggplot(d, aes(transect, CH4_best.flux, fill=transect)) +
  geom_boxplot(outlier.size=0.6,alpha=0.7,show.legend=FALSE) +
  scale_y_continuous(trans="asinh",breaks=c(0,1,5,25,100,400)) +
  labs(x=NULL,y=ch4_lab,title="CH4 by transect direction (raw)") + theme_b
save_both(p4a,"04a_CH4_by_transect_raw",w=6,h=5)
p4b <- ggplot(d, aes(tod, CH4_best.flux)) +
  geom_point(aes(colour=factor(height_cm)),alpha=0.7) +
  geom_smooth(method="loess",se=TRUE,colour="black",linewidth=0.7) +
  scale_y_continuous(trans="asinh",breaks=c(0,1,5,25,100,400)) +
  labs(x="time of day (h)",y=ch4_lab,colour="height (cm)",title="CH4 vs time of day (raw)") + theme_b
save_both(p4b,"04b_CH4_by_timeofday_raw",w=7,h=5)

# model-adjusted: asinh(CH4) ~ height + tod + transect + (1|tree), controlling for height
if (have_lme4) {
  suppressMessages(library(lme4))
  dm <- d %>% mutate(aCH4=asinh(CH4_best.flux), hz=scale(height_cm)[,1], todz=scale(tod)[,1])
  m <- lmer(aCH4 ~ hz + todz + transect + (1|TreeID), data=dm)
  fe <- summary(m)$coefficients
  cf <- data.frame(term=rownames(fe), est=fe[,1], se=fe[,2]) %>% filter(term!="(Intercept)") %>%
    mutate(term=recode(term, hz="height (per SD)", todz="time of day (per SD)",
                       `transectN–S`="N–S vs E–W", `transectNE–SW`="NE–SW vs E–W", `transectNW–SE`="NW–SE vs E–W"),
           lo=est-1.96*se, hi=est+1.96*se)
  p4c <- ggplot(cf, aes(est, reorder(term,est))) +
    geom_vline(xintercept=0,linetype="dashed",colour="grey60") +
    geom_errorbarh(aes(xmin=lo,xmax=hi),height=0.2,colour="grey50") + geom_point(size=2.6,colour="#1f4e79") +
    labs(x="effect on asinh(CH4 flux)  (±95% CI)", y=NULL,
         title="Adjusted effects on CH4 (mixed model, tree random effect)",
         subtitle="height + time-of-day + transect; transect & ToD controlled for height & tree") +
    theme_b + theme(plot.subtitle=element_text(size=8,colour="grey30"))
  save_both(p4c,"04c_CH4_adjusted_effects",w=8,h=4.5)
  cat("\nMixed model fixed effects:\n"); print(round(fe,3))
}

# ===== bonus) CH4 vs CO2 =====================================================
p5 <- ggplot(d, aes(CO2_best.flux, CH4_best.flux, colour=factor(height_cm))) +
  geom_hline(yintercept=0,linetype="dashed",linewidth=0.3,colour="grey70") +
  geom_point(alpha=0.75,size=1.8) +
  scale_y_continuous(trans="asinh",breaks=c(0,1,5,25,100,400)) +
  labs(x=co2_lab,y=ch4_lab,colour="height (cm)",
       title="CH4 vs CO2 per measurement",subtitle="do higher-respiration stems emit more CH4?") +
  theme_b + theme(plot.subtitle=element_text(size=8,colour="grey30"))
save_both(p5,"05_CH4_vs_CO2",w=7,h=5.5)

message("\nSaved plots 01..05 to: ", plots_dir)
message("species with >=2 trees: ", length(multi), " | slopes computed: ", nrow(sl))
