# Title of Manuscript: Keeping up with the heat: Long-term dynamics and 
# plasticity of heat tolerance in a tropical plant community

# Author of the code: Georgia G. Hernández
# Last update: 2026/02/17
 
#-------------------------------------------------------------------------------
# This script is meant to include the graphs and analyses for the manuscript 
# about heat tolerance across years at La Selva, Costa Rica. 
# This Rscript includes the data of T50 from La Selva of years 2017-2023.

#-------------------------------------------------------------------------------
# Read packages

library(tidyverse)
library(tidytext)
library(lme4)
library(lubridate)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(gtsummary)
library(skimr)
library(effects)
library(lemon)
library(gt)
library(gtExtras)
library(ggpmisc)
library(ggrepel)

#-------------------------------------------------------------------------------
# Read Data

# Import T50 data & temperature
t50_sum <- readxl::read_xlsx("./data/t50_temp_summary.xlsx")

#-------------------------------------------------------------------------------
# Questions answered:

# Q1. Does heat tolerance of the Zingiberales community vary across years? 
# Q2. Do species differ in their inter-annual plasticity, and does this 
# difference alter their relative heat tolerance rankings from year to year? 
# Q3. How does heat tolerance vary from year to year in response to local 
# temperature changes? 

#-------------------------------------------------------------------------------
# Data wrangling

# Data needs some managing to exclude species that are not repeated through the 
# years, some changes in the names, and turning the columns into readable 
# variables for the model and graphs.

# Turn years and T50 as numeric
t50_sum$year <- as.numeric(as.character(t50_sum$year))
t50_sum$T50 <- as.numeric(as.character(t50_sum$T50))

# Filter out species added throughout the years that were not in the first data collection
zingiberales_traits <- t50_sum |> 
  dplyr::filter(!(species %in% c("musa_acuminata", "musa_paradisiaca",
                                 "goeppertia_micans_whitemorph", "curcuma_longa",
                                 "zingiber_officinale", "ctenanthe_palustris",
                                 "hellenia_speciosa", "goeppertia_leucostachys",
                                 "canna_indica", "kaempferia_rotunda", 
                                 "heliconia_psittacorum", "goeppertia_hammeli")))

# Change names of G. micans green morphotype to only G. micans
zingiberales_traits <- zingiberales_traits |> 
  dplyr::mutate(species = recode(species, 
                                 "goeppertia_micans_greenmorph" = "goeppertia_micans"))

# Centering year at the initial period (all years minus first year)
# To make intercept interpretable
zingiberales_traits <- zingiberales_traits |> 
  mutate(year_c = (year -  min(year)))

# Set the column "file_date" as a date fromat with lubridate
zingiberales_traits$file_date <- as.Date(as.character(zingiberales_traits$file_date),
                                         format = "%Y%m%d")

# I will use years as categorical because the effect of year is not linear.
zingiberales_traits <- zingiberales_traits |>
  mutate(year_cat = as.factor(year))

# Centering variables to make intercept more meaningful
zingiberales_traits <- zingiberales_traits |> 
  mutate(
    Tmax_30d_day = 
      max_max_air_temp_month_daylight - mean(max_max_air_temp_month_daylight),
    Tmax_15d_day = 
      max_max_air_temp_15days_daylight - mean(max_max_air_temp_15days_daylight),
    ATmax_15d_day = 
      avg_max_air_temp_15days_daylight - mean(avg_max_air_temp_15days_daylight),
    Tmax_3d_day = 
      max_max_air_temp_3days_daylight - mean(max_max_air_temp_3days_daylight),
    Tmax_y_day =
      air_temp_max_yesterday_daylight - mean(air_temp_max_yesterday_daylight),
    Tmax_c = 
      air_temperature_max_collectingtime - mean(air_temperature_max_collectingtime),
    Tmax_30d_dark = 
      max_max_air_temp_month_darkness - mean(max_max_air_temp_month_darkness),
    Tmax_15d_dark = 
      max_max_air_temp_15days_darkness - mean(max_max_air_temp_15days_darkness),
    ATmax_15d_dark = 
      avg_max_air_temp_15days_darkness - mean(avg_max_air_temp_15days_darkness),
    Tmax_3d_dark = 
      max_max_air_temp_3days_darkness - mean(max_max_air_temp_3days_darkness),
    Tmax_y_dark =
      air_temp_max_yesterday_darkness - mean(air_temp_max_yesterday_darkness)
  ) |> 
  dplyr::select(
    species, genus, family, T50, year, year_c, 
    year_cat, file_date,
    Tmax_30d_day, Tmax_15d_day, ATmax_15d_day, Tmax_3d_day,  
    Tmax_y_day, Tmax_c, Tmax_30d_dark, Tmax_15d_dark, 
    ATmax_15d_dark, Tmax_3d_dark, Tmax_y_dark)

zingiberales_traits <- zingiberales_traits |>
  mutate(year_cat = as.factor(year))

# ------------------------------------------------------------------------------
# Create colors and labels used throughout graphs

# Create a named vector to map families to specific colors
family_colors <- c(
  "cannaceae"     = "#B0D4E4",   
  "costaceae"     = "#3B75AF",   
  "heliconiaceae" = "#BDE293",   
  "marantaceae"   = "#5A9A41",   
  "musaceae"      = "#F4A09C",   
  "zingiberaceae" = "#C82525")   

# Labels for families in graphs
fam_label <- c("cannaceae"     = "Cannaceae",
               "costaceae"     = "Costaceae",
               "heliconiaceae" = "Heliconiaceae",
               "marantaceae"   = "Marantaceae",
               "musaceae"      = "Musaceae",
               "zingiberaceae" = "Zingiberaceae")

# Change species names to a format: from Goeppertia micans to G. micans
my_labeller <- function(str) {
  out <- strsplit(str, split = "_")
  out <- sapply(out, function(x) {
    paste0(toupper(substr(x[1], 1, 1)),
           ". ",
           x[2])
  })
  return(out)
}

zingiberales_traits4 <- zingiberales_traits |> 
  dplyr::select(species:year)

# Turn dates into numbers with the same distance, just to vizualize data better
zingiberales_traits4 <- zingiberales_traits4 |> 
  mutate(year = if_else(year == 2017, 0, year)) |> 
  mutate(year = if_else(year == 2021, 10, year)) |> 
  mutate(year = if_else(year == 2022, 20, year)) |> 
  mutate(year = if_else(year == 2023, 30, year))

# New facet names for families
fam_names <- c("Cannaceae", "Costaceae", "Heliconiaceae",
               "Marantaceae", "Musaceae", "Zingiberaceae")

names(fam_names) <- c("cannaceae", "costaceae", "heliconiaceae",
                      "marantaceae", "musaceae", "zingiberaceae")

#-------------------------------------------------------------------------------
# Analyses & Graphs by questions

# General Model used:
set.seed(2025)

mod1 <- 
  lmer(T50 ~ year_cat + 
         (1 | family) + #Random intercept and slope for family
         (1 + year_cat|species:family), # species is nested within family
       #"year_cat|" is for species nested within family to have an estimate for 
       # each year & be able to answer Q2
       data = zingiberales_traits,
       control = lmerControl(optimizer = "optimx", 
                             optCtrl = list(method = "nlminb"))) # here the 
       # optimization algorithm is changed because it had convergence issues

summary(mod1)

# Q1. --------------------------------------------------------------------------
# Q1. Does heat tolerance of the Zingiberales community vary across years? 

anova(mod1) 
fixef(mod1)
# T50 is significantly different over the years

# Confidence interval, aka average increase in T50 by year
tbl_regression(mod1)

ranova(mod1) 
# Pr(>Chisq) 2.215e-13  ***
# T50 of Zingiberales species is significantly different over the years

# Average increase of T50 by each species per year comparison:
# T50 increase from 2017 to 2021
cbind.data.frame(
  "species" = rownames(ranef(mod1)[["species:family"]]),
  "T50_inc" = fixef(mod1)[2] + ranef(mod1)[["species:family"]][, 2]
) |>
  arrange(T50_inc)

# T50 increase from 2017 to 2022
cbind.data.frame(
  "species" = rownames(ranef(mod1)[["species:family"]]),
  "T50_inc" = fixef(mod1)[3] + ranef(mod1)[["species:family"]][, 3]
) |>
  arrange(T50_inc)

# T50 increase from 2017 to 2023
cbind.data.frame(
  "species" = rownames(ranef(mod1)[["species:family"]]),
  "T50_inc" = fixef(mod1)[4] + ranef(mod1)[["species:family"]][, 4]
) |>
  arrange(T50_inc)

# The random intercepts by species of the model
dotplot.ranef.mer(ranef(mod1))


# Check model:
plot(mod1)
hist(residuals(mod1))
qqnorm(residuals(mod1))
abline(a = 0, b = 1, col = 2) 

# Graph 1.----------------------------------------------------------------------
# Prediction intervals 
conf_dt <- zingiberales_traits |> 
  dplyr::select(year_cat, species, family) |> 
  distinct()

pred_sp <- function(x){
  predict(x, 
          newdata = conf_dt, 
          re.form = NULL)
}

# Collapse bootstrap into median, 95% PI
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, 
                           function(x) as.numeric(quantile(x, 
                                                           probs = .5, 
                                                           na.rm = TRUE))),
               lwr = apply(merBoot$t, 2, 
                           function(x) as.numeric(quantile(x, 
                                                           probs = .025, 
                                                           na.rm = TRUE))),
               upr = apply(merBoot$t, 2, 
                           function(x) as.numeric(quantile(x, 
                                                           probs = .975, 
                                                           na.rm = TRUE))),
               avg = apply(merBoot$t, 2, 
                           function(x) as.numeric(mean(x, na.rm = TRUE))),
               se = apply(merBoot$t, 2, 
                          function(x) as.numeric(sd(x, na.rm = TRUE)))
    )
  )
}

boot_pred <- lme4::bootMer(mod1, pred_sp, 
                           nsim = 1000, 
                           use.u = TRUE, 
                           type = "semiparametric",
                           parallel = "multicore",
                           ncpus = 3)

pred_ci <- data.frame(confint(boot_pred)) |>
  setNames(c("lwr", "upr")) |> 
  bind_cols(conf_dt) |> 
  mutate(species = as.character(species)) |> 
  mutate(estimate = predict(mod1, 
                            newdata = conf_dt, 
                            re.form = NULL))

# Filter species in more than 1 year
aux_species <- zingiberales_traits4 |> 
  count(year, species) |> 
  group_by(species) |> 
  mutate(mult_years = n() > 1) |> 
  ungroup()

aux_species <- aux_species |>
  filter(!mult_years) |> 
  pull(species)

T50_filtered <- zingiberales_traits4 |>
  filter(!species %in% aux_species) |> 
  drop_na()

T50_filtered_2 <- T50_filtered |> 
  mutate(year = as.character(as.numeric(year)))

# Summary of the data to get Standard Error
year_sp_sum <- T50_filtered_2 |> 
  drop_na() |>
  group_by(species, year) |>
  summarise(n = n(),
            mean_T50 = mean(T50),
            sd_T50 = sd(T50)) |>
  mutate(se_T50 = sd_T50/sqrt(n)) 

N_year <- zingiberales_traits4 |> 
  group_by(year) |> 
  summarise(n = n(), 
            mean_sp = mean(T50), 
            se_sp = sd(T50)) |> 
  ungroup() |> 
  mutate(se_sp = se_sp/sqrt(n)) |> 
  mutate(year_cat = year)

N_year <- N_year |> 
  mutate(year_cat = if_else(year == 0, 2017, year_cat)) |> 
  mutate(year_cat = if_else(year == 10, 2021, year_cat)) |> 
  mutate(year_cat = if_else(year == 20, 2022, year_cat)) |> 
  mutate(year_cat = if_else(year == 30, 2023, year_cat))

pred_ci <- pred_ci |> 
  mutate(species = as.factor(species))

# Graph
mult_sites <- ggplot(data = pred_ci) +
  geom_vline(data = N_year,
             aes(xintercept = mean_sp)) +
  geom_rect(data = N_year,
            aes(ymin = 0, ymax = 36,
                xmin = mean_sp - qnorm(.975) * se_sp,
                xmax = mean_sp + qnorm(.975) * se_sp),
            fill = "grey",
            inherit.aes = FALSE,
            alpha = .8) +
  geom_errorbar(aes(x = estimate, 
                    xmin = lwr,
                    xmax = upr, 
                    y = reorder_within(species, -lwr, year_cat)),
                width = 0, 
                alpha = 0.9, 
                linewidth = 0.9) +
  theme_classic() +
  scale_color_manual(labels = c("0" = "2017",
                                "10" = "2021",
                                "20" = "2022",
                                "30" = "2023")) +
  theme(axis.text.y = element_text(size = 14, face = "italic"),
        axis.text.x = element_text(size = 24),
        axis.title = element_text(size = 24),
        strip.text = element_text(size = 24),
        panel.grid.major.y = element_line(linetype = 2, linewidth = 0.5),
        strip.background = element_blank()) +
  labs(y = "Species",
       x = expression("T"[50]~"(°C)")) +
  guides(color = guide_legend(title = "Years:")) +
  scale_y_reordered(labels = my_labeller) +
  facet_wrap( ~ year_cat, nrow = 1, scales = "free_y")

mult_sites

# Q2. --------------------------------------------------------------------------
# Q2. Do species differ in their inter-annual plasticity, and does this 
# difference alter their relative heat tolerance rankings from year to year? 

# For this analysis I will use the Friedman Rank Sum test. From which I will be 
# able to tell whether some species have their location parameter consistently 
# higher than other (across years).

# From the R documentation: "The **null hypothesis** is that apart from an 
# effect of blocks (years), the location parameter of y is the same in each of 
# the groups (species)".

# Select all species that are only repeated in all years
species_filt <- split(zingiberales_traits, 
                      f = zingiberales_traits$year) |> 
  lapply(\(x) unique(pull(x, species)))

species_filt <- Reduce(intersect, species_filt)

test_dt <- filter(zingiberales_traits,
                  species %in% species_filt) |> 
  group_by(species, year, family) |> 
  summarise(T50 = mean(T50, na.rm = TRUE)) |> 
  ungroup()

# Ranking Test by species
testing_frid <- with(test_dt,
                     friedman.test(y = T50, 
                                   groups = species,
                                   blocks = year))
testing_frid # it's significant, therefore the ranking is constant over the years

frid_post <- with(test_dt,
                  PMCMRplus::frdAllPairsConoverTest(y = T50, 
                                                    groups = species,
                                                    blocks = year,
                                                    p.adj = "bon"))

test_fam <- test_dt |> 
  group_by(family, year) |> 
  summarise(T50 = mean(T50, na.rm = TRUE)) |> 
  ungroup()

# Ranking Test by family
testing_fam <- with(test_fam,
                    friedman.test(y = T50, 
                                  groups = family,
                                  blocks = year))
testing_fam # it's significant, the family ranking is constant over the years

# Based on the Friedman test results, the null hypothesis is rejected. This 
# mean that species are consistently in the bottom or top throughout the years. 

# Graph 2.----------------------------------------------------------------------

# Existing model and bootstrap code
conf_dt <- zingiberales_traits |>
  dplyr::select(year_cat, species, family) |>
  distinct()

pred_sp <- function(x){
  predict(x,
          newdata = conf_dt,
          re.form = NULL)
}

REsim2 <- 
  function(merMod, n.sims = 200, oddsRatio = FALSE, seed = NULL) {
    stopifnot(inherits(merMod, "merMod"))
    if (!is.null(seed)) 
      set.seed(seed)
    else if (!exists(".Random.seed", envir = .GlobalEnv)) 
      runif(1)
    mysim <- arm::sim(merMod, n.sims = n.sims)
    reDims <- length(mysim@ranef)
    tmp.out <- vector("list", reDims)
    names(tmp.out) <- names(mysim@ranef)
    for (i in c(1:reDims)) {
      zed <- apply(mysim@ranef[[i]], c(2, 3), function(x) as.data.frame(x) |> 
                     dplyr::summarise_all(.funs = 
                                            list("mean" = mean, 
                                                 "median" = median,
                                                 "lwr" = ~quantile(., .025),
                                                 "upr" = ~quantile(., .975),
                                                 "sd" = sd)), 
                   simplify = FALSE)
      zed <- do.call(rbind, zed)
      zed$X1 <- rep(dimnames(mysim@ranef[[i]])[[2]], length(dimnames(mysim@ranef[[i]])[[3]]))
      zed$X2 <- rep(dimnames(mysim@ranef[[i]])[[3]], each = length(dimnames(mysim@ranef[[i]])[[2]]))
      tmp.out[[i]] <- zed
      rm(zed)
      tmp.out[[i]]$groupFctr <- names(tmp.out)[i]
      tmp.out[[i]]$X1 <- as.character(tmp.out[[i]]$X1)
      tmp.out[[i]]$X2 <- as.character(tmp.out[[i]]$X2)
    }
    dat <- do.call(rbind, tmp.out)
    dat$groupID <- dat$X1
    dat$X1 <- NULL
    dat$term <- dat$X2
    dat$X2 <- NULL
    dat <- dat[, c("groupFctr", "groupID", "term", "mean", "median", 
                   "sd", "lwr", "upr")]
    rownames(dat) <- NULL
    if (oddsRatio == TRUE) {
      dat$median <- exp(dat$median)
      dat$mean <- exp(dat$mean)
      dat$sd <- NA
      return(dat)
    }
    else {
      return(dat)
    }
  }

sim_eff <- REsim2(mod1, n.sims = 1000, seed = 2025)

sim_eff <- sim_eff |> 
  filter(groupFctr == "species:family",
         term != "(Intercept)") |>
  mutate(species = sapply(str_split(groupID, pattern = ":"), \(x) x[1]),
         family = sapply(str_split(groupID, pattern = ":"), \(x) x[2]),
         year = substr(term, nchar(term) - 3, nchar(term))) |>
  as_tibble()

# Run the bootstrap
boot_pred <- lme4::bootMer(mod1, pred_sp,
                           nsim = 100,
                           use.u = TRUE,
                           type = "semiparametric",
                           parallel = "multicore",
                           ncpus = 12)

# Convert the bootstrap output to a tidy data frame, where each row will be 
# one simulation for one species in one year
tidy_boot <- as.data.frame(t(boot_pred$t)) |>
  bind_cols(conf_dt) |>
  pivot_longer(cols = starts_with("V"),
               names_to = "simulation",
               values_to = "predicted_t50")

# Separate the baseline (2017) predictions
baseline_preds <- tidy_boot |>
  filter(year_cat == 2017) |>
  dplyr::select(species, simulation, baseline_t50 = predicted_t50)

# Calculate the change from baseline for each simulation
change_from_baseline <- tidy_boot |>
  filter(year_cat != 2017) |> 
  left_join(baseline_preds, by = c("species", "simulation")) |>
  mutate(change_t50 = predicted_t50 - baseline_t50)

# Summarize the distribution of changes to get the median and 95% intervals
summary_plot_data <- change_from_baseline |>
  group_by(year_cat, species, family) |>
  summarize(
    estimate = median(change_t50, na.rm = TRUE),
    lwr = quantile(change_t50, 0.025, na.rm = TRUE),
    upr = quantile(change_t50, 0.975, na.rm = TRUE),
    .groups = 'drop')

# Graph
plot_change <- ggplot(sim_eff,
                      aes(x = median, 
                          y = reorder_within(species, -median, year), 
                          color = family)) +
  geom_vline(xintercept = 0,  # Add the vertical red line at zero
             linetype = "dashed", 
             color = "red", 
             alpha = 1) +
  geom_errorbarh(aes(xmin = lwr, # Add the error bars (using the bootstrapped intervals)
                     xmax = upr), 
                 height = 0.2, alpha = 1) +
  geom_point(size = 2) +
  facet_wrap(~ year, 
             scales = "free_y") + 
  scale_y_reordered() + # This makes the reordering within facets work
  scale_color_manual(values = family_colors,
                     labels = c(fam_label)) +
  labs(x = expression("Estimated changes T"[50]~"(°C)"),
       y = "Species",
       color = "Family:") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 8, face = "italic")) +
  scale_y_discrete(labels = my_labeller)

plot_change

# Graph 3.----------------------------------------------------------------------

# Calculate the mean T50 for each species in the baseline year (2017)
baseline_t50 <- zingiberales_traits |> 
  filter(year == 2017) |> 
  group_by(species) |> 
  summarise(T50_2017 = mean(T50, na.rm = TRUE), # Calculate mean T50 for 2017
            .groups = "drop")

# Calculate delta
delta_t50_data <- zingiberales_traits |>
  filter(year %in% c(2021, 2022, 2023)) |> # Focus on response years
  group_by(species, year) |>
  summarise(Mean_T50_Year = mean(T50, na.rm = TRUE),
            .groups = 'drop') |>
  left_join(baseline_t50, by = "species") |> # Join with 2017 baseline
  mutate(Delta_T50 = Mean_T50_Year - T50_2017) # Calculate change from baseline

delta_t50_data <- delta_t50_data |>
  left_join(zingiberales_traits |>
              select(species, family) |>
              distinct(), by = "species")

# Get one point per species, by averaging the delta T50 across 2021-2023
species_plot_data <- delta_t50_data |>
  group_by(species, T50_2017) |> # Group by species and its 2017 T50
  summarise(Avg_Delta_T50 = mean(Delta_T50, 
                                 na.rm = TRUE), # Average delta across years
            .groups = 'drop') |>
  left_join(zingiberales_traits |>
              select(species, family) |>
              distinct(), by = "species")

# Graph
# Create a new column with the formatted species names for plotting
species_plot_data <- species_plot_data |>
  mutate(formatted_species = my_labeller(species))

delta_allspp <- ggplot(species_plot_data, 
                       aes(x = T50_2017, y = Avg_Delta_T50)) +
  geom_smooth(method = "lm", 
              fill = "#d3d3d3",
              color = "black", 
              linetype = "solid", 
              se = TRUE) +
  geom_point(aes(color = family), 
             size = 5, 
             alpha = 1) +
  geom_text_repel(aes(label = formatted_species),   # Add text labels for species 
                  size = 4.5, 
                  max.overlaps = 40, 
                  segment.color = 'grey50') +
  geom_hline(yintercept = 0, 
             linetype = "dotted",
             color = "gray40") + # Add a horizontal line at 0 for Delta T50 (no change)
  labs(x = expression("T"[50]~"from baseline in 2017 (°C)"),
       y = expression(Average ~~ Delta* T[50] ~~ "from baseline in 2017 (°C)"),
       color = "Family:") +
  theme_classic(base_size = 22) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.title = element_text(face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(.14, .22)) +
  scale_color_manual(values = family_colors,
                     labels = c(fam_label))
delta_allspp

# Model that explains the variation in delta and T50 
trade_off_model <- lm(Avg_Delta_T50 ~ T50_2017, data = species_plot_data)
summary(trade_off_model)
# For every 1°C increase in baseline T50, the average delta T509 decreases by 0.40

# Coefficients
coefficients(trade_off_model)

# Confidence intervals for coefficients
confint(trade_off_model)

# Q3. --------------------------------------------------------------------------
# Q3. How does heat tolerance vary from year to year in response to local 
# temperature changes? 

my_dt <-
  t50_sum |>
  select(species:year, matches("temp")) |>
  pivot_longer(cols =
                 air_temp_max_yesterday_darkness:air_temperature_median_collectingtime) |>
  mutate(statistic = case_when(
    grepl("^avg_|_mean_", name) ~ "Mean",
    grepl("^max_|_max_", name) ~ "Max",
    TRUE ~ NA_character_
  ),
  daylight = ifelse(grepl("darkness$", name), "Darkness", "Daylight"),
  collect  = case_when(
    grepl("_collectingtime", name) ~ "0 days",
    grepl("_yesterday_", name) ~ "1 day",
    grepl("_3days_", name) ~ "3 days",
    grepl("_15days_", name) ~ "15 days",
    grepl("_month_", name) ~ "30 days",
    TRUE ~ NA_character_
  )) |>
  mutate(statistic = factor(statistic, levels = c("Min", "Max", 
                                                  "Median", "Mean")),
         daylight  = factor(daylight, levels = c("Daylight", "Darkness")),
         collect   = factor(collect, levels = c("0 days", 
                                                "1 day",
                                                "3 days",
                                                "15 days",
                                                "30 days"))) |>
  filter(!is.na(statistic), !is.na(collect))

# Graph 4.----------------------------------------------------------------------
temp_dn <- ggplot(data = my_dt,
                  aes(x = value, y = T50,
                      color = daylight)) +
  geom_point(alpha = 0.3, 
             size = 2) +
  geom_smooth(method = "lm", 
              se = FALSE,
              lwd = 1.1) +
  stat_poly_eq(size = 5,
               label.y = c(0.98, 0.98),
               label.x = c(0.53, 0.05), 
               hjust = 0) + # To add R squared to graph
  facet_grid(statistic ~ collect) +
  theme_bw(base_size = 22) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line(color = "gray90", 
                                          linetype = "dashed"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        strip.background = element_rect(fill = "white", 
                                        color = "white")) +
  labs(y = expression("T"[50]~"(°C)"),
       x = "Temperature (°C)",
       color = NULL) +
  scale_color_manual(values = c("#ffc242", "#2b2d42"), 
                     labels = c("Daylight" = "Daytime",
                                "Darkness" = "Nighttime"))
temp_dn

# Based on the R squared, the variable of the mean temperature from the last 
# 15 days, and daytime, is the best predictor to be use in the other models


# Checking if R2 are working
my_dt |>
  filter(collect == "15 days", statistic == "Mean",
         daylight == "Daylight") |>
  lm(T50 ~ value, data = _) |>
  summary()


#-------------------------------------------------------------------------------

# Supplementary Materials

# Graph S1.---------------------------------------------------------------------
# Model prediction intervals 
time_eff <- ggemmeans(mod1, terms = "year_cat")

time_eff <- time_eff |>  
  as.data.frame() 

time_eff <- time_eff |> 
  mutate(year = x) |> 
  dplyr::select(year, predicted, conf.low, conf.high)

time_eff <- time_eff |> 
  mutate(year = as.numeric(as.character(year)))

time_eff 

time_eff <- time_eff |> 
  mutate(year = if_else(year == 2017, 0, year)) |> 
  mutate(year = if_else(year == 2021, 10, year)) |> 
  mutate(year = if_else(year == 2022, 20, year)) |> 
  mutate(year = if_else(year == 2023, 30, year))

# Graph
T50vsYear <- ggplot() +
  geom_jitter(data = zingiberales_traits4,
              aes(x = year,
                  y = T50),
              alpha = 0.4,
              col = "#9999994c",
              size = 4,
              width = 1) +
  geom_linerange(data = time_eff, 
                 aes(x = year, 
                     ymin = conf.low, 
                     ymax = conf.high),
                 linewidth = .9, 
                 color = "black", 
                 alpha = 0.9) +
  geom_point(data = time_eff,
             aes(x = year, 
                 y = predicted),
             size = 3, 
             color = "red", 
             alpha = 1) +
  scale_x_continuous(breaks = c(0, 10, 20, 30), 
                     labels = c("2017", "2021", "2022", "2023")) +
  scale_y_continuous(breaks = seq(45, 55, 2)) +
  theme_classic() +
  labs(x = "Years", 
       y = expression("T"[50]~"(°C)")) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 25))

T50vsYear


# Graph S2.---------------------------------------------------------------------

T50vsfam <- ggplot() +
  geom_jitter(data = zingiberales_traits4, 
              aes(x = year, 
                  y = T50), 
              alpha = 0.4,  
              size = 3, 
              width = 1, 
              color = "grey60") +
  facet_rep_wrap(~ family,
                 labeller = labeller(family = fam_names),
                 repeat.tick.labels = 'top') +
  scale_x_continuous(breaks = c(0, 10, 20, 30), 
                     labels = c("2017", "2021", "2022", "2023")) +
  theme_classic() +
  labs(x = "Years", 
       y = expression("T"[50]~"(°C)")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 16), 
        strip.background = element_blank(), 
        axis.line=element_line())

T50vsfam

# Graph S3.---------------------------------------------------------------------
# Get maximum Delta T50 across 2021-2023
species_plot_data_max <- delta_t50_data |> 
  group_by(species, T50_2017) |>  # Group by species and its 2017 T50
  summarise(Max_Delta_T50 = max(Delta_T50, 
                                na.rm = TRUE), # Average delta across years
            .groups = 'drop') |>
  left_join(zingiberales_traits |> 
              select(species, family) |> 
              distinct(),
            by = "species")

# Graph
species_plot_data_max <- species_plot_data_max |>
  mutate(formatted_species = my_labeller(species))

max_delta_allspp <- ggplot(species_plot_data_max, 
                           aes(x = T50_2017, y = Max_Delta_T50)) +
  geom_smooth(method = "lm", fill = "#d3d3d3",
              color = "black", linetype = "solid", se = TRUE) +
  geom_point(aes(color = family), 
             size = 5, alpha = 1) +
  geom_text_repel(aes(label = formatted_species),   # Add text labels for species 
                  size = 4.5, 
                  max.overlaps = 40, 
                  segment.color = 'grey50') +
  geom_hline(yintercept = 0, 
             linetype = "dotted",
             color = "gray40") + # Add a horizontal line at 0 for Delta T50 (no change)
  labs(x = expression("T"[50]~"from baseline in 2017 (°C)"),
       y = expression(Maximum ~~ Delta* T[50] ~~ "from baseline in 2017 (°C)"),
       color = "Family:") +
  theme_classic(base_size = 22) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.title = element_text(face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(.14, .22)) +
  scale_color_manual(values = family_colors,
                     labels = c(fam_label))
max_delta_allspp

# Model that explains the variation in delta and T50 
max_trade_off_model <- lm(Max_Delta_T50 ~ T50_2017, 
                          data = species_plot_data_max)
summary(max_trade_off_model)
# For every 1°C increase in baseline T50, the average delta T509 decreases by 0.44

# Coefficients:
coefficients(max_trade_off_model)

# Confidence intervals for coefficients:
confint(max_trade_off_model)


# Extra graphs -----------------------------------------------------------------
# 1. Predictions of T50 by species based on the model used ---------------------
pred_graph <- ggplot(data = pred_ci,
                     aes(x = year_cat,
                         y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr,
                    ymax = upr)) +
  facet_wrap(family ~ species) +
  theme_classic() +
  labs(x = "Year", 
       y = expression("T"[50]~"(°C)")) +
  theme_classic() +
  theme(legend.position = "none", 
        strip.background = element_rect(color = "white"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 25))

pred_graph

# 2. Temperature estimate using the best model, showing the species effect -----

mod1_4 <- 
  lmer(T50 ~ year_cat + 
         (1 | family) +
         (1 + year_cat|species:family) +   
         Tmax_y_day,
       data = zingiberales_traits,
       control = lmerControl(optimizer = "optimx", 
                             optCtrl = list(method = "nlminb")))
summary(mod1_4)

# Get the original mean of Tmax_y_day to uncentered in plot:
# mean(zingiberales_traits$air_temp_max_yesterday_daylight)
# 31.273

# Graph facet years, legend species
spp_temp_y <- zingiberales_traits |>
  mutate(fit = predict(mod1_4, re.form = NULL)) |>
  ggplot(data = _,
         aes(x = Tmax_y_day, 
             y = T50, color = species)) +
  geom_point(pch = 16) +
  geom_line(aes(y = fit)) +
  theme_bw(base_size = 20) +
  labs(x = expression("T"[MAX_Y]~"(°C)"),
       y = expression("T"[50]~"(°C)"),
       color = "Species:") +
  facet_wrap(.~ year_cat) +
  theme(strip.background = element_blank(),
        legend.text = element_text(size = 10, face = "italic")) +
  scale_y_continuous(limits = c(47, 57), 
                     breaks = seq(47, 57, by = 2)) +
  scale_x_continuous(labels = function(x) x + 31.3) +  
  scale_color_discrete(labels = my_labeller)

spp_temp_y

# 3. Average change by year vs heat tolerance ----------------------------------
delta_years <- ggplot(delta_t50_data, 
                      aes(x = T50_2017, y = Delta_T50)) +
  geom_smooth(method = "lm", 
              color = "black", 
              linetype = "dashed", 
              se = TRUE) +
  geom_point(aes(color = family), 
             size = 3, 
             alpha = 0.8) +
  geom_hline(yintercept = 0, 
             linetype = "dotted", 
             color = "gray50") +
  labs(x = expression("T"[50]~"from baseline in 2017 (°C)"),
       y = expression(Average ~~ Delta* T[50] ~~ "from baseline in 2017 (°C)"),
       color = "Family:") +
  facet_wrap(~ year) +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.title = element_text(face = "bold"),
        legend.position = "right", 
        strip.text = element_text(size = 20),
        strip.background = element_blank()) +
  scale_color_manual(values = family_colors,
                     labels = c(fam_label))

delta_years

# 4. Average change by year vs heat tolerance by family ------------------------ 
delta_fam <- ggplot(species_plot_data, aes(x = T50_2017, y = Avg_Delta_T50)) +
  geom_smooth(method = "lm", 
              color = "black", 
              linetype = "dashed", 
              se = TRUE) +
  geom_point(aes(color = family), 
             size = 4, 
             alpha = 1) +
  geom_hline(yintercept = 0, 
             linetype = "dotted", 
             color = "gray50") +
  labs(x = expression("T"[50]~"from baseline in 2017 (°C)"),
       y = expression(Average ~~ Delta* T[50] ~~ "from baseline in 2017 (°C)"),
       color = "Family:") +
  facet_wrap(~ family,
             labeller = labeller(family = c(fam_label))) +
  theme_bw(base_size = 20) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.position = "none",
        strip.text = element_text(size = 20)) +
  scale_color_manual(values = family_colors,
                     labels = c(fam_label))

delta_fam

### END ------------------------------------------------------------------------
