library(tidyverse)
library(here)
library(readxl)
library(patchwork)

# ----- constants -----
r <- 8.314 # GAS CONSTANT
t <- 298 # TEMPERATURE IN K
mw <- 18.015e-6 # MOLECULAR WEIGHT OF H2O 

# parameters for processing the osmotic pressure data
time_step <- 5
n_lag <- 8
cutoff_time <- 80

# ----- theme for the ggplots 
theme_clean <- theme(
  axis.line = element_line(colour = "black", linewidth = 0.5),
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "lightgrey", linetype = "dashed"),
  legend.key = element_rect(fill = "transparent"),
  text = element_text(size = 20)
)

# ----- functions -----

# this function takes the raw data and organizes it into a more tidy format
process_data <- function(path) {
  raw <- read.table(here(path), header = FALSE, fill = TRUE)
  
  time <- raw %>%
    filter(is.na(V6)) %>%
    select(1:5) %>%
    set_names(c("time", "radius", "shells", "x", "species"))
  
  time_new <- time[rep(seq_len(nrow(time)), each = max(time$shells)), ]
  
  shelldata <- raw %>%
    filter(!is.na(V6))
  
  cbind(time_new, shelldata) %>%
    as.data.frame() %>%
    mutate(
      radius = radius * 1e4,
      x = x * 1e4,
      V2 = V2 * 1e4,
      time = time / 60)
}

# this function calculates the average osmotic pressure for the line plots
lineplots <- function(df, flatten_after_max = FALSE, extend_to = NULL) {
  
  df_summary <- df %>%
    group_by(time) %>%
    summarize(
      wa = mean(V3, na.rm = TRUE),
      wa_sd = sd(V3, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(osmotic = (-r * t) / mw * log(wa) / 101325 * 0.1,
      
      # propagate SD to osmotic
      osmotic_sd = abs((-r * t) / mw / 101325 * 0.1 * (wa_sd / wa)))
  
  # for the ones that cannot be calculated to the end of the experiment time
  if (flatten_after_max) {
    max_time <- df_summary$time[which.max(df_summary$osmotic)]
    max_osmotic <- max(df_summary$osmotic)
    
    df_summary <- df_summary %>%
      mutate(osmotic = ifelse(time > max_time, max_osmotic, osmotic))
    
    if (!is.null(extend_to) && extend_to > max(df_summary$time)) {
      extra_times <- seq(max(df_summary$time) + 1, extend_to, by = 1)
      df_summary <- bind_rows(df_summary, tibble(
        time = extra_times,
        wa = NA_real_,
        wa_sd = NA_real_,
        osmotic = max_osmotic,
        osmotic_sd = 0
      ))
    }
  }
  
  df_summary
}

# calculates delta Pi over delta t
dpidt <- function(df, step) {
  df %>%
    mutate(time = round(time)) %>%
    filter(time %% step == 0) %>%
    group_by(time) %>%
    summarize(
      wa = mean(wa),
      osmotic = mean(osmotic),
      .groups = "drop") %>%
    mutate(slope = (lead(osmotic) - osmotic) / (lead(time) - time))
}

# function for fitting the model
fit_model <- function(df) {
  m <- lm(log ~ dpidt + 0, data = df)
  coefficients(m)
}

# plots the linear models
plot_decay <- function(df, color, show_y = TRUE) {
  ggplot(df, aes(dpidt, log)) +
    geom_point(aes(shape = Media), color = color, size = 5, alpha = 0.5) +
    geom_errorbar(aes(ymin = log - sd, ymax = log + sd), width = 2) +
    geom_smooth(method = "lm", formula = y ~ 0 + x,
                se = TRUE, linetype = "dashed",
                color = "black", alpha = 0.2) +
    scale_shape_manual(
      values = c("AS" = 17, "NaCl" = 2),
      labels = c("AS", "5.4M NaCl")) +
    labs(
      x = bquote(frac(Delta*Pi~(MPa), Delta*t)),
      y = if (show_y) bquote(log(C/C[0])) else NULL) +
    ylim(-6, 1) +
    theme_clean +
    theme(text = element_text(size = 22))
}

# plots the predicted bacteria concentrations
bacteria_predict <- function(df, step, coeff) {
  df %>%
    select(time, osmotic) %>%
    arrange(time) %>%
    mutate(target = round(time / step) * step) %>%
    group_by(target) %>%
    slice(which.min(abs(time - target))) %>%
    ungroup() %>%
    mutate(
      slope = lag(osmotic) - osmotic,
      decay = -(as.numeric(coeff) * slope),
      decay = replace_na(decay, 0),
      log = cumsum(decay)) %>%
    mutate(
      flat = log[which.min(abs(time - cutoff_time))],
      log = ifelse(time >= cutoff_time, flat, log)) %>%
    select(-flat)
}


# ----- data -----

files <- list(
  as50 = 'as50_ecoli.dat',
  as30 = 'as30_ecoli.dat',
  as70 = 'as70_default.dat',
  as50_10x = 'as5010x_2609.dat',
  as50_0.1x = 'as50_0.1x_2909.dat',
  as40 = 'as40_3010.dat',
  as80 = 'as80_2910.dat',
  # as30_sep = 'as30_sep.dat',
  # as50_sep = 'as50_200326.dat',
  # as50_default = 'as50_default.dat',
  pbs40 = 'pbs40_0303.dat',
  rh2 = 'rhdiff_3_new.dat',
  as60 = 'as60_1212.dat',
  pbs50 = 'pbs50_1212.dat',
  rh3 = 'rhslow_0702.dat'
  # nacl_sat = 'nacl_sat.dat'
)

combined <- map(files, process_data)

# runs the lineplot list for all of the data
lineplots_list <- imap(combined, ~{
  if (.y %in% c("as30", "as50")) {
    lineplots(.x, TRUE, 125)
  } else {
    lineplots(.x)}})

# calculates the derivatives 
derivatives <- map(lineplots_list, ~ dpidt(.x, time_step))

# ----- plots -----

# Supplemental figure for the osmotic pressure change
p_osm <- lineplots_list$as50_10x %>%
  ggplot(aes(time, osmotic)) +
  
  geom_ribbon(aes(ymin = osmotic - osmotic_sd,
                  ymax = osmotic + osmotic_sd),
              fill = "#0072B2", alpha = 0.2) +
  
  geom_line(color = "#0072B2", linewidth = 1.1) +
  
  labs(x = "time (min)", y = bquote(Pi~(MPa))) +
  xlim(10, 30) +
  theme_clean

deriv_plot <- bind_rows(derivatives, .id = "dataset") %>%
  ggplot(aes(time, slope, color = dataset)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "time (min)",
    y = bquote(frac(Delta*Pi, Delta*t))
  ) +
  theme_clean

print(p_osm)
print(deriv_plot)

# ----- decay + fits -----

# read in the decay data for e coli
decay <- read_excel(here('osmoticdata.xlsx'), sheet = 'ecoli') %>%
  filter(!is.na(Media)) %>%
  mutate(osm = op / 2.438e3)

# read in the decay data for s epidermidis
decay_sep <- read_excel(here('osmoticdata.xlsx'), sheet = 'sep') %>%
  filter(!is.na(Media)) %>%
  mutate(osm = op / 2.438e3)

# extract the coefficients for the models
coeff_main <- fit_model(decay)
coeff_sep  <- fit_model(decay_sep)

# generates Figure 2A, 2B
fitplot     <- plot_decay(decay, "darkred", TRUE)
fitplot_sep <- plot_decay(decay_sep, "#0045A0", FALSE)

combinedfits <- fitplot | fitplot_sep

print(fitplot)
print(fitplot_sep)
print(combinedfits)

# ----- predictions -----

# predict the concentration (note as50 is actually artificial saliva at 60%)
bacteria_as50 <- bacteria_predict(lineplots_list$as50, time_step, coeff_main)
bacteria_pbs40 <- bacteria_predict(lineplots_list$pbs40, time_step, coeff_sep)

# generate the plots for figure 3
p_osm_all <- bind_rows(lineplots_list, .id = "dataset") %>%
  ggplot(aes(time, osmotic, color = dataset, fill = dataset)) +
  
  geom_ribbon(aes(ymin = osmotic - osmotic_sd,
                  ymax = osmotic + osmotic_sd),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "time (min)",
    y = bquote(Pi~(MPa)),
    color = "Dataset",
    fill = "Dataset"
  ) +
  
  scale_x_continuous(expand = c(0, 0)) +
  xlim(0, 120) +
  theme_clean

p_osm_all
