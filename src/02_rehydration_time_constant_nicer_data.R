library(dplyr, warn.conflicts = F)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(lubridate, warn.conflicts = F)
library(bbmle, warn.conflicts = F)
library(patchwork)
library(reticulate)

JESSICA.CSV <- '~/Downloads/PSInet/SRER_LATR_pdd_Apr_June.csv'
MET.CSV <- '~/Downloads/PSInet/neon_atmdaily.csv'
REHYDRATION.CUTOFF <- 0.95

facet.theme <- theme_linedraw() +
  theme(strip.background = element_blank(),
    strip.text = element_text(color = 'black', face = 'bold'))

require(reticulate)
if (Sys.info()['nodename'] == 'Gullveig') {
  use_python('/usr/local/python-env/suntransit/bin/python')
} else {
}
suntransit <- import('suntransit')

# TODO Examine how the timing of min, max psi changes over the season


# Loading data #################################################################

df.met <- read.csv(MET.CSV) %>%
  # Roll back "midnight" to 10h00 local
  mutate(date.rel = as.Date(date) - dhours(10),
    DOY.rel = as.integer(format(date.rel, '%j')))

df <- read.csv(JESSICA.CSV) %>%
  mutate(date = as.Date(date),
    datetime = ymd_hms(dt)) %>%
  rowwise() %>%
  mutate(datetime = with_tz(datetime, tzone = 'America/Phoenix')) %>%
  mutate(hour = hour(datetime) + minute(datetime) / 60,
    DOY = as.integer(format(datetime, '%j')))

df %>%
  filter(DOY < 100) %>%
  # filter(DOY %% 10 == 0) %>%
ggplot(mapping = aes(x = hour, y = WP_m)) +
  geom_ribbon(aes(ymin = WP_m - 1 * WP_sd, ymax = WP_m + 1 * WP_sd),
    fill = 'lightblue', alpha = 0.5) +
  geom_line() +
  facet_wrap(~ DOY, scales = 'free_x', nrow = 2)


# Rehydration period: Direct calculation #######################################

# I want each day to start and end at 10h00 local, as this is (likely) before
#   the minimum daily water potential
df.emp <- df %>%
  select(datetime, date, DOY, hour, WP_m, WP_sd) %>%
  # Roll back "midnight" to 10h00 local
  mutate(datetime.rel = datetime - dhours(10),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
  filter(DOY.rel > 90) %>%
  arrange(datetime) %>%
  left_join(df.met, by = 'date')

# Diagnostics
# NOTE: It's apparent that rehydration must occur in less than 24 hours
#   every day, which makes sense because daytime transpiration will start again
#   within that window
df.emp %>%
  filter(DOY.rel == 91 | DOY.rel %% 10 == 0) %>%
  # filter(DOY.rel < 100) %>%
  arrange(datetime) %>%
  mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
ggplot(mapping = aes(x = hour.rel, y = WP_m)) +
  geom_ribbon(aes(ymin = WP_m - 1 * WP_sd, ymax = WP_m + 1 * WP_sd),
    fill = 'lightblue', alpha = 0.5) +
  geom_line() +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
    labels = function (x) { (x + 10) %% 24 }) +
  facet_wrap(~ DOY.rel, scales = 'free_x', nrow = 2) +
  labs(x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)') +
  facet.theme
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_example_time_series_relative_to_10h00_local.png')

# Compute daily min, max water potential and when that is achieved
df.emp.agg <- df.emp %>%
  group_by(DOY.rel) %>%
  summarize(min.psi = min(WP_m),
    max.psi = max(WP_m),
    when.min = (hour.rel[which.min(WP_m)] + 10) %% 24,
    when.max = (hour.rel[which.max(WP_m)] + 10) %% 24) %>%
  # Then, compute how long it took to recharge
  mutate(diff = (24 - when.min) + when.max)

df.emp.agg %>%
  select(when.min, when.max) %>%
  gather(key = group, value = hour) %>%
  mutate(group = if_else(str_detect(group, 'min'), 'Minimum', 'Maximum')) %>%
ggplot(mapping = aes(x = hour)) +
  geom_histogram(aes(fill = group), bins = 48) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = brewer.pal(5, 'RdBu')[c(5,1)]) +
  labs(fill = NULL, x = 'Hour of Day', y = 'Count of Days') +
  theme_minimal() +
  theme(legend.position = 'top',
    legend.margin = margin(0, 0, -0.1, 0, 'cm'))
ggsave(width = 5, height = 4, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_hour_of_day_histograms_min-max_psi.png')

df.emp.agg %>%
ggplot(mapping = aes(x = DOY.rel, y = diff)) +
  geom_line() +
  geom_smooth(fill = 'lightblue') +
  labs(x = 'Day of Year', y = 'Rehydration Time (hours)') +
  theme_minimal()
ggsave(width = 5, height = 4, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_series_empirical.png')


# Plotting rehydration period ##################################################

df.cutoff <- data.frame(cutoff = seq(0.9, 1.1, 0.01), not.rehydrated.count = NA)
for (cutoff in df.cutoff$cutoff) {
  df.cutoff[df.cutoff$cutoff == cutoff,2] <- with(df.emp.agg %>%
    mutate(rehydrated = (lag(max.psi, 1) / max.psi) >= cutoff) %>%
    group_by(rehydrated) %>%
    tally() %>%
    filter(!rehydrated), n)
}

ggplot(df.cutoff, aes(x = cutoff, y = not.rehydrated.count / 90)) +
  geom_bar(stat = 'identity') +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = latex2exp::TeX('Pct. of Pre-Dawn Max $\\psi$ Considered "Full" Rehydration       '),
    y = 'Days without "Full" Rehydration (%)') +
  theme_minimal()
ggsave(width = 4, height = 3, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_cutoff_barplot.png')

df.emp %>%
  left_join(
    # Compare to, e.g., 95% of previous day's maximum, as a kind of smoothing; note
    #   that numerator and denom. are flipped because of negative values
    mutate(df.emp.agg, rehydrated = (lag(max.psi, 1) / max.psi) >= REHYDRATION.CUTOFF),
    by = 'DOY.rel') %>%
  filter(DOY.rel > 91, DOY.rel < 102) %>%
  mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
  rename(`Rehydrated?` = rehydrated) %>%
ggplot(mapping = aes(x = hour.rel, y = WP_m)) +
  geom_ribbon(aes(ymin = WP_m - 1 * WP_sd, ymax = WP_m + 1 * WP_sd),
    fill = 'lightgray', alpha = 0.5) +
  geom_line(aes(color = `Rehydrated?`), linewidth = 0.6) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
    labels = function (x) { (x + 10) %% 24 }) +
  scale_color_manual(values = brewer.pal(5, 'RdBu')[c(1,5)]) +
  facet_wrap(~ DOY.rel, scales = 'free_x', nrow = 2) +
  labs(x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)') +
  facet.theme +
  theme(legend.position = 'top',
    legend.margin = margin(0, 0, -0.2, 0, 'cm'))
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_series_relative_to_10h00_local.png')

require(viridis)
# Compare to, e.g., 95% of previous day's maximum, as a kind of smoothing; note
#   that numerator and denom. are flipped because of negative values
mutate(df.emp.agg, rehydrated = (lag(max.psi, 1) / max.psi) >= REHYDRATION.CUTOFF) %>%
  filter(!is.na(rehydrated)) %>%
  rename(`Rehydrated?` = rehydrated) %>%
  mutate(diff.psi = min.psi - max.psi) %>%
ggplot(mapping = aes(x = DOY.rel, y = `Rehydrated?`)) +
  geom_point(aes(color = min.psi), shape = 1) +
  scale_color_gradientn(colors = magma(10, direction = -1)) +
  labs(x = 'Day of Year', color = latex2exp::TeX('min($\\psi$)')) +
  theme_dark()
ggsave(width = 5, height = 2, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_series_dot_plot.png')

# Compare to, e.g., 95% of previous day's maximum, as a kind of smoothing; note
#   that numerator and denom. are flipped because of negative values
mutate(df.emp.agg, rehydrated = (lag(max.psi, 1) / max.psi) >= REHYDRATION.CUTOFF) %>%
  filter(!is.na(rehydrated)) %>%
  rename(`Rehydrated?` = rehydrated) %>%
  mutate(diff.psi = min.psi - max.psi) %>%
  select(`Rehydrated?`, `Mid-day WP` = min.psi, `Pre-dawn WP` = max.psi,
    `Rehydration Time` = diff, `WP Gradient` = diff.psi) %>%
  gather(key = group, value = value, -`Rehydrated?`) %>%
ggplot(mapping = aes(x = `Rehydrated?`, y = value)) +
  geom_boxplot(fill = 'lightgray') +
  geom_jitter(width = 0.1, size = 0.8) +
  labs(y = NULL) +
  facet_wrap(~ group, scales = 'free_y', nrow = 1) +
  facet.theme +
  theme(panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())
ggsave(width = 6.5, height = 2.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_status_and_covariates.png')

# How far off are we from equilibrium each day?
df.emp.agg %>%
  filter(!is.na(max.psi)) %>%
  mutate(perc.diff = -100 * ((max.psi - lag(max.psi, 1)) / lag(max.psi, 1))) %>%
  left_join(df.met, by = 'DOY.rel') %>%
  mutate(rain = if_else(ppt_mm > 1, DOY.rel, NA)) %>%
ggplot(mapping = aes(x = DOY.rel, y = perc.diff)) +
  geom_vline(aes(xintercept = rain), color = 'darkblue', linetype = 'dotted') +
  geom_bar(aes(fill = as.character(sign(perc.diff))), stat = 'identity') +
  geom_hline(aes(yintercept = yint), color = '#333333', linetype = 'dashed',
    data = data.frame(yint = c(-5, 5))) +
  scale_x_continuous(expand = c(0, 0.5)) +
  scale_fill_manual(values = brewer.pal(5, 'RdBu')[c(1,5)], guide = 'none') +
  labs(y = latex2exp::TeX('Difference from Previous Pre-Dawn Maximum $\\psi$'),
    x = 'Day of Year') +
  theme_minimal()
ggsave(width = 5, height = 4, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_change_in_daily_max_WP_alt.png')


# Change in rehydration time? ##################################################

df.emp.agg %>%
  select(DOY.rel, when.min, when.max) %>%
  gather(key = group, value = hour, -DOY.rel) %>%
  mutate(group = if_else(str_detect(group, 'min'), 'Minimum', 'Maximum')) %>%
  mutate(date = as.Date(sprintf('2023-%03d', DOY.rel), '%Y-%j')) %>%
  group_by(DOY.rel) %>%
  # America/Phoenix is GMT-07
  mutate(sunrise = (unlist(suntransit$sunrise_sunset(
    c(33, -112), dt = date[1]))[1] - 7) %% 24) %>%
  mutate(sunset = (unlist(suntransit$sunrise_sunset(
    c(33, -112), dt = date[1]))[2] - 7) %% 24) %>%
ggplot(mapping = aes(x = DOY.rel, y = hour)) +
  geom_ribbon(aes(ymin = sunrise, ymax = sunset), fill = 'gold', alpha = 0.5) +
  geom_line(aes(color = group)) +
  scale_color_manual(values = brewer.pal(5, 'RdBu')[c(5,1)]) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(color = NULL, x = 'Day of Year', y = 'Hour (Local Time)') +
  theme_minimal() +
  theme(legend.position = 'top', legend.margin = margin(0, 0, -0.3, 0, 'cm'))
ggsave(width = 6, height = 4, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_hour_of_day_min_max_with_daylight_hours.png')

# Example
# (unlist(suntransit$sunrise_sunset(c(33, -112), dt = as.Date('2015-03-31'))) - 7) %% 24



# Rehydration period: Smoothing ################################################

require(broom)
require(purrr)
smoother <- function (df0, window.size) {
  # Expects the columns: DOY, hour, value
  ff.b <- rep(1 / window.size, window.size)

  df0 %>%
    arrange(DOY, hour) %>%
    tidyr::nest(data = -c(DOY)) %>%
    mutate(smoothed = purrr::map(data,
      ~ signal::filtfilt(ff.b, c(1),
        c(rep(.x$value[1], window.size), .x$value,
          rep(.x$value[n()], window.size))))) %>%
    tidyr::unnest(smoothed) %>%
    group_by(DOY) %>%
    mutate(hour = seq(-window.size * 0.5, 23.5 + (window.size * 0.5), 0.5)) %>%
    filter(hour >= 0, hour <= 23.5) %>%
    left_join(df0, by = c('DOY', 'hour')) %>%
    # NOTE: The actual day starts at 10h00 local
    mutate(hour.rel = (hour + 10) %% 24) %>%
    select(-data)
}

df.tmp.6 <- df.emp %>%
  select(datetime.rel, DOY.rel, hour.rel, value = WP_m) %>%
  rename(DOY = DOY.rel, hour = hour.rel) %>%
  group_by(DOY) %>%
  filter(n() == 48) %>%
  smoother(window.size = 6) %>%
  group_by(DOY) %>%
  mutate(fd = smoothed - lag(smoothed, 1),
    fd.perc = -100 * (fd / smoothed),
    smoothed.rel = abs(100 * ((max(smoothed) - smoothed) / max(smoothed))), # Relative to max value
    # eq.max = smoothed == max(smoothed)) %>%
    eq.max = (hour > 8 & fd.perc > 0 & fd.perc < 1 & smoothed.rel < 3)) %>%
  arrange(DOY, hour)

df.tmp.12 <- df.emp %>%
  select(datetime.rel, DOY.rel, hour.rel, value = WP_m) %>%
  rename(DOY = DOY.rel, hour = hour.rel) %>%
  group_by(DOY) %>%
  filter(n() == 48) %>%
  smoother(window.size = 12) %>%
  group_by(DOY) %>%
  mutate(fd = smoothed - lag(smoothed, 1),
    fd.perc = -100 * (fd / smoothed),
    smoothed.rel = abs(100 * ((max(smoothed) - smoothed) / max(smoothed))), # Relative to max value
    eq.min = smoothed == min(smoothed),
    eq.max = smoothed == max(smoothed)) %>%
  arrange(DOY, hour)

df.tmp.6 %>%
  filter(DOY == 91 | DOY %% 10 == 0) %>%
  mutate(DOY = sprintf('DOY=%03d', DOY)) %>%
ggplot(mapping = aes(x = hour)) +
  geom_line(aes(y = value), color = 'black', linewidth = 0.6) +
  geom_line(aes(y = smoothed), color = 'red', linewidth = 0.6) +
  geom_vline(aes(xintercept = the.hour), color = 'darkblue', linetype = 'dashed',
    data = filter(df.tmp.6, eq.max, DOY == 91 | DOY %% 10 == 0) %>%
      summarize(the.hour = first(hour)) %>%
      mutate(DOY = sprintf('DOY=%03d', DOY))) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
    labels = function (x) { (x + 10) %% 24 }) +
  scale_color_manual(values = brewer.pal(5, 'RdBu')[c(1,5)]) +
  facet_wrap(~ DOY, scales = 'free_x', nrow = 2) +
  labs(x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)',
    subtitle = 'With a 6-hour low-pass moving window') +
  facet.theme +
  theme(legend.position = 'top',
    legend.margin = margin(0, 0, -0.2, 0, 'cm'))
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_example_time_series_smoothing_peak-ident_06hr.png')

df.tmp.12 %>%
  filter(DOY == 91 | DOY %% 10 == 0) %>%
  mutate(DOY = sprintf('DOY=%03d', DOY)) %>%
ggplot(mapping = aes(x = hour)) +
  geom_line(aes(y = value), color = 'black', linewidth = 0.6) +
  geom_line(aes(y = smoothed), color = 'red', linewidth = 0.6) +
  geom_vline(aes(xintercept = the.hour), color = 'darkblue', linetype = 'dashed',
    data = filter(df.tmp.12, eq.max, DOY == 91 | DOY %% 10 == 0) %>%
      summarize(the.hour = first(hour)) %>%
      mutate(DOY = sprintf('DOY=%03d', DOY))) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
    labels = function (x) { (x + 10) %% 24 }) +
  scale_color_manual(values = brewer.pal(5, 'RdBu')[c(1,5)]) +
  facet_wrap(~ DOY, scales = 'free_x', nrow = 2) +
  labs(x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)',
    subtitle = 'With a 12-hour low-pass moving window') +
  facet.theme +
  theme(legend.position = 'top',
    legend.margin = margin(0, 0, -0.2, 0, 'cm'))
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_example_time_series_smoothing_peak-ident_06hr.png')

df.smoothing <- df.tmp.12 %>%
  filter(eq.min | eq.max) %>%
  arrange(DOY, hour) %>%
  group_by(DOY) %>%
  filter(n() == 2) %>%
  mutate(diff = (hour.rel + 24) - lag(hour.rel, 1)) %>%
  select(datetime.rel, DOY.rel = DOY, hour.rel, smoothed, eq.min, eq.max, diff) %>%
  group_by(DOY.rel) %>%
  # NOTE: Arranging by eq.min necessary to get when.min, when.max
  arrange(DOY.rel, eq.min) %>%
  summarize(diff = median(diff, na.rm = T),
    when.max = hour.rel[1],
    when.min = hour.rel[2],
    max.psi = max(smoothed),
    min.psi = min(smoothed)) %>%
  # Compare to 95% of previous day's maximum, as a kind of smoothing; note
  #   that numerator and denom. are flipped because of negative values
  mutate(rehydrated = (lag(max.psi, 1) / max.psi) >= REHYDRATION.CUTOFF)

require(viridis)
df.smoothing %>%
  filter(!is.na(rehydrated)) %>%
  rename(`Rehydrated?` = rehydrated) %>%
  mutate(diff.psi = min.psi - max.psi) %>%
ggplot(mapping = aes(x = DOY.rel, y = `Rehydrated?`)) +
  geom_point(aes(color = min.psi), shape = 1) +
  scale_color_gradientn(colors = magma(10, direction = -1)) +
  labs(x = 'Day of Year', color = latex2exp::TeX('min($\\psi$)')) +
  theme_dark()

# How far off are we from equilibrium each day?
df.smoothing %>%
  filter(!is.na(max.psi)) %>%
  mutate(perc.diff = -100 * ((max.psi - lag(max.psi, 1)) / lag(max.psi, 1))) %>%
  left_join(df.met, by = 'DOY.rel') %>%
  mutate(rain = if_else(ppt_mm > 1, DOY.rel, NA)) %>%
ggplot(mapping = aes(x = DOY.rel, y = perc.diff)) +
  geom_vline(aes(xintercept = rain), color = 'darkblue', linetype = 'dotted') +
  geom_bar(aes(fill = as.character(sign(perc.diff))), stat = 'identity') +
  geom_hline(aes(yintercept = yint), color = '#333333', linetype = 'dashed',
    data = data.frame(yint = c(-5, 5))) +
  scale_x_continuous(expand = c(0, 0.5)) +
  scale_fill_manual(values = brewer.pal(5, 'RdBu')[c(1,5)], guide = 'none') +
  labs(y = latex2exp::TeX('Difference from Previous Pre-Dawn Maximum $\\psi$'),
    x = 'Day of Year') +
  theme_minimal()
ggsave(width = 5, height = 4, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_change_in_daily_max_WP_from_smoothing_alt.png')


# Model fitting framework ######################################################

rmse <- function (observed, predicted) {
  n <- length(predicted)
  # Make sure predicted and observed are the same length!
  stopifnot(n == length(observed))
  return((1/n) * sum(sqrt((predicted - observed)^2), na.rm = T))
}

model.exp <- function (t, pd = -2, tau = 4) {
  # Note this is 1 + ... instead of 1 - ... because WP is negative
  return(pd * (1 + exp(-t / tau)))
}

# Our objective function, calculating RMSE based on model predictions
rmse.model.exp <- function (params, df = df.out) {
  t <- df$t
  pd <- params[1]
  tau <- params[2]
  psi <- model.exp(t, pd, tau)
  return(rmse(df$WP_m, psi))
}

# The idealized function we're trying to fit
plot.function(model.exp, from = 0, to = 24, bty = 'n', xlab = 'Hour (Local Time)',
  ylab = 'Water Potential (MPa)', lwd = 2, xaxt = 'n')
axis(1, at = seq(0, 20, 5), labels = c(16, 20, 24, 4, 8))


# Fitting model parameters with optim ##########################################

# Sun sets at 19h00, rises after 05h00 local time, but stomata tend to close
#   early due to stress, so we'll start the clock at ??h00
HOUR.START <- 14
HOUR.END <- 7
df.out <- df %>%
  select(datetime, WP_m, DOY, hour) %>%
  arrange(datetime, hour) %>%
  filter(hour >= HOUR.START | hour <= HOUR.END) %>%
  group_by(DOY) %>%
  # Transform hours into a number of hours past HOUR.START
  mutate(t = if_else(hour < HOUR.START, hour + 24, hour) - HOUR.START) %>%
  # Create an identifier for each unique evening period
  # This is faster and can be done in a loop, but should be checked for
  #   generality
  mutate(datetime.rel = datetime - dhours(HOUR.START),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60)

# Diagnostics
df.out %>%
  filter(DOY.rel %% 10 == 0) %>%
ggplot(mapping = aes(x = datetime, y = WP_m)) +
  geom_line() +
  facet_wrap(~ DOY.rel, scales = 'free', ncol = 2)

df.out %>%
  filter(DOY.rel %% 10 == 0) %>%
ggplot(mapping = aes(x = hour.rel, y = WP_m)) +
  geom_line(aes(color = DOY.rel, group = DOY.rel)) +
  # Change X-axis so that the "day" starts at 10h00 local time
  scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
    labels = function (x) { (x + HOUR.START) %% 24 }) +
  labs(x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)',
    color = 'DOY')
ggsave(width = 5.5, height = 5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_curves_by_day_over_season.png')

# A data frame with the unique DOYs
df.doy <- df.out %>%
  group_by(DOY.rel) %>%
  summarize()
fit.params <- matrix(nrow = length(df.doy$DOY.rel), ncol = 3)
fit.scores <- matrix(nrow = nrow(fit.params), ncol = 1)
for (i in 1:nrow(fit.params)) {
  doy <- df.doy$DOY.rel[i]
  df.sub <- filter(df.out, DOY.rel == doy)
  if (nrow(df.sub) == 0) next() # Skip empty subsets
  wrap.rmse <- function (params) {
    rmse.model.exp(params, df = df.sub)
  }
  result <- optim(c(-2, 2), wrap.rmse,
    method = 'L-BFGS-B', lower = c(-9, 1), upper = c(0, 24))
    # control = list(pgtol = 0.01))
  fit.params[i,1] <- doy
  fit.params[i,2:3] <- result$par
  fit.scores[i,] <- wrap.rmse(result$par)
}

df.params <- bind_cols(
  as.data.frame(fit.params) %>%
    rename(DOY.rel = V1, psi0 = V2, tau = V3),
  as.data.frame(fit.scores) %>%
    rename(RMSE = V1))

df.params %>%
  left_join(df.met, by = 'DOY.rel') %>%
  mutate(rain = if_else(ppt_mm > 1, DOY.rel, NA)) %>%
ggplot(mapping = aes(x = DOY.rel, y = RMSE)) +
  geom_vline(aes(xintercept = rain), color = 'darkblue', linetype = 'dotted',
    linewidth = 1) +
  geom_bar(stat = 'identity') +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = 'Day of Year', y = 'Model RMSE (MPa)') +
  theme_minimal()
ggsave(width = 6.5, height = 3, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_RMSE_barplot_alt.png')

i <- 1
g1 <- NULL
for (doy in df.params$DOY.rel) {
  if (doy == 90 | doy %% 10 != 0) next()
  g0 <- df.out %>%
    filter(DOY.rel == doy) %>%
    mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
  ggplot(mapping = aes(x = hour.rel)) +
    geom_line(aes(y = WP_m), color = 'black', linewidth = 0.6) +
    # Change X-axis so that the "day" starts at 10h00 local time
    scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
      labels = function (x) { (x + 14) %% 24 }) +
    # scale_y_continuous(limits = c(-10, -2)) +
    labs(subtitle = sprintf('DOY=%03d', doy),
      x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)')

  g0 <- g0 + with(filter(df.params, DOY.rel == doy),
      geom_function(fun = model.exp, n = 1000, xlim = c(0, 16),
        args = list(pd = psi0, tau = tau), linetype = 'dashed', color = 'red'))

  if (i == 1) {
    g1 <- g0
  } else {
    g1 <- g1 + g0
  }
  i <- 1 + 1
}

g1 + plot_layout(nrow = 3, axis_titles = 'collect')
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_fits.png')

df.emp.agg %>%
  head

# 5-tau analysis
df.params %>%
  left_join(df.emp.agg, by = 'DOY.rel') %>%
  mutate(five.tau = tau * 5) %>%
  mutate(re.period = (when.max + 24) - when.min) %>%
  # One definition: Rehydration occurs if the 5*tau is less than the observed
  #   rehydration period
  mutate(rehydrated = five.tau < re.period) %>%
  filter(!is.na(rehydrated)) %>%
ggplot(mapping = aes(x = DOY.rel, y = rehydrated)) +
  geom_point(aes(color = re.period), shape = 1, size = 2) +
  scale_color_gradientn(colors = magma(10, direction = -1)) +
  labs(x = 'Day of Year', y = 'Rehydrated?',
    subtitle = latex2exp::TeX('Based on $5\\tau$ less than argmax($\\Psi$) - argmin($\\Psi$)'),
    color = 'Rehydration\nTime (hrs)') +
  theme_dark() +
  theme(plot.background = element_blank(),
    legend.background = element_blank())
ggsave(width = 5, height = 2.4, dpi = 172, bg = 'transparent',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_series_5tau.png')

df.params %>%
  left_join(df.smoothing, by = 'DOY.rel') %>%
  mutate(five.tau = tau * 5) %>%
  mutate(re.period = (when.max + 24) - when.min) %>%
  # One definition: Rehydration occurs if the 5*tau is less than the observed
  #   rehydration period
  mutate(rehydrated = five.tau < re.period) %>%
  filter(!is.na(rehydrated)) %>%
ggplot(mapping = aes(x = DOY.rel, y = rehydrated)) +
  geom_point(aes(color = re.period), shape = 1, size = 2) +
  scale_color_gradientn(colors = magma(10, direction = -1)) +
  labs(x = 'Day of Year', y = 'Rehydrated?',
    subtitle = latex2exp::TeX('Based on $5\\tau$ less than argmax($\\Psi$) - argmin($\\Psi$) (Smoothing)'),
    color = 'Rehydration\nTime (hrs)') +
  theme_dark() +
  theme(plot.background = element_blank(),
    legend.background = element_blank())
ggsave(width = 5, height = 2.4, dpi = 172, bg = 'transparent',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_series_5tau_smoothed.png')

df.params %>%
  left_join(df.emp.agg, by = 'DOY.rel') %>%
  mutate(five.tau = tau * 5) %>%
  mutate(date = as.Date(sprintf('2023-%03d', DOY.rel), '%Y-%j')) %>%
  group_by(DOY.rel) %>%
  # America/Phoenix is GMT-07
  mutate(sunrise = (unlist(suntransit$sunrise_sunset(
    c(33, -112), dt = date[1]))[1] - 7) %% 24) %>%
  mutate(sunset = (unlist(suntransit$sunrise_sunset(
    c(33, -112), dt = date[1]))[2] - 7) %% 24) %>%
  # One definition: Rehydration occurs if the 5*tau is less than the observed
  #   rehydration period
  mutate(re.time = (when.min + five.tau) %% 24,
    rehydrated = re.time < sunrise) %>%
  filter(!is.na(rehydrated)) %>%
ggplot(mapping = aes(x = DOY.rel, y = rehydrated)) +
  geom_point(aes(color = re.time), shape = 1, size = 2) +
  scale_color_gradientn(colors = magma(10, direction = -1)) +
  labs(x = 'Day of Year', y = 'Rehydrated?',
    subtitle = latex2exp::TeX('Based on (argmin($\\Psi$) + $5\\tau$) earlier than sunrise'),
    color = 'Rehydration\nTime (hrs)') +
  theme_dark() +
  theme(plot.background = element_blank(),
    legend.background = element_blank())
ggsave(width = 5, height = 2.4, dpi = 172, bg = 'transparent',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_series_5tau_less_than_sunrise.png')
