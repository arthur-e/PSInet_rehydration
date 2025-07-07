library(dplyr, warn.conflicts = F)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(lubridate, warn.conflicts = F)
library(bbmle, warn.conflicts = F)

JESSICA.CSV <- '~/Downloads/SRER_LATR_pdd_Apr_June.csv'

facet.theme <- theme_linedraw() +
  theme(strip.background = element_blank(),
    strip.text = element_text(color = 'black', face = 'bold'))


# Loading data #################################################################

df <- read.csv(JESSICA.CSV) %>%
  mutate(datetime = ymd_hms(dt)) %>%
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

# Sun sets at 19h00, rises after 05h00 local time
df.out <- df %>%
  select(datetime, WP_m, DOY, hour) %>%
  arrange(datetime, hour) %>%
  filter(hour >= 19 | hour <= 6) %>%
  group_by(DOY) %>%
  # Transform hours into a number of hours past 19h00
  mutate(t = if_else(hour < 19, hour + 24, hour) - 19) %>%
  # Create an identifier for each unique evening period
  # This is faster and can be done in a loop, but should be checked for
  #   generality
  mutate(datetime.rel = datetime - dhours(19),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60)

# No longer necessary because of datetime.rel
# Create an identifier for each unique evening period
# df.out <- df.out %>%
#   mutate(t.int = t - lag(t, 1),
#     t.int = if_else(is.na(t.int), 0, t.int))
# df.out$period <- 1
# j <- 1
# for (i in 2:nrow(df.out)) {
#   if (df.out[i,]$t.int < 0) {
#     j = j + 1
#   }
#   df.out[i,'period'] <- j
# }

df.out %>%
  filter(DOY.rel %% 10 == 0) %>%
ggplot(mapping = aes(x = datetime, y = WP_m)) +
  geom_line() +
  facet_wrap(~ DOY.rel, scales = 'free', ncol = 2)

df.out %>%
  filter(DOY.rel %% 10 == 0) %>%
ggplot(mapping = aes(x = hour.rel, y = WP_m)) +
  geom_line(aes(color = DOY.rel, group = DOY.rel))


# Empirical determination of periods, extremes #################################

# I want each day to start and end at 10h00 local, as this is (likely) before
#   the minimum daily water potential
df.emp <- df %>%
  select(datetime, DOY, hour, WP_m, WP_sd) %>%
  # Roll back "midnight" to 10h00 local
  mutate(datetime.rel = datetime - dhours(10),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
  filter(DOY.rel > 90) %>%
  arrange(datetime)

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
  # Compare to 95% of previous day's maximum, as a kind of smoothing; note
  #   that numerator and denom. are flipped because of negative values
  mutate(rehydrated = (lag(max.psi, 1) / max.psi) >= 0.95) %>%
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

# Diagnostic
df.emp %>%
  left_join(df.emp.agg, by = 'DOY.rel') %>%
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
df.emp.agg %>%
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

df.emp.agg %>%
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


# Empirical determination, with smoothing ######################################

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
  # filter(DOY == 91 | DOY %% 10 == 0) %>%
  filter(DOY < 100) %>%
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
  summarize(diff = median(diff, na.rm = T),
    max.psi = max(smoothed),
    min.psi = min(smoothed)) %>%
  # Compare to 95% of previous day's maximum, as a kind of smoothing; note
  #   that numerator and denom. are flipped because of negative values
  mutate(rehydrated = (lag(max.psi, 1) / max.psi) >= 0.95)

df.smoothing %>%
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



# Identifying rehydration by discrete derivative ###############################

df.emp %>%
  group_by(DOY.rel) %>%
  mutate(WP_fd = WP_m - lag(WP_m, 1),
    WP_fd_perc = -100 * (WP_fd / WP_m)) %>%
  mutate(eq = if_else(WP_fd_perc > 0,
    if_else(WP_fd_perc <= 1, TRUE, FALSE), FALSE)) %>%
  View


# Model fitting framework ######################################################

rmse <- function (observed, predicted) {
  n <- length(predicted)
  # Make sure predicted and observed are the same length!
  stopifnot(n == length(observed))
  return((1/n) * sum(sqrt((predicted - observed)^2), na.rm = T))
}

target.func <- function (t, pd = -2, tau = 4) {
  return(pd * (1 + exp(-t / tau)))
}

wrap.rmse <- function (pd, tau, df = df) {
  t <- df$t
  psi <- target.func(t, pd, tau)
  return(rmse(df$WP_m, psi))
}

plot.function(target.func, from = 0, to = 24, bty = 'n', xlab = 'Hour (Local Time)',
  ylab = 'Water Potential (MPa)', lwd = 2)


# Fitting model parameters with optim ##########################################

# df.out <-
# result <- optim(c(-2, 2), wrap.rmse,
#   method = 'L-BFGS-B', lower = c(-9, 1), upper = c(0, 24))

# require(broom)
# require(purrr)
# df %>%
#   nest(data = -DOY) %>%
#   mutate(fit = map(data, ~ lm(WP_m ~
