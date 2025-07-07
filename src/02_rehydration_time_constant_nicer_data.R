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
  filter(DOY %in% 110:119) %>%
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


# Empirical determination, more sophisticated ##################################

# In this section, we don't assume that rehydration can occur in less than 24
#   hours; i.e., we look forward 30 hours from 10h00 local

require(zoo)
df.emp$WP_m_min <- with( # NOTE: There is no rollmin() function
  mutate(df.emp, neg_WP_m = -WP_m), -(rollmax(neg_WP_m, k = 30, align = 'left', fill = NA)))
df.emp$WP_m_max <- with(
  mutate(df.emp, neg_WP_m = -WP_m), rollmax(WP_m, k = 30, align = 'left', fill = NA))

df.emp %>%
  # We used 30-hour windows starting at 10h00 local
  filter(hour == 10) %>%
  select(WP_m_max, WP_m_min) %>%
  gather(key = group, value = hour) %>%
  mutate(group = if_else(str_detect(group, 'min'), 'Minimum', 'Maximum')) %>%
ggplot(mapping = aes(x = hour)) +
  geom_density(aes(fill = group), alpha = 0.5, color = 'transparent') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = brewer.pal(5, 'RdBu')[c(5,1)]) +
  labs(fill = NULL, x = 'Water Potential (MPa)', y = 'Count of Days') +
  theme_minimal() +
  theme(legend.position = 'top',
    legend.margin = margin(0, 0, -0.1, 0, 'cm'))
ggsave(width = 5, height = 4, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_rehydration_time_min-max_water_potential_from_30hr_moving_window.png')

df.rehydration <- df.emp %>%
  group_by(DOY.rel) %>%
  filter(!is.na(WP_m_min)) %>%
  summarize(
    when.min = hour[which(WP_m == WP_m_min)],
    when.max = hour[which(WP_m == WP_m_max)],
    # We used 30-hour windows starting at 10h00 local
    min.psi = first(WP_m_min),
    max.psi = first(WP_m_max)) %>%
  # Compute whether recharge happened and how long it took
  mutate(max.psi.last = round(lag(max.psi, 1), 1),
    rehydration = round(max.psi, 1) >= max.psi.last,
    time.hours = (24 - when.min) + when.max) %>%
  filter(!is.na(rehydration))

df.emp %>%
  left_join(df.rehydration, by = 'DOY.rel') %>%
  filter(DOY.rel < 100) %>%
  select(datetime, datetime.rel, DOY.rel, hour, WP_m, min.psi, max.psi, rehydration) %>%
  mutate(hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
ggplot(mapping = aes(x = hour.rel, y = WP_m)) +
  geom_line(aes(color = rehydration, group = NULL)) +
  facet_wrap(~ DOY.rel)


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
