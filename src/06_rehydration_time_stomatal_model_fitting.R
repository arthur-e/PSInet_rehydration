library(dplyr, warn.conflicts = F)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(lubridate, warn.conflicts = F)
library(bbmle, warn.conflicts = F)
library(reticulate)
library(patchwork)
library(latex2exp)

JESSICA.CSV <- '~/Downloads/PSInet/SRER_LATR_pdd_Apr_June.csv'
MET.CSV <- '~/Downloads/PSInet/neon_atmdaily.csv'
SWGDN.CSV <- '/home/arthur.endsley/Workspace/NTSG/projects/Y2026_PSInet/data/202050725_Jessica_Creosote_site_met_data_MERRA_SWGDN.csv'

facet.theme <- theme_linedraw() +
  theme(strip.background = element_blank(),
    strip.text = element_text(color = 'black', face = 'bold'))

require(reticulate)
if (Sys.info()['nodename'] == 'Gullveig') {
  use_python('/usr/local/python-env/suntransit/bin/python')
} else {
}
suntransit <- import('suntransit')


# Loading data #################################################################

df.met <- read.csv(MET.CSV) %>%
  # Roll back "midnight" to 10h00 local
  mutate(datetime.rel = as.Date(date) - dhours(10),
    DOY.rel = as.integer(format(datetime.rel, '%j')))

df.par <- read.csv(SWGDN.CSV) %>%
  mutate(date = as.Date(as.character(date), '%Y%m%d')) %>%
  mutate(datetime = ymd_h(sprintf('%s %02d', format(date, '%Y-%m-%d'), hour)),
    datetime = with_tz(datetime, tzone = 'America/Phoenix')) %>%
  select(datetime, SWGDN = value) %>%
  # NOTE: About 45% of the incoming SWGDN is in the PAR spectrum
  mutate(PAR.Wm2 = 0.45 * SWGDN)

df <- read.csv(JESSICA.CSV) %>%
  mutate(date = as.Date(date),
    datetime = ymd_hms(dt)) %>%
  rowwise() %>%
  mutate(datetime = with_tz(datetime, tzone = 'America/Phoenix')) %>%
  mutate(hour = hour(datetime) + minute(datetime) / 60,
    DOY = as.integer(format(datetime, '%j'))) %>%
  left_join(df.par, by = 'datetime') %>%
  arrange(datetime) %>%
  # We can just copy the hourly PAR data to half-hourly steps because PAR
  #   is in units of power [W m-2]
  fill(PAR.Wm2, .direction = 'down') %>%
  # Now we convert to units of energy
  # NOTE Convert from [W m-2] to [MJ m-2]; 1 W == 1 J s-1,
  #   where the numerator is 1800 secs per half-hour and denom. is 1e6 J
  mutate(PAR.MJm2 = (1800 / 1e6) * PAR.Wm2)

# Sun sets at 19h00, rises after 05h00 local time, but stomata tend to close
#   early due to stress, so we'll start the clock at ??h00
HOUR.START <- 12
HOUR.END <- 8
df.out <- df %>%
  select(datetime, WP_m, WP_sd, PAR.MJm2, DOY, hour) %>%
  arrange(datetime, hour) %>%
  filter(hour >= HOUR.START | hour <= HOUR.END) %>%
  group_by(DOY) %>%
  # Transform hours into a number of hours past HOUR.START
  # mutate(t = if_else(hour < HOUR.START, hour + 24, hour) - HOUR.START) %>%
  # Create an identifier for each unique evening period
  # This is faster and can be done in a loop, but should be checked for
  #   generality
  mutate(datetime.rel = datetime - dhours(HOUR.START),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
  left_join(select(df.met, datetime.rel, Dmean), by = 'datetime.rel') %>%
  ungroup()

# Diagnostics
# df.out %>%
#   filter(DOY.rel %% 10 == 0) %>%
# ggplot(mapping = aes(x = datetime, y = WP_m)) +
#   geom_line() +
#   facet_wrap(~ DOY.rel, scales = 'free', ncol = 2)

# Interpolating met data
df.clean <- df.out %>%
  select(datetime, datetime.rel, DOY.rel, hour.rel, t, PAR.MJm2, Dmean, WP_m) %>%
  group_by(DOY.rel) %>%
  mutate(Dmean = median(Dmean, na.rm = T)) %>%
  fill(PAR.MJm2, .direction = 'downup')

# Filter each time series to *just* those hours *after* the minimum WP is reached
df.clean.min <- df.clean %>%
  group_by(DOY.rel) %>%
  filter(hour.rel >= hour.rel[which.min(WP_m)]) %>%
  # Transform hours into a number of hours past HOUR.START
  mutate(t = hour.rel - min(hour.rel))


# Tom's model (example) ########################################################

# Simulate the differential equation, forward in time
model.stomatal <- function (
    ts, vpd = 1000, irr = 0.5, psi0 = -3, tau = 4, psi.source = -2, k = 1, c1 = 1, c2 = 1) {
  wp <- numeric(length = length(ts)) # Output vector of WP
  # Initialize the WP time series with WP at first step
  wp[1] <- psi0
  # Now, each WP of each successive time step depends on the last
  for (t in 2:length(ts)) {
    dt <- ts[t] - ts[t-1]
    wp[t] <- wp[t-1] + (dt/tau) * (psi.source - wp[t-1] - (c1 * irr[t-1] * vpd[t-1]) / ((irr[t-1] + k) * (1 + c2 * vpd[t-1])))
  }
  return(wp)
}

# Confirmed we get the same result as Tom's script
result <- model.stomatal(1:10,
  vpd = c(30,20,20,10,5,5,5,3,3,3),
  irr = c(1000,900,700,400,50,0,0,0,0,0),
  tau = 5, psi0 = -2, psi.source = 0, c1 = 0.3, c2 = 0.1, k = 500)


# Model fitting framework ######################################################

rmse <- function (observed, predicted) {
  n <- length(predicted)
  # Make sure predicted and observed are the same length!
  stopifnot(n == length(observed))
  return((1/n) * sum(sqrt((predicted - observed)^2), na.rm = T))
}

rmse.model.stomatal <- function (params, df = df.clean.min) {
  ts <- df$t
  vpd <- df$Dmean
  irr <- df$PAR.MJm2
  psi0 <- df$WP_m[1]
  tau <- params[1]
  psi.source <- params[2]
  k <- params[3]
  c1 <- params[4]
  c2 <- params[5]
  pred <- model.stomatal(ts, vpd, irr, psi0, tau, psi.source, k, c1, c2)
  return(with(df, rmse(WP_m, pred)))
}

# Example
# with(filter(df.clean, DOY.rel == 110),
#   plot(model.stomatal(t, Dmean, PAR.MJm2, psi0 = WP_m[1]),
#     type = 'l', bty = 'n', ylab = 'Water Potential (MPa)', lwd = 2, xaxt = 'n'))
# abline(h = -2)
# axis(1, at = seq(0, 42, 6), labels = seq(0, 42, 6) / 2)


# Fitting the model ############################################################

# A data frame with the unique DOYs
df.doy <- df.clean.min %>%
  group_by(DOY.rel) %>%
  summarize()
fit.params <- matrix(nrow = length(df.doy$DOY.rel), ncol = 1 + 5)
fit.scores <- matrix(nrow = nrow(fit.params), ncol = 1)
for (i in 1:nrow(fit.params)) {
  doy <- df.doy$DOY.rel[i]
  df.sub <- filter(df.clean.min, DOY.rel == doy)
  if (nrow(df.sub) == 0) next() # Skip empty subsets
  wrap.rmse <- function (params) {
    rmse.model.stomatal(params, df = df.sub)
  }
  # Parameters are: tau, psi.source, k, c1, c2
  result <- optim(c(6, -1, 100, 0.1, 0.2), wrap.rmse, method = 'L-BFGS-B',
    lower = c(0, -6, 0, 0, 0.1), upper = c(24, 0, 1000, 1, 1))
    # control = list(pgtol = 0.01))
  fit.params[i,1] <- doy
  fit.params[i,2:6] <- result$par
  fit.scores[i,] <- wrap.rmse(result$par)
}

df.params <- bind_cols(
  as.data.frame(fit.params) %>%
    rename(DOY.rel = V1, tau = V2, psi.source = V3, k = V4, c1 = V5, c2 = V6),
  as.data.frame(fit.scores) %>%
    rename(RMSE = V1)) %>%
  # NOTE: We need dynamic inputs for this function
  right_join(df.clean.min, by = 'DOY.rel')


# Plotting results #############################################################

with(df.params, barplot(RMSE, names.arg = DOY.rel, xlab = 'Day of Year', ylab = 'RMSE (MPa)'))

i <- 1
g1 <- NULL
# for (doy in df.params$DOY.rel) {
for (doy in seq(110, 180, 10)) {
  g0 <- df.clean.min %>%
    filter(DOY.rel == doy) %>%
    mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
  ggplot(mapping = aes(x = hour.rel)) +
    geom_line(aes(y = WP_m), color = 'black', linewidth = 0.6) +
    # Change X-axis so that the "day" starts at 10h00 local time
    scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
      labels = function (x) { (x + 14) %% 24 }) +
    scale_y_continuous(limits = c(-9, -2)) +
    labs(subtitle = sprintf('DOY=%03d', doy),
      x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)')

  g0 <- g0 + with(filter(df.params, DOY.rel == doy),
      geom_function(fun = model.stomatal, n = max(t) - min(t),
        xlim = c(min(hour.rel), max(hour.rel)),
        # NOTE: Static inputs are indexed at their first value
        args = list(vpd = Dmean, irr = PAR.MJm2, psi0 = WP_m[1],
          tau = tau[1], psi.source = psi.source[1], k = k[1], c1 = c1[1], c2 = c2[1]),
        linetype = 'dashed', color = 'red'))
  g0 <- g0 + geom_hline(aes(yintercept = psi.source), color = 'blue',
    linetype = 'dashed', data = filter(df.params, DOY.rel == doy))

  if (i == 1) {
    g1 <- g0
  } else {
    g1 <- g1 + g0
  }
  i <- 1 + 1
}

g1 + plot_layout(nrow = 2, axis_titles = 'collect')
# ggsave(width = 6.5, height = 5, dpi = 172,
#   file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/Buckley_model_fits_examples2.png')


# Model inference ##############################################################

g2 <- df.params %>%
  mutate(date = as.Date(sprintf('2023-%03d', DOY.rel), '%Y-%j')) %>%
  group_by(DOY.rel, date) %>%
  filter(n() > 3) %>%
  summarize(tau = first(tau), psi.source = first(psi.source)) %>%
  gather(key = Parameter, value = value, -DOY.rel:-date) %>%
  left_join(select(df.met, DOY.rel, ppt_mm), by = 'DOY.rel') %>%
  mutate(precip = if_else(ppt_mm > 0, date, NA)) %>%
  mutate(Parameter = if_else(str_detect(Parameter, 'tau'), 'tau~paste("(hours)")', 'Psi[0]~paste("(MPa)")')) %>%
  group_by(DOY.rel) %>%
  # America/Phoenix is GMT-07
  mutate(sunrise = (unlist(suntransit$sunrise_sunset(
    c(33, -112), dt = date[1]))[1] - 7) %% 24) %>%
  mutate(sunset = (unlist(suntransit$sunrise_sunset(
    c(33, -112), dt = date[1]))[2] - 7) %% 24) %>%
  mutate(daylight = if_else(str_detect(Parameter, 'tau'), sunset - sunrise, NA),
    nighttime = if_else(str_detect(Parameter, 'tau'), (24 - sunset) + sunrise, NA)) %>%
ggplot(mapping = aes(x = date, y = value)) +
  # geom_area(aes(y = daylight), fill = '#ffcc00', alpha = 0.5) +
  geom_area(aes(y = nighttime), fill = 'azure3', alpha = 0.6) +
  geom_vline(aes(xintercept = precip), color = 'darkblue', linetype = 'dashed') +
  geom_line(aes(color = Parameter), linewidth = 0.6) +
  facet_wrap(~ Parameter, ncol = 1, strip.position = 'right',
    labeller = label_parsed, scales = 'free_y') +
  scale_color_manual(values = c('black', brewer.pal(5, 'RdBu')[1])) +
  scale_x_date(expand = c(0, 0), date_labels = '%b %d') +
  labs(x = NULL, y = NULL) +
  guides(color = 'none', fill = 'none') +
  facet.theme +
  theme(text = element_text(size = 14), panel.spacing = unit(1, 'lines'))
g2
# ggsave(width = 4, height = 5.5, dpi = 172,
#   file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/Buckley_model_parameter_time_series.png')
g2 + facet_wrap(~ Parameter, ncol = 2, strip.position = 'top',
    labeller = label_parsed, scales = 'free_y')
# ggsave(width = 6.5, height = 3.5, dpi = 172,
#   file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/Buckley_model_parameter_time_series_wide.png')


# TODO Add standard deviation of WP measurements
df.params %>%
  mutate(date = as.Date(sprintf('2023-%03d', DOY.rel), '%Y-%j')) %>%
  left_join(select(df.out, DOY.rel, hour.rel, WP_sd), by = c('DOY.rel', 'hour.rel')) %>%
  group_by(DOY.rel, date) %>%
  filter(n() > 3) %>%
  summarize(
    psi.source = first(psi.source),
    psi0 = WP_m[1],
    max.psi = max(WP_m)) %>%
  gather(key = Metric, value = value, -DOY.rel:-date) %>%
  mutate(Metric = if_else(str_detect(Metric, 'source'), 'psi[source]',
    if_else(str_detect(Metric, 'max'), 'psi[max]', 'psi[min]'))) %>%
  left_join(select(df.met, DOY.rel, ppt_mm), by = 'DOY.rel') %>%
  mutate(precip = if_else(ppt_mm > 0, date, NA)) %>%
ggplot(mapping = aes(x = date, y = value)) +
  geom_vline(aes(xintercept = precip), color = 'darkblue', linetype = 'dashed') +
  geom_line(aes(color = Metric)) +
  scale_x_date(expand = c(0, 0)) +
  scale_color_manual(values = c(brewer.pal(5, 'RdBu')[c(5,1)], 'black'),
    labels = parse(text = c('psi[max]', 'psi[min]', 'psi[source]'))) +
  guides(color = guide_legend('Water Potential')) +
  labs(x = NULL, y = 'Water Potential (MPa)') +
  theme_bw() +
  theme(text = element_text(size = 14),
    legend.position = 'top',
    legend.margin = margin(0, 0, -0.3, 0, 'cm'))
# ggsave(width = 6, height = 4.5, dpi = 172,
#   file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/Buckley_model_WP_source_sink_time_series.png')


df.params %>%
  mutate(date = as.Date(sprintf('2023-%03d', DOY.rel), '%Y-%j')) %>%
  group_by(DOY.rel, date) %>%
  filter(n() > 3) %>%
  summarize(psi.source = first(psi.source),
    psi0 = WP_m[1],
    max.psi = max(WP_m)) %>%
  mutate(diff = max.psi - psi.source) %>%
  left_join(select(df.met, DOY.rel, ppt_mm), by = 'DOY.rel') %>%
  mutate(precip = if_else(ppt_mm > 0, date, NA)) %>%
ggplot(mapping = aes(x = date, y = diff)) +
  geom_vline(aes(xintercept = precip), color = 'darkblue', linetype = 'dashed') +
  geom_bar(stat = 'identity') +
  scale_x_date(expand = c(0, 0)) +
  scale_color_manual(values = c(brewer.pal(5, 'RdBu')[c(5,1)], 'black'),
    labels = parse(text = c('psi[max]', 'psi[min]', 'psi[source]'))) +
  guides(color = guide_legend('Water Potential')) +
  labs(x = NULL, y = 'Difference in Water Potential (MPa)',
    subtitle = TeX('Difference: $\\Psi_{max} - \\Psi_{source}$')) +
  theme_bw() +
  theme(
    legend.position = 'top',
    legend.margin = margin(0, 0, -0.3, 0, 'cm'))
# ggsave(width = 6, height = 4, dpi = 172,
#   file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/Buckley_model_WP_diff_max_minus_source_time_series.png')


# How does tau very with WP? ###################################################

df.params %>%
  group_by(DOY.rel) %>%
  filter(n() > 3) %>%
  summarize(tau = first(tau),
    psi.source = first(psi.source),
    psi0 = WP_m[1],
    max.psi = max(WP_m),
    mean.psi = mean(WP_m)) %>%
  gather(key = Parameter, value = value, -DOY.rel:-tau) %>%
ggplot(mapping = aes(x = tau, y = value)) +
  geom_point() +
  facet_wrap(~ Parameter)
