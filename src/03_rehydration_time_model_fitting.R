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

# I want each day to start and end at 10h00 local, as this is (likely) before
#   the minimum daily water potential
df.emp <- df %>%
  select(datetime, date, DOY, hour, WP_m, WP_sd) %>%
  # Roll back "midnight" to 10h00 local
  mutate(datetime.rel = datetime - dhours(10),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
  filter(DOY.rel > 90) %>%
  arrange(datetime)

# Compute daily min, max water potential and when that is achieved
df.emp.agg <- df.emp %>%
  group_by(DOY.rel) %>%
  summarize(min.psi = min(WP_m),
    max.psi = max(WP_m),
    when.min = (hour.rel[which.min(WP_m)] + 10) %% 24,
    when.max = (hour.rel[which.max(WP_m)] + 10) %% 24) %>%
  # Then, compute how long it took to recharge
  mutate(diff = (24 - when.min) + when.max)

# Sun sets at 19h00, rises after 05h00 local time, but stomata tend to close
#   early due to stress, so we'll start the clock at ??h00
HOUR.START <- 12
HOUR.END <- 8
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
# df.out %>%
#   filter(DOY.rel %% 10 == 0) %>%
# ggplot(mapping = aes(x = datetime, y = WP_m)) +
#   geom_line() +
#   facet_wrap(~ DOY.rel, scales = 'free', ncol = 2)


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
model.mixed <- function (t, pd = -2, tau = 4, y.int = -4, slope = 0.5, f.linear = 0.5) {
  .exp <- model.exp(t, pd, tau)
  .linear <- y.int + (t * slope)
  return(((1 - f.linear) * .exp) + (f.linear * .linear))
}

# Our objective function, calculating RMSE based on exponential model predictions
rmse.model.exp <- function (params, df = df.out) {
  t <- df$t
  pd <- params[1]
  tau <- params[2]
  psi <- model.exp(t, pd, tau)
  return(rmse(df$WP_m, psi))
}

# Our objective function, calculating RMSE based on mixed model predictions
rmse.model.mixed <- function (params, df = df.out) {
  t <- df$t
  pd <- params[1]
  tau <- params[2]
  y.int <- params[3]
  slope <- params[4]
  f.linear <- params[5]
  psi <- model.mixed(t, pd, tau, y.int, slope, f.linear)
  return(rmse(df$WP_m, psi))
}

# The idealized function we're trying to fit
plot.function(function(t) {model.exp(t, tau = 2)}, from = 0, to = 24, bty = 'n', xlab = 'Hour (Local Time)',
  ylab = 'Water Potential (MPa)', lwd = 2, xaxt = 'n', col = 'blue')
plot.function(function(t) {model.exp(t, tau = 6)}, from = 0, to = 24, bty = 'n', xlab = 'Hour (Local Time)',
  ylab = 'Water Potential (MPa)', lwd = 2, xaxt = 'n', col = 'red', add = TRUE)
legend(16, -3.5, legend = c('tau=2', 'tau=6'),
  lwd = 2, col = c('blue', 'red'))
axis(1, at = seq(0, 20, 5), labels = c(16, 20, 24, 4, 8))

# The mixed model
plot.function(model.mixed, from = 0, to = 24, bty = 'n', xlab = 'Hour (Local Time)',
  ylab = 'Water Potential (MPa)', lwd = 2, xaxt = 'n')
axis(1, at = seq(0, 20, 5), labels = c(16, 20, 24, 4, 8))


# Fitting exponential model parameters with optim ##############################

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
# ggsave(width = 6.5, height = 3, dpi = 172, bg = 'white',
#   file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_RMSE_barplot_alt.png')

i <- 1
g1 <- NULL
# for (doy in df.params$DOY.rel) {
for (doy in 100:111) {
  g0 <- df.out %>%
    filter(DOY.rel == doy) %>%
    mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
  ggplot(mapping = aes(x = hour.rel)) +
    geom_line(aes(y = WP_m), color = 'black', linewidth = 0.6) +
    # Change X-axis so that the "day" starts at 10h00 local time
    scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
      labels = function (x) { (x + 14) %% 24 }) +
    scale_y_continuous(limits = c(-6, -2)) +
    labs(subtitle = sprintf('DOY=%03d', doy),
      x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)')

  g0 <- g0 + with(filter(df.params, DOY.rel == doy),
      geom_function(fun = model.exp, n = 1000, xlim = c(0, with(df.out, max(t))),
        args = list(pd = psi0, tau = tau), linetype = 'dashed', color = 'red'))
  g0 <- g0 + geom_hline(aes(yintercept = psi0), color = 'blue',
    linetype = 'dashed', data = filter(df.params, DOY.rel == doy))

  if (i == 1) {
    g1 <- g0
  } else {
    g1 <- g1 + g0
  }
  i <- 1 + 1
}

g1 + plot_layout(nrow = 3, axis_titles = 'collect')
# ggsave(width = 7, height = 4.5, dpi = 172,
#   file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_fits.png')

# It's apparent that there are days when a steady-state is reached (DOY 109-111)
#   but the exponential model doesn't reach a plateau; i.e., the predicted
#   psi0 isn't reached before WP starts to decline again
df.params %>% filter(DOY.rel %in% 100:111)

# We can quantify the steady state in a couple of ways...
df.steady.state <- df.params %>%
  mutate(steady.state.lag = psi0 >= lag(psi0, 1)) %>%
  left_join(df.out, by = 'DOY.rel') %>%
  group_by(DOY.rel) %>%
  arrange(DOY.rel, hour.rel) %>%
  summarize(steady.state.lag = first(steady.state.lag),
    steady.state = last(psi0) >= last(WP_m))

# Virtually every day, the estimated steady-state psi0 is higher than
#   yesterday's actual WP
with(df.steady.state, table(steady.state))

# But the estimated psi0 might be higher or lower than yesterday's estimate
with(df.steady.state, table(steady.state.lag))

# There's no pattern to the days where psi0 >= lag(psi, 1)
require(viridis)
df.steady.state %>%
ggplot(mapping = aes(x = DOY.rel, y = steady.state.lag)) +
  geom_point(shape = 1, size = 2) +
  scale_color_gradientn(colors = magma(10, direction = -1)) +
  labs(x = 'Day of Year', y = 'Rehydrated?',
    subtitle = latex2exp::TeX('Based on (argmin($\\Psi$) + $5\\tau$) earlier than sunrise'),
    color = 'Rehydration\nTime (hrs)') +
  theme_dark() +
  theme(plot.background = element_blank(),
    legend.background = element_blank())

# Nothing significant here, I think: As the rehydration curve becomes more
#   linear, the purely exponential model flattens early (i.e., tau -> 0),
#   and the right end of the exponential meets or undershoots that actual WP
df.out %>%
  left_join(df.params, by = 'DOY.rel') %>%
  mutate(psi.hat = model.exp(t, psi0, tau)) %>%
  group_by(DOY.rel) %>%
  filter(n() == 41) %>%
  summarize(diff.realized = first(psi0) - max(WP_m),
    diff.estimated = first(psi0) - max(psi.hat)) %>%
ggplot(mapping = aes(x = DOY.rel, y = diff.realized)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Day of Year', y = 'Difference (MPa)')


# Fitting mixed model parameters with optim ####################################

# A data frame with the unique DOYs
df.doy <- df.out %>%
  group_by(DOY.rel) %>%
  summarize()
fit.params <- matrix(nrow = length(df.doy$DOY.rel), ncol = 6)
fit.scores <- matrix(nrow = nrow(fit.params), ncol = 1)
for (i in 1:nrow(fit.params)) {
  doy <- df.doy$DOY.rel[i]
  df.sub <- filter(df.out, DOY.rel == doy)
  if (nrow(df.sub) == 0) next() # Skip empty subsets
  wrap.rmse <- function (params) {
    rmse.model.mixed(params, df = df.sub)
  }
  # Gradually increase the initial estimate of the linear fraction ("f.linear")
  initial.f.linear <- doy / 180
  result <- optim(c(-2, 2, -4, 0.5, initial.f.linear), wrap.rmse,
    method = 'L-BFGS-B', lower = c(-9, 1, -10, -100, 0), upper = c(0, 24, -1, 100, 1))
    # control = list(pgtol = 0.01))
  fit.params[i,1] <- doy
  fit.params[i,2:6] <- result$par
  fit.scores[i,] <- wrap.rmse(result$par)
}

df.params <- bind_cols(
  as.data.frame(fit.params) %>%
    rename(DOY.rel = V1, psi0 = V2, tau = V3, y.int = V4, slope = V5, f.linear = V6),
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
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_mixed_RMSE_barplot_forcing_linear_by_DOY.png')

df.params.good <- df.params %>%
  filter(DOY.rel > 142, DOY.rel != 181, RMSE < 0.2)
df.params.bad <- df.params %>%
  filter(DOY.rel > 142, DOY.rel != 181, RMSE > 0.5)

# Good results
i <- 1
g1 <- NULL
for (doy in df.params.good$DOY.rel[1:6]) {
  g0 <- df.out %>%
    filter(DOY.rel == doy) %>%
    mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
  ggplot(mapping = aes(x = hour.rel)) +
    geom_line(aes(y = WP_m), color = 'black', linewidth = 0.6) +
    # Change X-axis so that the "day" starts at 10h00 local time
    scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
      labels = function (x) { (x + HOUR.START) %% 24 }) +
    scale_y_continuous(expand = c(0, 0), limits = c(-9, -5)) +
    labs(subtitle = sprintf('DOY=%03d', doy),
      x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)')

  g0 <- g0 + with(filter(df.params.good, DOY.rel == doy),
      geom_function(fun = model.mixed, n = 1000, xlim = c(0, 17),
        args = list(pd = psi0, tau = tau, y.int = y.int, slope = slope, f.linear = f.linear),
        linetype = 'dashed', color = 'red'))

  if (i == 1) {
    g1 <- g0
  } else {
    g1 <- g1 + g0
  }
  i <- 1 + 1
}

g1 + plot_layout(nrow = 2, axis_titles = 'collect')
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_mixed_fits_low-RMSE.png')

# Bad results
i <- 1
g1 <- NULL
for (doy in sample(df.params.bad$DOY.rel, size = 6)) {
  g0 <- df.out %>%
    filter(DOY.rel == doy) %>%
    mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
  ggplot(mapping = aes(x = hour.rel)) +
    geom_line(aes(y = WP_m), color = 'black', linewidth = 0.6) +
    # Change X-axis so that the "day" starts at 10h00 local time
    scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
      labels = function (x) { (x + HOUR.START) %% 24 }) +
    scale_y_continuous(expand = c(0, 0), limits = c(-9, -5)) +
    labs(subtitle = sprintf('DOY=%03d', doy),
      x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)')

  g0 <- g0 + with(filter(df.params.bad, DOY.rel == doy),
      geom_function(fun = model.mixed, n = 1000, xlim = c(0, 17),
        args = list(pd = psi0, tau = tau, y.int = y.int, slope = slope, f.linear = f.linear),
        linetype = 'dashed', color = 'red'))

  if (i == 1) {
    g1 <- g0
  } else {
    g1 <- g1 + g0
  }
  i <- 1 + 1
}

g1 + plot_layout(nrow = 2, axis_titles = 'collect')
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_mixed_fits_high-RMSE.png')


# Fitting mixed model, sweeping the mixing ratio ###############################

# A data frame with the unique DOYs
df.doy <- df.out %>%
  group_by(DOY.rel) %>%
  summarize()
fit.params <- matrix(nrow = length(df.doy$DOY.rel), ncol = 6)
fit.scores <- matrix(nrow = nrow(fit.params), ncol = 1)

# Here, we do an inner loop over different mixing ratios, selecting for the
#   minimum RMSE
for (i in 1:nrow(fit.params)) {
  doy <- df.doy$DOY.rel[i]
  df.sub <- filter(df.out, DOY.rel == doy)
  if (nrow(df.sub) == 0) next() # Skip empty subsets

  # Gradually increase the initial estimate of the linear fraction ("f.linear")
  .scores <- matrix(nrow = 21, ncol = 1)
  .params <- matrix(nrow = 21, ncol = 5)
  j <- 1
  for (m.ratio in seq(0, 1, 0.05)) {
    wrap.rmse <- function (params) {
      params <- c(params, m.ratio)
      rmse.model.mixed(params, df = df.sub)
    }
    result <- optim(c(-2, 2, -4, 0.5), wrap.rmse,
      method = 'L-BFGS-B', lower = c(-9, 1, -10, -100), upper = c(0, 24, -1, 100))
    .params[j,1:5] <- c(result$par, m.ratio)
    .scores[j,] <- wrap.rmse(result$par)
    j <- j + 1
  }

  fit.params[i,1] <- doy
  fit.params[i,2:6] <- .params[which.min(.scores),]
  fit.scores[i,] <- .scores[which.min(.scores),]
}

df.params2 <- bind_cols(
  as.data.frame(fit.params) %>%
    rename(DOY.rel = V1, psi0 = V2, tau = V3, y.int = V4, slope = V5, f.linear = V6),
  as.data.frame(fit.scores) %>%
    rename(RMSE = V1))

df.params2 %>%
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
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_mixed_RMSE_barplot_forcing_linear_by_DOY_minimize-RMSE.png')

i <- 1
g1 <- NULL
for (doy in c(c(177, 179),
    sample(filter(df.params2, DOY.rel != 181)$DOY.rel, size = 7))) {
  g0 <- df.out %>%
    filter(DOY.rel == doy) %>%
    mutate(DOY.rel = sprintf('DOY=%03d', DOY.rel)) %>%
  ggplot(mapping = aes(x = hour.rel)) +
    geom_line(aes(y = WP_m), color = 'black', linewidth = 0.6) +
    # Change X-axis so that the "day" starts at 10h00 local time
    scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
      labels = function (x) { (x + HOUR.START) %% 24 }) +
    # scale_y_continuous(expand = c(0, 0), limits = c(-9, -5)) +
    labs(subtitle = sprintf('DOY=%03d', doy),
      x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)')

  g0 <- g0 + with(filter(df.params2, DOY.rel == doy),
      geom_function(fun = model.mixed, n = 1000, xlim = c(0, 17),
        args = list(pd = psi0, tau = tau, y.int = y.int, slope = slope, f.linear = f.linear),
        linetype = 'dashed', color = 'red'))

  if (i == 1) {
    g1 <- g0
  } else {
    g1 <- g1 + g0
  }
  i <- 1 + 1
}

g1 + plot_layout(nrow = 3, axis_titles = 'collect')
ggsave(width = 7, height = 4.5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_mixed_fits_minimize-RMSE.png')

df.params2 %>%
  select(DOY.rel, `Fraction Linear` = f.linear, tau) %>%
  gather(key = Parameter, value = value, -DOY.rel) %>%
  left_join(df.met, by = 'DOY.rel') %>%
  mutate(rained = if_else(ppt_mm > 0, DOY.rel, NA)) %>%
ggplot(mapping = aes(x = DOY.rel, y = value)) +
  # geom_bar(stat = 'identity') +
  geom_line(color = 'darkred') +
  geom_vline(aes(xintercept = rained), color = 'darkblue', linetype = 'dashed') +
  facet_wrap(~ Parameter, scale = 'free_y', ncol = 2) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = 'Day of Year', y = 'Parameter Value') +
  facet.theme
ggsave(width = 7, height = 3, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_RC_circuit_model_mixed_fits_minimizing_RMSE_showing_parameters.png')


# Five-tau analysis ############################################################

require(viridis)
df.params2 %>%
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
