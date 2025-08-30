library(dplyr, warn.conflicts = F)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(lubridate, warn.conflicts = F)
library(patchwork)
library(latex2exp)

JESSICA.CSV <- '~/Downloads/PSInet/all_auto_joined.csv'
MET.CSV <- '/home/arthur.endsley/Workspace/NTSG/projects/Y2026_PSInet/data/PSInet_DB_select_auto_sites_met_data.csv.gz'

facet.theme <- theme_linedraw() +
  theme(strip.background = element_blank(),
    strip.text = element_text(color = 'black', face = 'bold'))

df <- read.csv(JESSICA.CSV) %>%
  select(dataset.name = dataset_name, indiv.id = individual_id, plot.id = plot_id,
    sensor.id = sensor_id, date, time, timezone, latitude = latitude_wgs84,
    longitude = longitude_wgs84, water_potential_mean, water_potential_n, genus,
    specific_epithet, organ, canopy_position) %>%
  filter(!is.na(timezone)) %>%
  filter(str_detect(timezone, 'America')) %>%
  # Because rowwise() is SO SLOW, we have to do this the hard way
  # rowwise() %>%
  # mutate(datetime = ymd_hms(sprintf('%s %s', date, time)),
  #   datetime = with_tz(datetime, tzone = timezone)) %>%
  mutate(datetime = ymd_hms(sprintf('%s %s', date, time)),
    hour = hour(datetime) + minute(datetime) / 60)

df.tzones <- df %>%
  group_by(dataset.name) %>%
  summarize(tzone = first(timezone))

# Unfortunately, none of the met data is for these sites
df.met <- read.csv(MET.CSV) %>%
  rename(dataset.name = site_id) %>%
  mutate(datetime.UTC = ymd_hms(sprintf('%s %s:00', date_UTC, time_UTC))) %>%
  left_join(df.tzones, by = 'dataset.name') %>%
  # NOTE: Phoenix is on standard time all the time; for simplicity, convert
  #   all times to a single GMT offset
  # NOTE: We're doing it this way because rowwise() is too slow
  mutate(datetime = if_else(str_detect(tzone, 'Phoenix|Denver'), datetime.UTC - dhours(7),
    if_else(str_detect(tzone, 'Los_Angeles'), datetime.UTC - dhours(8),
      if_else(str_detect(tzone, 'New_York'), datetime.UTC - dhours(5), NA))),
    hour = hour(datetime) + minute(datetime) / 60) %>%
  select(-contains('UTC'), -tzone) %>%
  select(dataset.name, datetime, hour, everything()) %>%
  # Compute AVP (Gates 1980, Biophysical Ecology, p.311),
  #   SVP (August-Roche-Magnus), and then VPD
  mutate(avp = (QV2M * PS) / (0.622 + (0.379 * QV2M)),
    svp = 610.7 * exp((17.38 * (T2M - 273.15)) / (239 + (T2M - 273.15))),
    VPD = svp - avp) %>%
  mutate(VPD = if_else(VPD < 0, 0, VPD),
    PAR.MJm2 = 0.45 * (0.0036 * (24 / 1) * SWGDN)) %>%
  select(-avp, -svp, -PS, -QV2M, -T2M, -SWGDN)

# They are ALL "stem" measurements (makes sense, stem psychrometers)
# At each site, they are ALL same "canopy_position"
# Each dataset is only one genus
# Each dataset is only one species
# Datasets may have multiple plots
# Plots may have multiple individuals
# Individuals may have multiple sensors

df.agg <- df %>%
  # Aggregate across sensors and individuals
  group_by(dataset.name, plot.id, genus, datetime, hour) %>%
  summarize(WP = mean(water_potential_mean, na.rm = T)) %>%
  filter(WP > -10) %>%
  mutate(DOY = as.integer(format(datetime, '%j'))) %>%
  group_by(dataset.name, plot.id, DOY) %>%
  filter(n() >= 48) %>%
  mutate(site.id = sprintf('%s %s (%s)', dataset.name, plot.id, genus)) %>%
  ungroup()

# Taking it as granted that the "time" is already local because the drawdown
#   minima look reasonable
df.agg %>%
ggplot(mapping = aes(x = hour, y = WP)) +
  geom_line(aes(group = DOY), alpha = 0.1) +
  geom_smooth() +
  facet_wrap(~ site.id, scales = 'free_y') +
  theme(legend.position = 'top') +
  labs(x = 'Hour (Local Time, MST)')

# Sun sets at 19h00, rises after 05h00 local time, but stomata tend to close
#   early due to stress, so we'll start the clock at ??h00
HOUR.START <- 10
HOUR.END <- 8
df.agg.filt <- df.agg %>%
  filter(dataset.name != 'Kno_1', site.id != 'Ker_5 Holter 1') %>%
  left_join(select(df.met, -hour), by = c('dataset.name', 'datetime')) %>%
  select(site.id, datetime, DOY, hour, WP, VPD, PAR.MJm2) %>%
  arrange(datetime, hour) %>%
  group_by(site.id, DOY) %>%
  # Create an identifier for each unique evening period
  # This is faster and can be done in a loop, but should be checked for
  #   generality
  mutate(datetime.rel = datetime - dhours(HOUR.START),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
  arrange(site.id, datetime.rel, DOY.rel, hour.rel) %>%
  # Interpolate from hourly to half-hourly
  fill(VPD:PAR.MJm2, .direction = 'downup') %>%
  ungroup()

# Diagnostic
df.agg.filt %>%
ggplot(mapping = aes(x = hour.rel, y = WP)) +
  geom_line(aes(group = DOY.rel, color = DOY.rel)) +
  facet_wrap(~ site.id, scales = 'free_y')

# Filter each time series to *just* those hours *after* the minimum WP is reached
df.modeling <- df.agg.filt %>%
  group_by(site.id, DOY.rel) %>%
  # Start at minimum WP
  filter(hour.rel >= hour.rel[which.min(WP)]) %>%
  # End at shortly (3 hours) after maximum WP
  filter(hour.rel <= hour.rel[which.max(WP)] + 3) %>%
  # Transform hours into a number of hours past HOUR.START
  mutate(t = hour.rel - min(hour.rel)) %>%
  filter(!is.na(var(WP)), var(WP) > 0)

# df.modeling %>%
# ggplot(mapping = aes(x = hour.rel, y = WP)) +
#   geom_line(aes(group = DOY.rel, color = DOY.rel)) +
#   facet_wrap(~ site.id, scales = 'free_y')

# df.modeling %>%
#   filter(site.id == 'Ker_3 1 (Quercus)') %>%
# ggplot(mapping = aes(x = hour.rel, y = WP)) +
#   geom_line() +
#   facet_wrap(~ DOY.rel, scales = 'free_y')


# Modeling framework ###########################################################

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

rmse <- function (observed, predicted) {
  n <- length(predicted)
  # Make sure predicted and observed are the same length!
  stopifnot(n == length(observed))
  return((1/n) * sum(sqrt((predicted - observed)^2), na.rm = T))
}

rmse.model.stomatal <- function (params, df = df.clean.min) {
  ts <- df$t
  vpd <- df$VPD
  irr <- df$PAR.MJm2
  psi0 <- df$WP[1]
  tau <- params[1]
  psi.source <- params[2]
  k <- params[3]
  c1 <- params[4]
  c2 <- params[5]
  pred <- model.stomatal(ts, vpd, irr, psi0, tau, psi.source, k, c1, c2)
  return(with(df, rmse(WP, pred)))
}


# Fitting the model ############################################################

# A data frame with the unique DOYs
df.doy <- df.modeling %>%
  group_by(site.id, DOY.rel) %>%
  summarize()
fit.params <- matrix(nrow = nrow(df.doy), ncol = 2 + 5)
fit.scores <- matrix(nrow = nrow(fit.params), ncol = 1)
for (i in 1:nrow(fit.params)) {
  doy <- df.doy$DOY.rel[i]
  site.name <- df.doy$site.id[i]
  df.sub <- filter(df.modeling, site.id == site.name, DOY.rel == doy)
  if (nrow(df.sub) == 0) next() # Skip empty subsets
  wrap.rmse <- function (params) {
    rmse.model.stomatal(params, df = df.sub)
  }
  # Parameters are: tau, psi.source, k, c1, c2
  result <- optim(c(6, -1, 100, 0.1, 0.2), wrap.rmse, method = 'L-BFGS-B',
    lower = c(1, -6, 1, 0, 0.1), upper = c(24, 0, 1000, 1, 1))
    # control = list(pgtol = 0.01))
  fit.params[i,1] <- site.name
  fit.params[i,2] <- doy
  fit.params[i,3:7] <- result$par
  fit.scores[i,] <- wrap.rmse(result$par)
}

df.params0 <- bind_cols(
  as.data.frame(fit.params) %>%
    rename(site.id = V1, DOY.rel = V2, tau = V3, psi.source = V4, k = V5, c1 = V6, c2 = V7) %>%
    mutate(DOY.rel = as.integer(DOY.rel),
      across(tau:c2, .fns = as.numeric)),
  as.data.frame(fit.scores) %>%
    rename(RMSE = V1))

df.params <- df.params0 %>%
  # NOTE: We need dynamic inputs for this function
  right_join(df.modeling, by = c('site.id', 'DOY.rel'))

df.params.filt <- df.params %>%
  group_by(site.id, DOY.rel) %>%
  filter(n() >= 24) %>%
  summarize(site.id = first(site.id), DOY.rel = first(DOY.rel))


# Plotting results #############################################################

N.PLOTS <- 15
idx <- as.integer(runif(N.PLOTS, 1, nrow(df.params.filt)))
g1 <- NULL
for (i in idx) {
  site.name <- df.params.filt[i,]$site.id
  doy <- df.params.filt[i,]$DOY.rel
  g0 <- df.modeling %>%
    filter(site.id == site.name, DOY.rel == doy) %>%
    mutate(DOY.rel = sprintf('%s DOY=%03d', site.name, DOY.rel)) %>%
  ggplot(mapping = aes(x = hour.rel)) +
    geom_line(aes(y = WP), color = 'black', linewidth = 0.6) +
    # Change X-axis so that the "day" starts at 10h00 local time
    # scale_x_continuous(expand = c(0, 0), breaks = seq(2, 24, 6),
    #   labels = function (x) { (x + 14) %% 24 }) +
    # scale_y_continuous(limits = c(-9, -2)) +
    labs(subtitle = sprintf('%s DOY=%03d', site.name, doy),
      x = 'Hour of Day (Local Time)', y = 'Water Potential (MPa)')

  g0 <- g0 + with(filter(df.params, site.id == site.name, DOY.rel == doy),
      geom_function(fun = model.stomatal, n = max(t) - min(t),
        xlim = c(min(hour.rel), max(hour.rel)),
        # NOTE: Static inputs are indexed at their first value
        args = list(vpd = VPD, irr = PAR.MJm2, psi0 = WP[1],
          tau = tau[1], psi.source = psi.source[1], k = k[1], c1 = c1[1], c2 = c2[1]),
        linetype = 'dashed', color = 'red'))
  g0 <- g0 + geom_hline(aes(yintercept = psi.source), color = 'blue',
    linetype = 'dashed', data = filter(df.params, site.id == site.name, DOY.rel == doy))

  if (i == idx[1]) {
    g1 <- g0
  } else {
    g1 <- g1 + g0
  }
  i <- 1 + 1
}

g1 + plot_layout(nrow = 3, axis_titles = 'collect')


# Model inference ##############################################################

df.params0 %>%
ggplot(mapping = aes(x = site.id, y = RMSE)) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal()

# TODO Going to need datasets with longer records, in terms of multiple days

g2 <- df.params0 %>%
  # NOTE: Every study (so far) takes place in less than one year
  mutate(date = as.Date(sprintf('2023-%03d', DOY.rel), '%Y-%j')) %>%
  select(site.id, DOY.rel, date, tau) %>%
ggplot(mapping = aes(x = date, y = tau)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ site.id, scales = 'free_y') +
  scale_color_manual(values = c('black', brewer.pal(5, 'RdBu')[1])) +
  scale_x_date(expand = c(0, 0), date_labels = '%b %d') +
  labs(x = NULL, y = NULL) +
  guides(color = 'none', fill = 'none') +
  facet.theme +
  theme(text = element_text(size = 14), panel.spacing = unit(1, 'lines'))
g2
