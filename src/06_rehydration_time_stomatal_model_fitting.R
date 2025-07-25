library(dplyr, warn.conflicts = F)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(lubridate, warn.conflicts = F)
library(bbmle, warn.conflicts = F)
library(patchwork)

JESSICA.CSV <- '~/Downloads/PSInet/SRER_LATR_pdd_Apr_June.csv'
MET.CSV <- '~/Downloads/PSInet/neon_atmdaily.csv'
SWGDN.CSV <- '/home/arthur.endsley/Workspace/NTSG/projects/Y2026_PSInet/data/202050725_Jessica_Creosote_site_met_data_MERRA_SWGDN.csv'

facet.theme <- theme_linedraw() +
  theme(strip.background = element_blank(),
    strip.text = element_text(color = 'black', face = 'bold'))


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
  select(datetime, WP_m, PAR.MJm2, DOY, hour) %>%
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
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
  left_join(select(df.met, datetime.rel, Dmean), by = 'datetime.rel')

# Diagnostics
df.out %>%
  filter(DOY.rel %% 10 == 0) %>%
ggplot(mapping = aes(x = datetime, y = WP_m)) +
  geom_line() +
  facet_wrap(~ DOY.rel, scales = 'free', ncol = 2)


# Model fitting framework ######################################################

rmse <- function (observed, predicted) {
  n <- length(predicted)
  # Make sure predicted and observed are the same length!
  stopifnot(n == length(observed))
  return((1/n) * sum(sqrt((predicted - observed)^2), na.rm = T))
}

model.stomatal <- function (
    ts, vpd = 1000, irr = 0.5, psi0 = -3, tau = 4, psi.source = -2, k = 1, c1 = 1, c2 = 1) {
  wp <- numeric(length = length(ts))
  wp[1] <- psi0 + (1/tau) * (psi.source - psi0 - (c1 * irr * vpd) / ((irr + k) * (1 + c2 * vpd)))
  for (t in ts[2:length(ts)]) {
    wp[t] <- wp[t-1] + (1/tau) * (psi.source - wp[t-1] - (c1 * irr * vpd) / ((irr + k) * (1 + c2 * vpd)))
  }
  return(wp)
}

plot(model.stomatal(1:24), type = 'l', bty = 'n', ylab = 'Water Potential (MPa)',
  lwd = 2, xaxt = 'n', ylim = c(-3, -2))
abline(h = -2)
axis(1, at = seq(0, 20, 5), labels = c(16, 20, 24, 4, 8))
