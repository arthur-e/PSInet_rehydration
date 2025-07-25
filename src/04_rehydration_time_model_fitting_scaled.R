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

# Normalize and then scale the data
df.emp.scaled <- df.emp %>%
  select(DOY, hour, DOY.rel, hour.rel, datetime, WP_m, WP_sd) %>%
  left_join(mutate(df.emp.agg,
      min.psi.prev = lag(min.psi, 1),
      max.psi.prev = lag(max.psi, 1)),
    by = 'DOY.rel') %>%
  group_by(DOY.rel) %>%
  arrange(DOY.rel, hour.rel) %>%
  mutate(WP.norm = (WP_m - WP_m[1] - (min.psi.prev - WP_m[1])),
    WP.scaled = -10 + (WP.norm / max(WP.norm))) %>%
  ungroup()

df.emp.scaled %>%
  select(DOY.rel, max.psi.prev, hour.rel, WP_m, WP.scaled, WP.norm) %>%
  gather(key = field, value = value, -DOY.rel:-hour.rel) %>%
  filter(DOY.rel %in% c(100, 101, 140, 141, 160, 161)) %>%
ggplot(mapping = aes(x = hour.rel, y = value)) +
  geom_line() +
  facet_grid(field ~ DOY.rel, scale = 'free_y') +
  labs(x = 'Hour of Day (Relative)', y = 'Water Potential (MPa)') +
  facet.theme
ggsave(width = 7.5, height = 5, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_time_series_data_scaled_example.png')

# Sun sets at 19h00, rises after 05h00 local time, but stomata tend to close
#   early due to stress, so we'll start the clock at ??h00
HOUR.START <- 12
HOUR.END <- 8
df.out <- df.emp.scaled %>%
  select(datetime, WP_m, WP.norm, WP.scaled, DOY, hour) %>%
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
