library(dplyr, warn.conflicts = F)
library(tidyr)
library(ggplot2)
library(lubridate)
library(stringr)

JESSICA.CSV <- '~/Downloads/PSInet/all_auto_joined.csv'

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
  mutate(datetime = ymd_hms(sprintf('%s %s', date, time))) %>%
  mutate(hour = hour(datetime) + minute(datetime) / 60)

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
  mutate(site.id = sprintf('%s %s', dataset.name, plot.id)) %>%
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
df.out <- df.agg %>%
  filter(dataset.name != 'Kno_1', site.id != 'Ker_5 Holter 1') %>%
  select(site.id, datetime, WP, DOY, hour) %>%
  arrange(datetime, hour) %>%
  filter(hour >= HOUR.START | hour <= HOUR.END) %>%
  group_by(site.id, DOY) %>%
  # Transform hours into a number of hours past HOUR.START
  mutate(t = if_else(hour < HOUR.START, hour + 24, hour) - HOUR.START) %>%
  # Create an identifier for each unique evening period
  # This is faster and can be done in a loop, but should be checked for
  #   generality
  mutate(datetime.rel = datetime - dhours(HOUR.START),
    DOY.rel = as.integer(format(datetime.rel, '%j')),
    hour.rel = hour(datetime.rel) + minute(datetime.rel) / 60) %>%
  ungroup()

df.out %>%
ggplot(mapping = aes(x = hour.rel, y = WP)) +
  geom_line(aes(group = DOY.rel, color = DOY.rel)) +
  facet_wrap(~ site.id, scales = 'free_y')

# Need to compute both the timing of minimum WP but also the t index
# df.clean.min <- df.out %>%
#   group_by(site.id, DOY.rel) %>%
#   filter(hour.rel >= hour.rel[which.min(WP_m)])

