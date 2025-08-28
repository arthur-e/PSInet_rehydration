library(dplyr, warn.conflicts = F)
library(tidyr)
library(ggplot2)
library(lubridate)

JESSICA.CSV <- '~/Downloads/PSInet/all_auto_joined.csv'

df <- read.csv(JESSICA.CSV) %>%
  select(indiv.id = individual_id, plot.id = plot_id, sensor_id, date, time,
    timezone, latitude = latitude_wgs84, longitude = longitude_wgs84,
    water_potential_mean, water_potential_n, genus, specific_epithet, organ,
    canopy_position) %>%
  # Necessary to do this rowwise because tzone argument can't use a column value;
  #   ACTUALLY, still doesn't work
  # rowwise() %>%
  # mutate(datetime = ymd_hms(sprintf('%s %s', date, time), tzone = timezone)) %>%
  mutate(datetime = ymd_hms(sprintf('%s %s', date, time))) %>%
  mutate(hour = hour(datetime) + minute(datetime) / 60)

with(df, table(genus, specific_epithet))

with(df %>%
  filter(genus == 'Juniperus', specific_epithet == 'osteosperma'), table(indiv.id, plot.id))

with(df %>%
  filter(genus == 'Juniperus', specific_epithet == 'osteosperma') %>%
  filter(indiv.id == 6), table(sensor_id))

# TODO We need to convert hour of day to local time
df %>%
  filter(genus == 'Juniperus', specific_epithet == 'osteosperma') %>%
  filter(water_potential_mean > -10) %>%
  mutate(DOY = as.integer(format(datetime, '%j'))) %>%
  filter(DOY %in% 150:180) %>%
  filter(indiv.id == 6) %>%
ggplot(mapping = aes(x = hour, y = water_potential_mean)) +
  geom_line(aes(color = sensor_id)) +
  facet_wrap(~ DOY, scales = 'free_x') +
  theme(legend.position = 'top') +
  labs(x = 'Hour (Local Time, MST)')
