library(dplyr, warn.conflicts = F)
library(tidyr)
library(ggplot2)
library(lubridate)

JESSICA.CSV <- '~/Downloads/all_auto_joined.csv'

df <- read.csv(JESSICA.CSV) %>%
  select(indiv.id = individual_id, plot.id = plot_id, sensor_id, date, time,
    timezone, latitude = latitude_wgs84, longitude = longitude_wgs84,
    water_potential_mean, water_potential_n, genus, specific_epithet, organ,
    canopy_position) %>%
  # mutate(date = as.Date(date)) %>%
  mutate(date = ymd_hms(sprintf('%s %s', date, time)))

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
  group_by(indiv.id, genus, specific_epithet, date) %>%
  summarize(water_potential_mean = mean(water_potential_mean)) %>%
  mutate(DOY = as.integer(format(date, '%j'))) %>%
  filter(DOY %in% 150:180) %>%
  filter(indiv.id == 6) %>%
  View
  mutate(hour = hour(date)) %>%
ggplot(mapping = aes(x = hour, y = water_potential_mean)) +
  geom_line() +
  facet_wrap(~ DOY, scales = 'free_x') +
  theme(legend.position = 'top')
