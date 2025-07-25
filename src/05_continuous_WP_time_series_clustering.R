library(dplyr, warn.conflicts = F)
library(tidyr)
library(stringr)
library(ggplot2)
library(proxy, warn.conflicts = F)
library(dtw, warn.conflicts = F)
library(dtwclust, warn.conflicts = F)
library(lubridate, warn.conflicts = F)

JESSICA.CSV <- '~/Downloads/PSInet/SRER_LATR_pdd_Apr_June.csv'
MET.CSV <- '~/Downloads/PSInet/neon_atmdaily.csv'

facet.theme <- theme_linedraw() +
  theme(strip.background = element_blank(),
    strip.text = element_text(color = 'black', face = 'bold'))

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
    WP.scaled = 1 + (WP.norm / max(WP.norm))) %>%
  ungroup()

# Sun sets at 19h00, rises after 05h00 local time, but stomata tend to close
#   early due to stress, so we'll start the clock at ??h00
HOUR.START <- 12
HOUR.END <- 8
df.emp.byrow <- df.emp.scaled %>%
  group_by(DOY.rel) %>%
  filter(n() == 48, DOY.rel > 91) %>%
  filter(hour >= HOUR.START | hour <= HOUR.END) %>%
  arrange(DOY.rel, hour.rel) %>%
  head
  select(DOY.rel, WP.scaled) %>%
  ungroup()

# Create a list structure, for compatibility with tsclust()
i <- 1
wp.list <- list()
for (doy in unique(df.emp.byrow$DOY.rel)) {
  wp.list[[i]] <- with(filter(df.emp.byrow, DOY.rel == doy), WP.scaled)
  i <- i + 1
}
plot(wp.list[[1]], type = 'l')
plot(wp.list[[60]], type = 'l')


# Dynamic time warping: Dynamic barycentric averaging ##########################

# This next set of outputs uses dynamic time warping proper.

# With DTW barycentric averaging (DBA)
dtw.dba.5 <- tsclust(wp.list, type = 'part', k = 5L,
  seed = 99, distance = 'dtw_basic', centroid = 'dba')
table(dtw.dba.5@cluster)
plot(dtw.dba.5, type = 'centroids',
  clus = 1:5, linetype = 'solid', size = 1)


# Dynamic time warping: Soft DTW Centroids #####################################

# Here, we go beyond averaging and try a few different ways to define centroids
dtw.sdtw.3 <- tsclust(wp.list, type = 'part', k = 3L,
  seed = 99, distance = 'dtw_basic', centroid = 'sdtw_cent')
table(dtw.sdtw.3@cluster)
plot(dtw.sdtw.3, type = 'centroids',
  clus = 1:3, linetype = 'solid', size = 1)

dtw.sdtw.5 <- tsclust(wp.list, type = 'part', k = 5L,
  seed = 99, distance = 'dtw_basic', centroid = 'sdtw_cent')
table(dtw.sdtw.5@cluster)
plot(dtw.sdtw.5, type = 'centroids',
  clus = 1:5, linetype = 'solid', size = 1)


# Dynamic time warping: PAM ####################################################

# Here, we try partitioning around medoids
dtw.pam.3 <- tsclust(wp.list, type = 'part', k = 3L,
  seed = 99, distance = 'dtw_basic', centroid = 'pam')
table(dtw.pam.3@cluster)
plot(dtw.pam.3, type = 'centroids',
  clus = 1:3, linetype = 'solid', size = 1)
data.frame(DOY = unique(df.emp.byrow$DOY.rel), Cluster = dtw.pam.3@cluster) %>%
ggplot(mapping = aes(DOY, Cluster)) +
  geom_point()

# N=4 shows that we're just starting to fit noise
dtw.pam.4 <- tsclust(wp.list, type = 'part', k = 4L,
  seed = 99, distance = 'dtw_basic', centroid = 'pam')
table(dtw.pam.4@cluster)
plot(dtw.pam.4, type = 'centroids',
  clus = 1:4, linetype = 'solid', size = 1) +
  # Change X-axis so that the "day" starts at 10h00 local time
  scale_x_continuous(expand = c(0, 0), breaks = seq(2, 48, 6),
    labels = function (x) { ((x / 2) + HOUR.START) %% 24 }) +
  labs(title = 'DTW with Partitioning around Medoids',
    x = 'Hour of Day (Local Time)', y = 'Scaled WP') +
  facet.theme
ggsave(width = 6, height = 4, dpi = 172,
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_DTW_with_PAM_and_scaled_WP_N4_centroids.png')

data.frame(DOY = unique(df.emp.byrow$DOY.rel), Cluster = dtw.pam.4@cluster) %>%
  filter(Cluster <= 3) %>%
ggplot(mapping = aes(DOY, Cluster)) +
  geom_point() +
  scale_y_continuous(breaks = 1:3) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank())
ggsave(width = 4, height = 2.5, dpi = 172, bg = 'white',
  file = '~/Workspace/NTSG/projects/Y2026_PSInet/outputs/rehydration_time_disequilibrium/20250707_DTW_with_PAM_and_scaled_WP_N4_cluster_time_series.png')

dtw.pam.5 <- tsclust(wp.list, type = 'part', k = 5L,
  seed = 99, distance = 'dtw_basic', centroid = 'pam')
table(dtw.pam.5@cluster)
plot(dtw.pam.5, type = 'centroids',
  clus = 1:5, linetype = 'solid', size = 1)
data.frame(DOY = unique(df.emp.byrow$DOY.rel), Cluster = dtw.pam.5@cluster) %>%
ggplot(mapping = aes(DOY, Cluster)) +
  geom_point()
