library(dplyr, warn.conflicts = F)
library(tidyr)
library(ggplot2)
library(lubridate)
library(bbmle)

JESSICA.CSV <- '~/Downloads/SRER_LATR_pdd_Apr_June.csv'

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

df.out <-
result <- optim(c(-2, 2), wrap.rmse,
  method = 'L-BFGS-B', lower = c(-9, 1), upper = c(0, 24))

require(broom)
require(purrr)
df %>%
  nest(data = -DOY) %>%
  mutate(fit = map(data, ~ lm(WP_m ~
