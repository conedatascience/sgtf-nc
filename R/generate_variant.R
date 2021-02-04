library(tidyverse)
library(data.table)


# pull data ---------------------------------------------------------------

var_concern <- "https://raw.githubusercontent.com/nicholasdavies/newcovid/master/fitting_data/sgtf-2021-01-18.csv"
dat <- data.table::fread(var_concern)
setDT(dat)
dat[,perc:=sgtf/(sgtf+other)]


# visualize data ----------------------------------------------------------

dat %>%
  ggplot(aes(date,perc,colour = nhs_name))+
  geom_line()+
  labs(
    title = "B1.1.7 Percentage of Variants Circulating",
    colour = "NHS Region",
    y = "Percent",
    x = NULL
  )+
  scale_y_continuous(labels = scales::percent)


# format for fitting ------------------------------------------------------

dat_2 <- dat %>%
  group_by(nhs_name) %>%
  mutate(t = 1:n()) %>%
  as.data.table() %>%
  .[,n:=sgtf+other] %>%
  as.data.frame() %>%
  filter(perc<1&perc>0)

library(rstanarm)

fit_new_variant <- stan_betareg(
  perc ~ t + nhs_name, data = dat_2,link = "logit",cores = 2,
  algorithm = "sampling" # just for speed of example
)
summary(fit_new_variant)
(fit_coef <- broom.mixed::tidyMCMC(fit_new_variant))

fitted_int <- fit_coef$estimate[fit_coef$term=="(Intercept)"]
fitted_slope <- fit_coef$estimate[fit_coef$term=="t"]
# (Intercept)
# -4.648421
# t
# 0.07249007
# simulating --------------------------------------------------------------

days <- 365

pred <- exp(-exp(-(fitted_int +fitted_slope*1:days)))
pred <-arm::invlogit(fitted_int +fitted_slope*1:days)

dat_2 %>%
  ggplot(aes(t,perc))+
  geom_line(aes(colour = nhs_name))+
  geom_line(data = data.frame(pred = pred,
                              t = 1:days),
            aes(t, pred), inherit.aes = FALSE, lty = "dashed")

first_detected <- as.Date("2021-01-29")

first_import <- first_detected - 6

detected_overlay <- data.frame(pred = pred,
                               t = 1:days)
readr::write_csv(fit_coef, here::here("output", "new-variant-coefficients.csv"))
readr::write_csv(detected_overlay, here::here("output", "fitted-uk-variant.csv"))
readr::write_csv(dat_2, here::here("data-raw", "uk-data.csv"))

new_variant_impact <-data.frame(pred = pred,
                                t = 1:days) %>%
  mutate(date = seq.Date(as.Date(first_import), by = "day", length.out = days)) %>%
  select(date, t, perc_variant = pred) %>%
  as.data.table() %>%
  .[,beta_multipler:=(1.5*perc_variant)+1*(1-perc_variant)]

write_csv(new_variant_impact, here::here("data-raw", "new-variant.csv"))
