# Purpose: Simulate Range of Values for NC and Guilford County

library(tidyverse)
library(data.table)
source("R/run_ode.R")

# pull base data ----------------------------------------------------------

raw_rt <- data.table::fread(here::here("data-raw", "rt-raw.csv"))

## Population
pop_nc <- nccovid::nc_population$july_2020[nccovid::nc_population$county=="STATE"]
pop_guilford <- nccovid::nc_population$july_2020[nccovid::nc_population$county=="Guilford"]

## Recovered
incubation_time <- 5
reporting_delay <- 2
incubation_offset <- incubation_time + reporting_delay
underascertainment <- 2.5
nc_introduction <- as.Date("2021-01-23")
guilford_introduction <- as.Date("2021-01-28")

cases_raw <- nccovid::pull_covid_summary()
cases_county_raw <- nccovid::get_covid_state(select_county = "Guilford")

recovered_nc <-sum(cases_raw[date<=nc_introduction-incubation_offset]$daily_cases)*underascertainment
recovered_guilford <- sum(cases_county_raw[date<=nc_introduction-incubation_offset]$cases_daily) *underascertainment

infected_days <- 8

infected_nc <- cases_raw[date==nc_introduction-incubation_offset]$daily_cases * infected_days
infected_guilford <- cases_county_raw[date==nc_introduction-incubation_offset]$cases_daily * infected_days

exposed_nc <- infected_nc *.5
exposed_guilford <- infected_guilford *.5

## vaccines
vax_raw <- nccovid::get_vaccinations()

vax_nc <- vax_raw[date<=nc_introduction-incubation_offset][date==max(date)][,sum(dose_1,na.rm=TRUE)]
vax_guilford <- vax_raw[county=="Guilford"][date<=nc_introduction-incubation_offset][date==max(date)][,sum(dose_1,na.rm=TRUE)]

vax_nc_rate <- round(147575/7) # Week of Feb 2 state total dose 1 allotments
vax_guilford_rate <- round(5775/7) # Week of Feb 2 dose 1 allotsments, no mass

## R estimates

reff_nc <- raw_rt[county=="North Carolina" & date==nc_introduction-incubation_offset]
reff_guilford <- raw_rt[county=="Guilford" & date==guilford_introduction-incubation_offset]

## Susceptibles
nc_s = pop_nc - vax_nc - recovered_nc - exposed_nc - infected_nc
guilford_s = pop_guilford - vax_guilford - recovered_guilford - exposed_guilford - infected_guilford


simulation_grid_nc <-data.frame(region = "North Carolina",
                             s_in = nc_s,
                             e_in = exposed_nc,
                             i_in = infected_nc,
                             r_in = recovered_nc,
                             r_vax_in = vax_nc,
                             #tranmission_increase = 1.5,
                             vax_rate_in= vax_nc_rate,
                             VE = .95,
                             start_date  = nc_introduction - incubation_offset) %>%
  crossing(R_eff = c(reff_nc$bottom, reff_nc$median, reff_nc$top),
           transmission_increase = c(1,1.5))

simulation_grid_guilford <-data.frame(region = "Guilford",
                                s_in = guilford_s,
                                e_in = exposed_guilford,
                                i_in = infected_guilford,
                                r_in = recovered_guilford,
                                vax_rate_in= vax_guilford_rate,
                                VE = .95,
                                start_date  = guilford_introduction - incubation_offset) %>%
  crossing(R_eff = c(reff_nc$bottom, reff_nc$median, reff_nc$top),
           transmission_increase = c(1,1.5))


guilford_sim <- pmap(simulation_grid_guilford, run_ode)
nc_sim <- pmap(simulation_grid_nc, run_ode)


# extract reff ------------------------------------------------------------
combined_wide_reff <- guilford_sim %>%
  map("summary_data") %>%
  bind_rows() %>%
  mutate(scenario = ifelse(R_eff_in == max(R_eff_in), "Upper",
                           ifelse(R_eff_in == min(R_eff_in),"Lower", "Median"))) %>%
  bind_rows({

    nc_sim %>%
      map("summary_data") %>%
      bind_rows() %>%
      mutate(scenario = ifelse(R_eff_in == max(R_eff_in), "Upper",
                               ifelse(R_eff_in == min(R_eff_in),"Lower", "Median")))

  }) %>%
  select(date, reff_evolution, region, scenario, transmission_increase) %>%
  pivot_wider(names_from = scenario, values_from = reff_evolution)


# save outputs ------------------------------------------------------------

write_rds(list(guilford_sim= guilford_sim,
               nc_sim = nc_sim,
               combined_wide_reff = combined_wide_reff),
          here::here("output", "simulated-evolutions.rds"))


# combined_wide_reff%>%
#   filter(!is.na(Median)) %>%
#   mutate(variant = ifelse(transmission_increase==1.5, "With B1.1.7", "No B1.1.7")) %>%
#   ggplot(aes(date, Median, colour = as.factor(variant)))+
#   geom_line()+
#   geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey80", alpha = .2)+
#   geom_hline(yintercept = 1)+
#   labs(
#     title = "Guilford County Reproduction Number"
#   )+
#   theme_minimal()+
#   facet_wrap(~region)

