# Purpose: Simulate Range of Values for NC and Guilford County

library(tidyverse)
library(data.table)


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

vax_nc_rate <- round(150000/7)
vax_guilford_rate <- round(4975/7)

## R estimates

reff_nc <- raw_rt[county=="North Carolina" & date==nc_introduction-incubation_offset]
reff_guilford <- raw_rt[county=="Guilford" & date==guilford_introduction-incubation_offset]

simulation_grid <-











