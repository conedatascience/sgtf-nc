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

























