library(data.table)
library(ggplot2)
library(dplyr)

odin_model <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S[]) <- max(S[i] - n_SE[i] - n_vax_removed,0)
  update(E[]) <- E[i] + n_SE[i] - n_EI[i]
  update(I[]) <- I[i] + n_EI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  update(RVax[]) <- RVax[i] + min(n_vax_removed,S[i])
  update(Reff[]) <- S[i]/N[i]*transmission_multiplier*beta/gamma
  ## Individual probabilities of transition:
  p_SE[] <- 1 - exp(-beta * I[i] / N[i]*transmission_multiplier)
  p_EI <- 1 - exp(-delta)
  p_IR <- 1 - exp(-gamma)
  n_vax_removed <- rbinom(vax_rate, .95)
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE[] <- rbinom(S[i], p_SE[i])
  n_EI[] <- rbinom(E[i], p_EI)
  n_IR[] <- rbinom(I[i], p_IR)

  ## Total population size
  N[] <- S[i] + E[i] + I[i] + R[i] + RVax[i]

  ## Initial states:
  initial(S[]) <- S_ini
  initial(E[]) <- E_ini
  initial(I[]) <- I_ini
  initial(R[]) <- R_ini
  initial(RVax[]) <- R_vax_ini

  initial(Reff[]) <- beta/gamma

  transmission_multiplier <- transmission*1/(1 + exp(-(-4.648421 +0.07249007 *time)))+(1-1/(1 + exp(-(-4.648421 +0.07249007 *time))))
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  E_ini <- user(2)
  I_ini <- user(1)
  R_ini <- user(0)
  R_vax_ini <- user(0)
  Reff <- user(1)
  beta <- user(0.2)
  gamma <- user(0.1)
  delta <- user(.3)
  vax_rate <- user(0)
  transmission <- user(1.5)

  ## Number of replicates
  nsim <- user(100)
  dim(N) <- nsim
  dim(S) <- nsim
  dim(E) <- nsim
  dim(I) <- nsim
  dim(R) <- nsim
  dim(Reff) <- nsim
  dim(RVax) <- nsim
  dim(p_SE) <- nsim
  dim(n_SE) <- nsim
  dim(n_EI) <- nsim
  dim(n_IR) <- nsim
  update(time) <- time + 1
  initial(time) <- 0
  })


pop <- nccovid::nc_population[nccovid::nc_population$county=="STATE"]$july_2020

cases <- nccovid::pull_covid_summary()

observed_cases <- cases[date<as.Date("2021-01-23")][,.(cases = sum(daily_cases))]

current_infection <- tail(cases,1)$daily_cases*10

likely_cases <- round(observed_cases*2.5)

likely_exposed <- round(current_infection*.25)

vax <- nccovid::get_vaccinations()[date==max(date)]

vaccinated <- vax[,sum(dose_1)]

s_in <- pop - likely_exposed - current_infection - likely_cases$cases-vaccinated
R_eff <- 1.00
gamma_in = 1/10
beta_in <- R_eff/(s_in/pop)*gamma_in
x <- odin_model(S_ini = s_in,
                E_ini = likely_exposed,
                I_ini = current_infection,
                R_ini = likely_cases$cases,
                R_vax_ini = vaccinated,
                beta = beta_in,
                delta = .2,
                gamma = .1,
                transmission = 1.5,
                vax_rate = 21428)
.1/.1
x_res <- x$run(0:150)

x_res <- as.data.table(x_res)

date_box <- data.table(time = 0:150,
                       date =seq.Date(as.Date("2021-01-23"), by = "day",
                                      length.out = 151))

x_long <- melt(x_res, id.vars = c("step", "time"))
x_long[,compartment := stringr::str_extract(variable, "[:alpha:]+")]
x_long[,iteration := stringr::str_extract(variable, "[:digit:]+?")]
x_long[,id:=rep(1:100,times = 906)]
x_long <- merge(x_long, date_box, by = "time", all.x = TRUE)
x_syn <- dcast(x_long, formula = time + date + id~compartment , value.var = "value")
x_syn <- x_syn[order(time)][,new_cases:=shift(S,n = 1)-S, by = "id"]
x_syn[,perc_S:= S/pop]
theme_set(theme_bw())
x_syn %>%
  mutate(date = as.Date(date)) %>%
  ggplot(aes(date,I, group = as.factor(id)))+
  geom_line(alpha = .02)

x_syn %>%
  filter(time!=0) %>%
  filter(id %in% 1:100) %>%
  ggplot(aes(time,new_cases/4))+
  geom_line(alpha = .05, aes( group = id))+
  geom_smooth()

x_syn %>%
  filter(time!=0) %>%
  ggplot(aes(date,Reff))+
  geom_line(alpha = .05, aes( group = id))+
  geom_smooth()+
  geom_hline(yintercept = 1, lty = 2, colour = "orange")
