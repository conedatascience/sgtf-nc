#' Run Stochastic SEIR Model + Vaccinations with Variant Effect
#'
#'@param R_eff a double, >0, the effective reproduction number
#'@param delta_in a double, the
#'@param gamma_in
#'
run_ode <- function(R_eff = 1,
                    delta_in = 1/6,
                    gamma_in=1/10,
                    s_in=1000,
                    e_in=2,
                    i_in=1,
                    r_in=0,
                    r_vax_in=0,
                    pop = s_in+e_in + i_in+ r_in+ r_vax,
                    tranmission_increase=1,
                    vax_rate_in=0,
                    start_date = "2021-01-23",
                    horizon = 150,
                    VE = .95){
  require(data.table)
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
    n_vax_removed <- rbinom(vax_rate, VE)
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

    transmission_multiplier <- transmission*exp(-exp(-(-4.648421 +0.07249007 *time)))+(1-exp(-exp(-(-4.648421 +0.07249007 *time))))
    ## User defined parameters - default in parentheses:
    S_ini <- user(1000)
    E_ini <- user(2)
    I_ini <- user(1)
    R_ini <- user(0)
    R_vax_ini <- user(0)
    R_eff <- user(1)
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

  beta_in <- R_eff/(s_in/pop)*gamma_in
  x <- odin_model(S_ini = s_in,
                  E_ini = e_in,
                  I_ini = i_in,
                  R_ini = r_in,
                  R_vax_ini = r_vax_in,
                  beta = beta_in,
                  delta = delta_in,
                  gamma = gamma_in,
                  transmission = tranmission_increase,
                  vax_rate = vax_rate_in)

  x_res <- x$run(0:horizon)

  x_res <- data.table::as.data.table(x_res)

  date_box <- data.table::data.table(time = 0:horizon,
                         date =seq.Date(as.Date(start_date), by = "day",
                                        length.out = (horizon+1)))

  x_long <- melt(x_res, id.vars = c("step", "time"))
  x_long[,compartment := stringr::str_extract(variable, "[:alpha:]+")]
  x_long[,iteration := stringr::str_extract(variable, "[:digit:]+?")]
  x_long[,id:=rep(1:100,times = nrow(date_box)*6)]
  x_long <- merge(x_long, date_box, by = "time", all.x = TRUE)
  x_syn <- dcast(x_long, formula = time + date + id~compartment , value.var = "value")
  x_syn <- x_syn[order(time)][,new_cases:=shift(S,n = 1)-S, by = "id"]
  x_syn <- x_syn[,perc_S:= S/pop]

  list(x_syn = x_syn, x_long = x_long)
}
