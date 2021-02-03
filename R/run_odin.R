odin_model <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S[]) <- S[i] - n_SE[i] - vax_rate
  update(E[]) <- E[i] + n_SE[i] - n_EI[i]
  update(I[]) <- I[i] + n_EI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i] + vax_rate

  ## Individual probabilities of transition:
  p_SE[] <- 1 - exp(-beta * I[i] / N[i]*transmission_multiplier)
  p_EI <- 1 - exp(-delta)
  p_IR <- 1 - exp(-gamma)

  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE[] <- rbinom(S[i], p_SE[i])
  n_EI[] <- rbinom(E[i], p_EI)
  n_IR[] <- rbinom(I[i], p_IR)

  ## Total population size
  N[] <- S[i] + E[i] + I[i] + R[i]

  ## Initial states:
  initial(S[]) <- S_ini
  initial(E[]) <- E_ini
  initial(I[]) <- I_ini
  initial(R[]) <- R_ini

  transmission_multiplier <- transmission*exp(-exp(-(-4.648421 +0.07249007 *time)))+(1-exp(-exp(-(-4.648421 +0.07249007 *time))))
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  E_ini <- user(2)
  I_ini <- user(1)
  R_ini <- user(0)
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
  dim(p_SE) <- nsim
  dim(n_SE) <- nsim
  dim(n_EI) <- nsim
  dim(n_IR) <- nsim
  update(time) <- time + 1
  initial(time) <- 0
  })

x <- odin_model(transmission = 4, vax_rate = 1)
x_res <- x$run(0:100)
sir_col_transp <- paste0(sir_col, "66")
sir_col <- c("#8c8cd9", "#cc0044", "#999966", "#A20BE2")
matplot(x_res[, 1], x_res[, -1], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = rep(sir_col_transp, each = 100), lty = 1)
x_res[, 4]

