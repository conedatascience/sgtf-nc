library(tidyverse)
library(data.table)

rt_history <- fread(here::here("data-raw", "rt-raw.csv"))
sgtf <- fread(here::here("data-raw", "new-variant.csv"))

guilford <- rt_history[county=="Guilford"]

guilford_sgtf <- merge(guilford, sgtf, by = "date", all.x = TRUE)

guilford_sgtf[,new_r:=median*beta_multipler]
guilford_sgtf[,date_n:=as.numeric(date)]

plot(new_r~date, data = guilford_sgtf, ylim = c(0,1.2))
abline(h = 1)

fit<- mgcv::gam(median~s(date_n, k = 8), data = guilford_sgtf)
