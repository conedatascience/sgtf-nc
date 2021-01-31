# Pulling Rt values
#

rt <- nccovid::pull_estimates()

data.table::fwrite(rt, here::here("data-raw", "rt-raw.csv"))
