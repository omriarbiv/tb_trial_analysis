library(parallel)

set.seed(4)
nsim <- 1e5
ncores <- detectCores() - 1

# TDR power analysis ----------------------------------------------------------
tdr_data <- read.csv("tdr.csv")

# Find rate of adverse events for patients receiving SDR
sdr_meta_rate <- with(tdr_data,
                      meta::metarate(
                        event = sdr_event,
                        time = sdr_days * n_sdr,
                        studlab = study,
                        method = "Inverse"
                  ))

mean_sdr_rate <- exp(sdr_meta_rate$TE.random)

# Double rate of adverse events for patients receiving TDR
mean_tdr_rate <- mean_sdr_rate * 2

calc_power_tdr <- function(nsim, n_sdr, n_tdr, tdr_days, 
                           days_post_tdr, sdr_days, 
                           sdr_rate = mean_sdr_rate,
                           tdr_rate = mean_tdr_rate) {
  
  # Simulate nsim total events expected in HRZE arm
  sdr_sim <- replicate(nsim, sum(rpois(n_sdr, sdr_rate * sdr_days)))
  
  # Simulate nsim total events expected in TDR arm
  tdr_sim <- replicate(nsim, sum(rpois(n_tdr, tdr_rate * tdr_days) +
                                   rpois(n_tdr, sdr_rate * days_post_tdr)))
  
  # Calculate incident-rate ratio in each simulation
  irr <- lapply(1:nsim, \(.x) {
    epitools::rateratio(c(sdr_sim[.x], tdr_sim[.x],
                          n_sdr * sdr_days, n_tdr * (tdr_days + days_post_tdr)
    ))
  })
  
  irr_p <- sapply(1:nsim, \(.x) irr[[.x]]$p.value[2, 1])
  irr_val <- sapply(1:nsim, \(.x) irr[[.x]]$measure[2, 1])
  
  # Power of IRR with p value < 0.05 that found TDR had more events than SDR
  power <- sum(irr_p < 0.05 & irr_val > 1) / nsim
  return(power)
}

t1 <- Sys.time()
tdr_data$power <- with(tdr_data,
                       mclapply(1:nrow(tdr_data), \(.x) {
                         calc_power_tdr(nsim, n_sdr[.x], n_tdr[.x],
                                        tdr_days[.x], 
                                        after_tdr_sdr_days[.x],
                                        sdr_days[.x])
                       }, mc.cores = ncores)) |> 
                      simplify2array()
Sys.time() - t1

# Find max power in each study comparison
tdr_by_study <- split(tdr_data, tdr_data$study)
tdr_max_power <- lapply(tdr_by_study, \(x) max(x$power))

# Make summary table
tdr_power <- data.frame(
  study = unique(tdr_data$study), 
  power = unsplit(tdr_max_power, unique(tdr_data$study))
)

# PaMZ power analysis ---------------------------------------------------------

pamz_data <- read.csv("pamz.csv")

rhze_meta_rate <- with(pamz_data,
                       meta::metarate(
                         event = rhze_event,
                         time = rhze_days * rhze_n,
                         studlab = study,
                         method = "Inverse"
                       ))

mean_rhze_rate <- exp(rhze_meta_rate$TE.random)
mean_pamz_rate <- mean_rhze_rate * 2

# This function is slightly different than `calc_power_tdr` because unlike 
# the TDR studies, participants randomized to PaMZ were not given RHZE afterwards 
# (whereas many of the studies of TDR gave TDR for some time than SDR afterwards).

calc_power_pamz <- function(nsim, rhze_n, pamz_n, 
                            rhze_days, pamz_days,
                            rhze_rate = mean_rhze_rate,
                            pamz_rate = mean_pamz_rate) {
  rhze_sim <- replicate(nsim, sum(rpois(rhze_n, rhze_rate * rhze_days)))
  pamz_sim <- replicate(nsim, sum(rpois(pamz_n, pamz_rate * pamz_days)))
  
  irr <- lapply(1:nsim, \(.x) {
    epitools::rateratio(c(rhze_sim[.x], pamz_sim[.x],
                          rhze_n * rhze_days, pamz_n * pamz_days))
  })
  
  irr_p <- sapply(1:nsim, \(.x) irr[[.x]]$p.value[2, 1])
  irr_val <- sapply(1:nsim, \(.x) irr[[.x]]$measure[2, 1])
  power <- sum(irr_p < 0.05 & irr_val > 1) / nsim
  return(power)
}

pamz_data$power <- with(pamz_data,
                       mclapply(1:nrow(pamz_data), \(.x) {
                         calc_power_pamz(nsim, rhze_n[.x], pamz_n[.x],
                                         rhze_days[.x], pamz_days[.x])
                       }, mc.cores = ncores)) |> 
                    simplify2array()

# Find max power in each study comparison
pamz_by_study <- split(pamz_data, pamz_data$study)
pamz_max_power <- lapply(pamz_by_study, \(x) max(x$power))

# Make summary table
pamz_power <- data.frame(
  study = unique(pamz_data$study), 
  power = unsplit(pamz_max_power, unique(pamz_data$study))
)
