library(parallel)

set.seed(4)
nsim <- 1e6
irr_scenario <- c(1.5, 2, 3)
ncores <- 4

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
  power <- (sum(irr_p < 0.05 & irr_val > 1) / nsim) * 100
  return(power)
}


for (irr in irr_scenario) {
	tdr_data[[paste0("power", irr)]] <- with(
		tdr_data, mclapply(1:nrow(tdr_data), \(.x) {
				 calc_power_tdr(nsim, n_sdr[.x], n_tdr[.x], tdr_days[.x], after_tdr_sdr_days[.x],
												sdr_days[.x], sdr_rate = mean_sdr_rate, 
												tdr_rate = mean_sdr_rate*irr)
						 }, mc.cores = ncores, mc.set.seed = TRUE)) |> 
						simplify2array()
}

# Find max power in each study comparison
tdr_max_power <- replicate(3, NA, F)
names(tdr_max_power) <- paste0("power", irr_scenario)

tdr_by_study <- split(tdr_data, tdr_data$study)
 
for(irr in irr_scenario) { 
	s <- paste0("power", irr)
	tdr_max_power[[s]] <- sapply(tdr_by_study, \(x) max(x[[s]]))
}

# Make summary table
tdr_power <- data.frame(
  study = unique(tdr_data$study), 
  power1.5 = unsplit(tdr_max_power[["power1.5"]], unique(tdr_data$study)),
  power2 = unsplit(tdr_max_power[["power2"]], unique(tdr_data$study)),
  power3 = unsplit(tdr_max_power[["power3"]], unique(tdr_data$study))
)

write.csv(tdr_power, "tdr_power_table.csv")

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
  power <- (sum(irr_p < 0.05 & irr_val > 1) / nsim) * 100
  return(power)
}

for (irr in irr_scenario) {
	pamz_data[[paste0("power", irr)]] <- with(
		pamz_data, mclapply(1:nrow(pamz_data), \(.x) {
                        calc_power_pamz(nsim, rhze_n[.x], pamz_n[.x],
                                        rhze_days[.x], pamz_days[.x],
                                        rhze_rate = mean_rhze_rate, 
																				pamz_rate = mean_rhze_rate*irr)
                       }, mc.cores = ncores, mc.set.seed = TRUE)) |> 
                    simplify2array()
}

# Find max power in each study comparison
pamz_max_power <- replicate(3, NA, F)
names(pamz_max_power) <- paste0("power", irr_scenario)

pamz_by_study <- split(pamz_data, pamz_data$study)
 
for(irr in irr_scenario) { 
	s <- paste0("power", irr)
	pamz_max_power[[s]] <- sapply(pamz_by_study, \(x) max(x[[s]]))
}

# Make summary table
pamz_power <- data.frame(
  study = unique(pamz_data$study), 
  power1.5 = unsplit(pamz_max_power[["power1.5"]], unique(pamz_data$study)),
  power2 = unsplit(pamz_max_power[["power2"]], unique(pamz_data$study)),
  power3 = unsplit(pamz_max_power[["power3"]], unique(pamz_data$study))
)

write.csv(pamz_power, "pamz_power.csv")
