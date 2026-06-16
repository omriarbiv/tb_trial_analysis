d <- read.csv("75 percentile.csv")

# Process
dp <- d[d$adult_study, ]
dp <- dp[!dp$study == "paton", ]
dp <- dp[!is.na(dp$study), ]
dp <- dp[order(dp$year, dp$study_name, decreasing = T), ]

dp$calc_75 <- qnorm(0.75, dp$mean, dp$sd)
dp$comb_75 <- ifelse(!is.na(dp$iqr_high), dp$iqr_high, dp$calc_75)

# Adding star for calculated numbers
dp$study_name[!is.na(dp$calc_75)] <- 
  paste0(dp$study_name[!is.na(dp$calc_75)], "*") 

# Fix names (adding et al.,)
dp$study_name <- sub("^([[:alpha:]-]+) ([0-9]{4})", "\\1 et al., \\2", 
                     dp$study_name)

# Fix names (change exp to experimental)
dp$study_name <- sub("exp", "experimental", dp$study_name)

# Removing study with ranges
dp2 <- dp[!is.na(dp$comb_75), ]

# Add colours (closest to global age [55])
# Generate 110 colours with 55 in the middle
col_range <- c(hcl.colors(55), hcl.colors(55, rev = T))
# age_range <- 1:(55*2)
dp2$col <- col_range[round(dp2$comb_75)]

pdf("75-2.pdf", 8, 8)
with(dp2, 
     dotchart(comb_75,
              labels = study_name,
              xlab = "75th percentile of participant age",
              groups = factor(regimen),
              xlim = c(30, 81), bg = "grey10",
              xaxp = c(30, 70, 5),
              lwd = 0.1,
              pt.cex = 1.4,
              lcolor = "grey50",
              frame.plot = F))
abline(v = 55, col = "grey10", lwd = 2)
abline(v = seq(40, 70, by = 10), lty = 2)
dev.off()