library(epitrix)
library(EpiEstim)

tt <- seq(0, 30, 0.1)
r <- 0.1
xmax <- max(tt)
ymax <- max(exp(r*tt))


data("Flu2009")
R_flu <- epitrix::r2R0(r, Flu2009$si_distr)
mean_GT_flu <- sum(Flu2009$si_distr * (seq_along(Flu2009$si_distr)-1))

data("Measles1861")
R_measles <- epitrix::r2R0(r, Measles1861$si_distr)
mean_GT_measles <- sum(Measles1861$si_distr * (seq_along(Measles1861$si_distr)-1))

## Alpha variant example with uncertainty on generation time
# Volz et al. Nature 2021 - central assumption is assuming a gamma-distributed
# generation time with mean 6.4 days and CV of 0.6612.
# Alternative uses a GT which is twice shorter on average.
# Davies et al. GT of either 5.5 days (and SD of 2.1 days)
# or mean 3.6 days (and SD of 3.1 days)

## Volz version
#mean_GT_alpha_long <- 6.4
#mean_GT_alpha_short <- .5 * mean_GT_alpha_long
#cv_GT_alpha <- 0.6612
#GT_distr_alpha_long <- EpiEstim::discr_si(tt, mu = mean_GT_alpha_long,
#                    sigma = mean_GT_alpha_long * cv_GT_alpha)
#GT_distr_alpha_short <- EpiEstim::discr_si(tt, mu = mean_GT_alpha_short,
#                    sigma = mean_GT_alpha_short * cv_GT_alpha)
## Davies version
GT_distr_alpha_long <- EpiEstim::discr_si(tt, mu = 5.5, sigma = 2.1)
GT_distr_alpha_short <- EpiEstim::discr_si(tt, mu = 3.6, sigma = 3.1)

# daily growth rate between 0 and 6% https://assets.publishing.service.gov.uk/media/600ac5538fa8f56553673aaf/s1011-spi-m-o-consensus-statement.pdf
r_alpha <- 0.03
R_alpha_long <- epitrix::r2R0(r_alpha, GT_distr_alpha_long)
R_alpha_short <- epitrix::r2R0(r_alpha, GT_distr_alpha_short)

pdf("figs/R_versus_r_flu_measles.pdf", width = 10, height = 5)

################################################################################
### Create panel layout
################################################################################

par(mfrow = c(1, 2),
    las = 0, xaxs="i", yaxs="i", mar = c(5, 5, 2, 2))

################################################################################
### Flu-like natural history
################################################################################

## exponential growth curve obtained with flu-like natural history
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax),
     main = paste0("R = ", signif(R_flu, 3), "; Flu-like natural history"))
mtext("Time", 1, 3)
mtext("Incidence", 2, 3)
abline(v = mean_GT_flu*seq_len(30), lty = 2, col = "grey")

################################################################################
### same r but Measles-like natural history
################################################################################

## exponential growth curve obtained with measles-like natural history
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax),
     main = paste0("R = ", signif(R_measles, 3), "; Measles-like natural history"))
mtext("Time", 1, 3)
mtext("Incidence", 2, 3)
abline(v = mean_GT_measles*seq_len(30), lty = 2, col = "grey")

dev.off()


pdf("figs/R_versus_r_covid_alpha_variant.pdf", width = 10, height = 5)

################################################################################
### Create panel layout
################################################################################

par(mfrow = c(1, 2),
    las = 0, xaxs="i", yaxs="i", mar = c(5, 5, 2, 2))

################################################################################
### long GT
################################################################################

## exponential growth curve obtained with long GT
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax),
     main = paste0("R = ", signif(R_alpha_long, 2), "; Long generation time"))
mtext("Time", 1, 3)
mtext("Incidence", 2, 3)
abline(v = mean_GT_alpha_long*seq_len(30), lty = 2, col = "grey")

################################################################################
### short GT
################################################################################

## exponential growth curve obtained with measles-like natural history
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax),
     main = paste0("R = ", signif(R_alpha_short, 2), "; Short generation time"))
mtext("Time", 1, 3)
mtext("Incidence", 2, 3)
abline(v = mean_GT_alpha_short*seq_len(30), lty = 2, col = "grey")

dev.off()

################################################################################
### Numbers
################################################################################

# Flu
epitrix::R02AR(R0 = R_flu)
1 - 1/R_flu

# Measles
epitrix::R02AR(R0 = R_measles)
1 - 1/R_measles


# Alpha long GT
epitrix::R02AR(R0 = R_alpha_long)
1 - 1/R_alpha_long

# Alpha short GT
epitrix::R02AR(R0 = R_alpha_short)
1 - 1/R_alpha_short
