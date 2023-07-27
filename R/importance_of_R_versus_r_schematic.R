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

pdf("figs/R_versus_r.pdf", width = 10, height = 5)

################################################################################
### Create panel layout
################################################################################

par(mfrow = c(1, 2),
    las = 0, xaxs="i", yaxs="i", mar = c(5, 5, 2, 2))

################################################################################
### R illustration, baseline panel
################################################################################

## exponential growth curve obtained with flu-like natural history
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax),
     main = paste0("R = ", signif(R_flu, 2), "; Flu-like natural history"))
mtext("Time", 1, 3)
mtext("Incidence", 2, 3)
abline(v = mean_GT_flu*seq_len(30), lty = 2, col = "grey")

## exponential growth curve obtained with measles-like natural history
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax),
     main = paste0("R = ", signif(R_measles, 2), "; Measles-like natural history"))
mtext("Time", 1, 3)
mtext("Incidence", 2, 3)
abline(v = mean_GT_measles*seq_len(30), lty = 2, col = "grey")

dev.off()
