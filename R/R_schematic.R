#install.packages("plotrix")             # Install plotrix package
#install.packages("igraph")
library("plotrix")                      # Load plotrix package
library(igraph)
iArrows <- igraph:::igraph.Arrows

source("R/functions.R")

par(las = 0, xaxs="i", yaxs="i", mar = c(7, 7, 2, 2))

epsilon <- 0.1
r <- 0.2
tt <- seq(0, 8, 0.01)
xmax <- max(tt)
ymax <- max(exp(r*tt))


################################################################################
### R illustration, panel 1
################################################################################

## exponential growth curve
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax))
mtext("Time", 1, 5)
mtext("Incidence", 2, 3)
text(0.9*xmax, 0.95*ymax, expression(I[t]~"="~I[0]~e^rt))

#### little men ####

## primary case
col <- "grey"
lwd <- 2
draw_man(x0 = 0, y0 = 0, radius = 0.2, col = col, lwd = lwd)

## secondary cases
gt <- 5.4
draw_man(x0 = gt, y0 = 0, radius = 0.2, col = col, lwd = lwd)
draw_man(x0 = gt, y0 = 1, radius = 0.2, col = col, lwd = lwd)
draw_man(x0 = gt, y0 = 2, radius = 0.2, col = col, lwd = lwd)

#### generation time ####
par(xpd = TRUE)
segments(0, -1, gt, -1, lty = 2)
segments(gt - radius, -1 - 0.5*radius, gt, -1)
segments(gt - radius, -1 + 0.5*radius, gt, -1)
text(gt/2, -1.25, "GT")
par(xpd = FALSE)

#### Reproduction number ####
iArrows(radius, 1+epsilon, gt + radius, 3+epsilon,
        h.lwd=1, sh.lwd=1, sh.col="black", h.lty = 1, sh.lty = 2,
        curve=0.75, width=1.75, size=1.25)
text(gt/2, 4, "x R")
