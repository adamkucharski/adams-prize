#install.packages("plotrix")             # Install plotrix package
#install.packages("igraph")
library("plotrix")                      # Load plotrix package
library(igraph)
iArrows <- igraph:::igraph.Arrows

source("R/functions.R")

par(las = 0, xaxs="i", yaxs="i", mar = c(7, 7, 2, 2))

epsilon <- 0.1
r <- 0.25
tt <- seq(0, 8, 0.01)
xmax <- max(tt)
ymax <- max(exp(r*tt))


################################################################################
### Create panel layout
################################################################################

par(mfrow = c(2, 2))

################################################################################
### R illustration, baseline panel
################################################################################

## exponential growth curve
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax))
mtext("Time", 1, 5)
mtext("Incidence", 2, 3)
text(0.9*xmax, 0.95*ymax, expression(I[t]~"="~I[0]~e^rt))

#### little men ####

## primary case
col <- "black"
lwd <- 3
radius <- 0.15
draw_man(x0 = 0, y0 = 0, radius = radius, col = col, lwd = lwd)

## secondary cases
gt <- 5.45
R <- 4
for(i in seq(0, R-1))
{
  draw_man(x0 = gt, y0 = i, radius = radius, col = col, lwd = lwd)
}
#### generation time ####
par(xpd = TRUE)
ygt <- -.8
segments(0, ygt, gt, ygt, lty = 2)
segments(gt - radius, ygt - 0.5*radius, gt, ygt)
segments(gt - radius, ygt + 0.5*radius, gt, ygt)
text(gt/2, ygt-.3, "GT")
par(xpd = FALSE)

#### Reproduction number ####
iArrows(radius, 1+epsilon, gt + radius, R+epsilon,
        h.lwd=1, sh.lwd=1, sh.col="black", h.lty = 1, sh.lty = 2,
        curve=0.75, width=1.75, size=1.25)
text(gt/2, R+.5, "x R")

################################################################################
### R illustration, under-reporting
################################################################################

p_report <- 0.5

## exponential growth curve
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax), col = "grey",axes = FALSE)
axis(side = 1)
yaxis <- seq(0, 16, 2)
axis(side = 2, at = yaxis, labels = yaxis*2)
lines(tt, exp(r*tt)*p_report)
mtext("Time", 1, 5)
mtext("Incidence", 2, 3)
text(0.9*xmax, 0.95*ymax, expression(I[0]~e^rt))
text(0.9*xmax, 0.95*ymax*p_report, expression(rho~I[0]~e^rt))

#### little men ####

## primary case
col <- "black"
col2 <- "grey"
lwd <- 2
radius <- 0.1
draw_man(x0 = 0, y0 = 0, radius = radius, height = .5, col = col, lwd = lwd)
draw_man(x0 = 0, y0 = .5, radius = radius, height = .5, col = col2, lwd = lwd)

## secondary cases
gt <- 5.5
R <- 4
for(i in seq(0, R-1))
{
  draw_man(x0 = gt, y0 = i/2, radius = radius, height = .5, col = col, lwd = lwd)
}
for(i in seq(R, 2*R-1))
{
  draw_man(x0 = gt, y0 = i/2, radius = radius, height = .5, col = col2, lwd = lwd)
}
#### generation time ####
par(xpd = TRUE)
ygt <- -.8
segments(0, ygt, gt, ygt, lty = 2)
segments(gt - radius, ygt - 0.5*radius, gt, ygt)
segments(gt - radius, ygt + 0.5*radius, gt, ygt)
text(gt/2, ygt-.3, "GT")
par(xpd = FALSE)

#### Reproduction number ####
iArrows(radius, 1+epsilon, gt + radius, R+epsilon,
        h.lwd=1, sh.lwd=1, sh.col="black", h.lty = 1, sh.lty = 2,
        curve=0.75, width=1.75, size=1.25)
text(gt/2, R+.5, "x R")

