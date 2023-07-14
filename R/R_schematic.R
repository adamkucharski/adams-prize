#install.packages("plotrix")
#install.packages("igraph")
#install.packages("epitrix")
library(plotrix)
library(igraph)
library(epitrix)
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
arrow_end <- 0.15
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
ygt <- -1.2
segments(0, ygt, gt+radius, ygt, lty = 2)
segments(gt - radius, ygt - arrow_end, gt+radius, ygt)
segments(gt - radius, ygt + arrow_end, gt+radius, ygt)
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
ygt <- -1.2
segments(0, ygt, gt+radius, ygt, lty = 2)
segments(gt - radius, ygt - arrow_end, gt+radius, ygt)
segments(gt - radius, ygt + arrow_end, gt+radius, ygt)
text(gt/2, ygt-.3, "GT")
par(xpd = FALSE)

#### Reproduction number ####
iArrows(radius, 1+epsilon, gt + radius, R+epsilon,
        h.lwd=1, sh.lwd=1, sh.col="black", h.lty = 1, sh.lty = 2,
        curve=0.75, width=1.75, size=1.25)
text(gt/2, R+.5, "x R")


################################################################################
### R illustration, GT is not the same for all
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
arrow_end <- 0.15
draw_man(x0 = 0, y0 = 0, radius = radius, col = col, lwd = lwd)

## secondary cases
gt <- 5.45
R <- 4
x_jit <- c(-0., 1.1, -0.4, 0)
for(i in seq(0, R-1))
{
  draw_man(x0 = gt + x_jit[i+1], y0 = i, radius = radius, col = col, lwd = lwd)
}

#### generation time ####
par(xpd = TRUE)
ygt <- -1.2
segments(0, ygt, gt+radius, ygt, lty = 2)
segments(gt - radius, ygt - arrow_end, gt+radius, ygt)
segments(gt - radius, ygt + arrow_end, gt+radius, ygt)
text(gt/2, ygt-.3, "GT")
par(xpd = FALSE)

#### Reproduction number ####
iArrows(radius, 1+epsilon, gt + radius, R+epsilon,
        h.lwd=1, sh.lwd=1, sh.col="black", h.lty = 1, sh.lty = 2,
        curve=0.75, width=1.75, size=1.25)
text(gt/2, R+.5, "x R")

### for a similar r and mean GT, R decreases as the variance in GT increases.
epitrix::r2R0(r, c(0, 0, 0, 0, 0, 0, 1))
epitrix::r2R0(r, c(0, 0, 0, 0, 0, .33, .33, .33))
epitrix::r2R0(r, c(0, 0, 0, 0, 0.2, .2, .2, .2, .2))
### OR: in other words, for the same R and same mean GT, r increases as the variance in GT increases
epitrix::r2R0(r, c(0, 0, 0, 0, 0, 0, 1))
epitrix::r2R0(r*1.01427, c(0, 0, 0, 0, 0, .33, .33, .33))
epitrix::r2R0(r*1.045, c(0, 0, 0, 0, 0.2, .2, .2, .2, .2))

### if you add variance, you have more cases coming in earlier, which leads to r increasing
### this is more emphasised with distributions that have more weight on small values
### e.g. exp GT
epitrix::r2R0(r, c(0, 0, 0, 0, 0, 0, 1))
epitrix::r2R0(r, c(0, 0, 0, 0, 0.2, .2, .2, .2, .2))
epitrix::r2R0(r, distcrete::distcrete("exp", 1, rate = 1/6.5, w = 0))
#sum((0:100) * distcrete::distcrete("exp", 1, rate = 1/6.5, w = 0)$d(0:100))
epitrix::r2R0(r, distcrete::distcrete("gamma", 1, shape = 4, scale = 1.5, w = 0))
#sum((0:100) * distcrete::distcrete("gamma", 1, shape = 4.35, scale = 1.5, w = 0)$d(0:100))

### imagine you have a dirac GT.
### if you want to create the same shape of cases over 3 days, you need the distribution over these
### 3 days to be exponentially distr (like the cases) which means this corresponds to a larger mean GT
### to preserve the same mean GT, you'd then need to bring the GT to a lower value, or increaase the r.



