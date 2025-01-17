#install.packages("plotrix")             # Install plotrix package
#install.packages("igraph")
#install.packages("scales")
library("plotrix")                      # Load plotrix package
library(igraph)
library(scales)
iArrows <- igraph:::igraph.Arrows

source("R/functions.R")

epsilon <- 0.1
r <- 0.25
tt <- seq(0, 8, 0.01)
xmax <- max(tt)
ymax <- max(exp(r*tt))


pdf("figs/R_estim.pdf", width = 10, height = 10)

################################################################################
### Create panel layout
################################################################################

par(mfrow = c(2, 2),
    las = 0, xaxs="i", yaxs="i", mar = c(7, 7, 2, 2))

################################################################################
### R illustration, baseline panel
################################################################################

## exponential growth curve
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax))
mtext("Time", 1, 5)
mtext("Incidence", 2, 3)
text(0.9*xmax, 0.95*ymax, expression(I[t]~"="~I[0]~e^rt))

# GT line
par(xpd = TRUE)
ygt <- -1
segments(gt + radius, ygt, gt + radius, 1.1*ymax, col = "lightgrey")
par(xpd = FALSE)

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

## panel name
par(xpd = TRUE)
text(-1, 8.5, "A", cex = 2)
par(xpd = FALSE)

################################################################################
### R illustration, GT is not the same for all
################################################################################

## exponential growth curve
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n", col = "grey",
     xlim = range(tt), ylim = c(0, 1.1*ymax))

#### generation time option 2####
par(xpd = TRUE)
ygt <- 0
xpolygon <- c(gt +range(c(x_jit, x_jit + 2*radius)))
polygon(x = c(xpolygon, rev(xpolygon)),
        y = c(ygt, ygt, 1.1*ymax, 1.1*ymax), col = scales::alpha("lightgrey", 0.7),
        border = NA)
par(xpd = FALSE)

## redraw exponential growth curve on top
lines(tt, exp(r*tt), col = "grey")
mtext("Time", 1, 5)
mtext("Incidence", 2, 3)

lines(tt, exp(r*1.14*tt))
#text(0.8*xmax, 0.98*ymax, expression(I[t]~"="~I[0]~e^rt))

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
x_jit <- c(.6, .15, -0.15, -.6)
for(i in seq(0, R-1))
{
  draw_man(x0 = gt + x_jit[i+1], y0 = i, radius = radius, col = col, lwd = lwd)
}

#### generation time option 1 ####
#par(xpd = TRUE)
#ygt <- -1.2
# option 1
# polygon(x = c(gt +range(x_jit), rev(gt +range(x_jit))),
#         y = c(ygt, ygt, 0.65*ygt, 0.65*ygt), col = scales::alpha("grey", 0.5),
#         border = NA)
# text(gt, ygt-.3, "GT")
# par(xpd = FALSE)

#### Reproduction number ####
# iArrows(radius, 1+epsilon, gt +min(x_jit) + radius, R+epsilon,
#         h.lwd=1, sh.lwd=1, sh.col="black", h.lty = 1, sh.lty = 2,
#         curve=0.75, width=1.75, size=1.25)
# text(gt/2, R+.5, "x R")


## panel name
par(xpd = TRUE)
text(-1, 8.5, "B", cex = 2)
par(xpd = FALSE)

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
#text(0.9*xmax, 0.95*ymax, expression(I[0]~e^rt))
#text(0.9*xmax, 0.95*ymax*p_report, expression(rho~I[0]~e^rt))

# # GT line
# par(xpd = TRUE)
# ygt <- 0
# segments(gt + radius, ygt, gt + radius, 1.1*ymax, col = "lightgrey")
# par(xpd = FALSE)

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
# #### generation time ####
# par(xpd = TRUE)
# ygt <- -1.2
# segments(0, ygt, gt+radius, ygt, lty = 2)
# segments(gt - radius, ygt - arrow_end, gt+radius, ygt)
# segments(gt - radius, ygt + arrow_end, gt+radius, ygt)
# text(gt/2, ygt-.3, "GT")
# par(xpd = FALSE)
#
# #### Reproduction number ####
iArrows(radius, 1+epsilon, gt + radius, R+epsilon,
         h.lwd=1, sh.lwd=1, sh.col= col2, h.lty = 1, sh.lty = 2,
         curve=0.75, width=1.75, size=1.25)
text(gt/2, R+.5, "x R", col = col2)

iArrows(radius, .5+epsilon, gt + radius, R/2+epsilon,
        h.lwd=1, sh.lwd=1, sh.col= col, h.lty = 1, sh.lty = 2,
        curve=0.2, width=1.75, size=1.25)
text(gt*5/6, R/2+.25, "x R", col = col)


## panel name
par(xpd = TRUE)
text(-1, 8.5, "C", cex = 2)
par(xpd = FALSE)

################################################################################
### R illustration, delays in reporting
################################################################################

p_not_report_same_day <- 0.25


## exponential growth curve
plot(tt, exp(r*tt), type = "l", xlab = "", ylab = "", bty = "n",
     xlim = range(tt), ylim = c(0, 1.1*ymax), col = "grey")
mtext("Time", 1, 5)
mtext("Incidence", 2, 3)
#text(0.9*xmax, 0.95*ymax, expression(I[t]~"="~I[0]~e^rt))

lines(tt, exp(r*0.8*tt))

# # GT line
# par(xpd = TRUE)
# ygt <- 0
# segments(gt + radius, ygt, gt + radius, 1.1*ymax, col = "lightgrey")
# par(xpd = FALSE)

#### little men ####

## primary case
col <- "black"
col2 <- "grey"
lwd <- 3
radius <- 0.15
arrow_end <- 0.15
draw_man(x0 = 0, y0 = 0, radius = radius, col = col, lwd = lwd)

## secondary cases
gt <- 5.45
R <- 4
for(i in seq(0, R-2))
{
  draw_man(x0 = gt, y0 = i, radius = radius, col = col, lwd = lwd)
}
i <- R-1
draw_man(x0 = gt, y0 = i, radius = radius, col = col2, lwd = lwd)
#
# #### generation time ####
# par(xpd = TRUE)
# ygt <- -1.2
# segments(0, ygt, gt+radius, ygt, lty = 2)
# segments(gt - radius, ygt - arrow_end, gt+radius, ygt)
# segments(gt - radius, ygt + arrow_end, gt+radius, ygt)
# text(gt/2, ygt-.3, "GT")
# par(xpd = FALSE)
#
# #### Reproduction number ####
iArrows(radius, 1+epsilon, gt + radius, R+epsilon,
        h.lwd=1, sh.lwd=1, sh.col= col2, h.lty = 1, sh.lty = 2,
        curve=0.75, width=1.75, size=1.25)
text(gt/2, R+.5, "x R", col = col2)


## panel name
par(xpd = TRUE)
text(-1, 8.5, "D", cex = 2)
par(xpd = FALSE)

dev.off()


