library(EpiEstim)
library(epitrix)
library(scales)

incidence <- readRDS("outputs/data_outcomes2.rds")

# round and start on first day with data

incidence <- round(incidence)
incidence <- incidence[-seq(1, 3)]

####################################
# assume a GT or SI
####################################

mean_SI <- 6
sd_SI <- 2

distr_SI <- discr_si(seq(0, 21), mu = mean_SI, sigma = sd_SI)

####################################
### estimation with WL
####################################

y <- log(incidence+1)
x <- seq_along(y)

plot(x, y, type = "l")
dt <- 7
t1 <- seq(1, 1+dt)
t2 <- seq(2*dt+1, 6*dt)
t3 <- seq(7*dt+1, 8*dt)
t4 <- seq(9*dt+1, 11*dt)

lm1 <- lm(y[t1]~x[t1])
lm2 <- lm(y[t2]~x[t2])
lm3 <- lm(y[t3]~x[t3])
lm4 <- lm(y[t4]~x[t4])

r1 <- summary(lm1)$coefficients[2, ]
r2 <- summary(lm2)$coefficients[2, ]
r3 <- summary(lm3)$coefficients[2, ]
r4 <- summary(lm4)$coefficients[2, ]


pred.lm1 <- predict(lm1, interval = "confidence")
lines(x[t1], pred.lm1[,"fit"], col = "red")
lines(x[t1], pred.lm1[,"lwr"], col = "red", lty = 2)
lines(x[t1], pred.lm1[,"upr"], col = "red", lty = 2)

pred.lm2 <- predict(lm2, interval = "confidence")
lines(x[t2], pred.lm2[,"fit"], col = "red")
lines(x[t2], pred.lm2[,"lwr"], col = "red", lty = 2)
lines(x[t2], pred.lm2[,"upr"], col = "red", lty = 2)

pred.lm3 <- predict(lm3, interval = "confidence")
lines(x[t3], pred.lm3[,"fit"], col = "red")
lines(x[t3], pred.lm3[,"lwr"], col = "red", lty = 2)
lines(x[t3], pred.lm3[,"upr"], col = "red", lty = 2)

pred.lm4 <- predict(lm4, interval = "confidence")
lines(x[t4], pred.lm4[,"fit"], col = "red")
lines(x[t4], pred.lm4[,"lwr"], col = "red", lty = 2)
lines(x[t4], pred.lm4[,"upr"], col = "red", lty = 2)

R1 <- lm2R0_sample(lm1, w = distr_SI)
R2 <- lm2R0_sample(lm2, w = distr_SI)
R3 <- lm2R0_sample(lm3, w = distr_SI)
R4 <- lm2R0_sample(lm4, w = distr_SI)

quantile(R1, c(.5, .025, .975))
quantile(R2, c(.5, .025, .975))
quantile(R3, c(.5, .025, .975))
quantile(R4, c(.5, .025, .975))

R_WL <- data.frame(t = x,
                   R = rep(NA, length(x)),
                   R_low = rep(NA, length(x)),
                   R_up = rep(NA, length(x)))

R_WL$R[t1] <- mean(R1)
R_WL$R_low[t1] <- quantile(R1, .025)
R_WL$R_up[t1] <- quantile(R1, .975)

R_WL$R[t2] <- mean(R2)
R_WL$R_low[t2] <- quantile(R2, .025)
R_WL$R_up[t2] <- quantile(R2, .975)

R_WL$R[t3] <- mean(R3)
R_WL$R_low[t3] <- quantile(R3, .025)
R_WL$R_up[t3] <- quantile(R3, .975)

R_WL$R[t4] <- mean(R4)
R_WL$R_low[t4] <- quantile(R4, .025)
R_WL$R_up[t4] <- quantile(R4, .975)


R_WL_1 <- R_WL[t1,]
R_WL_2 <- R_WL[t2,]
R_WL_3 <- R_WL[t3,]
R_WL_4 <- R_WL[t4,]

####################################
### estimation with Cori et al.
####################################

config <- make_config(incid = incidence,
                      si_distr = distr_SI,
                      t_start = seq_along(x)[-1],
                      t_end = seq_along(x)[-1])

R_cori <- estimate_R(incid = incidence,
                method = "non_parametric_si",
                config = config)

####################################
### estimation with WT
####################################

config$n_sim <- 100

R_WT <- wallinga_teunis(incid = incidence,
                          method = "non_parametric_si",
                          config = config)

####################################
### summary
####################################

## TODO: generate a plot which compares the different estimates
#plot(R_cori, what = "R")
#plot(R_WT, what = "R")
# somehow plot R_WL

plot(R_cori$R$t_end, R_cori$R$`Mean(R)`, type = "l",
     xlab = "Time",
     ylab = "R estimates",
     col = "blue",
     ylim = c(0, 5))#ylim = c(0, 25))
polygon(c(R_cori$R$t_end, rev(R_cori$R$t_end)),
        c(R_cori$R$`Quantile.0.025(R)`, rev(R_cori$R$`Quantile.0.975(R)`)),
        col = alpha("blue", .2), border = NA)
abline(h = 1, col = "grey", lty = 2)
lines(R_WT$R$t_end, R_WT$R$`Mean(R)`, type = "l", col = "turquoise")
polygon(c(R_WT$R$t_end, rev(R_WT$R$t_end)),
        c(R_WT$R$`Quantile.0.025(R)`, rev(R_WT$R$`Quantile.0.975(R)`)),
        col = alpha("turquoise", .2), border = NA)
lines(R_WL$t, R_WL$R, type = "l", col = "violet")
polygon(c(R_WL_1$t, rev(R_WL_1$t)),
        c(R_WL_1$R_low, rev(R_WL_1$R_up)),
        col = alpha("violet", .2), border = NA)
polygon(c(R_WL_2$t, rev(R_WL_2$t)),
        c(R_WL_2$R_low, rev(R_WL_2$R_up)),
        col = alpha("violet", .2), border = NA)
polygon(c(R_WL_3$t, rev(R_WL_3$t)),
        c(R_WL_3$R_low, rev(R_WL_3$R_up)),
        col = alpha("violet", .2), border = NA)
polygon(c(R_WL_4$t, rev(R_WL_4$t)),
        c(R_WL_4$R_low, rev(R_WL_4$R_up)),
        col = alpha("violet", .2), border = NA)
legend("topright", c("Cori", "Wallinga and Teunis", "Wallinga and Lipsitch"),
       col = c("blue", "turquoise", "violet"), lty = 1)




