
library(EpiEstim)
library(epitrix)
library(scales)

col_EpiEstim <- "blue"
col_WT <- "turquoise"
col_WL <- "violet"

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
### estimation with WL on weekly time windows
####################################

par(mfrow = c(2, 1))

y <- log(incidence+1)
x <- seq_along(y)

plot(x, y,
     xlim = c(0, 85), ylim = c(0, 7),
     xlab = "Days",
     ylab = "Log daily incidence")
dt <- 7
times <- lapply(seq(1, length(y)-5, dt), function(i) seq(i, i + dt - 1))
lm_all <- lapply(seq_along(times), function(i) lm(y[times[[i]]]~x[times[[i]]]))

r_all <- t(sapply(seq_along(times), function(i) summary(lm_all[[i]])$coefficients[2, ]))
R_WL_all <- list()

R_WL <- data.frame(t = x,
                   R = rep(NA, length(x)),
                   R_low = rep(NA, length(x)),
                   R_up = rep(NA, length(x)))

for(i in seq_along(times)) {
  pred.lm <- predict(lm_all[[i]], interval = "confidence")
  lines(x[times[[i]]], pred.lm[,"fit"], col = col_WL)
  polygon(c(x[times[[i]]], rev(x[times[[i]]])),
          c(pred.lm[,"lwr"], rev(pred.lm[,"upr"])),
          col = alpha(col_WL, .2), border = NA)
  R_WL_all[[i]] <- lm2R0_sample(lm_all[[i]], w = distr_SI)
  R_WL$R[times[[i]]] <- mean(R_WL_all[[i]])
  R_WL$R_low[times[[i]]] <- quantile(R_WL_all[[i]], 0.025)
  R_WL$R_up[times[[i]]] <- quantile(R_WL_all[[i]], 0.975)
}

R_WL <- R_WL[-which(is.na(R_WL$R)), ]

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

plot(R_cori$R$t_end, R_cori$R$`Mean(R)`, type = "l",
     xlim = c(0, 85), ylim = c(0, 5),
     xlab = "Days",
     ylab = "R estimates",
     col = col_EpiEstim)
polygon(c(R_cori$R$t_end, rev(R_cori$R$t_end)),
        c(R_cori$R$`Quantile.0.025(R)`, rev(R_cori$R$`Quantile.0.975(R)`)),
        col = alpha(col_EpiEstim, .2), border = NA)
abline(h = 1, col = "grey", lty = 2)
lines(R_WT$R$t_end, R_WT$R$`Mean(R)`, type = "l", col = col_WT)
polygon(c(R_WT$R$t_end, rev(R_WT$R$t_end)),
        c(R_WT$R$`Quantile.0.025(R)`, rev(R_WT$R$`Quantile.0.975(R)`)),
        col = alpha(col_WT, .2), border = NA)
lines(R_WL$t, R_WL$R, type = "l", col = col_WL)
polygon(c(R_WL$t, rev(R_WL$t)),
        c(R_WL$R_low, rev(R_WL$R_up)),
        col = alpha(col_WL, .2), border = NA)
legend("topright", c("EpiEstim", "Wallinga and Teunis", "Wallinga and Lipsitch"),
       col = c(col_EpiEstim, col_WT, col_WL), lty = 1)

