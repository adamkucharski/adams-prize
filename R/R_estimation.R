####################################
# Rt estimation and deconvolution on simulated data
# Collated to run in single script
####################################

####################################
# load libraries
####################################

library(EpiEstim)
library(epitrix)
library(scales)
library(MASS)
library(tidyverse)
library(EpiNow2)

# Generate simulated data -------------------------------------------------

# Simulate infection dynamics
growth_infections <- c(rep(0.1,40),rep(-1,5),rep(0.05,40))

initial_infections <- 20
data_infections <- initial_infections
for(ii in 2:length(growth_infections)){
  data_infections <- c(data_infections,data_infections[ii-1]*exp(growth_infections[ii]))
}

#data_infections <- c(seq(220,1000,20),rev(seq(200,1000,200)),rep(200,40)) # linear example

x_infections <- 1:length(data_infections) # time labels
n_inf <- length(data_infections) # number of days to consider

# Set delay function pmf
mean_p <- 10
scale_p <- 1
# alternative value to have no variability in incub:
# scale_p <- 0.001
shift_p <- 0

# Plot delay function
max_days <- 30
#p_by_day <- function(x){dgamma(x,shape=mean_p/scale_p,scale=scale_p)}
p_by_day <- function(x){pgamma(x+1,shape=mean_p/scale_p,scale=scale_p) - pgamma(x,shape=mean_p/scale_p,scale=scale_p)}
#plot(1:max_days,p_by_day(1:max_days))
tmp <- epitrix::gamma_shapescale2mucv(shape = mean_p/scale_p, scale = scale_p)
tmp$sd <- tmp$mu * tmp$cv

# Define transition matrix to construct outcome data
f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
n_delay_days <- 50 # maximum delay period to consider

for(ii in 1:n_inf){
  i_max <- min(ii+n_delay_days-1,n_inf)
  j_max <- min(n_inf-ii+1,n_delay_days)
  f_matrix[ii:i_max,ii] <- p_by_day(0:(j_max-1)) # fill matrix entries

}

# Run simulation function to generate expected delayed outcomes from original infections
data_outcomes <- f_matrix %*% data_infections

# Run again but with slight and heavy noise on delayed outcome observations
set.seed(5)
data_outcomes_slight <- round(data_outcomes * rlnorm(length(data_outcomes),0,0.0001) ) # add some noise
data_outcomes2 <- round(data_outcomes * rlnorm(length(data_outcomes),0,0.1) ) # add some noise


# Rt inference ------------------------------------------------------------

col_EpiEstim <- "blue"
col_WT <- "turquoise"
col_WL <- "violet"

# round and start on first day with data
incidence <- data_outcomes2
incidence_valid <- round(incidence)>0
incidence <- incidence[incidence_valid]
n_days_removed <- min(which(incidence_valid)) - 1

# generate matching indices
x_incidence <- x_infections[incidence_valid]

####################################
# assume a GT or SI
####################################

mean_SI <- 6
sd_SI <- 2
max_SI <- 25

# EpiEstim discretisation:
#distr_SI <- discr_si(seq(0, max_SI), mu = mean_SI, sigma = sd_SI)
#mean_SI_EpiEstim_discr <- sum(distr_SI*(seq(0, max_SI)))

# EpiNow2 / standard discretisation:
serial_interval_covid <-
  dist_spec(
    mean = mean_SI,
    sd = sd_SI,
    max = max_SI,
    distribution = "gamma"
  )
mean_SI_epinow <- sum(serial_interval_covid$np_pmf*(seq_len(max_SI)))
var_SI_epinow <- sum(serial_interval_covid$np_pmf*(seq_len(max_SI))^2)
distr_SI <- c(0, serial_interval_covid$np_pmf)

####################################
### estimation with WL on weekly time windows
####################################

y <- log(incidence+1)
x <- seq_along(y)

dt <- 7
times <- lapply(seq(1, length(y)-6, dt), function(i) seq(i, i + dt - 1))
lm_all <- lapply(seq_along(times), function(i) lm(y[times[[i]]]~x[times[[i]]]))

r_all <- t(sapply(seq_along(times), function(i) summary(lm_all[[i]])$coefficients[2, ]))
R_WL_all <- list()

R_WL <- data.frame(t = x,
                   R = rep(NA, length(x)),
                   R_low = rep(NA, length(x)),
                   R_up = rep(NA, length(x)))

for(i in seq_along(times)) {
  pred.lm <- predict(lm_all[[i]], interval = "confidence")
  #lines(x[times[[i]]], pred.lm[,"fit"], col = col_WL)
  #polygon(c(x[times[[i]]], rev(x[times[[i]]])),
  #        c(pred.lm[,"lwr"], rev(pred.lm[,"upr"])),
  #        col = alpha(col_WL, .2), border = NA)
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
### estimation with EpiNow2
####################################

# Make data.frame of cases for EpiNow2 - need to define dates
min_date <- as.Date("2020-03-01")
case_data <- data.frame(date = min_date+x_incidence, confirm = incidence)

### Now the SI distr used is harmonised between EpiNow2 and EpiEstim
plot(1:max_SI, serial_interval_covid$np_pmf)
lines(0:max_SI, distr_SI)

# Delay infection to outcome
# Uses: epiparameter::convert_params_to_summary_stats(distribution = "gamma", shape = mean_p/scale_p, scale = scale_p)

max_incub <- 25
incubation_time_covid <- dist_spec(
  mean = mean_p, #from simulation model
  sd = tmp$sd, #from parameterisation above
  max = max_incub,
  distribution = "gamma"
)
mean_discr_incub <- sum(incubation_time_covid$np_pmf*(seq_len(max_incub)))
var_discr_incub <- sum(incubation_time_covid$np_pmf*(seq_len(max_incub))^2)

## NOTE this is slow
R_epinow <- epinow(
  # cases
  reported_cases = case_data,
  # delays
  generation_time = generation_time_opts(serial_interval_covid),
  # delays
  delays = delay_opts(incubation_time_covid),
  # computation
  stan = stan_opts(
    cores = 4, samples = 1000, chains = 3,
    control = list(adapt_delta = 0.99)
  )
)

# Check output: base::plot(epinow_estimates)


####################################
### summary plots for Rt
####################################

pdf(paste("outputs/R_plot.pdf",sep=""),width=6,height=5)

cex_legend <- .75
par(mfrow=c(2,2),mgp=c(2.5,0.7,0),mar = c(3.5,3.5,1,1))
letter_x <- 1
# Plot incidence
# plot(x, y,
#      xlim = c(0, 85), ylim = c(0, 7),
#      xlab = "Days",
#      ylab = "Log daily incidence")

plot(x_infections,data_infections,ylim=c(0,1.5e3),xlab="time",ylab="events",type="l",lwd=2,
     xlim = c(0, 85))
lines(x_infections,data_outcomes2,col="red",lwd=1,lty=1)
legend("topright", c("infections", "outcomes"),
       col = c("black", "red"), lty = 1, cex = cex_legend)
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1

# Plot on log scale
plot(x_infections,data_infections,ylim=c(1,3e4),xlab="time",ylab="events (log scale)",type="l",lwd=2,log="y",
     xlim = c(0, 85))
lines(x_infections,data_outcomes2,col="red",lwd=1,lty=1)
legend("topright", c("infections", "outcomes"),
       col = c("black", "red"), lty = 1, cex = cex_legend)
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1


# Plot Rt estimates by outcome data
plot(R_cori$R$t_end + n_days_removed, R_cori$R$`Mean(R)`, type = "l",
     xlim = c(0, 85), ylim = c(0, 5),
     xlab = "Days",
     ylab = "R estimates",
     col = col_EpiEstim)
polygon(c(R_cori$R$t_end + n_days_removed, rev(R_cori$R$t_end)),
        c(R_cori$R$`Quantile.0.025(R)`, rev(R_cori$R$`Quantile.0.975(R)`)),
        col = alpha(col_EpiEstim, .2), border = NA)
abline(h = 1, col = "grey", lty = 2)
lines(R_WT$R$t_end + n_days_removed, R_WT$R$`Mean(R)`, type = "l", col = col_WT)
polygon(c(R_WT$R$t_end + n_days_removed, rev(R_WT$R$t_end)),
        c(R_WT$R$`Quantile.0.025(R)`, rev(R_WT$R$`Quantile.0.975(R)`)),
        col = alpha(col_WT, .2), border = NA)
lines(R_WL$t + n_days_removed, R_WL$R, type = "l", col = col_WL)
polygon(c(R_WL$t + n_days_removed, rev(R_WL$t)),
        c(R_WL$R_low, rev(R_WL$R_up)),
        col = alpha(col_WL, .2), border = NA)
legend("topright", c("EpiEstim", "Wallinga and Teunis", "Wallinga and Lipsitch"),
       col = c(col_EpiEstim, col_WT, col_WL), lty = 1, cex = cex_legend)
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1

# Plot Rt estimates by infection date via EpiNow2 and simple shift for EpiNow
R_estimates <- R_epinow$estimates$summarised |> filter(variable=="R")
x_numeric <- as.numeric(R_estimates$date-min_date) # convert to numerical values

shift_epiestim <- mean_discr_incub # mean_p # shift by half the delay period? ANNE TO CHECK
## should be ok although technically the mean of the discrete distribution is slightly different
## not sure if it's sum(incubation_time_covid$np_pmf*(0:24)) or sum(incubation_time_covid$np_pmf*(1:25))
## howver this is not, I suspect, why the curves are shifted.
## I think this is because I used daily estimates for EpiEstim and EpiNow2 must be doing some smoothing
## seems to be in the gp = gp_opts() arguments that you change it in EpiNow2 but not sure.


plot(x_numeric,R_estimates$median,yaxs="i",ylab="R estimates",ylim=c(0,5),xlab="days",type="l",col="darkorange",
     xlim = c(0, 85))
polygon(c(x_numeric,rev(x_numeric)),c(R_estimates$lower_90,rev(R_estimates$upper_90)),
        col=rgb(1,0.5,0,0.2),border=NA)
abline(h = 1, col = "grey", lty = 2)
lines(R_cori$R$t_end + n_days_removed - shift_epiestim, R_cori$R$`Mean(R)`, col = col_EpiEstim) # Plot shifted EpiNow
polygon(c(R_cori$R$t_end+n_days_removed-shift_epiestim, rev(R_cori$R$t_end+n_days_removed-shift_epiestim)),
        c(R_cori$R$`Quantile.0.025(R)`, rev(R_cori$R$`Quantile.0.975(R)`)),
        col = alpha(col_EpiEstim, .2), border = NA)
legend("topright", c("EpiEstim shifted", "EpiNow2"),
       col = c(col_EpiEstim, "darkorange"), lty = 1, cex = cex_legend)
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1

dev.off()

# Analysis and plots for deconvolution ------------------------------------

# Define inversion function
invert_f <- ginv(f_matrix)

# Function to estimate infection incidence from delayed outcomes
estimate_infections_mat <- function(delayed_outcomes){

  # Define transition matrix and calculate entries
  n_inf <- length(delayed_outcomes)
  f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
  n_delay_days <- 30 # maximum incubation period to consider

  for(ii in 1:n_inf){
    i_max <- min(ii+n_delay_days-1,n_inf)
    j_max <- min(n_inf-ii+1,n_delay_days)
    f_matrix[ii:i_max,ii] <- p_by_day(0:(j_max-1))
  }

  # Calculate Moore-Penrose pseudoinverse matrix
  invert_f <- ginv(f_matrix)

  # Apply inversion matrix
  output <- invert_f %*% delayed_outcomes
  output

}

# Run simulation function to generate delayed outcomes from original infections
data_outcomes <- f_matrix %*% data_infections

pdf(paste("outputs/convolution_plot.pdf",sep=""),width=4,height=6)
cex_legend <- .75
par(mfrow=c(3,1),mgp=c(2.5,0.7,0),mar = c(3.5,3.5,1,1))
letter_x <- 1
ymax <- 1.3e3

# Plot original incidence vs inference based on deconvolution and shift
plot(x_infections,data_infections,ylim=c(0,ymax),xlab="time",ylab="daily events")
lines(x_infections,data_outcomes,col="red",lwd=1,lty=1)
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1

# Plot simple deconvolution
inc1 <- estimate_infections_mat(data_outcomes)
lines(inc1,col="blue",lwd=2)

# Plot simple shift
lines(x_infections-mean_p,data_outcomes,col="orange",lwd=2)

# Add caption
legend("topright", c("true infections","delayed outcomes","deconvolved outcomes","shifted outcomes"),
       col = c("black", "red","blue","orange"), lty = 1, cex = cex_legend)

# Run again but with noise on delayed outcome observations
#data_outcomes2 <- data_outcomes * rlnorm(length(data_outcomes),0,0.0001) # add some noise

# Plot noisy deconvolution
plot(data_infections,yaxs="i",ylab="daily events",ylim=c(0,ymax),xlab="days")

inc2 <- estimate_infections_mat(data_outcomes_slight)
lines(inc2,col="blue",lwd=1)
lines(data_outcomes_slight,col="red",lwd=2)

legend("topright", c("true infections","\'noisy\' outcomes","deconvolved outcomes"),
       col = c("black","red","blue"), lty = 1, cex = cex_legend, bg = "white")
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1

# EpiNow2 inference

# Extract estimates
infection_estimates <- R_epinow$estimates$summarised |> filter(variable=="infections")
x_numeric <- as.numeric(infection_estimates$date-min_date) # convert to numerical values

x_obs <- as.numeric(R_epinow$estimates$observations$date-min_date)
y_obs <- as.numeric(R_epinow$estimates$observations$confirm)

plot(x_infections,data_infections,yaxs="i",ylab="daily events",ylim=c(0,ymax),xlab="days")
lines(data_outcomes2,col="red",lwd=2)
polygon(c(x_numeric,rev(x_numeric)),c(infection_estimates$lower_90,rev(infection_estimates$upper_90)),
        col=rgb(1,0.5,0,0.2),border=NA)
lines(x_numeric,infection_estimates$median,lwd=2,col="darkorange")

legend("topright", c("true infections","noisy outcomes","EpiNow2 estimated infections"),
       col = c("black","red","darkorange"), lty = 1, cex = cex_legend)
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1


# dev.copy(png,paste0("outputs/convolution_plot.png"), units="cm",width=15,height=15,res=150)
# dev.off()
# dev.copy(pdf,paste("outputs/convolution_plot.pdf",sep=""),width=4,height=6)
dev.off()


