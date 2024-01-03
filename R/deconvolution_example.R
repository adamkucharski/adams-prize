# - - - - - - - - - - - - - - - - - - - - - - -
# Deconvolution of simulated infection data
# - - - - - - - - - - - - - - - - - - - - - - -

# Load libraries
library(MASS)
library(tidyverse)

# setwd("~/Documents/GitHub/adams-prize")

# Simulate infection dynamics
xx <- 0:1000
data_infections <- c(seq(220,1000,20),rev(seq(200,1000,200)),rep(200,40))

n_inf <- length(data_infections) # number of days to consider
data_infections <- data_infections #* rlnorm(n_inf,0,0.2) # add some noise

# Set delay function pmf
#p_by_day <- #epiparameter::epidist("SARS_CoV_2_wildtype","incubation")$pmf
mean_p <- 10
scale_p <- 1
shift_p <- 0

# Plot delay function
max_days <- 30
p_by_day <- function(x){dgamma(x,shape=mean_p/scale_p,scale=scale_p)}
plot(1:max_days,p_by_day(1:max_days))

# Define transition matrix to construct outcome data
f_matrix <- matrix(0,nrow=n_inf,ncol=n_inf)
n_delay_days <- 50 # maximum delay period to consider

for(ii in 1:n_inf){
  i_max <- min(ii+n_delay_days-1,n_inf)
  j_max <- min(n_inf-ii+1,n_delay_days)

  f_matrix[ii:i_max,ii] <- p_by_day(0:(j_max-1)) # fill matrix entries

}


# Quick simulation --------------------------------------------------------

# Simulate outcomes
letter_x <- 1

data_outcomes <- f_matrix %*% data_infections

par(mfrow=c(1,1),mgp=c(2,0.7,0),mar = c(3,3,1,1))
ymax <- 1.1e3

# Plot original incidence
plot(data_infections,yaxs="i",ylab="daily incidence (%)",ylim=c(0,ymax),xlab="days")
points(data_outcomes,col="red")
title(main=LETTERS[letter_x],adj=0);letter_x <- letter_x+1

# Run deconvolution methods -----------------------------------------------


# Inversion function
invert_f <- ginv(f_matrix)

# Function to estimate infection incidence from delayed outcomes
estimate_infections <- function(delayed_outcomes){

  # Define transition matrix -
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
#data_outcomes <- data_outcomes * rlnorm(length(data_outcomes),0,0.05) # add some noise

par(mfrow=c(1,2),mgp=c(2,0.7,0),mar = c(3,3,1,1))
letter_x <- 1

# Plot original incidence vs inference based on deconvolution and shift
plot(data_infections,yaxs="i",ylab="daily cases",ylim=c(0,ymax),xlab="days")
lines(data_outcomes,col="red",lwd=2)
title(main=LETTERS[1],adj=0);letter_x <- letter_x+1

# Plot deconvolution
inc1 <- estimate_infections(data_outcomes)
lines(inc1,col="blue",lwd=2)

# Plot simple shift
lines(tail(data_outcomes,-mean_p),col="orange",lwd=2)

# Add caption
x_shift <- 50; y_shift <- 1000; t_size <- 0.7
graphics::text(labels="true infections",x=x_shift,y=y_shift,adj=0,col="black",cex=t_size)
graphics::text(labels="delayed outcomes",x=x_shift,y=(y_shift-50),adj=0,col="red",cex=t_size)
graphics::text(labels="deconvolved outcomes",x=x_shift,y=(y_shift-100),adj=0,col="blue",cex=t_size)
graphics::text(labels="shifted outcomes",x=x_shift,y=(y_shift-150),adj=0,col="orange",cex=t_size)


# Run again but with noise on delayed outcome observations
data_outcomes2 <- data_outcomes * rlnorm(length(data_outcomes),0,0.0001) # add some noise

# Plot original incidence
plot(data_infections,yaxs="i",ylab="daily cases",ylim=c(0,ymax),xlab="days")

inc2 <- estimate_infections(data_outcomes2)
lines(inc2,col="blue",lwd=1)
lines(data_outcomes2,col="red",lwd=2)

graphics::text(labels="noisy outcomes",x=x_shift,y=(y_shift-50),adj=0,col="red",cex=t_size)
graphics::text(labels="deconvolved outcomes",x=x_shift,y=(y_shift-100),adj=0,col="blue",cex=t_size)


title(main=LETTERS[2],adj=0);letter_x <- letter_x+1

dev.copy(png,paste0("outputs/convolution_plot.png"),units="cm",width=20,height=10,res=150)
dev.off()


saveRDS(data_outcomes, "outputs/data_outcomes.rds")
saveRDS(data_outcomes2, "outputs/data_outcomes2.rds")
