R0 <- seq(0,15,0.1)
hit <- (R0 - 1)/R0

#ve <- seq(0, 0.9, 0.1)
ve <- seq(0.2, .8, 0.2)
hit_ve <- sapply(ve, function(e) (R0 - 1)/(R0*e))

colnames(hit_ve) <- ve
rownames(hit_ve) <- R0

palette <- viridisLite::viridis(11)
#palette <- viridisLite::magma(11)
#palette <- viridisLite::turbo(11)

plot(R0, hit, type = "l",
     xlab = expression(R[0]),
     ylab = "Herd immunity threshold",
     ylim = c(0, 1), xaxs="i", yaxs="i",
     col = palette[1],
     lwd = 2)

for(i in seq(ncol(hit_ve))) {
  lines(R0, hit_ve[,i],
        col = palette[11-i], lwd = 2)
}

legend("bottomright",
       paste("VE =",c(ve, 1)),
       col = c(palette[11-seq(ncol(hit_ve))], palette[1]), lty = 1, lwd = 2)
