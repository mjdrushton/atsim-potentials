

coul <- function(r, qi, qj){
	return ((qi * qj)/(4.0*pi*0.0055264*r))
}

buck <- function(r, A, rho, C){
	return (A*exp(-r/rho) - C/r^6)
}

morse <- function(r, D, gamma, r_star){
	return (D*(exp(-2*gamma*(r-r_star)))-2*exp(-gamma*(r-r_star)))
}

pot <- function(r, qi, qj, A, rho, C, D, gamma, r_star){
	return (coul(r, qi,qj) + buck(r, A, rho, C) + morse(r, D, gamma, r_star))
}


A <- 693.650934
rho <- 0.327022
C <- 0.0

D <- 0.577190
gamma <- 1.650000
r_star <- 2.369000

q_O <- -1.2
q_U <- 2.4

UO_pot <- function(r){
	return(pot(r, q_O, q_U, A, rho, C, D, gamma, r_star))
}



xlim <- c(0, 4)
ylim <- c(-25, 30)

pdf("analytical_plot.pdf", fonts = "Helvetica")
r <- seq(0.5, 4, length.out = 200)
r_coarse <- seq(0.98, 4, length.out = 10)

UO <- UO_pot(r)
UO_coarse <- UO_pot(r_coarse)

plot(r, UO, xlim = xlim, ylim = ylim, type = "l", xlab = NA, ylab = NA)
points(r_coarse, UO_coarse, col = 'red', pch=16)
abline(h=0, lty="dashed")
abline(v=r_coarse, col = "lightblue")
dev.off()

pdf("table_plot.pdf")
plot(r_coarse, UO_coarse, xlim = xlim, ylim = ylim, type = "l", xlab = NA, ylab = NA)
points(r_coarse, UO_coarse, col = 'red', pch=16)
abline(h=0, lty="dashed")
dev.off()