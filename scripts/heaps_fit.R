#!/usr/bin/env Rscript

require(tidyverse)

objectFun <- function(p, x, y){
    y.hat <- p[1] * x^(p[2]) + p[3]
    #J <- sqrt(sum((y - y.hat)^2))/length(x)
    J <- sum((y - y.hat)^2) # sum of squared error
    return(J)
}

# ENABLE command line arguments
args <- commandArgs(TRUE)
x <- read.delim(args[1])

fit <- with(x, optim(c(0,0,0), objectFun, gr = NULL, nth.genome, base.pairs/max(base.pairs), method = "L-BFGS-B", lower = c(-100, -100), upper = c(100, 100)))

print(fit)

print(min(x$base.pairs))
print(max(x$base.pairs))
print(fit$par[1]*max(x$base.pairs))
print(fit$par[2]*max(x$base.pairs))
print(fit$par[3]*max(x$base.pairs))

z <- max(x$base.pairs)
m <- z/1e9

f <- function(x) { fit$par[1] * x^fit$par[2] + fit$par[3] }
n <- max(x$nth.genome)
print(z * (f(n) - f(n-1)))
print(z * (f(2) - f(1)))
#print(f(n) - f(n-1))
pdf(NULL)
ggplot(x, aes(x=nth.genome, y=base.pairs/1e9)) + geom_point(alpha=I(1/10)) + stat_function(fun=function(x) (fit$par[1] * x^fit$par[2] + fit$par[3]) * m) + scale_y_continuous("observed pangenome size (Gbp)") + scale_x_continuous(paste("Nth included genome (", max(x$permutation)+1 ," permutations) with gamma=", round(fit$par[2], digits=3), sep = "")) + expand_limits(x = 0, y = 0)
ggsave(args[2], height=5, width=9)
