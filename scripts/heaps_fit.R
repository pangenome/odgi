#!/usr/bin/env Rscript

require(tidyverse)

objectFun <- function(p, x, y){                                                                                                                      
  y.hat <- p[1] * x^(p[2])
  J <- sqrt(sum((y - y.hat)^2))/length(x)
  return(J)
}

# ENABLE command line arguments
args <- commandArgs(TRUE)
x <- read.delim(args[1])

fit <- with(x, optim(c(0,0), objectFun, gr = NULL, nth.genome, base.pairs/max(base.pairs), method = "L-BFGS-B", lower = c(-100, -100), upper = c(100, 100)))

print(fit)

ggplot(x, aes(x=nth.genome, y=base.pairs/max(base.pairs))) + geom_point() + stat_function(fun=function(x) fit$par[1] * x^fit$par[2])
ggsave(args[2])
