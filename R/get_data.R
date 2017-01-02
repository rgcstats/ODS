set.seed(94114)
xu <- rgamma(n=100,shape=3,scale=1)/3
yu <- rnorm(n=100) + 5 + 0.5*xu
pi.u <- sqrt( 0.1 + 0.9*(yu-mean(yu))^2/var(yu) )
pi.u <- pi.u * 10 / sum(pi.u)
s.indicators <- 1*(runif(100)<=pi.u)
s <- c(1:100)[s.indicators]
ys <- yu[s] ; xs <- xu[s]
yr <- yu[-s] ; xr <- xu[-s]
pi.s <- pi.u[s]
population.example <- data.frame(x=xu,y=yu,pi=pi.u,s=s.indicators)
save(population.example,file="data/population_example.Rdata")
