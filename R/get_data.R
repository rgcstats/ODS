set.seed(94114)
N <- 1000
n0 <- 200
xu <- rgamma(n=N,shape=3,scale=1)/3
#yu <- rnorm(n=100) + 5 + 0.5*xu
x1u <- rnorm(n=N)
x2u <- 1*(runif(N)<=0.4)
yu <- 1 + 0.5*x1u - 0.5*x2u + rnorm(N)
pi.u <- sqrt( 0.1 + 0.9*(yu-mean(yu))^2/var(yu) )
pi.u <- pi.u * n0 / sum(pi.u)
s.indicators <- 1*(runif(N)<=pi.u)
s <- c(1:N)[s.indicators==1]
pi.s <- pi.u[s]
population.example <- data.frame(x=xu,x1=x1u,x2=x2u,y=yu,pi=pi.u,s=s.indicators)
sample.example <- population.example[s,]
save(population.example,sample.example,file="data/population_example.rdata")

