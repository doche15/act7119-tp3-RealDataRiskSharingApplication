lam <- 1:4

alp <- 5:8
bet <- lam

set.seed(2000)
n <- 1e6
reaX1 <- replicate(n, {sum(rgamma(rpois(1, lam[1]),
                                               alp[1], bet[1]))})

reaX2 <- replicate(n, {sum(rgamma(rpois(1, lam[2]),
                                     alp[2], bet[2]))})

reaX3 <- replicate(n, {sum(rgamma(rpois(1, lam[3]),
                                     alp[3], bet[3]))})

reaX4 <- replicate(n, {sum(rgamma(rpois(1, lam[4]),
                                     alp[4], bet[4]))})

realisations <- cbind(reaX1, reaX2, reaX3, reaX4)

reaecm2 <- (lam[2] * alp[2] / bet[2]) / sum(lam * alp / bet) * rowSums(realisations)

mean(reaecm2)
hist(reacm2)