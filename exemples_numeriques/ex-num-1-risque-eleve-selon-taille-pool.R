###
### Travail 3 du cours ACT-7119
###
### Analyser l'impact de la taille du pool sur un risque artificiel élevé.
###
##

params_data <- read.csv("params_belgian.csv")
summary(params_data)
primes_pures <- params_data$lambda * (params_data$alpha / params_data$beta)

# Selon les valeurs maximales des paramètres, on crée un risque élevé.
params_risque_eleve <- c(1, 3, 0.00025)
prime_pure_risque_eleve <- params_risque_eleve[1] *
    (params_risque_eleve[2] / params_risque_eleve[3])

c(prime_pure_risque_eleve = prime_pure_risque_eleve, prime_pure_moyenne = mean(primes_pures))


# Calculs de partage de risque
source("algorithms/cond_mean_rsharing_FFT_Blier.R")
source("algorithms/mean_prop_rsharing.R")
source("algorithms/linear_reg_rsharing.R")

total_pool_size <- nrow(params_data)

set.seed(11)
pool10 <- rbind(params_risque_eleve,
                params_data[sample(seq(total_pool_size), 9),])
sum(pool10$lambda * (pool10$alpha / pool10$beta))

pool100 <- rbind(params_risque_eleve,
                 params_data[sample(seq(total_pool_size), 99),])
sum(pool100$lambda * (pool100$alpha / pool100$beta))

pool1000 <- rbind(params_risque_eleve,
                  params_data[sample(seq(total_pool_size), 999),])
sum(pool1000$lambda * (pool1000$alpha / pool1000$beta))

pool5000 <- rbind(params_risque_eleve,
                   params_data[sample(seq(total_pool_size), 4999),])
sum(pool5000$lambda * (pool5000$alpha / pool5000$beta))


pool10_contrib <- compute_cm(2^15, 1, pool10, "unbiased")
pool100_contrib <- compute_cm(2^18, 1, pool100, "unbiased")
pool1000_contrib <- compute_cm(2^18, 1, pool1000, "unbiased")
pool10000_contrib <- compute_cm(2^18, 5, pool10000, "unbiased")

