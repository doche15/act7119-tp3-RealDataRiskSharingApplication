###
### Travail 3 du cours ACT-7119
###
### Analyser l'impact de la taille du pool sur un risque artificiel élevé.
###
##


library(ggplot2)


params_data <- read.csv("params_belgian.csv")
summary(params_data)
primes_pures <- params_data$lambda * (params_data$alpha / params_data$beta)

# Selon les valeurs maximales des paramètres, on crée un risque élevé.
params_risque_eleve <- c(0.5, 3, 0.0005)
prime_pure_risque_eleve <- params_risque_eleve[1] *
    (params_risque_eleve[2] / params_risque_eleve[3])

c(prime_pure_risque_eleve = prime_pure_risque_eleve, prime_pure_moyenne = mean(primes_pures))


# Calculs de partage de risque
subpool1 <- read.csv("donnees-sous-pools/sous-pool-1-k3--2025-04-04.csv")
subpool1 <- subpool1[,c(2, 3, 4)]

c(prime_pure_risque_eleve = prime_pure_risque_eleve,
  prime_pure_moyenne = mean(subpool1$lambda * (subpool1$alpha / subpool1$beta)))

source("algorithms/cond_mean_rsharing_FFT_Blier.R")
source("algorithms/mean_prop_rsharing.R")
source("algorithms/linear_reg_rsharing.R")

total_pool_size <- nrow(subpool1)

set.seed(11)
pool500 <- rbind(params_risque_eleve,
                 subpool1[sample(seq(total_pool_size), 499),])
ES_p500 <- sum(pool500$lambda * (pool500$alpha / pool500$beta))

# pool750 <- rbind(params_risque_eleve,
#                  subpool1[sample(seq(total_pool_size), 749),])
# sum(pool750$lambda * (pool750$alpha / pool750$beta))
#
# pool1000 <- rbind(params_risque_eleve,
#                   subpool1[sample(seq(total_pool_size), 999),])
# sum(pool1000$lambda * (pool1000$alpha / pool1000$beta))


# Règle cm
nfft <- 2^16
h <- 5
support_S <- 0:(nfft - 1) * h

POOL_SIZE_1 <- 500
# POOL_SIZE_2 <- 750
# POOL_SIZE_3 <- 1000

pool500_contrib <- compute_cm(nfft, h, pool500, "unbiased") # temps de calculs : 56.5 secondes.
# pool750_contrib <- compute_cm(nfft, h, pool750, "unbiased") # temps de calculs : 88 secondes.
# pool1000_contrib <- compute_cm(nfft, h, pool1000, "unbiased") # temps de calculs : 117.8 secondes

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, pool500_contrib$fs)
sum(pool500_contrib$fs)

# plot(0:(nfft-1) * h, pool750_contrib$fs)
# sum(pool750_contrib$fs)
#
# plot(0:(nfft-1) * h, pool1000_contrib$fs)
# sum(pool1000_contrib$fs)

# 25 000 à 75 000


results_contrib_pool500 <- data.frame(s = support_S[(10000:60000) / h],
                                      unif_rsharing = support_S[(10000:60000) / h] / POOL_SIZE_1,
                                      cm_rsharing = pool500_contrib$contrib[[1]][(10000:60000) / h])

results_contrib_pool500_joined <- rbind(data.frame(s = results_contrib_pool500$s,
                                                   contribution = results_contrib_pool500$unif_rsharing,
                                                   methode = rep("Uniforme", nrow(results_contrib_pool500))),
                                        data.frame(s = results_contrib_pool500$s,
                                                   contribution = results_contrib_pool500$cm_rsharing,
                                                   methode = rep("Par espérance conditionnelle", nrow(results_contrib_pool500))))


# Référence :
# https://stackoverflow.com/questions/14794599/how-to-change-line-width-in-ggplot
fig1 <- ggplot(results_contrib_pool500_joined) +
    geom_line(aes(s / 1000, contribution, col = methode), linewidth=1.6) +
    geom_vline(xintercept = ES_p500 / 1000, linetype=2) +
    annotate("text", x = 32000 / 1000, y = 32000, label = "E[S]", size=10) +
    xlab("s (en milliers)") +
    ylab("Contribution") +
    xlim(10000 / 1000, 60000 / 1000) +
    ylim(0, 40000) +
    theme_classic() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=24),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))
ggsave(fig1, file="risque-eleve-1.pdf", width=10, height=8)


