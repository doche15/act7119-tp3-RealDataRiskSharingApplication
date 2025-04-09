###
### Travail 3 du cours ACT-7119
###
### Analyser l'impact de la taille du pool sur la convergence de l'espérance 
### conditionnelle vers l'espérance marginale
###
### On veut observer le comportement des 3 premiers risques; ce sont les mêmes
### dans chaque pool
##

library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
display.brewer.all(n=5, type="seq", exact.n=TRUE,colorblindFriendly=F)

source("algorithms/cond_mean_rsharing_FFT_Blier.R")

summary(params_data)

# trouve les trois risques représentatifs ----
# Importation des données de sous-pools
subpool1 <- read.csv("donnees-sous-pools/sous-pool-1-k3--2025-04-04.csv")
subpool2 <- read.csv("donnees-sous-pools/sous-pool-2-k3--2025-04-04.csv")
subpool3 <- read.csv("donnees-sous-pools/sous-pool-3-k3--2025-04-04.csv")


# Constantes
SUBPOOL1_SIZE <- nrow(subpool1)
SUBPOOL2_SIZE <- nrow(subpool2)
SUBPOOL3_SIZE <- nrow(subpool3)

RISK_ECH_SIZE <- 1000
PURE_PREMIUM_RATE_RISK_CHOICE <- 0.05


# Échantillonnage des 1 000 risques des sous-pools
set.seed(11)

subpool1_ech <- subpool1[sample(seq(SUBPOOL1_SIZE), RISK_ECH_SIZE),]
subpool2_ech <- subpool2[sample(seq(SUBPOOL2_SIZE), RISK_ECH_SIZE),]
subpool3_ech <- subpool3[sample(seq(SUBPOOL3_SIZE), RISK_ECH_SIZE),]

prime_pure_moyennes_ech <-c(mean(subpool1_ech$prime_pure),
                            mean(subpool2_ech$prime_pure),
                            mean(subpool3_ech$prime_pure))


str(subpool1_ech)

# Choix du risque typique analysé dans chaque échantillon de sous-pool
id_risk_ech_subpool1 <- subpool1_ech$id[sample(
  which(
    (1 - PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[1] <= subpool1_ech$prime_pure &
      subpool1_ech$prime_pure <= (1 + PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[1]
  ),
  1
)]

id_risk_ech_subpool2 <- subpool2_ech$id[sample(
  which(
    (1 - PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[2] <= subpool2_ech$prime_pure &
      subpool2_ech$prime_pure <= (1 + PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[2]
  ),
  1
)]

id_risk_ech_subpool3 <- subpool3_ech$id[sample(
  which(
    (1 - PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[3] <= subpool3_ech$prime_pure &
      subpool3_ech$prime_pure <= (1 + PURE_PREMIUM_RATE_RISK_CHOICE) * prime_pure_moyennes_ech[3]
  ),
  1
)]

# création des pools ----
n <- c(3, 50, 100, 1000) # taille des pools étudiés

set.seed(21)
pool_nmax <- rbind(params_data[c(id_risk_ech_subpool1,
                                 id_risk_ech_subpool2,
                                 id_risk_ech_subpool3),],
                   params_data[sample(seq(nrow(params_data))[-c(id_risk_ech_subpool1,
                                                                id_risk_ech_subpool2,
                                                                id_risk_ech_subpool3)], 
                                      max(n) - n[1]),])

# calcul des espérances conditionnelles ----
kmax <- 2^18
h <- 1

cm_id1 <- compute_cm(kmax = kmax,
                     h = h,
                     params = pool_nmax[1,],
                     discr_method = 'unbiased')

cm_id2 <- compute_cm(kmax = kmax,
                     h = h,
                     params = pool_nmax[2,],
                     discr_method = 'unbiased')

cm_id3 <- compute_cm(kmax = kmax,
                     h = h,
                     params = pool_nmax[3,],
                     discr_method = 'unbiased')

cm_n3 <- compute_cm(kmax = kmax,
                    h = h,
                    params = pool_nmax[seq(n[1]),],
                    discr_method = 'unbiased')

cm_n50 <- compute_cm(kmax = kmax,
                    h = h,
                    params = pool_nmax[seq(n[2]),],
                    discr_method = 'unbiased')

cm_n100 <- compute_cm(kmax = kmax,
                      h = h,
                      params = pool_nmax[seq(n[3]),],
                      discr_method = 'unbiased')

cm_n1000 <- compute_cm(kmax = kmax,
                      h = h,
                      params = pool_nmax[seq(n[4]),],
                      discr_method = 'unbiased')

# verif
val_range_n1 = seq(1.8e4)
plot(cm_id1$fs[val_range_n1])
sum(cm_id1$fs[val_range_n1])

plot(cm_id2$fs[val_range_n1])
sum(cm_id2$fs[val_range_n1])

plot(cm_id3$fs[val_range_n1])
sum(cm_id3$fs[val_range_n1])

val_range_n3 = seq(2e4)
plot(cm_n3$fs[val_range_n3])
sum(cm_n3$fs[val_range_n3])

val_range_n50 = seq(5e4)
plot(cm_n50$fs[val_range_n50])
sum(cm_n50$fs[val_range_n50])

val_range_n100 = seq(5e4)
plot(cm_n100$fs[val_range_n100])
sum(cm_n100$fs[val_range_n100])

val_range_n1000 = seq(4e4, 2e5)
plot(cm_n1000$fs[val_range_n1000])
sum(cm_n1000$fs[val_range_n1000])

# paramètres graphique ----
hjust = -0.25
xlimits = c(0, 280)

wrapper_graph_x = function(vec) c(xlimits[1], vec, xlimits[2])
wrapper_graph_y = function(vec) c(0, vec, 1)

# distribution de S ----
# densité
fs1_id1 <- cm_id1$fs[val_range_n1]
sum(fs1_id1)
fs1_id2 <- cm_id2$fs[val_range_n1]
sum(fs1_id2)
fs1_id3 <- cm_id3$fs[val_range_n1]
sum(fs1_id3)
fs3 <- cm_n3$fs[val_range_n3]
sum(fs3)
fs50 <- cm_n50$fs[val_range_n50]
sum(fs50)
fs100 <- cm_n100$fs[val_range_n100]
sum(fs100)
fs1000 <- cm_n1000$fs[val_range_n1000]
sum(fs1000)

# répartition
Fs1_id1 <- cumsum(fs1_id1)
Fs1_id1[length(Fs1_id1)]
Fs1_id2 <- cumsum(fs1_id2)
Fs1_id2[length(Fs1_id2)]
Fs1_id3 <- cumsum(fs1_id3)
Fs1_id3[length(Fs1_id3)]
Fs3 <- cumsum(fs3)
Fs3[length(Fs3)]
Fs50 <- cumsum(fs50)
Fs50[length(Fs50)]
Fs100 <- cumsum(fs100)
Fs100[length(Fs100)]
Fs1000 <- cumsum(fs1000)
Fs1000[length(Fs1000)]

# analyse du risque 1 dans chacun des pools ----
# val_range_n3_insta = seq(1.7e4) # évite instabilité graphique
# contributions
contrib_id1_risk1 <- unlist(cm_id1$contrib)[val_range_n1]
contrib_n3_risk1 <- unlist(cm_n3$contrib[1])[val_range_n3]
contrib_n50_risk1 <- unlist(cm_n50$contrib[1])[val_range_n50]
contrib_n100_risk1 <- unlist(cm_n100$contrib[1])[val_range_n100]
contrib_n1000_risk1 <- unlist(cm_n1000$contrib[1])[val_range_n1000]

# prime pure
sum(contrib_id1_risk1 * fs1_id1)
sum(contrib_n3_risk1 * fs3)
sum(contrib_n50_risk1 * fs50)
sum(contrib_n100_risk1 * fs100)
sum(contrib_n1000_risk1 * fs1000)
(purepremium1 <- pool_nmax$lambda[1] * pool_nmax$alpha[1] / pool_nmax$beta[1])

# graphique
df_risk1_n1 <- data.frame(x = contrib_id1_risk1,
                         y = Fs1_id1,
                         group = 'n = 1')

df_risk1_n3 <- data.frame(x = contrib_n3_risk1, 
                         y = Fs3,
                         group = 'n = 3')

df_risk1_n50 <- data.frame(x = wrapper_graph_x(contrib_n50_risk1), 
                         y = wrapper_graph_y(Fs50),
                         group = 'n = 50')

df_risk1_n100 <- data.frame(x = wrapper_graph_x(contrib_n100_risk1), 
                           y = wrapper_graph_y(Fs100),
                           group = 'n = 100')

df_risk1_n1000 <- data.frame(x = wrapper_graph_x(contrib_n1000_risk1), 
                            y = wrapper_graph_y(Fs1000),
                            group = 'n = 1000')

df1_tot <- rbind(df_risk1_n1,
                df_risk1_n3,
                df_risk1_n50,
                df_risk1_n100,
                df_risk1_n1000)

df1_tot$group = factor(df1_tot$group, levels = c('n = 1', 'n = 3', 'n = 50', 'n = 100', 'n = 1000'))

plot_risk1 = ggplot(df1_tot, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  labs(color = "Taille du pool",
       subtitle = paste0("Prime pure = ", round(purepremium1, 2)),
       y = expression("P" * "[" * E(X[1] ~ "|" ~ S) <= x * "]"),
       title = substitute(paste("Risque 1 : ", lambda[1] == lambda_val, ', ',
                                   alpha[1] == alpha_val, ', ', beta[1] == beta_val),
                             list(lambda_val = round(pool_nmax$lambda[1], 3),
                                  alpha_val = round(pool_nmax$alpha[1], 3),
                                  beta_val = round(pool_nmax$beta[1], 4)))) +
  xlim(xlimits[1], xlimits[2]) +
  ylim(0, 1) +
  geom_vline(xintercept = purepremium1, linetype = 'dashed') +
  theme_bw() +
  scale_color_brewer(palette = "YlOrRd") +
  theme(legend.position = 'bottom')

# analyse du risque 2 dans chacun des pools ----
# val_range_n3 = seq(2e4) # évite instabilité graphique
# contributions
contrib_id2_risk2 <- unlist(cm_id2$contrib)[val_range_n1]
contrib_n3_risk2 <- unlist(cm_n3$contrib[2])[val_range_n3]
contrib_n50_risk2 <- unlist(cm_n50$contrib[2])[val_range_n50]
contrib_n100_risk2 <- unlist(cm_n100$contrib[2])[val_range_n100]
contrib_n1000_risk2 <- unlist(cm_n1000$contrib[2])[val_range_n1000]

# prime pure
sum(contrib_id2_risk2 * fs1_id2)
sum(contrib_n3_risk2 * fs3)
sum(contrib_n50_risk2 * fs50)
sum(contrib_n100_risk2 * fs100)
sum(contrib_n1000_risk2 * fs1000)
(purepremium2 <- pool_nmax$lambda[2] * pool_nmax$alpha[2] / pool_nmax$beta[2])

# graphique
df_risk2_n1 <- data.frame(x = contrib_id2_risk2, 
                         y = Fs1_id2,
                         group = 'n = 1')

df_risk2_n3 <- data.frame(x = contrib_n3_risk2, 
                         y = Fs3,
                         group = 'n = 3')

df_risk2_n50 <- data.frame(x = wrapper_graph_x(contrib_n50_risk2), 
                         y = wrapper_graph_y(Fs50),
                         group = 'n = 50')

df_risk2_n100 <- data.frame(x = wrapper_graph_x(contrib_n100_risk2), 
                           y = wrapper_graph_y(Fs100),
                           group = 'n = 100')

df_risk2_n1000 <- data.frame(x = wrapper_graph_x(contrib_n1000_risk2), 
                            y = wrapper_graph_y(Fs1000),
                            group = 'n = 1000')

df_tot2 <- rbind(df_risk2_n1,
                df_risk2_n3,
                df_risk2_n50,
                df_risk2_n100,
                df_risk2_n1000)

df_tot2$group = factor(df_tot2$group, levels = c('n = 1', 'n = 3', 'n = 50', 'n = 100', 'n = 1000'))

plot_risk2 = ggplot(df_tot2, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  labs(color = "Taille du pool",
       subtitle = paste0("Prime pure = ", round(purepremium2, 2)),
       y = expression("P" * "[" * E(X[2] ~ "|" ~ S) <= x * "]"),
       title = substitute(paste("Risque 2 : ", lambda[2] == lambda_val, ', ',
                                   alpha[2] == alpha_val, ', ', beta[2] == beta_val),
                             list(lambda_val = round(pool_nmax$lambda[2], 3),
                                  alpha_val = round(pool_nmax$alpha[2], 3),
                                  beta_val = round(pool_nmax$beta[2], 4)))) +
  xlim(xlimits[1], xlimits[2]) +
  ylim(0, 1) +
  geom_vline(xintercept = purepremium2, linetype = 'dashed') +
  theme_bw() +
  scale_color_brewer(palette = "YlOrRd") +
  theme(legend.position = 'bottom')

# analyse du risque 3 dans chacun des pools ----
# val_range_n3 = seq(1.6e4) # évite instabilité graphique
# contributions
contrib_id3_risk3 <- unlist(cm_id3$contrib)[val_range_n1]
contrib_n3_risk3 <- unlist(cm_n3$contrib[3])[val_range_n3]
contrib_n50_risk3 <- unlist(cm_n50$contrib[3])[val_range_n50]
contrib_n100_risk3 <- unlist(cm_n100$contrib[3])[val_range_n100]
contrib_n1000_risk3 <- unlist(cm_n1000$contrib[3])[val_range_n1000]

# prime pure
sum(contrib_id3_risk3 * fs1_id3)
sum(contrib_n3_risk3 * fs3)
sum(contrib_n50_risk3 * fs50)
sum(contrib_n100_risk3 * fs100)
sum(contrib_n1000_risk3 * fs1000)
(purepremium3 <- pool_nmax$lambda[3] * pool_nmax$alpha[3] / pool_nmax$beta[3])

# graphique
df_risk3_n1 <- data.frame(x = contrib_id3_risk3, 
                         y = Fs1_id3,
                         group = 'n = 1')

df_risk3_n3 <- data.frame(x = contrib_n3_risk3, 
                         y = Fs3,
                         group = 'n = 3')

df_risk3_n50 <- data.frame(x = contrib_n50_risk3, 
                         y = Fs50,
                         group = 'n = 50')

df_risk3_n100 <- data.frame(x = contrib_n100_risk3, 
                           y = Fs100,
                           group = 'n = 100')

df_risk3_n1000 <- data.frame(x = wrapper_graph_x(contrib_n1000_risk3), 
                            y = wrapper_graph_y(Fs1000),
                            group = 'n = 1000')

df3_tot <- rbind(df_risk3_n1,
                df_risk3_n3,
                df_risk3_n50,
                df_risk3_n100,
                df_risk3_n1000)

df3_tot$group = factor(df3_tot$group, levels = c('n = 1', 'n = 3', 'n = 50', 'n = 100', 'n = 1000'))

plot_risk3 = ggplot(df3_tot, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  labs(color = "Taille du pool",
       subtitle = paste0("Prime pure = ", round(purepremium3, 2)),
       y = expression("P" * "[" * E(X[3] ~ "|" ~ S) <= x * "]"),
       title = substitute(paste("Risque 3 : ", lambda[3] == lambda_val, ', ',
                                   alpha[3] == alpha_val, ', ', beta[3] == beta_val),
                             list(lambda_val = round(pool_nmax$lambda[3], 3),
                                  alpha_val = round(pool_nmax$alpha[3], 3),
                                  beta_val = round(pool_nmax$beta[3], 4)))) +
  xlim(xlimits[1], xlimits[2]) +
  ylim(0, 1) +
  geom_vline(xintercept = purepremium3, linetype = 'dashed') +
  theme_bw() +
  scale_color_brewer(palette = "YlOrRd") +
  theme(legend.position = 'bottom')

# enregistre graphiques ----
# title = "Fonction de répartition des espérances conditionnelles"
combined_plot = ggarrange(plot_risk1,
                          plot_risk2,
                          plot_risk3, 
                          ncol = 3,
                          common.legend = TRUE,
                          legend = 'bottom')

ggsave("risk_all.pdf", plot = combined_plot, width = 13, height = 4, dpi = 300)





