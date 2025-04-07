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
library(RColorBrewer)

source("algorithms/cond_mean_rsharing_FFT_Blier.R")

summary(params_data)

# création des pools ----
n <- c(3, 100, 1000) # taille des pools étudiés

set.seed(21)
pool_n3 <- params_data[sample(nrow(params_data), n[1]),]

set.seed(21)
pool_n100 <- params_data[sample(nrow(params_data), n[2]),]

set.seed(21)
pool_n1000 <- params_data[sample(nrow(params_data), n[3]),]

# calcul des espérances conditionnelles ----
kmax <- 2^18
h <- 1

cm_n3 <- compute_cm(kmax = kmax,
                    h = h,
                    params = pool_n3,
                    discr_method = 'unbiased')

cm_n100 <- compute_cm(kmax = kmax,
                      h = h,
                      params = pool_n100,
                      discr_method = 'unbiased')

cm_n1000 <- compute_cm(kmax = kmax,
                       h = h,
                       params = pool_n1000,
                       discr_method = 'unbiased')

# verif
val_range_n3 = seq(2e4)
plot(cm_n3$fs[val_range_n3])
sum(cm_n3$fs[val_range_n3])

val_range_n100 = seq(5e4)
plot(cm_n100$fs[val_range_n100])
sum(cm_n100$fs[val_range_n100])

val_range_n1000 = seq(4e4, 2e5)
plot(cm_n1000$fs[val_range_n1000])
sum(cm_n1000$fs[val_range_n1000])

# paramètres graphique ----
hjust = -0.3
xlimits = c(0, 200)

# analyse du risque 1 dans chacun des pools ----
risk <- 1

val_range_n3 = seq(1.7e4) # évite instabilité graphique
# contributions
contrib_n3_risk <- unlist(cm_n3$contrib[risk])[val_range_n3]
contrib_n100_risk <- unlist(cm_n100$contrib[risk])[val_range_n100]
contrib_n1000_risk <- unlist(cm_n1000$contrib[risk])[val_range_n1000]

# densité
fs3 <- cm_n3$fs[val_range_n3]
sum(fs3)
fs100 <- cm_n100$fs[val_range_n100]
sum(fs100)
fs1000 <- cm_n1000$fs[val_range_n1000]
sum(fs1000)

# prime pure
sum(contrib_n3_risk * fs3)
sum(contrib_n100_risk * fs100)
sum(contrib_n1000_risk * fs1000)
purepremium <- pool_n3$lambda[risk] * pool_n3$alpha[risk] / pool_n3$beta[risk]

# répartition
Fs3 <- cumsum(fs3)
Fs3[length(Fs3)]
Fs100 <- cumsum(fs100)
Fs100[length(Fs100)]
Fs1000 <- cumsum(fs1000)
Fs1000[length(Fs1000)]

# graphique
df_risk_n3 <- data.frame(x = contrib_n3_risk, 
                         y = Fs3,
                         group = 'n = 3')

df_risk_n100 <- data.frame(x = c(0, contrib_n100_risk, 200), 
                           y = c(0, Fs100, 1),
                           group = 'n = 100')

df_risk_n1000 <- data.frame(x = c(0, contrib_n1000_risk, 200), 
                            y = c(0, Fs1000, 1),
                            group = 'n = 1000')

df_tot <- rbind(df_risk_n3,
                df_risk_n100,
                df_risk_n1000)

df_tot$group = factor(df_tot$group, levels = c('n = 3', 'n = 100', 'n = 1000'))

plot_risk1 = ggplot(df_tot, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  labs(color = "Taille du pool",
       y = expression("P" * "[" * E(X[1] ~ "|" ~ S) <= x * "]"),
       title = substitute(paste("Risque 1 : ", lambda[1] == lambda_val, ', ',
                                   alpha[1] == alpha_val, ', ', beta[1] == beta_val),
                             list(lambda_val = round(pool_n3$lambda[risk], 3),
                                  alpha_val = round(pool_n3$alpha[risk], 3),
                                  beta_val = round(pool_n3$beta[risk], 4)))) +
  xlim(xlimits[1], xlimits[2]) +
  ylim(0, 1) +
  geom_vline(xintercept = purepremium, linetype = 'dashed') +
  theme_bw() +
  annotate("text",
           x = purepremium,
           y = 0.4,
           label = paste0("Prime pure = ", round(purepremium, 2)),
           hjust = hjust,
           size = 5) +
  theme(legend.position = 'bottom')

# analyse du risque 2 dans chacun des pools ----
risk <- 2

val_range_n3 = seq(2e4) # évite instabilité graphique
# contributions
contrib_n3_risk <- unlist(cm_n3$contrib[risk])[val_range_n3]
contrib_n100_risk <- unlist(cm_n100$contrib[risk])[val_range_n100]
contrib_n1000_risk <- unlist(cm_n1000$contrib[risk])[val_range_n1000]

# densité
fs3 <- cm_n3$fs[val_range_n3]
sum(fs3)
fs100 <- cm_n100$fs[val_range_n100]
sum(fs100)
fs1000 <- cm_n1000$fs[val_range_n1000]
sum(fs1000)

# prime pure
sum(contrib_n3_risk * fs3)
sum(contrib_n100_risk * fs100)
sum(contrib_n1000_risk * fs1000)
purepremium <- pool_n3$lambda[risk] * pool_n3$alpha[risk] / pool_n3$beta[risk]

# répartition
Fs3 <- cumsum(fs3)
Fs3[length(Fs3)]
Fs100 <- cumsum(fs100)
Fs100[length(Fs100)]
Fs1000 <- cumsum(fs1000)
Fs1000[length(Fs1000)]

# graphique
df_risk_n3 <- data.frame(x = contrib_n3_risk, 
                         y = Fs3,
                         group = 'n = 3')

df_risk_n100 <- data.frame(x = c(0, contrib_n100_risk, 200), 
                           y = c(0, Fs100, 1),
                           group = 'n = 100')

df_risk_n1000 <- data.frame(x = c(0, contrib_n1000_risk, 200), 
                            y = c(0, Fs1000, 1),
                            group = 'n = 1000')

df_tot <- rbind(df_risk_n3,
                df_risk_n100,
                df_risk_n1000)

df_tot$group = factor(df_tot$group, levels = c('n = 3', 'n = 100', 'n = 1000'))

plot_risk2 = ggplot(df_tot, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  labs(color = "Taille du pool",
       y = expression("P" * "[" * E(X[2] ~ "|" ~ S) <= x * "]"),
       title = substitute(paste("Risque 2 : ", lambda[2] == lambda_val, ', ',
                                   alpha[2] == alpha_val, ', ', beta[2] == beta_val),
                             list(lambda_val = round(pool_n3$lambda[risk], 3),
                                  alpha_val = round(pool_n3$alpha[risk], 3),
                                  beta_val = round(pool_n3$beta[risk], 4)))) +
  xlim(xlimits[1], xlimits[2]) +
  ylim(0, 1) +
  geom_vline(xintercept = purepremium, linetype = 'dashed') +
  theme_bw() +
  annotate("text",
           x = purepremium,
           y = 0.4,
           label = paste0("Prime pure = ", round(purepremium, 2)),
           hjust = hjust,
           size = 5) +
  theme(legend.position = 'bottom')

# analyse du risque 3 dans chacun des pools ----
risk <- 3
val_range_n3 = seq(1.6e4) # évite instabilité graphique

# contributions
contrib_n3_risk <- unlist(cm_n3$contrib[risk])[val_range_n3]
contrib_n100_risk <- unlist(cm_n100$contrib[risk])[val_range_n100]
contrib_n1000_risk <- unlist(cm_n1000$contrib[risk])[val_range_n1000]

# densité
fs3 <- cm_n3$fs[val_range_n3]
sum(fs3)
fs100 <- cm_n100$fs[val_range_n100]
sum(fs100)
fs1000 <- cm_n1000$fs[val_range_n1000]
sum(fs1000)

# prime pure
sum(contrib_n3_risk * fs3)
sum(contrib_n100_risk * fs100)
sum(contrib_n1000_risk * fs1000)
purepremium <- pool_n3$lambda[risk] * pool_n3$alpha[risk] / pool_n3$beta[risk]

# répartition
Fs3 <- cumsum(fs3)
Fs3[length(Fs3)]
Fs100 <- cumsum(fs100)
Fs100[length(Fs100)]
Fs1000 <- cumsum(fs1000)
Fs1000[length(Fs1000)]

# graphique
df_risk_n3 <- data.frame(x = contrib_n3_risk, 
                         y = Fs3,
                         group = 'n = 3')


test1 = c(0, Fs100, 1)
test2 = c(0, contrib_n100_risk, 200)
df_risk_n100 <- data.frame(x = test2, 
                           y = test1,
                           group = 'n = 100')

test3 = c(0, Fs1000, 1)
test4 = c(0, contrib_n1000_risk, 200)
df_risk_n1000 <- data.frame(x = test4, 
                            y = test3,
                            group = 'n = 1000')

df_tot <- rbind(df_risk_n3,
                df_risk_n100,
                df_risk_n1000)

df_tot$group = factor(df_tot$group, levels = c('n = 3', 'n = 100', 'n = 1000'))

plot_risk3 = ggplot(df_tot, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  labs(color = "Taille du pool",
       y = expression("P" * "[" * E(X[3] ~ "|" ~ S) <= x * "]"),
       title = substitute(paste("Risque 3 : ", lambda[3] == lambda_val, ', ',
                                   alpha[3] == alpha_val, ', ', beta[3] == beta_val),
                             list(lambda_val = round(pool_n3$lambda[risk], 3),
                                  alpha_val = round(pool_n3$alpha[risk], 3),
                                  beta_val = round(pool_n3$beta[risk], 4)))) +
  xlim(xlimits[1], xlimits[2]) +
  ylim(0, 1) +
  geom_vline(xintercept = purepremium, linetype = 'dashed') +
  theme_bw() +
  annotate("text",
           x = purepremium,
           y = 0.4,
           label = paste0("Prime pure = ", round(purepremium, 2)),
           hjust = hjust,
           size = 5) +
  theme(legend.position = 'bottom')

# légende ----
legend <- get_legend(ggplot(df_tot, aes(x = x, y = y, color = group)) +
                       geom_line(linewidth = 1) +
                       labs(color = "Taille du pool",
                            y = expression("P" * "[" * E(X[3] ~ "|" ~ S) <= x * "]"),
                            title = substitute(paste("Risque 3 : ", lambda[3] == lambda_val, ', ',
                                                     alpha[3] == alpha_val, ', ', beta[3] == beta_val),
                                               list(lambda_val = round(pool_n3$lambda[risk], 3),
                                                    alpha_val = round(pool_n3$alpha[risk], 3),
                                                    beta_val = round(pool_n3$beta[risk], 4)))) +
                       xlim(xlimits[1], xlimits[2]) +
                       ylim(0, 1) +
                       geom_vline(xintercept = purepremium, linetype = 'dashed') +
                       theme_bw() +
                       annotate("text",
                                x = purepremium,
                                y = 0.4,
                                label = paste0("Prime pure = ", round(purepremium, 2)),
                                hjust = hjust,
                                size = 5))

legend_plot <- plot_grid(legend)

# enregistre graphiques ----
# title = "Fonction de répartition des espérances conditionnelles"
ggsave("risk1.png", plot = plot_risk1, width = 5, height = 4, dpi = 300)
ggsave("risk2.png", plot = plot_risk2, width = 5, height = 4, dpi = 300)
ggsave("risk3.png", plot = plot_risk3, width = 5, height = 4, dpi = 300)
ggsave("legend.png", plot = legend_plot, width = 3, height = 1, dpi = 300)
