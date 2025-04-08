###
### Travail 3, ACT-7119
### Exemple numérique : Profils de risque et sous-pools
###
##


# Importation des paquetages
library(ggplot2)
library(ggpubr)


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


# Tableau 1
summary(subpool1_ech)
summary(subpool2_ech)
summary(subpool3_ech)


# Tableau 2
risk_ech_subpool1 <- subpool1_ech[subpool1_ech$id == id_risk_ech_subpool1,]
risk_ech_subpool2 <- subpool2_ech[subpool2_ech$id == id_risk_ech_subpool2,]
risk_ech_subpool3 <- subpool3_ech[subpool3_ech$id == id_risk_ech_subpool3,]


# Partage de risque
## Restructuration des échantillons de sous-pools pour que le risque
## analysé de chaque sous-pool soit à la ligne 1 des échantillons
subpool1_ech <- rbind(
    risk_ech_subpool1,
    subpool1_ech[!(subpool1_ech$id == id_risk_ech_subpool1),]
)
nrow(subpool1_ech)


subpool2_ech <- rbind(
    risk_ech_subpool2,
    subpool2_ech[!(subpool2_ech$id == id_risk_ech_subpool2),]
)
nrow(subpool2_ech)


subpool3_ech <- rbind(
    risk_ech_subpool3,
    subpool3_ech[!(subpool3_ech$id == id_risk_ech_subpool3),]
)
nrow(subpool3_ech)


## Partage de risque
source("algorithms/uniform_rsharing.R")
source("algorithms/cond_mean_rsharing_FFT_Blier.R")
source("algorithms/mean_prop_rsharing.R")
source("algorithms/linear_reg_rsharing.R")

#load("exemples_numeriques/results_contrib_rsharing_souspools.RData")
# Échantillon du sous-pool 1
nfft <- 2^17
h <- 5
support_S <- (c(0, seq(nfft - 1)) * h)

ES_sp1 <- sum(subpool1_ech$lambda * subpool1_ech$alpha / subpool1_ech$beta)
EX_sp1_r1 <- subpool1_ech$lambda[1] * subpool1_ech$alpha[1] / subpool1_ech$beta[1]
VarS_sp1 <- sum((subpool1_ech$lambda * subpool1_ech$alpha / (subpool1_ech$beta^2)) * (subpool1_ech$alpha + 1))
VarX_sp1_r1 <- (subpool1_ech$lambda[1] * subpool1_ech$alpha[1] / (subpool1_ech$beta[1]^2)) * (subpool1_ech$alpha[1] + 1)

unif_rsharing_sp1_r1 <- support_S / RISK_ECH_SIZE
mean_prop_rsharing_sp1_r1 <- (EX_sp1_r1 / ES_sp1) * support_S
linear_reg_rsharing_sp1_r1 <- EX_sp1_r1 + (VarX_sp1_r1 / VarS_sp1) * (support_S - ES_sp1)
cm_rsharing_sp1 <- compute_cm(nfft, h, subpool1_ech, "unbiased")
# Temps : 241.2 secondes, 2025-04-07, avec "unbiased".
cm_rsharing_sp1_r1 <- cm_rsharing_sp1$contrib[[1]]

results_contrib_rsharing_sp1 <- data.frame(s = support_S,
                                           unif_rsharing = unif_rsharing_sp1_r1,
                                           mean_prop_rsharing = mean_prop_rsharing_sp1_r1,
                                           linear_reg_rsharing = linear_reg_rsharing_sp1_r1,
                                           cm_rsharing = cm_rsharing_sp1_r1)

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, cm_rsharing_sp1$fs)
sum(cm_rsharing_sp1$fs)
# Entre environ 30 000 et 90 000.

cm_rsharing_sp1_r1[30100 / h]
cm_rsharing_sp1_r1[89900 / h]

results_contrib_rsharing_sp1_reduced <- results_contrib_rsharing_sp1[(30100:130900) / h,]
#results_contrib_rsharing_sp1_reduced <- results_contrib_rsharing_sp1[(30100:120000) / h,]
# results_contrib_rsharing_sp1_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp1_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp1_reduced$unif_rsharing,
#                                                                 methode = rep("Uniforme", nrow(results_contrib_rsharing_sp1_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp1_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp1_reduced$mean_prop_rsharing,
#                                                                 methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp1_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp1_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp1_reduced$linear_reg_rsharing,
#                                                                 methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp1_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp1_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp1_reduced$cm_rsharing,
#                                                                 methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp1_reduced))))


results_contrib_rsharing_sp1_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp1_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp1_reduced$mean_prop_rsharing - results_contrib_rsharing_sp1_reduced$unif_rsharing) / results_contrib_rsharing_sp1_reduced$unif_rsharing,
                                                                methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp1_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp1_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp1_reduced$linear_reg_rsharing - results_contrib_rsharing_sp1_reduced$unif_rsharing) / results_contrib_rsharing_sp1_reduced$unif_rsharing,
                                                                methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp1_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp1_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp1_reduced$cm_rsharing - results_contrib_rsharing_sp1_reduced$unif_rsharing) / results_contrib_rsharing_sp1_reduced$unif_rsharing,
                                                                methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp1_reduced))))


# Référence :
# https://stackoverflow.com/questions/14794599/how-to-change-line-width-in-ggplot
# fig1 <- ggplot(results_contrib_rsharing_sp1_reduced_joined) +
#     geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
#     xlab("s (en milliers)") +
#     ylab("Contribution") +
#     xlim(30000 / 1000, 90000 / 1000) +
#     ylim(20, 100) +
#     theme_classic() +
#     theme(axis.text = element_text(size=20),
#           axis.title = element_text(size=20),
#           legend.text = element_text(size=14),
# 	    legend.title = element_text(size=16),
#           legend.position = "bottom") +
#     guides(col = guide_legend(title = "Méthode"))


# Échantillon du sous-pool 2
ES_sp2 <- sum(subpool2_ech$lambda * subpool2_ech$alpha / subpool2_ech$beta)
EX_sp2_r1 <- subpool2_ech$lambda[1] * subpool2_ech$alpha[1] / subpool2_ech$beta[1]
VarS_sp2 <- sum((subpool2_ech$lambda * subpool2_ech$alpha / (subpool2_ech$beta^2)) * (subpool2_ech$alpha + 1))
VarX_sp2_r1 <- (subpool2_ech$lambda[1] * subpool2_ech$alpha[1] / (subpool2_ech$beta[1]^2)) * (subpool2_ech$alpha[1] + 1)

unif_rsharing_sp2_r1 <- support_S / RISK_ECH_SIZE
mean_prop_rsharing_sp2_r1 <- (EX_sp2_r1 / ES_sp2) * support_S
linear_reg_rsharing_sp2_r1 <- EX_sp2_r1 + (VarX_sp2_r1 / VarS_sp2) * (support_S - ES_sp2)
cm_rsharing_sp2 <- compute_cm(nfft, h, subpool2_ech, "unbiased")
# Temps : 252.5 secondes, 2025-04-07, avec "unbiased".

cm_rsharing_sp2_r1 <- cm_rsharing_sp2$contrib[[1]]

results_contrib_rsharing_sp2 <- data.frame(s = support_S,
                                           unif_rsharing = unif_rsharing_sp2_r1,
                                           mean_prop_rsharing = mean_prop_rsharing_sp2_r1,
                                           linear_reg_rsharing = linear_reg_rsharing_sp2_r1,
                                           cm_rsharing = cm_rsharing_sp2_r1)

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, cm_rsharing_sp2$fs)
sum(cm_rsharing_sp2$fs)
# Entre environ 50 000 et 120 000.

cm_rsharing_sp2_r1[50000 / h]
cm_rsharing_sp2_r1[120000 / h]


results_contrib_rsharing_sp2_reduced <- results_contrib_rsharing_sp2[(30100:200900) / h,]
#results_contrib_rsharing_sp2_reduced <- results_contrib_rsharing_sp2[(50000:120000) / h,]

# results_contrib_rsharing_sp2_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp2_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp2_reduced$unif_rsharing,
#                                                                 methode = rep("Uniforme", nrow(results_contrib_rsharing_sp2_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp2_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp2_reduced$mean_prop_rsharing,
#                                                                 methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp2_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp2_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp2_reduced$linear_reg_rsharing,
#                                                                 methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp2_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp2_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp2_reduced$cm_rsharing,
#                                                                 methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp2_reduced))))


results_contrib_rsharing_sp2_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp2_reduced$mean_prop_rsharing - results_contrib_rsharing_sp2_reduced$unif_rsharing) / results_contrib_rsharing_sp2_reduced$unif_rsharing,
                                                                methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp2_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp2_reduced$linear_reg_rsharing - results_contrib_rsharing_sp2_reduced$unif_rsharing) / results_contrib_rsharing_sp2_reduced$unif_rsharing,
                                                                methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp2_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp2_reduced$cm_rsharing - results_contrib_rsharing_sp2_reduced$unif_rsharing) / results_contrib_rsharing_sp2_reduced$unif_rsharing,
                                                                methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp2_reduced))))


j <-  data.frame(s = results_contrib_rsharing_sp2_reduced$s,
                 contribution = (results_contrib_rsharing_sp2_reduced$linear_reg_rsharing - results_contrib_rsharing_sp2_reduced$unif_rsharing) / results_contrib_rsharing_sp2_reduced$unif_rsharing,
                 methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp2_reduced)))

(j$contribution[10000] - j$contribution[10]) / (j$s[10000] - j$s[10])

(j$contribution[20000] - j$contribution[10000]) / (j$s[20000] - j$s[10000])


# Référence :
# https://stackoverflow.com/questions/14794599/how-to-change-line-width-in-ggplot
#ggplot(results_contrib_rsharing_sp2_reduced_joined) +
#    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
#    xlab("s (en milliers)") +
#    ylab("Contribution") +
#    xlim(50000 / 1000, 120000 / 1000) +
#    ylim(50, 130) +
#    theme_classic() +
#    theme(axis.text = element_text(size=12),
#          axis.title = element_text(size=12),
#          legend.text = element_text(size=9),
#          legend.position = "bottom") +
#    guides(col = guide_legend(title = "Méthode"))


# Échantillon du sous-pool 3
ES_sp3 <- sum(subpool3_ech$lambda * subpool3_ech$alpha / subpool3_ech$beta)
EX_sp3_r1 <- subpool3_ech$lambda[1] * subpool3_ech$alpha[1] / subpool3_ech$beta[1]
VarS_sp3 <- sum((subpool3_ech$lambda * subpool3_ech$alpha / (subpool3_ech$beta^2)) * (subpool3_ech$alpha + 1))
VarX_sp3_r1 <- (subpool3_ech$lambda[1] * subpool3_ech$alpha[1] / (subpool3_ech$beta[1]^2)) * (subpool3_ech$alpha[1] + 1)

unif_rsharing_sp3_r1 <- support_S / RISK_ECH_SIZE
mean_prop_rsharing_sp3_r1 <- (EX_sp3_r1 / ES_sp3) * support_S
linear_reg_rsharing_sp3_r1 <- EX_sp3_r1 + (VarX_sp3_r1 / VarS_sp3) * (support_S - ES_sp3)
cm_rsharing_sp3 <- compute_cm(nfft, h, subpool3_ech, "unbiased")
# Temps : 250.8 secondes, 2025-04-04, avec "unbiased".
cm_rsharing_sp3_r1 <- cm_rsharing_sp3$contrib[[1]]

results_contrib_rsharing_sp3 <- data.frame(s = support_S,
                                           unif_rsharing = unif_rsharing_sp3_r1,
                                           mean_prop_rsharing = mean_prop_rsharing_sp3_r1,
                                           linear_reg_rsharing = linear_reg_rsharing_sp3_r1,
                                           cm_rsharing = cm_rsharing_sp3_r1)

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, cm_rsharing_sp3$fs)
sum(cm_rsharing_sp3$fs)
# Entre environ 110 000 et 220 000.

cm_rsharing_sp3_r1[110000 / h]
cm_rsharing_sp3_r1[220000 / h]


results_contrib_rsharing_sp3_reduced <- results_contrib_rsharing_sp3[(100100:219900) / h,]
#results_contrib_rsharing_sp3_reduced <- results_contrib_rsharing_sp3[(110000:219900) / h,]

# results_contrib_rsharing_sp3_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$unif_rsharing,
#                                                                 methode = rep("Uniforme", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$mean_prop_rsharing,
#                                                                 methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$linear_reg_rsharing,
#                                                                 methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$cm_rsharing,
#                                                                 methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp3_reduced))))


results_contrib_rsharing_sp3_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp3_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp3_reduced$mean_prop_rsharing - results_contrib_rsharing_sp3_reduced$unif_rsharing) / results_contrib_rsharing_sp3_reduced$unif_rsharing,
                                                                methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp3_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp3_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp3_reduced$linear_reg_rsharing - results_contrib_rsharing_sp3_reduced$unif_rsharing) / results_contrib_rsharing_sp3_reduced$unif_rsharing,
                                                                methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp3_reduced))),
                                                     data.frame(s = results_contrib_rsharing_sp3_reduced$s,
                                                                contribution = (results_contrib_rsharing_sp3_reduced$cm_rsharing - results_contrib_rsharing_sp3_reduced$unif_rsharing) / results_contrib_rsharing_sp3_reduced$unif_rsharing,
                                                                methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp3_reduced))))


# save(results_contrib_rsharing_sp1,
#      results_contrib_rsharing_sp2,
#      results_contrib_rsharing_sp3, file = "exemples_numeriques/results_contrib_rsharing_souspools.RData")

# Référence :
# https://stackoverflow.com/questions/14794599/how-to-change-line-width-in-ggplot
#ggplot(results_contrib_rsharing_sp3_reduced_joined) +
#    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
#    xlab("s (en milliers)") +
#    ylab("Contribution") +
#    xlim(110000 / 1000, 220000 / 1000) +
#    ylim(50, 300) +
#    theme_classic() +
#    theme(axis.text = element_text(size=12),
#          axis.title = element_text(size=12),
#          legend.text = element_text(size=9),
#          legend.position = "bottom") +
#    guides(col = guide_legend(title = "Méthode"))


fig1 <- ggplot(results_contrib_rsharing_sp1_reduced_joined) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept = ES_sp1 / 1000, linetype=2) +
    annotate("text", x = 54000 / 1000, y = -0.01, label = "E[S]", size=10) +
    xlab("s (en milliers)") +
    ylab("Différence par rapport à la règle uniforme (%)") +
    title("Sous-pool 1") +
    xlim(30000 / 1000, 220000 / 1000) +
    theme_classic() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=24),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          title = element_text(size=22),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))
ggsave(fig1, file="profil-risque-1.pdf", width=10, height=8)


fig2 <- ggplot(results_contrib_rsharing_sp2_reduced_joined) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept = ES_sp2 / 1000, linetype=2) +
    annotate("text", x = 95000 / 1000, y = 0.1, label = "E[S]", size=10) +
    xlab("s (en milliers)") +
    ylab("Différence par rapport à la règle uniforme (%)") +
    title("Sous-pool 2") +
    xlim(30000 / 1000, 220000 / 1000) +
    theme_classic() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16),
          title = element_text(size=22),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))
ggsave(fig2, file="profil-risque-2-scale-commune.pdf", width=10, height=8)


fig3 <- ggplot(results_contrib_rsharing_sp3_reduced_joined) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept = ES_sp3 / 1000, linetype=2) +
    annotate("text", x = 173000 / 1000, y = 0.06, label = "E[S]", size=10) +
    xlab("s (en milliers)") +
    ylab("Différence par rapport à la règle uniforme (%)") +
    title("Sous-pool 3") +
    xlim(30000 / 1000, 220000 / 1000) +
    theme_classic() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16),
          title = element_text(size=22),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))
ggsave(fig3, file="profil-risque-3.pdf", width=10, height=8)


# Exportation finale des figures
# obj_ggarrange <- ggarrange(fig1, fig2, fig3,
#                            nrow = 3, common.legend = TRUE,
#                            labels = c("Sous-pool 1", "Sous-pool 2",
#                                       "Sous-pool 3"),
#                            legend = "bottom",
#                            font.label = list(size=20))

obj_ggarrange <- ggarrange(fig1, fig2, fig3,
                           nrow = 3, common.legend = TRUE,
                           legend = "bottom",
                           labels = c("Sous-pool 1",
                                      "Sous-pool 2",
                                      "Sous-pool 3"),
                           hjust = -1,
                           font.label = list(size=20))


ggsave(obj_ggarrange, file="profils-de-risque-ex-num-1.pdf",
       width = 10, height = 24)




# Tests supplémentaires
# (1) Mélange de sous-pools
# (2) Variation de n


# (1)

# risk_ech_subpool1 <- subpool1_ech[subpool1_ech$id == id_risk_ech_subpool1,]
# risk_ech_subpool2 <- subpool2_ech[subpool2_ech$id == id_risk_ech_subpool2,]
# risk_ech_subpool3 <- subpool3_ech[subpool3_ech$id == id_risk_ech_subpool3,]
set.seed(11)
subpool123_ech <- rbind(
    risk_ech_subpool1,
    risk_ech_subpool2,
    risk_ech_subpool3,
    subpool1[sample(seq(SUBPOOL1_SIZE), 499),],
    subpool2[sample(seq(SUBPOOL2_SIZE), 499),],
    subpool3[sample(seq(SUBPOOL3_SIZE), 499),]
)

subpool123_ech_SIZE <- nrow(subpool123_ech)


# Échantillon du sous-pool 1
nfft <- 2^17
h <- 5
support_S <- (c(0, seq(nfft - 1)) * h)

ES_pool_melange <- sum(subpool123_ech$lambda * subpool123_ech$alpha / subpool123_ech$beta)
EX_pool_melange_r1 <- subpool123_ech$lambda[1] * subpool123_ech$alpha[1] / subpool123_ech$beta[1]
EX_pool_melange_r2 <- subpool123_ech$lambda[2] * subpool123_ech$alpha[2] / subpool123_ech$beta[2]
EX_pool_melange_r3 <- subpool123_ech$lambda[3] * subpool123_ech$alpha[3] / subpool123_ech$beta[3]

VarS_pool_melange <- sum((subpool123_ech$lambda * subpool123_ech$alpha / (subpool123_ech$beta^2)) * (subpool123_ech$alpha + 1))

VarX_pool_melange_r1 <- (subpool123_ech$lambda[1] * subpool123_ech$alpha[1] / (subpool123_ech$beta[1]^2)) * (subpool123_ech$alpha[1] + 1)
VarX_pool_melange_r2 <- (subpool123_ech$lambda[2] * subpool123_ech$alpha[2] / (subpool123_ech$beta[2]^2)) * (subpool123_ech$alpha[2] + 1)
VarX_pool_melange_r3 <- (subpool123_ech$lambda[3] * subpool123_ech$alpha[3] / (subpool123_ech$beta[3]^2)) * (subpool123_ech$alpha[3] + 1)

unif_rsharing_pool_melange <- support_S / subpool123_ech_SIZE

mean_prop_rsharing_pool_melange_r1 <- (EX_pool_melange_r1 / ES_pool_melange) * support_S
mean_prop_rsharing_pool_melange_r2 <- (EX_pool_melange_r2 / ES_pool_melange) * support_S
mean_prop_rsharing_pool_melange_r3 <- (EX_pool_melange_r3 / ES_pool_melange) * support_S

linear_reg_rsharing_pool_melange_r1 <- EX_pool_melange_r1 + (VarX_pool_melange_r1 / VarS_pool_melange) * (support_S - ES_pool_melange)
linear_reg_rsharing_pool_melange_r2 <- EX_pool_melange_r2 + (VarX_pool_melange_r2 / VarS_pool_melange) * (support_S - ES_pool_melange)
linear_reg_rsharing_pool_melange_r3 <- EX_pool_melange_r3 + (VarX_pool_melange_r3 / VarS_pool_melange) * (support_S - ES_pool_melange)




cm_rsharing_pool_melange <- compute_cm(nfft, h, subpool123_ech, "unbiased")
# Temps : 367.4 secondes, 2025-04-07, avec "unbiased".
cm_rsharing_pool_melange_r1 <- cm_rsharing_pool_melange$contrib[[1]]
cm_rsharing_pool_melange_r2 <- cm_rsharing_pool_melange$contrib[[2]]
cm_rsharing_pool_melange_r3 <- cm_rsharing_pool_melange$contrib[[3]]

# results_contrib_rsharing_sp1 <- data.frame(s = support_S,
#                                            unif_rsharing = unif_rsharing_sp1_r1,
#                                            mean_prop_rsharing = mean_prop_rsharing_sp1_r1,
#                                            linear_reg_rsharing = linear_reg_rsharing_sp1_r1,
#                                            cm_rsharing = cm_rsharing_sp1_r1)

# À quelles valeurs de "S" présenter le graphique?
plot(0:(nfft-1) * h, cm_rsharing_pool_melange$fs)
sum(cm_rsharing_pool_melange$fs)

results_contrib_rsharing_reduced_joined_sp123_r1 <- rbind(data.frame(s = support_S[(100000:200000) / h],
                                                                     contribution = (mean_prop_rsharing_pool_melange_r1[(100000:200000) / h] - unif_rsharing_pool_melange[(100000:200000) / h]) / unif_rsharing_pool_melange[(100000:200000) / h],
                                                                     methode = rep("Au prorata des espérances", length(support_S[(100000:200000) / h]))),
                                                          data.frame(s = support_S[(100000:200000) / h],
                                                                     contribution = (linear_reg_rsharing_pool_melange_r1[(100000:200000) / h] - unif_rsharing_pool_melange[(100000:200000) / h]) / unif_rsharing_pool_melange[(100000:200000) / h],
                                                                     methode = rep("Régression linéaire", length(support_S[(100000:200000) / h]))),
                                                          data.frame(s = support_S[(100000:200000) / h],
                                                                     contribution = (cm_rsharing_pool_melange_r1[(100000:200000) / h] - unif_rsharing_pool_melange[(100000:200000) / h]) / unif_rsharing_pool_melange[(100000:200000) / h],
                                                                     methode = rep("Espérance conditionnelle", length(support_S[(100000:200000) / h]))))

results_contrib_rsharing_reduced_joined_sp123_r2 <- rbind(data.frame(s = support_S[(100000:200000) / h],
                                                                     contribution = (mean_prop_rsharing_pool_melange_r2[(100000:200000) / h] - unif_rsharing_pool_melange[(100000:200000) / h]) / unif_rsharing_pool_melange[(100000:200000) / h],
                                                                     methode = rep("Au prorata des espérances", length(support_S[(100000:200000) / h]))),
                                                          data.frame(s = support_S[(100000:200000) / h],
                                                                     contribution = (linear_reg_rsharing_pool_melange_r2[(100000:200000) / h] - unif_rsharing_pool_melange[(100000:200000) / h]) / unif_rsharing_pool_melange[(100000:200000) / h],
                                                                     methode = rep("Régression linéaire", length(support_S[(100000:200000) / h]))),
                                                          data.frame(s = support_S[(100000:200000) / h],
                                                                     contribution = (cm_rsharing_pool_melange_r2[(100000:200000) / h] - unif_rsharing_pool_melange[(100000:200000) / h]) / unif_rsharing_pool_melange[(100000:200000) / h],
                                                                     methode = rep("Espérance conditionnelle", length(support_S[(100000:200000) / h]))))

results_contrib_rsharing_reduced_joined_sp123_r3 <- rbind(data.frame(s = support_S,
                                                                     contribution = (mean_prop_rsharing_pool_melange_r3 - unif_rsharing_pool_melange) / unif_rsharing_pool_melange,
                                                                     methode = rep("Au prorata des espérances", length(support_S))),
                                                          data.frame(s = support_S,
                                                                     contribution = (linear_reg_rsharing_pool_melange_r3 - unif_rsharing_pool_melange) / unif_rsharing_pool_melange,
                                                                     methode = rep("Régression linéaire", length(support_S))),
                                                          data.frame(s = support_S,
                                                                     contribution = (cm_rsharing_pool_melange_r3- unif_rsharing_pool_melange) / unif_rsharing_pool_melange,
                                                                     methode = rep("Espérance conditionnelle", length(support_S))))



fig1 <- ggplot(results_contrib_rsharing_reduced_joined_sp123_r1) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept = ES_pool_melange / 1000, linetype=2) +
    annotate("text", x = 160000 / 1000, y = -0.425, label = "E[S]", size=10) +
    xlab("s (en milliers)") +
    ylab("Différence par rapport à la règle uniforme (%)") +
    title("Risque 1") +
    xlim(100000 / 1000, 200000 / 1000) +
    ylim(-0.6, -0.4) +
    theme_classic() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=24),
          legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          title = element_text(size=22),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))


fig2 <- ggplot(results_contrib_rsharing_reduced_joined_sp123_r2) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept = ES_pool_melange / 1000, linetype=2) +
    annotate("text", x = 160000 / 1000, y = 0.03, label = "E[S]", size=10) +
    xlab("s (en milliers)") +
    ylab("Différence par rapport à la règle uniforme (%)") +
    title("Risque 2") +
    xlim(100000 / 1000, 200000 / 1000) +
    ylim(-0.15, 0.05) +
    theme_classic() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16),
          title = element_text(size=22),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))


fig3 <- ggplot(results_contrib_rsharing_reduced_joined_sp123_r3) +
    geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept = ES_pool_melange / 1000, linetype=2) +
    annotate("text", x = 160000 / 1000, y = 0.73, label = "E[S]", size=10) +
    xlab("s (en milliers)") +
    ylab("Différence par rapport à la règle uniforme (%)") +
    title("Risque 3") +
    xlim(100000 / 1000, 200000 / 1000) +
    ylim(0.6, 0.75) +
    theme_classic() +
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16),
          title = element_text(size=22),
          legend.position = "bottom") +
    guides(col = guide_legend(title = "Méthode"))


# Exportation finale des figures
obj_ggarrange <- ggarrange(fig1, fig2, fig3,
                           nrow = 3, common.legend = TRUE,
                           legend = "bottom",
                           labels = c("Risque 1",
                                      "Risque 2",
                                      "Risque 3"),
                           hjust = -1,
                           font.label = list(size=20))


ggsave(obj_ggarrange, file="profils-de-risque-ex-num-1-melange-sous-pools.pdf",
       width = 10, height = 24)




# Validation
# results_contrib_rsharing_sp3_reduced <- results_contrib_rsharing_sp3[(80000:250000) / h,]
#
# results_contrib_rsharing_sp3_reduced_joined <- rbind(data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$unif_rsharing,
#                                                                 methode = rep("Uniforme", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$mean_prop_rsharing,
#                                                                 methode = rep("Au prorata des espérances", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$linear_reg_rsharing,
#                                                                 methode = rep("Régression linéaire", nrow(results_contrib_rsharing_sp3_reduced))),
#                                                      data.frame(s = results_contrib_rsharing_sp3_reduced$s,
#                                                                 contribution = results_contrib_rsharing_sp3_reduced$cm_rsharing,
#                                                                 methode = rep("Espérance conditionnelle", nrow(results_contrib_rsharing_sp3_reduced))))
# ggplot(results_contrib_rsharing_sp3_reduced_joined) +
#     geom_line(aes(s / 1000, contribution, col = methode), size=1.4) +
#     xlab("s (en milliers)") +
#     ylab("Contribution") +
#     xlim(80000 / 1000, 250000 / 1000) +
#     ylim(50, 300) +
#     theme_classic() +
#     theme(axis.text = element_text(size=12),
#           axis.title = element_text(size=12),
#           legend.text = element_text(size=9),
#           legend.position = "bottom") +
#     guides(col = guide_legend(title = "Méthode"))
