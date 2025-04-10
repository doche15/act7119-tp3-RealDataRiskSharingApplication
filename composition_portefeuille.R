### ACT-7119: Modèles de risque avec dépendance et mesures de risque
### TP3 : Partage de risque
## Extraction des profils de risque
library(ggplot2)
library(maps)
library(RColorBrewer)

# Importation des paramètres et des données
params_data <- read.csv("params_belgian.csv")
source("preprocessing/BelgianMTPL_P.R")

# Fonction de calcul de l'espérance conditionnelle
source("algorithms/cond_mean_rsharing_FFT_Blier.R")

## Analyse du portefeuille
summary(params_data)
lambdas <- params_data$lambda
alphas <- params_data$alpha
betas <- params_data$beta

compute_cm(2^15, 0.5, params_data[1:1000, ], "unbiased")

hist(lambdas, freq = FALSE)
hist(alphas, freq = FALSE)
hist(betas, freq = FALSE)

purepremiums <- lambdas * alphas / betas
hist(purepremiums, freq = FALSE)

df_params <- cbind(test, alpha = alphas, beta = betas, lambda = lambdas,
                   premium = purepremiums)

ggplot(df_params, aes(x = purepremiums, y = after_stat(density))) +
    geom_histogram(fill = "skyblue", binwidth = 5) +
    geom_vline(xintercept = 101.3983, color = "black", linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = 101.3983, y = 0.0063, label = "Prime moyenne: 101.3983", hjust = -0.1, color = "black") +
    theme_minimal() +
    labs(x = "Prime pure (euros)", y = "Densité")
ggsave("repart_primes.pdf")

summary(df_params)

ggplot(df_params, aes(x = long, y = lat, size = premium,
                      col = premium)) + geom_point(alpha = 0.2)

mapBelgium <- borders("world", regions="Belgium", colour="gray50", fill="white")

ggplot(df_params, aes(x = long, y = lat,
                      col = lambda)) + mapBelgium + geom_point(alpha = 1) +
    scale_color_distiller(palette = "OrRd", direction=1)

ggplot(df_params, aes(x = bm, y = lambda)) + geom_point(alpha = 0.2)


aggregate(premium ~ coverage, data = df_params, FUN = mean)
aggregate(premium ~ sex, data = df_params, FUN = mean)
aggregate(premium ~ use, data = df_params, FUN = mean)
aggregate(premium ~ fuel, data = df_params, FUN = mean)
aggregate(premium ~ ageph, data = df_params, FUN = mean) # intéressant
ggplot(aggregate(premium ~ ageph, data = df_params, FUN = mean), aes(x = ageph, y = premium)) +
  geom_point() # belle quadratique

aggregate(lambda ~ coverage, data = df_params, FUN = mean)
aggregate(lambda ~ sex, data = df_params, FUN = mean)
aggregate(lambda ~ use, data = df_params, FUN = mean)
aggregate(lambda ~ fuel, data = df_params, FUN = mean)
aggregate(lambda ~ ageph, data = df_params, FUN = mean) # intéressant
ggplot(aggregate(lambda ~ ageph, data = df_params, FUN = mean), aes(x = ageph, y = lambda)) +
  geom_point() # similaire

aggregate(alpha ~ coverage, data = df_params, FUN = mean)
aggregate(alpha ~ sex, data = df_params, FUN = mean)
aggregate(alpha ~ use, data = df_params, FUN = mean)
aggregate(alpha ~ fuel, data = df_params, FUN = mean)
aggregate(alpha ~ ageph, data = df_params, FUN = mean) # intéressant
ggplot(aggregate(alpha ~ ageph, data = df_params, FUN = mean), aes(x = ageph, y = alpha)) +
  geom_point()


