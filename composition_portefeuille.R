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


compute_cm(2^15, 0.5, params_data[1:1000, ], "unbiased")

## Analyse du portefeuille
summary(params_data)
lambdas <- params_data$lambda
alphas <- params_data$alpha
betas <- params_data$beta

hist(lambdas, freq = FALSE)
hist(alphas, freq = FALSE)
hist(betas, freq = FALSE)

purepremiums <- lambdas * alphas / betas
hist(purepremiums, freq = FALSE)

df_params <- cbind(test, alpha = alphas, beta = betas, lambda = lambdas,
                   premium = purepremiums)
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

aggregate(lambda ~ coverage, data = df_params, FUN = mean)
aggregate(lambda ~ sex, data = df_params, FUN = mean)
aggregate(lambda ~ use, data = df_params, FUN = mean)
aggregate(lambda ~ fuel, data = df_params, FUN = mean)
aggregate(lambda ~ ageph, data = df_params, FUN = mean) # intéressant


aggregate(alpha ~ coverage, data = df_params, FUN = mean)
aggregate(alpha ~ sex, data = df_params, FUN = mean)
aggregate(alpha ~ use, data = df_params, FUN = mean)
aggregate(alpha ~ fuel, data = df_params, FUN = mean)
aggregate(alpha ~ ageph, data = df_params, FUN = mean) # intéressant

