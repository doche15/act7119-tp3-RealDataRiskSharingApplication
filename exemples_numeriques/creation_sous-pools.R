###
### Travail 3 du cours ACT-7119
###
### Création des sous-pools
### Utilisation de l'Algorithme des k-moyennes sur les variables de latitude
### et de longitude du jeu de données.
###
##


# %%%%%%%%%%%%%%%%%%%%%
# Création de sous-pools.
# Par quantile selon le paramètre de fréquence.
# %%%%%%%%%%%%%%%%%%%%%
source("preprocessing/BelgianMTPL_P.R")
# test
# Jeu de données "beMTPL97" de CASdatasets selon ma compréhension.
str(test)
params_data <- read.csv("params_belgian.csv")
primes_pures <- params_data$lambda * (params_data$alpha / params_data$beta)

data_subpools <- data.frame(
    id = seq(length(params_data$lambda)),
    lambda = params_data$lambda,
    alpha = params_data$alpha,
    beta = params_data$beta,
    prime_pure = primes_pures,
    ageph = test$ageph,
    sex = test$sex,
    power = test$power,
    long = test$long,
    lat = test$lat
)

str(data_subpools)

n_subpools <- 3
quantile_sep <- c(0, seq(n_subpools) * (1 / n_subpools))

for (i in seq(n_subpools))
{
    # Création des sous-pools selon la prime pure ("lambda")
    quantile_inf <- quantile(data_subpools$prime_pure, quantile_sep[i])
    quantile_sup <- quantile(data_subpools$prime_pure, quantile_sep[i + 1])

    temp_data_subpool <- NULL
    if (i == n_subpools)
    {
        temp_data_subpool <- data_subpools[quantile_inf <= data_subpools$prime_pure & data_subpools$prime_pure <= quantile_sup,]
    }
    else
    {
        temp_data_subpool <- data_subpools[quantile_inf <= data_subpools$prime_pure & data_subpools$prime_pure < quantile_sup,]
    }

    write.csv(temp_data_subpool, paste0("donnees-sous-pools/sous-pool-", i, "-k", n_subpools, "--", Sys.Date(), ".csv"), row.names = FALSE)
}

subpool1 <- read.csv("donnees-sous-pools/sous-pool-1-k3--2025-04-04.csv")
subpool2 <- read.csv("donnees-sous-pools/sous-pool-2-k3--2025-04-04.csv")
subpool3 <- read.csv("donnees-sous-pools/sous-pool-3-k3--2025-04-04.csv")

summary(subpool1$prime_pure) # valid, OK
summary(subpool2$prime_pure) # valid, OK
summary(subpool3$prime_pure) # valid, OK

nrow(subpool1) + nrow(subpool2) + nrow(subpool3) # valid, OK.


summary(subpool1)
table(subpool1$sex) / nrow(subpool1)

summary(subpool2)
table(subpool2$sex) / nrow(subpool2)

summary(subpool3)
table(subpool3$sex) / nrow(subpool3)


# %%%%%%%%%%%%%%%%%%%%%
# Avec k-moyennes.
# %%%%%%%%%%%%%%%%%%%%%
#
# source("preprocessing/BelgianMTPL_P.R")
# # Jeu de données : "test"
#
#
# library(ggplot2)
# set.seed(11)
#
# # Selon les notes de cours du cours ACT-4114.
# ngroups <- 3 # fixé pour avoir un nombre fixe de sous-pools.
# long_lat_data <- data.frame(long = test$long, lat = test$lat)
#
# scaled_data <- scale(long_lat_data)
#
# kmeans_obj <- kmeans(scaled_data, ngroups, nstart = 10)
#
# long_lat_data_with_groups <- long_lat_data %>%
#     mutate(groups = kmeans_obj$cluster)
#
# ggplot(long_lat_data_with_groups) +
#     geom_point(aes(lat, long, col = groups))
#
#
# # Avec prime pure? Je ne sais pas si cela est rigoureux
# # méthodologiquement-parlant.
# params_data <- read.csv("params_belgian.csv")
# primes_pures <- params_data$lambda * (params_data$alpha / params_data$beta)
#
# ngroups <- 3 # fixé pour avoir un nombre fixe de sous-pools.
# long_lat_pprem_data <- data.frame(long = test$long,
#                                   lat = test$lat,
#                                   prime_pure = primes_pures)
#
# scaled_data <- scale(long_lat_pprem_data)
#
# kmeans_obj <- kmeans(scaled_data, ngroups, nstart = 10)
#
# data_with_groups <- long_lat_pprem_data %>%
#     mutate(pool = kmeans_obj$cluster)
#
# ggplot(data_with_groups) +
#     geom_point(aes(lat, long, size = prime_pure, col = pool))
#
# aggregate(prime_pure ~ pool, data = data_with_groups, FUN = mean)
#
#
# write.csv(data_with_groups, "sous-pools-k-3.csv", row.names = FALSE)
