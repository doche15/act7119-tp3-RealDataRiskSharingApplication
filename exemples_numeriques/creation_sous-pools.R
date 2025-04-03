###
### Travail 3 du cours ACT-7119
###
### Création des sous-pools
### Utilisation de l'Algorithme des k-moyennes sur les variables de latitude
### et de longitude du jeu de données.
###
##


source("preprocessing/BelgianMTPL_P.R")
# Jeu de données : "test"


library(ggplot2)
set.seed(11)

# Selon les notes de cours du cours ACT-4114.
ngroups <- 3 # fixé pour avoir un nombre fixe de sous-pools.
long_lat_data <- data.frame(long = test$long, lat = test$lat)

scaled_data <- scale(long_lat_data)

kmeans_obj <- kmeans(scaled_data, ngroups, nstart = 10)

long_lat_data_with_groups <- long_lat_data %>%
    mutate(groups = kmeans_obj$cluster)

ggplot(long_lat_data_with_groups) +
    geom_point(aes(lat, long, col = groups))


# Avec prime pure? Je ne sais pas si cela est rigoureux
# méthodologiquement-parlant.
params_data <- read.csv("params_belgian.csv")
primes_pures <- params_data$lambda * (params_data$alpha / params_data$beta)

ngroups <- 3 # fixé pour avoir un nombre fixe de sous-pools.
long_lat_pprem_data <- data.frame(long = test$long,
                                  lat = test$lat,
                                  prime_pure = primes_pures)

scaled_data <- scale(long_lat_pprem_data)

kmeans_obj <- kmeans(scaled_data, ngroups, nstart = 10)

data_with_groups <- long_lat_pprem_data %>%
    mutate(pool = kmeans_obj$cluster)

ggplot(data_with_groups) +
    geom_point(aes(lat, long, size = prime_pure, col = pool))

aggregate(prime_pure ~ pool, data = data_with_groups, FUN = mean)


write.csv(data_with_groups, "sous-pools-k-3.csv", row.names = FALSE)
