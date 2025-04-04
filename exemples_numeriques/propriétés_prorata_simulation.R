### ACT-7119: Modèles de risque avec dépendance et mesures de risque
### TP3 : Partage de risque
###
### Règle de partage de risque au prorata
###
### Démontration des propriétés par simulation
##

source("algorithms/mean_prop_rsharing.R")

params_data <- read.csv("params_belgian.csv")
summary(params_data)

sim_pois_gamma <- numeric(nrow(params_data)) # initialise vecteur des simulations

for (i in nrow(params_data)){
  
  set.seed(232) # point d'ancrage
  
  sim_pois_gamma[i] <- simulate_PoisComp_gamma(params_data$lambda[i],
                                               params_data$alpha[i],
                                               params_data$beta[i])
  
}

s <- sum(sim_pois_gamma) # perte totale

prop_rsharing <- compute_mean_prop_rsharing(s, params_data) # contributions

# conservation ----
## reshuffling ----
set.seed(21) # point d'ancrage
id <- sample(nrow(params_data), nrow(params_data))

params_data_shuffled <- params_data[id,] # shuffle

prop_rsharing_reshuffling <- compute_mean_prop_rsharing(s, params_data_shuffled)

sum(round(prop_rsharing[id], 8) == round(prop_rsharing_reshuffling, 8)) # contributions inchangées

## normalization ----
# crée un risque nul artificiellement (en évitant une division par zéro; beta = 1)
params_data_normalization <- rbind(c(0, 0, 1), # ne respecte pas le support des paramètres, mais R le traite bien comme une constante nulle
                                   params_data)

sim_pois_gamma_norm <- c(0, sim_pois_gamma) # ajoute la constante nulle

s <- sum(sim_pois_gamma_norm) # perte totale

prop_rsharing_normalization <- compute_mean_prop_rsharing(s, params_data_normalization)

prop_rsharing_normalization[1] # contribution nulle
sum(round(prop_rsharing_normalization[-1], 8) == round(prop_rsharing, 8)) # autres contributions inchangées

## translativity ----
# comme on staisfait au reshuffling, c'est suffisant de tester la translativity
# sur un seul participant (on le fait pour i = 1)
constante <- 10

params_data_translativity <- cbind(params_data,
                                   cst = c(constante, rep(0, nrow(params_data) - 1)))

sim_pois_gamma_trans <- sim_pois_gamma + params_data_translativity$cst # ajoute la constante

s <- sum(sim_pois_gamma_trans) # perte totale (augmente de "constante = 10")

prop_rsharing_translativity <- compute_mean_prop_rsharing_translativity(s, params_data_translativity)

round(prop_rsharing_translativity[1], 8) == round(prop_rsharing[1] + constante, 8) # on ne satisfait pas la translativité

## positive homogeneity ----
facteur <- 1.5

params_data_homogeneity <- cbind(params_data,
                                 fact = c(facteur, rep(1, nrow(params_data) - 1)))

sim_pois_gamma_homogeneity <- sim_pois_gamma * params_data_homogeneity$fact # ajoute le facteur

s <- sum(sim_pois_gamma_homogeneity) # perte totale

prop_rsharing_homogeneity <- compute_mean_prop_rsharing_homogeneity(s, params_data_homogeneity)

round(prop_rsharing_homogeneity[1], 8) == round(prop_rsharing[1] * facteur, 8) # même chose

## constancy ----
# crée une perte certaine artificiellement (on le fait pour i = 1)
perte <- 10

# mélange de normalization et translativity
params_data_constancy <- cbind(rbind(c(0, 0, 1), params_data), # ajoute risque nul
                               cst = c(perte, rep(0, nrow(params_data)))) # ajoute perte certaine

sim_pois_gamma_trans <- sim_pois_gamma_norm + params_data_constancy$cst # ajoute la perte certaine

s <- sum(sim_pois_gamma_trans) # perte totale (augmente de "perte = 10")

prop_rsharing_constancy <- compute_mean_prop_rsharing_translativity(s, params_data_constancy)

# on ne satisfait pas la constancy
round(prop_rsharing_constancy[1], 8) == perte
sum(round(prop_rsharing_constancy[-1], 8) == round(prop_rsharing, 8))

## no-ripoff ----
# Dans notre scénario, on a :
# F_X^{-1} (1) = \infty car X ~ PoiComp Gamma, donc
sum(prop_rsharing <= Inf)

## actuarial fairness ----
# à revoir si la méthodologie est correcte
n <- 10000

prop_rsharing_list <- sim_pois_gamma_list <- list()

for (years in seq(n)){
  
  sim_pois_gamma <- numeric(nrow(params_data)) # initialise vecteur des simulations
  
  set.seed(years) # point d'ancrage
  
  for (i in nrow(params_data)){
    
    sim_pois_gamma[i] <- simulate_PoisComp_gamma(params_data$lambda[i],
                                                 params_data$alpha[i],
                                                 params_data$beta[i])
    
  }
  
  s <- sum(sim_pois_gamma) # perte totale
  
  sim_pois_gamma_list[[years]] <- sim_pois_gamma
  prop_rsharing_list[[years]] <- compute_mean_prop_rsharing(s, params_data)
  
}

moy_si <- moy_contrib <- numeric(nrow(params_data))
for (indiv in nrow(params_data)){
  
  moy_si[indiv] <- mean(sapply(sim_pois_gamma_list, function(x) x[indiv]))
  moy_contrib[indiv] <- mean(sapply(prop_rsharing_list, function(x) x[indiv]))
  
}

# improvement ----
# à suivre...


