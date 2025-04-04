### ACT-7119: Modèles de risque avec dépendance et mesures de risque
### TP3 : Partage de risque
### Calcul de la règle de partage de risque au prorata des espérances de
### [Denuit et al., 2022].
###
### Par simulation.
##


simulate_PoisComp_gamma <- function(p_lambda, p_alpha, p_beta)
{
    freq <- rpois(1, p_lambda)
    ifelse(freq == 0, 0, sum(rgamma(freq, p_alpha, p_beta)))
}


compute_mean_prop_rsharing <- function(s, params_data)
{
    ES <- sum(params_data$lambda * params_data$alpha / params_data$beta)
    EXi <- params_data$lambda * params_data$alpha / params_data$beta

    s * (EXi / ES)
}


compute_mean_prop_rsharing_by_simulation <- function(nsim, params_data)
{
    set.seed(11)
    pool_size <- length(params_data$lambda)
    pool_realizations_increment_Sn <- numeric(nsim)

    # Barre de progression
    progress_bar <- txtProgressBar(1, pool_size, 1, "*", style = 3)
    for (i in seq(pool_size))
    {
        # On ne conserve pas les réalisations des risques individuels, car
        # ellkes ne sont pas utilisées. On souhaite seulement conserver les
        # "nsim" réalisations de S_n.
        pool_realizations_increment_Sn <- pool_realizations_increment_Sn +
            replicate(
                nsim,
                simulate_PoisComp_gamma(
                    params_data$lambda[i],
                    params_data$alpha[i],
                    params_data$beta[i]
                )
            )

        setTxtProgressBar(progress_bar, i)
    }

    ES <- sum(params_data$lambda * params_data$alpha / params_data$beta)
    linear_factors <- matrix((params_data$lambda * params_data$alpha / params_data$beta) / ES,
                             nrow = pool_size)

    Sn_realizations <- matrix(pool_realizations_increment_Sn, ncol = nsim)

    # Dimension de la matrice résultante des nsim contributions pour chaque
    # risque : pool_size x nsim
    list(contributions = linear_factors %*% Sn_realizations, Sn_realizations = Sn_realizations)
}

compute_mean_prop_rsharing_translativity <- function(s, params_data)
{
  ES <- sum(params_data$lambda * params_data$alpha / params_data$beta)
  EXi <- params_data$lambda * params_data$alpha / params_data$beta + params_data$cst
  
  s * (EXi / ES)
}

compute_mean_prop_rsharing_homogeneity <- function(s, params_data)
{
  ES <- sum(params_data$lambda * params_data$alpha / params_data$beta)
  EXi <- params_data$lambda * params_data$alpha / params_data$beta * params_data$fact
  
  s * (EXi / ES)
}


# Test 1
# params_data <- read.csv("params_belgian.csv")
# s <- 1000
# xx <- compute_mean_prop_rsharing(s, params_data[1:1000,])
#
# length(xx)
#
# c(sum(xx), s) # OK


# Test 2
# params_data <- read.csv("params_belgian.csv")
# xx <- compute_mean_prop_rsharing(10000, params_data[1:1000,])
#
# dim(xx$contributions)
#
# mean_contributions <- rowMeans(xx$contributions)
# ES <- sum(params_data$lambda[1:1000] * params_data$alpha[1:1000] / params_data$beta[1:1000])
# c(sum(mean_contributions), ES) # OK

