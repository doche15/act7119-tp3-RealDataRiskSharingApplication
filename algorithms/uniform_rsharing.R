### ACT-7119: Modèles de risque avec dépendance et mesures de risque
### TP3 : Partage de risque
### Calcul de la règle de partage de risque unniforme [Denuit et al., 2022].
###
##


compute_unif_rsharing <- function(s, params_data)
{
    pool_size <- length(params_data$lambda)

    rep(s / pool_size, pool_size)
}