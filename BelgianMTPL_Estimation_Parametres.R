#### ACT-7119 - Travail pratique 3
### Partage de risque
### Exemple numérique sur un portefeuille de risques hétérogènes indépendants
### avec fréquence Poisson et sévérité Gamma

### Modèles entraînés dans le projet de recherche sur le boosting probabiliste
### de Chevalier et Côté (2025)

### https://arxiv.org/abs/2412.14916

library(reticulate)
library(dplyr)
library(catboost)

source("preprocessing/BelgianMTPL_P.R")

load("models/rebal_ratio_belgian_gamma.RData")
load("models/rebal_ratio_belgian_pois.RData")

### Fréquence Poisson avc CatBoost
load("models/catb_pois.RData")

test_pool <- catboost.load_pool(data=test[ , c(-1, -2)],
                                baseline = as.matrix(log(test$Exposure)))

lambda_pois <- catboost.predict(mod.catboost_pois_final, test_pool,
                                        prediction_type = "Exponent") * rebal_ratio_belgianmtpl[5]

### Sévérité Gamma avec NGBoost

py$test_x_oh <- r_to_py(test_x[, -1])
py_run_string("import ngboost
from ngboost import NGBRegressor
from ngboost.distns import Gamma
import numpy as np
from sklearn.tree import DecisionTreeRegressor

file = os.path.join('models',  'ngb_gamma.joblib')
print(file)
model = load(file)
preds_alpha = model.pred_dist(test_x_oh, init_score =  np.zeros((len(test_x_oh), 2))).alpha
preds_beta = model.pred_dist(test_x_oh, init_score =  np.zeros((len(test_x_oh), 2))).beta
y_ngb = model.predict(test_x_oh, init_score =  np.zeros((len(test_x_oh), 2)))")

preds_gamma <- data.frame("alpha" = py$preds_alpha * rebal_ratio_belgianmtpl[6],
                          "beta" = py$preds_beta)


write.csv(data.frame("lambda" = lambda_pois,
           "alpha" = preds_gamma$alpha,
           "beta" = preds_gamma$beta), row.names = FALSE, file = "params_belgian.csv")
