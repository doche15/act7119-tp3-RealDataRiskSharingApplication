#### From point to probabilistic gradient boosting for claim
#### frequency and severity prediction.

### Adapted for the work in ACT-7119

library(fastDummies)
library(dplyr)

dataset_final <- readRDS("mtpl_data.rds")

seed <- 20231103
set.seed(seed)

# train-test: 85 %-15 %
nrow(dataset_final) * 0.15

test_id <- sample(1:nrow(dataset_final), size = 24500, replace = FALSE)
dat_train <- dataset_final[-test_id, ]
dat_test <- dataset_final[test_id, ]

# Remove line index
row.names(dat_train) <- NULL
row.names(dat_test) <- NULL

# Training data
train <- dat_train[, c(-1, -3, -4, -6, -16)]
colnames(train)[2] <- "Total.Claim.Amount" # will only be used to assess tariff structure,
colnames(train)[1] <- "Exposure"
# not for severity modelling (where we would
# use average.)

train$coverage <- as.factor(train$coverage)
train$sex <- as.factor(train$sex)
train$fuel <- as.factor(train$fuel)
train$use <- as.factor(train$use)
train$fleet <- as.factor(train$fleet)

# Test data
test <- dat_test[, c(-1, -3, -4, -6, -16)]
colnames(test)[2] <- "Total.Claim.Amount" # same comment
colnames(test)[1] <- "Exposure"

test$coverage <- as.factor(test$coverage)
test$sex <- as.factor(test$sex)
test$fuel <- as.factor(test$fuel)
test$use <- as.factor(test$use)
test$fleet <- as.factor(test$fleet)

# transformation for implementations requiring one-hot encoding
full_data <- rbind(train, test)
id_train <- 1:nrow(train)

categorielles <- c(3, 5, 9, 10, 11)

tmp <- full_data[, categorielles]
tmp1 <- dummy_cols(tmp, remove_first_dummy = TRUE)
tmp1 <- tmp1[, -(1:length(categorielles))]
d <- data.frame(full_data[, -categorielles], tmp1)
m <- as.matrix(d)
variablereponse <- 1
train_y <- m[id_train, variablereponse]
train_x <- m[id_train, -variablereponse]

test_y <- m[-id_train, variablereponse]
test_x <- m[-id_train, -variablereponse]

noms_variables <- colnames(train_x)