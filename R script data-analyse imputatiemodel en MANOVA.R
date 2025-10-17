# R script for the analysis of OCT and cognition data in first-episode psychosis patients

--------------------------------------------------------------------------------

## Imputation model to deal with missing data

# Packages
library(mice)
library(naniar)
library(dplyr)

# Check datatypes
str(OCT_data)

# Check missing data
colSums(is.na(OCT_data)) 

# Set up the variables with missing data for imputation
impute_vars <- c("MT_scan1", "RNFL_scan1", "CD_ratio_scan1", "cup_volume_scan1")

# Set predictive mean matching (pmm) for the four retinal variables
meth <- rep("", ncol(OCT_data))
names(meth) <- colnames(OCT_data)
meth[impute_vars] <- "pmm"

# Create a predictor matrix for the four retinal variables and the two covariates
pred <- matrix(0, ncol = ncol(OCT_data), nrow = ncol(OCT_data))
colnames(pred) <- colnames(OCT_data)
rownames(pred) <- colnames(OCT_data)

# Puts 1 for the predictors and covariates
predictor_cols <- c(impute_vars, "age_FEP", "sex")
for (r in impute_vars) {
  pred[r, predictor_cols] <- 1
}

# Ensure that a retinal variable does not use itself to predict values for missing data
diag(pred) <- 0

# Create and run imputation model
imp <- mice(OCT_data,
            m = 20,
            maxit = 50,
            seed = 123,
            method = meth,
            predictorMatrix = pred)

# Visualize imputations
plot(imp)

## Check imputation model

# Create a list of the imputations
imputations <- complete(imp, "all")

# Create a new variables containing all variables you want to check
vars <- c("MT_scan1", "RNFL_scan1", "CD_ratio_scan1", "cup_volume_scan1")

# Create an overview of the means and standard deviations
stats <- lapply(vars, function(v) {
  # Observed data
  obs <- OCT_data[[v]]
  obs_mean <- mean(obs, na.rm=TRUE)
  obs_sd <- sd(obs, na.rm=TRUE)
  
  # Imputed data over all imputations
  imp_vals <- unlist(lapply(imputations, function(x) x[[v]]))
  imp_mean <- mean(imp_vals, na.rm=TRUE)
  imp_sd <- sd(imp_vals, na.rm=TRUE)
  
  data.frame(
    variable = v,
    observed_mean = obs_mean,
    observed_sd = obs_sd,
    imputed_mean = imp_mean,
    imputed_sd = imp_sd
  )
})

stats_df <- bind_rows(stats)
print(stats_df)

# Process all imputations in dataset and make it complete
OCT_data_complete <- complete(imp, 1)

# Fit the imputation model
fit_MT <- with(imp, lm(MT_scan1 ~ Group + age_FEP + sex))
fit_RNFL <- with(imp, lm(RNFL_scan1 ~ Group + age_FEP + sex))
fit_cup_volume <- with(imp, lm(cup_volume_scan1 ~ Group + age_FEP + sex))
fit_CD_ratio <- with(imp, lm(CD_ratio_scan1 ~ Group + age_FEP + sex))

# Pool the data of the imputation model
pooled_MT <- pool(fit_MT)
pooled_RNFL <- pool(fit_RNFL)
pooled_cup_volume <- pool(fit_cup_volume)
pooled_CD_ratio <- pool(fit_CD_ratio)

--------------------------------------------------------------------------------

## Permutation MANCOVA

# Packages
library(ggplot2)
library(MVN)
library(vegan)

## Assumption checks
# Assumption of normality - univariate
by(OCT_data_complete$MT_scan1, OCT_data_complete$Group, shapiro.test)
by(OCT_data_complete$RNFL_scan1, OCT_data_complete$Group, shapiro.test)
by(OCT_data_complete$CD_ratio_scan1, OCT_data_complete$Group, shapiro.test)
by(OCT_data_complete$cup_volume_scan1, OCT_data_complete$Group, shapiro.test)

ggplot(OCT_data_complete, aes(sample = MT_scan1)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~Group)

ggplot(OCT_data_complete, aes(sample = RNFL_scan1)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~Group)

ggplot(OCT_data_complete, aes(sample = cup_volume_scan1)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~Group)

ggplot(OCT_data_complete, aes(sample = CD_ratio_scan1)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~Group)

# Assumption of normality - multivariate
# Select only the dependent variables
DVs <- OCT_data_complete[, c("MT_scan1", "RNFL_scan1", "CD_ratio_scan1", "cup_volume_scan1")]
DVs <- colnames(DVs)
df_mvn <- OCT_data_complete[, c(DVs, "Group"), drop = FALSE]

# Check multivariate normality (Mardia)
mvn(data = df_mvn, subset = "Group")

## Permutation MANCOVA

# Create a matrix of the dependent variables
DVs_matrix <- as.matrix(OCT_data_complete[, DVs])

# Run the permutation MANCOVA
set.seed(123)
perm_mancova <- adonis2(DVs_matrix ~ Group + age_FEP + sex,
                        data = OCT_data_complete, 
                        permutations = 999,
                        method = "euclidean")

perm_mancova

## Permutation MANCOVA for the different variables (age, sex, group)

# Define dependent variables
DVs <- c("MT_scan1", "RNFL_scan1", "CD_ratio_scan1", "cup_volume_scan1")
DVs_matrix <- as.matrix(OCT_data_complete[, DVs])

# Define covariates
factors <- c("Group", "age_FEP", "sex")

# Loop over these covariates and run separate permutation MANCOVAs
results_list <- list()
set.seed(123)
for(fac in factors){
  formula <- as.formula(paste("DVs_matrix ~", fac))
  perm_test <- adonis2(
    formula,
    data = OCT_data_complete,
    permutations = 999,
    method = "euclidean")

# Add factor name
perm_test$Factor <- fac
results_list[[fac]] <- perm_test}

# Combine results into one table
results_table <- do.call(rbind, lapply(results_list, function(x){
       data.frame(
          Factor = x$Factor[1],
          Df = x$Df[1],
          SumOfSqs = x$SumOfSqs[1],
          R2 = x$R2[1],
          F = x$F[1],
          p_value = x$`Pr(>F)`[1] )
  }))

 print(results_table)
 
--------------------------------------------------------------------------------

## Linear regression analysis - assumption checks
   
# Packages
library(lmtest)
library(car)

# Filter dataset on Group = 1
OCT_group1 <- subset(OCT_data_complete, Group == 1)

# Model MT_scan1
model_mt <- lm(BACS_scan1 ~ MT_scan1 + age_FEP + sex + ap_scan1 + smoking_status,
               data = OCT_group1)

## Check assumptions
# Assumption of linearity
plot(OCT_group1$MT_scan1, OCT_group1$BACS_scan1,
     xlab = "MT_scan1", ylab = "BACS_scan1")
abline(lm(BACS_scan1 ~ MT_scan1, data = OCT_group1), col = "red")

# Homoscedasticity (visual)
plot(fitted(model_mt), resid(model_mt),
     +      xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")

# Homoscedasticity (formal)
bptest(model_mt)

# Normality of residuals (visual)
qqnorm(resid(model_mt))
qqline(resid(model_mt), col = "red")

# Normality of residuals - Shapiro-Wilk test (formal)
shapiro.test(resid(model_mt))

# Multicollinearity
vif(model_mt)


## Model RNFL_scan1
model_rnfl <- lm(BACS_scan1 ~ RNFL_scan1 + age_FEP + sex + ap_scan1 + smoking_status,
                 data = OCT_group1)

## Check assumptions
# Assumption of linearity
plot(OCT_group1$RNFL_scan1, OCT_group1$BACS_scan1,
     xlab = "RNFL_scan1", ylab = "BACS_scan1")
abline(lm(BACS_scan1 ~ RNFL_scan1, data = OCT_group1), col = "red")

# Homoscedasticity (visual)
plot(fitted(model_rnfl), resid(model_rnfl),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted - RNFL_scan1")
abline(h = 0, col = "red")

# Homoscedasticity (formal)
bptest(model_rnfl)

# Normality residuals (visual)
qqnorm(resid(model_rnfl), main = "QQ-plot - RNFL_scan1")
qqline(resid(model_rnfl), col = "red")

# Normality residuals - Shapiro-Wilk test (formal)
shapiro.test(resid(model_rnfl))

# Multicollinearity
vif(model_rnfl)


## Model CD_ratio_scan1
model_cd <- lm(BACS_scan1 ~ CD_ratio_scan1 + age_FEP + sex + ap_scan1 + smoking_status,
               data = OCT_group1)

## Check assumptions
# Assumption of linearity
plot(OCT_group1$CD_ratio_scan1, OCT_group1$BACS_scan1,
     xlab = "CD_ratio_scan1", ylab = "BACS_scan1")
abline(lm(BACS_scan1 ~ CD_ratio_scan1, data = OCT_group1), col = "red")

# Homoscedasticity (visual)
plot(fitted(model_cd), resid(model_cd),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted - CD_ratio_scan1")
abline(h = 0, col = "red")

# Homoscedasticity (formal)
bptest(model_cd)

# Normality residuals (visual)
qqnorm(resid(model_cd), main = "QQ-plot - CD_ratio_scan1")
qqline(resid(model_cd), col = "red")

# Normality residuals - Shapiro-Wilk test (formal)
shapiro.test(resid(model_cd))

# Multicollinearity
vif(model_cd)


## Model cup_volume_scan1
model_cup <- lm(BACS_scan1 ~ cup_volume_scan1 + age_FEP + sex + ap_scan1 + smoking_status,
                data = OCT_group1)

## Check assumptions
# Assumption of linearity
plot(OCT_group1$cup_volume_scan1, OCT_group1$BACS_scan1,
     xlab = "cup_volume_scan1", ylab = "BACS_scan1")
abline(lm(BACS_scan1 ~ cup_volume_scan1, data = OCT_group1), col = "red")

# Homoscedasticity (visual)
plot(fitted(model_cup), resid(model_cup),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted - cup_volume_scan1")
abline(h = 0, col = "red")

# Homoscedasticity (formal)
bptest(model_cup)

# Normality residuals (visual)
qqnorm(resid(model_cup), main = "QQ-plot - cup_volume_scan1")
qqline(resid(model_cup), col = "red")

# Normality residuals - Shapiro-Wilk test (formal)
shapiro.test(resid(model_cup))

# Multicollinearity
vif(model_cup)

--------------------------------------------------------------------------------

## Linear regression analysis
# Packages
library(lmtest)
library(sandwich)
library(dplyr)
library(broom)

# Creat a variable containing all predictors
predictors <- c("MT_scan1", "RNFL_scan1", "CD_ratio_scan1", "cup_volume_scan1")

# Create an empty list for the results
results_list <- list()

# Loop over each predictor
for (pred in predictors) {
  
  # Formula
  f <- as.formula(paste("BACS_scan1 ~", pred, "+ age_FEP + sex + ap_scan1 + smoking_status"))
  
  # Model
  model <- lm(f, data = OCT_group1)
  
  # Robust coefficients
  coefs <- coeftest(model, vcov = vcovHC(model, type = "HC3"))
  
  # 95%-CI
  ci <- confint(model)
  
  # Select predictor-row 
  row_pred <- coefs[2, ]  # 2nd row  = the predictor
  
  # Combine results into the empty list made previously
  results_list[[pred]] <- data.frame(
    Predictor = pred,
    Estimate = row_pred[1],
    Std_Error = row_pred[2],
    t_value = row_pred[3],
    p_value = row_pred[4],
    CI_lower = ci[2, 1],
    CI_upper = ci[2, 2]
  )
}

# Combine results to create a table
results_table <- bind_rows(results_list)

# Apply FDR-correction
results_table$p_FDR <- p.adjust(results_table$p_value, method = "BH")

# Print results
print(results_table)
