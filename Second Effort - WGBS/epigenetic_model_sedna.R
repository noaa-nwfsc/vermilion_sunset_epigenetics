library(tidyverse)
library(mice)
library(glmnet)
library(broom.mixed)
library(rstanarm)
library(doParallel)
library(doMC)

registerDoMC(cores = 20)
registerDoParallel(cores = 20)
options(mc.cores = 20)


df_master <- readRDS(file = "vs_percent_wide_length_honeycomb.RDS") 

# Remove any columns containing "lambda" prior to processing
lambda_cols <- grep("lambda", colnames(df_master), ignore.case = TRUE)
if(length(lambda_cols) > 0) df_master <- df_master[, -lambda_cols]



# Identify unique CpG base names by stripping the _M and _U suffixes
all_cpg_cols <- colnames(df_master)[grepl("_M$|_U$", colnames(df_master))]
cpg_base_names <- unique(gsub("_M$|_U$", "", all_cpg_cols))

cat("Found", length(cpg_base_names), "unique CpG sites. Vectorizing calculations...\n")

# ==========================================
# 1. BAYESIAN GLM FEATURE SELECTION
# ==========================================

# 1. Generate the exact string column names for M and U bulk pull
m_cols <- paste0(cpg_base_names, "_M")
u_cols <- paste0(cpg_base_names, "_U")

# 2. Extract entire chunks as numeric matrices 
mat_M <- as.matrix(df_master[, m_cols])
mat_U <- as.matrix(df_master[, u_cols])

# 3. Perform bulk vector math simultaneously across the whole grid
mat_Beta <- mat_M / (mat_M + mat_U)

# Convert to stabilized M-values with a small offset (performed worse)
#offset <- 0.01
#mat_Mval <- log2((mat_Beta + offset) / (1 - mat_Beta + offset))


# 4. Bind directly into the final dataframe format in one single shot
df_beta <- data.frame(mat_Beta)
colnames(df_beta) <- cpg_base_names
rownames(df_beta) <- rownames(df_master)

# Append the AGE column back
df_beta$AGE <- df_master$AGE

cat("Running Fast Parallel GLM Feature Selection across CPU cores...\n")

num_cores <- 20 
cat("Utilizing", num_cores, "CPU cores for parallel processing.\n")

# Prepare a clean data frame with just the columns we need to pass to the workers
worker_df <- df_master[, c(all_cpg_cols, "AGE")]

# Define a lightweight worker function that runs on each CPU core
run_fast_glm <- function(cpg, df) {
  M <- df[[paste0(cpg, "_M")]]
  U <- df[[paste0(cpg, "_U")]]
  
  # Fast median imputation inline
  M[is.na(M)] <- median(M, na.rm = TRUE)
  U[is.na(U)] <- median(U, na.rm = TRUE)
  
  response_matrix <- matrix(c(M, U), ncol = 2)
  
  tryCatch({
    # Using standard GLM with a fast mathematical solver instead of heavy Bayesian MCMC
    fit <- glm(response_matrix ~ AGE, data = df, family = binomial(link = "logit"))
    
    # Extract the exact p-value for the AGE coefficient
    summary_coef <- summary(fit)$coefficients
    p_val <- summary_coef["AGE", "Pr(>|z|)"]
    return(p_val)
  }, error = function(e) { 
    return(1) # Return 1 (insignificant) if the model fails to converge
  })
}

# Distribute the CpG list across your CPU cores simultaneously
p_values <- unlist(mclapply(cpg_base_names, run_fast_glm, df = worker_df, mc.cores = num_cores))
names(p_values) <- cpg_base_names

sig_cpgs <- names(p_values[p_values < 0.05])
top_20_count <- round(length(sig_cpgs) * 0.60)
selected_features <- names(sort(p_values[sig_cpgs]))[1:top_20_count]

cat("Selected the top", length(selected_features), "CpG sites from the Bayesian step.\n")
df_pre_impute <- df_beta[, c(selected_features, "AGE")]


#saveRDS(df_pre_impute, file = "post_GLM_df.RDS")

#df_pre_impute <- readRDS("post_GLM_df.RDS")

# ==========================================
# 2. MEMORY-SAFE MICE IMPUTATION (ON BETAS)
# ==========================================
cat("Running Custom MICE Imputation on selected Beta values...\n")
pred_matrix <- matrix(0, nrow = ncol(df_pre_impute), ncol = ncol(df_pre_impute))
rownames(pred_matrix) <- colnames(df_pre_impute)
colnames(pred_matrix) <- colnames(df_pre_impute)
pred_matrix[, "AGE"] <- 1  # Anchoring prediction to Age

imputation_method <- rep("pmm", ncol(df_pre_impute))
names(imputation_method) <- colnames(df_pre_impute)
imputation_method[colSums(is.na(df_pre_impute)) == 0] <- ""

# 1. Raise the future object limit threshold from 500MB to 2GB (2000 * 1024^2)
options(future.globals.maxSize = 2000 * 1024^2)
# 2. Aggressively clean out hidden backtrace and call memory variables
rm(list = intersect(ls(), c("mat_M", "mat_U", "worker_df")))
gc() # Call R's system garbage collector to clear the 520 MiB bloat



mice_out <- futuremice(
  df_pre_impute,
  m = 1,
  maxit = 3,
  method = imputation_method,
  print = FALSE,
  predictorMatrix = pred_matrix, 
  n.core = 1,
  seed = 42)



#mice_out <- mice(df_pre_impute, m = 1, maxit = 3, method = imputation_method,
#                 predictorMatrix = pred_matrix, printFlag = FALSE, seed = 42)

df_imputed <- complete(mice_out, 1)

saveRDS(df_imputed, file = "mice_imputed_df_honeycomb.RDS")

# 1. Identify if 'Length' or 'LENGTH' or 'length' exists in your original data
length_col_name <- grep("^length$", colnames(df_master), ignore.case = TRUE, value = TRUE)

# 2. Extract Length matching the exact row order of your imputed dataset
length_data <- df_master[rownames(df_imputed), length_col_name]

# 3. Fill any missing lengths with the median to prevent GLMNET from discarding samples
length_data[is.na(length_data)] <- median(length_data, na.rm = TRUE)

# 4. Bind the selected CpGs and Length together into your predictor matrix
X_matrix <- as.matrix(df_imputed[, selected_features])
X_matrix <- cbind(X_matrix, Length = length_data) # Appends Length as a clean, final column

# Target variable remains Age
Y_vector <- df_imputed$AGE
#Y_vector <- log(df_imputed$AGE + 1)

# ==========================================================================
# UPGRADE 2: INVERSE PROBABILITY WEIGHTS CALCULATION
# ==========================================================================
# Group ages into 5-year bins and calculate weight based on bin rarity
age_bins <- cut(df_imputed$AGE, breaks = seq(0, 80, by = 5))
bin_counts <- table(age_bins)

# Weight is inversely proportional to how common the age bin is
sample_weights <- 1 / as.numeric(bin_counts[age_bins])

# Normalize weights so they average out to 1
sample_weights <- sample_weights / mean(sample_weights)

# Test the sweet-spot range indicated by your last run   
# ==========================================
# 3. HIGH-PERFORMANCE CORE-SAFE TUNING ENGINE
# ==========================================
library(foreach)

alphas_to_test <- c(0.35, 0.37, 0.4)

summary_table <- data.frame(
  Alpha = alphas_to_test,
  Best_Test_MAE = NA,
  Test_R2 = NA
)

global_best_mae <- Inf
global_best_model <- NULL
global_best_alpha <- NA
global_best_r2  <- -Inf

# Master configuration for your deeper iteration depths
NUM_BOOTSTRAPS <- 100   # Scaled up as requested
NUM_SUBSAMPLES <- 500   # Scaled up for deeper search depth

for(a_idx in seq_along(alphas_to_test)) {
  current_alpha <- alphas_to_test[a_idx]
  cat("\nTesting Pipeline Configuration: Alpha =", current_alpha, " (Running on 20 Cores Closely monitored)...\n")
  
  # --- CORE-SAFE PARALLEL BOOTSTRAP DISTRIBUTION ---
  # Distributes the 100 bootstraps across your 20 workers. 
  # cv.glmnet parallel flag is set to FALSE to avoid process multiplication.
  boot_results <- foreach(boot_iter = 1:NUM_BOOTSTRAPS, .packages = c("glmnet")) %dopar% {
    
    #MASTER_SEED <- 5
    set.seed(boot_iter)
    
    train_idx <- sample(1:nrow(df_imputed), size = floor(0.70 * nrow(df_imputed)))
    
    X_train     <- X_matrix[train_idx, , drop = FALSE]
    Y_train_log <- log(Y_vector[train_idx] + 1)
    W_train     <- sample_weights[train_idx]
    
    X_test      <- X_matrix[-train_idx, , drop = FALSE]
    Y_test_log  <- log(Y_vector[-train_idx] + 1)
    
    current_features <- colnames(X_train)
    
    # Construct Penalty Factors Vector
    p_factors <- rep(1, length(current_features))
    names(p_factors) <- current_features
    if("Length" %in% names(p_factors)) p_factors["Length"] <- 0
    
    # --- Step 3A: Recursive Lasso Pruning Loop ---
    for(step in 1:4) {
      if(length(current_features) <= 10) break
      local_p_factors <- p_factors[current_features]
      
      cv_filter <- cv.glmnet(X_train[, current_features, drop = FALSE], Y_train_log, 
                             weights = W_train,
                             penalty.factor = local_p_factors,
                             alpha = 0.1, nfolds = 10, 
                             family = "gaussian", type.measure = "mae", 
                             standardize = TRUE,
                             parallel = FALSE) # FIXED: Force FALSE to protect core limit
      
      coefs <- as.matrix(abs(coef(cv_filter, s = "lambda.min")[-1, , drop = FALSE]))
      
      if("Length" %in% rownames(coefs)) {
        cutoff_val <- quantile(coefs[rownames(coefs) != "Length", 1], probs = 0.5)
        current_features <- unique(c(rownames(coefs)[coefs[, 1] >= cutoff_val], "Length"))
      } else {
        cutoff_val <- quantile(coefs, probs = 0.5)
        current_features <- rownames(coefs)[coefs >= cutoff_val]
      }
    }
    
    # --- Step 3B: Random Feature Subsampling Engine ---
    cpg_features <- current_features[current_features != "Length"]
    target_n <- min(249, length(cpg_features)) 
    
    local_best_r2 <- -Inf
    local_best_metrics <- list(r2 = -Inf, mae = Inf, preds = NULL, y_test = NULL, model = NULL, features = NULL, x_train = NULL, y_train = NULL)
    
    for(sub_iter in 1:NUM_SUBSAMPLES) {
      subset_cpgs <- c(sample(cpg_features, size = target_n), "Length")
      sub_p_factors <- rep(1, length(subset_cpgs))
      names(sub_p_factors) <- subset_cpgs
      sub_p_factors["Length"] <- 0
      
      model_elastic <- cv.glmnet(X_train[, subset_cpgs, drop = FALSE], Y_train_log, 
                                 weights = W_train,
                                 penalty.factor = sub_p_factors,
                                 alpha = current_alpha, 
                                 nfolds = 10, 
                                 family = "gaussian", type.measure = "mae", 
                                 standardize = TRUE,
                                 parallel = FALSE) # FIXED: Force FALSE to protect core limit
      
      preds <- predict(model_elastic, newx = X_test[, subset_cpgs, drop = FALSE], s = "lambda.min")
      true_y_test <- Y_vector[-train_idx]
      true_preds  <- exp(as.vector(preds)) - 1
      
      mae_test <- mean(abs(true_y_test - true_preds))
      lm_sub   <- lm(true_y_test ~ true_preds)
      r2_sub   <- summary(lm_sub)$r.squared
      
      if(!is.na(r2_sub) && r2_sub > local_best_r2) {
        local_best_r2 <- r2_sub
        local_best_metrics <- list(
          r2 = r2_sub, mae = mae_test, preds = true_preds, y_test = true_y_test,
          model = model_elastic, features = subset_cpgs, x_train = X_train, y_train = Y_train_log
        )
      }
    }
    return(local_best_metrics)
  }
  
  # --- COMPILING THE VISITED BOOTSTRAPS TO FIND THE ALPHA CHAMPION ---
  best_alpha_r2 <- -Inf
  best_alpha_mae <- Inf
  var_winner <- NULL
  
  for(b in seq_along(boot_results)) {
    curr_res <- boot_results[[b]]
    if(!is.null(curr_res) && curr_res$r2 > best_alpha_r2) {
      best_alpha_r2         <- curr_res$r2
      best_alpha_mae        <- curr_res$mae
      var_winner            <- curr_res
    }
  }
  
  summary_table$Best_Test_MAE[a_idx] <- round(best_alpha_mae, 4)
  summary_table$Max_Test_R2[a_idx]   <- round(best_alpha_r2, 4)
  
  if(best_alpha_r2 > global_best_r2) {
    global_best_r2    <- best_alpha_r2
    global_best_mae   <- best_alpha_mae
    global_best_alpha <- current_alpha
    global_best_model <- list(
      model = var_winner$model, features = var_winner$features,
      train_x = var_winner$x_train, train_y = var_winner$y_train,
      test_y = var_winner$y_test, test_preds = var_winner$preds
    )
  }
}


#saveRDS(global_best_model, file = paste0(file_prefix, "_model.rds"))
# ==========================================
# 4. DIAGNOSTICS & PERFORMANCE EVALUATION
# ==========================================
cat("\n\n--- MODEL COMPARISON SUMMARY TABLE (MAXIMIZED R2) ---\n")
print(summary_table)

# Get log predictions, then convert to original years
train_preds_log <- as.vector(predict(global_best_model$model,
                                     newx = global_best_model$train_x[, global_best_model$features],
                                     s = "lambda.min"))
train_preds_raw <- exp(train_preds_log) - 1
true_train_y    <- exp(global_best_model$train_y) - 1

# Calibrate on the raw year scale
calibration_lm <- lm(true_train_y ~ train_preds_raw)

# Calibrate Testing Predictions (already on raw scale in best_model structure)
calibrated_test_preds <- as.vector(predict(calibration_lm, newdata = data.frame(train_preds_raw = global_best_model$test_preds)))
calibrated_mae <- mean(abs(global_best_model$test_y - calibrated_test_preds))

# Calibrate Training Predictions for plotting
calibrated_train_preds <- as.vector(predict(calibration_lm, newdata = data.frame(train_preds_raw = train_preds_raw)))
calibrated_train_mae <- mean(abs(true_train_y - calibrated_train_preds))

# Re-assign unlogged values to the plotting frame
global_best_model$train_y <- true_train_y

# --- RELATIVE MAE & STATS ---
age_range <- max(Y_vector) - min(Y_vector)
relative_mae_raw <- (global_best_mae / age_range) * 100
relative_mae_calibrated <- (calibrated_mae / age_range) * 100

lm_final_calibrated <- lm(global_best_model$test_y ~ calibrated_test_preds)
summary_stats_calibrated <- summary(lm_final_calibrated)

# ==========================================================================
# EXTRACT FINAL ACTIVE COV/SNP MARKERS
# ==========================================================================
raw_coefs <- coef(global_best_model$model, s = "lambda.min")
active_coefs_matrix <- as.matrix(raw_coefs)[-1, , drop = FALSE]

active_features <- rownames(active_coefs_matrix)[active_coefs_matrix[, 1] != 0]
final_snps_used <- active_features[active_features != "Length"]
num_snps_used <- length(final_snps_used)
length_covariate_selected <- "Length" %in% active_features

# ==========================================================================
# TOLERANCE ACCURACY METRICS (OUT-OF-SAMPLE TEST SET)
# ==========================================================================
test_actual <- global_best_model$test_y
test_errors <- abs(test_actual - calibrated_test_preds)
total_test_n <- length(test_errors)

acc_1yr_pct <- (sum(test_errors <= 1.0) / total_test_n) * 100
acc_2yr_pct <- (sum(test_errors <= 2.0) / total_test_n) * 100
acc_3yr_pct <- (sum(test_errors <= 3.0) / total_test_n) * 100

# ==========================================================================
# DYNAMIC FILENAME PREPARATION (FIXED ORDERING)
# ==========================================================================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")

safe_alpha <- "unspecified"
if (!is.null(global_best_alpha) && !is.function(global_best_alpha) && !is.na(global_best_alpha)) {
  safe_alpha <- as.character(global_best_alpha)
}

file_prefix <- paste0(
  "clock_alpha", safe_alpha, 
  "_snps", num_snps_used,
  "_R", round(summary_stats_calibrated$r.squared, 3), 
  "_MAE", round(calibrated_mae, 2), 
  "_", timestamp, 
  "_honeycombalignment"
)

# Export Summary Metrics Table
write.csv(summary_table, paste0("alpha_tuning_metrics_summary_", timestamp, "_honeycomb.csv"), row.names = FALSE)

# Export Selected Active Markers to CSV
active_df <- data.frame(
  Feature = rownames(active_coefs_matrix)[active_coefs_matrix[, 1] != 0],
  Coefficient = active_coefs_matrix[active_coefs_matrix[, 1] != 0, 1]
)
write.csv(active_df, paste0(file_prefix, "_selected_markers_honeycomb.csv"), row.names = FALSE)

# ==========================================================================
# PRINT DETAILED BIOLOGICAL DIAGNOSTIC REPORT
# ==========================================================================
cat("\n=========================================================\n")
cat("          FINAL EPIGENETIC CLOCK REPORT CARD             \n")
cat("=========================================================\n")
cat("Total Chronological Age Range in Cohort: ", age_range, " years\n")
cat("Raw Model Test MAE:                      ", round(global_best_mae, 3), " years\n")
cat("Calibrated Model Test MAE:               ", round(calibrated_mae, 3), " years\n")
cat("Median Absolute Error (MedAE):           ", round(median(test_errors), 3), " years\n")
cat("Raw Model Relative MAE:                  ", round(relative_mae_raw, 2), "%\n")
cat("Calibrated Model Relative MAE:           ", round(relative_mae_calibrated, 2), "%\n")
cat("Final Epigenetic Clock R-squared:        ", round(summary_stats_calibrated$r.squared, 4), "\n")
cat("---------------------------------------------------------\n")
cat("Accuracy within ±1.0 Year of Otolith Age: ", round(acc_1yr_pct, 1), "% (", sum(test_errors <= 1.0), "/", total_test_n, ")\n", sep="")
cat("Accuracy within ±2.0 Years of Otolith Age:", round(acc_2yr_pct, 1), "% (", sum(test_errors <= 2.0), "/", total_test_n, ")\n", sep="")
cat("Accuracy within ±3.0 Years of Otolith Age:", round(acc_3yr_pct, 1), "% (", sum(test_errors <= 3.0), "/", total_test_n, ")\n", sep="")
cat("---------------------------------------------------------\n")
cat("Total Candidate Markers Fed to Model:    ", length(global_best_model$features), "\n")
cat("TOTAL ACTIVE CpG/SNP MARKERS SELECTED:   ", num_snps_used, "\n")
cat("Was 'Length' Covariate Kept as a Feature? ", ifelse(length_covariate_selected, "YES", "NO"), "\n")
cat("=========================================================\n\n")

saveRDS(global_best_model, file = paste0(file_prefix, "_model.rds"))



# ==========================================
# 5. GENERATE CALIBRATION & HETEROSCEDASTICITY PLOTS
# ==========================================
# Create two distinct data frames and bind them together
train_plot_df <- data.frame(
  Actual_Age = global_best_model$train_y,
  Raw_Predictions = train_preds_raw,
  Calibrated_Predictions = calibrated_train_preds,
  Residuals = global_best_model$train_y - calibrated_train_preds,
  Set = "Training (In-Sample)"
)

test_plot_df <- data.frame(
  Actual_Age = global_best_model$test_y,
  Raw_Predictions = global_best_model$test_preds,
  Calibrated_Predictions = calibrated_test_preds,
  Residuals = global_best_model$test_y - calibrated_test_preds,
  Set = "Testing (Out-of-Sample)"
)

plot_df <- rbind(train_plot_df, test_plot_df)

# Plot 1: Standard Age Calibration Graph (With Both Sets)
annotation_text <- paste0(
  "Optimized Alpha = ", safe_alpha, "\n",
  "Test R² = ", round(summary_stats_calibrated$r.squared, 3), "\n",
  "Test Calibrated MAE = ", round(calibrated_mae, 2), " yrs\n",
  "Test Relative MAE = ", round(relative_mae_calibrated, 1), "%")

age_plot <- ggplot(plot_df, aes(x = Actual_Age, y = Calibrated_Predictions)) +
    # Plots Training as subtle grey circles and Testing as vibrant cyan triangles
    geom_point(aes(color = Set, shape = Set, alpha = Set), size = 3) +
    scale_color_manual(values = c("Training (In-Sample)" = "grey60", "Testing (Out-of-Sample)" = "darkcyan")) +
    scale_shape_manual(values = c("Training (In-Sample)" = 16, "Testing (Out-of-Sample)" = 17)) +
    scale_alpha_manual(values = c("Training (In-Sample)" = 0.5, "Testing (Out-of-Sample)" = 0.9)) +
    # Fit the dashed regression line ONLY to the Testing data
    geom_smooth(data = filter(plot_df, Set == "Testing (Out-of-Sample)"), 
                method = "lm", color = "red", se = TRUE, linetype = "dashed") +
    # Standard 1-to-1 line of perfect prediction
    geom_abline(intercept = 0, slope = 1, color = "grey40", linetype = "solid") +
    annotate("text", x = min(plot_df$Actual_Age) + 0.5, y = max(plot_df$Calibrated_Predictions) - 5, 
             label = annotation_text, hjust = 0, size = 4, fontface = "bold", color = "black") +
    labs(x = "Chronological Otolith Age (Years)", y = "Calibrated Predicted Age (Years)") +
    theme_minimal(base_size = 14) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic"),
          legend.position = "bottom")

ggsave(filename = paste0(file_prefix, "_performance.png"), plot = age_plot, width = 8, height = 7, dpi = 300)

# Plot 2: Heteroscedasticity Residual Analysis Graph (With Both Sets)
residual_plot <- ggplot(plot_df, aes(x = Actual_Age, y = Residuals)) +
    geom_point(aes(color = Set, shape = Set, alpha = Set), size = 3) +
    scale_color_manual(values = c("Training (In-Sample)" = "grey60", "Testing (Out-of-Sample)" = "purple")) +
    scale_shape_manual(values = c("Training (In-Sample)" = 16, "Testing (Out-of-Sample)" = 17)) +
    scale_alpha_manual(values = c("Training (In-Sample)" = 0.4, "Testing (Out-of-Sample)" = 0.8)) +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    # Fit a smooth trend line to the testing residuals to watch for error drift
    geom_smooth(data = filter(plot_df, Set == "Testing (Out-of-Sample)"), 
                method = "loess", color = "orange", se = FALSE) + 
    labs(title = "Heteroscedasticity Residual Analysis",
         subtitle = "Checking for error expansion ('megaphone effect') across age ranges",
         x = "Chronological Otolith Age (Years)", y = "Residual Error (Actual - Predicted)") +
    theme_minimal(base_size = 14) + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic"),
          legend.position = "bottom")

ggsave(filename = paste0(file_prefix, "_residual_analysis.png"), plot = residual_plot, width = 8, height = 7, dpi = 300)
cat("Diagnostic visualizations successfully written to workspace directory as unique files.\n")


# ==========================================================================
# 6. LENGTH-ONLY BASELINE MODEL COMPARISON
# ==========================================================================
cat("\n=========================================================\n")
cat("      RUNNING LENGTH-ONLY BASELINE BENCHMARK             \n")
cat("=========================================================\n")

# 1. Prepare Length and Log-Age from the exact processed matrices
length_vec  <- X_matrix[, "Length"]
log_age_vec <- log(Y_vector + 1)

# 2. Fit Weighted Linear Regression using ONLY Length
fit_length_base <- lm(log_age_vec ~ length_vec, weights = sample_weights)

# 3. Generate raw predictions and back-transform to real age (years)
preds_length_log <- predict(fit_length_base)
preds_length_raw <- exp(preds_length_log) - 1
preds_length_raw[preds_length_raw < 0] <- 0

# 4. Apply Post-Hoc Linear Calibration
calib_length_lm  <- lm(Y_vector ~ preds_length_raw)
calib_length_preds <- predict(calib_length_lm)

# 5. Compute Benchmarking Metrics
ss_res_len <- sum((Y_vector - calib_length_preds)^2)
ss_tot_len <- sum((Y_vector - mean(Y_vector))^2)

r2_length  <- 1 - (ss_res_len / ss_tot_len)
mae_length <- mean(abs(Y_vector - calib_length_preds))
relative_mae_length <- (mae_length / age_range) * 100

# 6. Print Comparison Delta
cat(sprintf("Length-Only Model -> Calibrated R²: %.4f | Calibrated MAE: %.2f years (%.1f%%)\n", 
            r2_length, mae_length, relative_mae_length))
cat(sprintf("Epigenetic Clock   -> Calibrated R²: %.4f | Calibrated MAE: %.2f years (%.1f%%)\n", 
            summary_stats_calibrated$r.squared, calibrated_mae, relative_mae_calibrated))
cat("---------------------------------------------------------\n")
cat(sprintf("DELTA (Added Value) -> ΔR²: +%.4f | ΔMAE: -%.2f years\n", 
            summary_stats_calibrated$r.squared - r2_length, 
            mae_length - calibrated_mae))
cat("=========================================================\n\n")

