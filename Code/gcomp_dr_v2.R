
#Additional functions to organize output from g-comp and DR
############################################################
## Helper function:
## Construct Direct, Indirect, Total, and Overall effects
## from the estimated average potential outcomes (APOs).
##
## Input:
##   apo: data frame with columns alpha, estimand, estimate
##   y1_name: label for \bar y(1, alpha)
##   y0_name: label for \bar y(0, alpha)
##   ym_name: label for \bar y(alpha)
##
## Output:
##   data frame with columns:
##   estimation, alpha0, alpha1, type
############################################################
make_effects_from_apo <- function(apo,
                                  y1_name,
                                  y0_name,
                                  ym_name) {
  
  ## Sort alphas from largest to smallest so output matches the desired order
  alpha_vals <- sort(unique(apo$alpha), decreasing = TRUE)
  
  ## -----------------------------
  ## Direct effects
  ## DE(alpha) = ybar(1, alpha) - ybar(0, alpha)
  ## For direct effects, alpha0 = alpha1 = alpha
  ## -----------------------------
  de <- do.call(rbind, lapply(alpha_vals, function(a) {
    y1 <- apo$estimate[apo$alpha == a & apo$estimand == y1_name]
    y0 <- apo$estimate[apo$alpha == a & apo$estimand == y0_name]
    
    data.frame(
      estimation = y1 - y0,
      alpha0 = a,
      alpha1 = a,
      type = "Direct"
    )
  }))
  
  ## Create all pairwise alpha comparisons:
  ## (0.75, 0.50), (0.75, 0.25), (0.50, 0.25), etc.
  pairs <- combn(alpha_vals, 2, simplify = FALSE)
  
  ## -----------------------------
  ## Indirect effects
  ## IE(alpha0, alpha1) = ybar(0, alpha0) - ybar(0, alpha1)
  ## -----------------------------
  ie <- do.call(rbind, lapply(pairs, function(pr) {
    a0 <- pr[1]
    a1 <- pr[2]
    
    y0_a0 <- apo$estimate[apo$alpha == a0 & apo$estimand == y0_name]
    y0_a1 <- apo$estimate[apo$alpha == a1 & apo$estimand == y0_name]
    
    data.frame(
      estimation = y0_a0 - y0_a1,
      alpha0 = a0,
      alpha1 = a1,
      type = "Indirect"
    )
  }))
  
  ## -----------------------------
  ## Total effects
  ## TE(alpha0, alpha1) = ybar(1, alpha0) - ybar(0, alpha1)
  ## -----------------------------
  te <- do.call(rbind, lapply(pairs, function(pr) {
    a0 <- pr[1]
    a1 <- pr[2]
    
    y1_a0 <- apo$estimate[apo$alpha == a0 & apo$estimand == y1_name]
    y0_a1 <- apo$estimate[apo$alpha == a1 & apo$estimand == y0_name]
    
    data.frame(
      estimation = y1_a0 - y0_a1,
      alpha0 = a0,
      alpha1 = a1,
      type = "Total"
    )
  }))
  
  ## -----------------------------
  ## Overall effects
  ## OE(alpha0, alpha1) = ybar(alpha0) - ybar(alpha1)
  ## -----------------------------
  oe <- do.call(rbind, lapply(pairs, function(pr) {
    a0 <- pr[1]
    a1 <- pr[2]
    
    ym_a0 <- apo$estimate[apo$alpha == a0 & apo$estimand == ym_name]
    ym_a1 <- apo$estimate[apo$alpha == a1 & apo$estimand == ym_name]
    
    data.frame(
      estimation = ym_a0 - ym_a1,
      alpha0 = a0,
      alpha1 = a1,
      type = "Overall"
    )
  }))
  
  ## Stack the four effect types together in the desired conceptual order
  rbind(de, ie, te, oe)
}


############################################################
## Helper function:
## Reorder the final results table to match the screenshot:
##
## Direct
## Var DE
## Indirect
## Var IE
## Total
## Var TE
## Overall
## Var OE
##
## Within each type, sort by descending alpha0 then alpha1.
############################################################
reorder_effect_table <- function(df) {
  type_order <- c("Direct", "Var DE",
                  "Indirect", "Var IE",
                  "Total", "Var TE",
                  "Overall", "Var OE")
  
  type_rank <- match(df$type, type_order)
  
  ord <- order(type_rank, -df$alpha0, -df$alpha1)
  df[ord, , drop = FALSE]
}

############################################################
## G-computation estimator with bootstrap confidence intervals
##
## What this function does:
##   1. Fit an outcome regression model
##   2. Estimate ybar(1, alpha), ybar(0, alpha), ybar(alpha)
##   3. Convert Averge Potential Outcomes (APOs) into DE/IE/TE/OE
##   4. Bootstrap by resampling individuals with replacement,
##      carrying along each sampled individual's summary
##      neighbor information
##   5. Return point estimates, bootstrap variances,
##      and percentile confidence intervals
############################################################
gcomp_with_bootstrap <- function(data,
                                 alpha = c(0.75, 0.50, 0.25),
                                 z_covariates = c("var1", "var2"),
                                 zn_covariates = c("avg_var1", "avg_var2"),
                                 outcome_name = "outcome",
                                 treatment_name = "treatment",
                                 s_name = "na_a",
                                 d_name = "na",
                                 B = 500,
                                 ci_level = 0.95,
                                 seed = 123) {
  
  ##########################################################
  ## Inner function:
  ## Fit the g-computation estimator one time on a dataset
  ## and return the DE/IE/TE/OE effect table.
  ##########################################################
  fit_once <- function(dat) {
    
    ## Build the outcome model formula:
    ## outcome ~ A + S + A:S + Z_i + Z_Ni
    q_terms <- c(treatment_name,
                 s_name,
                 paste0(treatment_name, ":", s_name),
                 z_covariates,
                 zn_covariates)
    
    q_formula <- as.formula(
      paste(outcome_name, "~", paste(q_terms, collapse = " + "))
    )
    
    ## Fit a logistic regression outcome model
    ## Assumes binary outcome
    Qfit <- glm(q_formula, family = binomial(), data = dat)
    
    ## Number of individuals in the dataset
    n <- nrow(dat)
    
    ########################################################
    ## Estimate ybar(a, alpha):
    ## Average over individuals of the model-based predicted
    ## outcome under fixed a and stochastic S_i ~ Binomial(d_i, alpha)
    ########################################################
    est_ybar_a_alpha <- function(a, alpha_val) {
      mu_i <- numeric(n)
      
      for (i in seq_len(n)) {
        ## Degree for individual i
        d_i <- dat[[d_name]][i]
        
        ## All possible treated-neighbor counts for unit i
        s_vals <- 0:d_i
        
        ## Policy probability pi(s_i; alpha)
        p_s <- dbinom(s_vals, size = d_i, prob = alpha_val)
        
        ## Build prediction data holding covariates fixed at unit i's values
        ## and setting A_i = a, S_i = s over all possible s values
        newdat <- data.frame(tmpA = rep(a, length(s_vals)),
                             tmpS = s_vals)
        names(newdat) <- c(treatment_name, s_name)
        
        for (v in z_covariates) newdat[[v]] <- dat[[v]][i]
        for (v in zn_covariates) newdat[[v]] <- dat[[v]][i]
        
        ## Predict E(Y_i | A_i=a, S_i=s, Z_i, Z_Ni)
        q_pred <- predict(Qfit, newdata = newdat, type = "response")
        
        ## Average predicted outcomes over S_i under policy alpha
        mu_i[i] <- sum(q_pred * p_s)
      }
      
      ## Average over individuals
      mean(mu_i)
    }
    
    ########################################################
    ## Estimate ybar(alpha):
    ## Average over both A_i and S_i under the joint policy:
    ##   A_i ~ Bernoulli(alpha)
    ##   S_i ~ Binomial(d_i, alpha)
    ########################################################
    est_ybar_alpha <- function(alpha_val) {
      mu_i <- numeric(n)
      
      for (i in seq_len(n)) {
        d_i <- dat[[d_name]][i]
        s_vals <- 0:d_i
        p_s <- dbinom(s_vals, size = d_i, prob = alpha_val)
        
        ## Prediction data for A_i = 0 across all s
        newdat0 <- data.frame(tmpA = rep(0, length(s_vals)),
                              tmpS = s_vals)
        
        ## Prediction data for A_i = 1 across all s
        newdat1 <- data.frame(tmpA = rep(1, length(s_vals)),
                              tmpS = s_vals)
        
        names(newdat0) <- c(treatment_name, s_name)
        names(newdat1) <- c(treatment_name, s_name)
        
        ## Hold covariates fixed at the observed values for individual i
        for (v in z_covariates) {
          newdat0[[v]] <- dat[[v]][i]
          newdat1[[v]] <- dat[[v]][i]
        }
        for (v in zn_covariates) {
          newdat0[[v]] <- dat[[v]][i]
          newdat1[[v]] <- dat[[v]][i]
        }
        
        ## Predicted outcomes under A=0 and A=1
        q0 <- predict(Qfit, newdata = newdat0, type = "response")
        q1 <- predict(Qfit, newdata = newdat1, type = "response")
        
        ## Average over A_i and S_i under the policy alpha
        mu_i[i] <- sum(((1 - alpha_val) * q0 + alpha_val * q1) * p_s)
      }
      
      mean(mu_i)
    }
    
    ## Compute APOs for each alpha value
    apo <- do.call(rbind, lapply(alpha, function(a0) {
      data.frame(
        alpha = a0,
        estimand = c("ybar(1,alpha)", "ybar(0,alpha)", "ybar(alpha)"),
        estimate = c(est_ybar_a_alpha(1, a0),
                     est_ybar_a_alpha(0, a0),
                     est_ybar_alpha(a0))
      )
    }))
    
    ## Convert APOs to causal effects
    make_effects_from_apo(
      apo = apo,
      y1_name = "ybar(1,alpha)",
      y0_name = "ybar(0,alpha)",
      ym_name = "ybar(alpha)"
    )
  }
  
  ## Set seed for reproducibility
  set.seed(seed)
  
  ## Compute point estimates from the original data
  point_eff <- fit_once(data)
  
  ## Number of individuals and number of effect rows
  n <- nrow(data)
  n_eff <- nrow(point_eff)
  
  ## Storage matrix for bootstrap draws:
  ## rows = bootstrap samples, columns = effect rows
  boot_mat <- matrix(NA_real_, nrow = B, ncol = n_eff)
  
  ## ----------------------------------------
  ## Bootstrap loop
  ## Resample rows with replacement.
  ## Each sampled row carries along its observed
  ## summary neighbor information.
  ## ----------------------------------------
  # for (b in seq_len(B)) {
  #   idx <- sample.int(n, size = n, replace = TRUE)
  #   d_b <- data[idx, , drop = FALSE]
  #   
  #   fit_b <- try(fit_once(d_b), silent = TRUE)
  #   
  #   if (!inherits(fit_b, "try-error")) {
  #     boot_mat[b, ] <- fit_b$estimation
  #   }
  # }
  
  b <- 1
  attempt <- 0
  max_attempts <- B * 20
  
  while (b <= B && attempt < max_attempts) {
    attempt <- attempt + 1
    
    idx <- sample.int(n, size = n, replace = TRUE)
    d_b <- data[idx, , drop = FALSE]
    
    warned <- FALSE
    
    fit_b <- tryCatch(
      {
        withCallingHandlers(
          fit_once(d_b),
          warning = function(w) {
            if (grepl("fitted probabilities numerically 0 or 1 occurred",
                      conditionMessage(w))) {
              warned <<- TRUE
              invokeRestart("muffleWarning")
            }
          }
        )
      },
      error = function(e) NULL
    )
    
    if (!warned && !is.null(fit_b)) {
      boot_mat[b, ] <- fit_b$estimation
      b <- b + 1
    }
  }
  
  ## Compute percentile CI bounds
  alpha_tail <- (1 - ci_level) / 2
  lo <- alpha_tail
  hi <- 1 - alpha_tail
  
  ## Bootstrap variances and percentile confidence limits
  var_vec <- apply(boot_mat, 2, var, na.rm = TRUE)
  lci_vec <- apply(boot_mat, 2, quantile, probs = lo, na.rm = TRUE)
  uci_vec <- apply(boot_mat, 2, quantile, probs = hi, na.rm = TRUE)
  
  ## Create variance rows corresponding to each point-estimate row
  var_rows <- point_eff
  var_rows$estimation <- var_vec
  
  ## Rename effect type to variance label
  var_rows$type[point_eff$type == "Direct"] <- "Var DE"
  var_rows$type[point_eff$type == "Indirect"] <- "Var IE"
  var_rows$type[point_eff$type == "Total"] <- "Var TE"
  var_rows$type[point_eff$type == "Overall"] <- "Var OE"
  
  ## Main rows get CI columns
  out_main <- point_eff
  out_main$lower <- lci_vec
  out_main$upper <- uci_vec
  
  ## Variance rows do not get CI columns
  out_var <- var_rows
  out_var$lower <- NA_real_
  out_var$upper <- NA_real_
  
  ## Stack estimates and variances together
  out <- rbind(out_main, out_var)
  
  ## Keep only requested columns
  out <- out[, c("estimation", "alpha0", "alpha1", "type", "lower", "upper")]
  
  ## Reorder rows to match desired table layout
  out <- reorder_effect_table(out)
  
  ## Clean row names
  rownames(out) <- NULL
  
  list(
    results = out,
    point_estimates = point_eff,
    bootstrap_draws = boot_mat
  )
}

## Later: Add  M-estimation (use geex)
############################################################
## Augmented IPW2 estimator with bootstrap (douby-robust)
##
## What this function does:
##   1. Fit:
##      - outcome model Q(A,S,Z,Z_N)
##      - individual treatment model g_A(A | Z)
##      - neighbor treatment model g_S(S | A,Z,Z_N)
##   2. Compute the DR APOs
##   3. Convert APOs into DE/IE/TE/OE
##   4. Bootstrap by resampling individuals with replacement
##   5. Return point estimates, bootstrap variances,
##      and percentile confidence intervals
############################################################
dr_with_bootstrap <- function(data,
                              alpha = c(0.75, 0.50, 0.25),
                              z_covariates = c("var1", "var2"),
                              zn_covariates = c("avg_var1", "avg_var2"),
                              outcome_name = "outcome",
                              treatment_name = "treatment",
                              s_name = "na_a",
                              d_name = "na",
                              nots_name = "notna_a",
                              weight_bound = 1e-6,
                              B = 500,
                              ci_level = 0.95,
                              seed = 123) {
  
  ##########################################################
  ## Inner function:
  ## Fit the DR estimator once on a dataset
  ## and return the DE/IE/TE/OE effect table.
  ##########################################################
  fit_once_dr <- function(dat) {
    
    ## ----------------------------------------
    ## Outcome model Q(A,S,Z,Z_N)
    ## ----------------------------------------
    q_terms <- c(treatment_name,
                 s_name,
                 paste0(treatment_name, ":", s_name),
                 z_covariates,
                 zn_covariates)
    
    q_formula <- as.formula(
      paste(outcome_name, "~", paste(q_terms, collapse = " + "))
    )
    Qfit <- glm(q_formula, family = binomial(), data = dat)
    
    ## ----------------------------------------
    ## Individual treatment model g_A(A | Z)
    ## ----------------------------------------
    gA_formula <- as.formula(
      paste(treatment_name, "~", paste(z_covariates, collapse = " + "))
    )
    gA_fit <- glm(gA_formula, family = binomial(), data = dat)
    
    ## ----------------------------------------
    ## Neighbor treatment model g_S(S | A, Z, Z_N)
    ## Fit as a grouped binomial regression
    ## using cbind(S, d-S)
    ## ----------------------------------------
    gS_formula <- as.formula(
      paste0("cbind(", s_name, ", ", nots_name, ") ~ ",
             treatment_name, " + ",
             paste(c(z_covariates, zn_covariates), collapse = " + "))
    )
    gS_fit <- glm(gS_formula, family = binomial(), data = dat)
    
    ## Pull observed values
    n <- nrow(dat)
    A_obs <- dat[[treatment_name]]
    S_obs <- dat[[s_name]]
    d_obs <- dat[[d_name]]
    Y_obs <- dat[[outcome_name]]
    
    ## ----------------------------------------
    ## Compute fitted probabilities for observed A_i
    ## and bound them away from 0 and 1
    ## ----------------------------------------
    pA1 <- predict(gA_fit, newdata = dat, type = "response")
    pA1 <- pmin(pmax(pA1, weight_bound), 1 - weight_bound)
    
    gA_obs <- ifelse(A_obs == 1, pA1, 1 - pA1)
    gA_obs <- pmax(gA_obs, weight_bound)
    
    ## ----------------------------------------
    ## Compute fitted probabilities for observed S_i
    ## under the grouped binomial model
    ## ----------------------------------------
    pS_obs <- predict(gS_fit, newdata = dat, type = "response")
    #Do we really need this step?
    pS_obs <- pmin(pmax(pS_obs, weight_bound), 1 - weight_bound)
    
    gS_obs <- dbinom(S_obs, size = d_obs, prob = pS_obs)
    gS_obs <- pmax(gS_obs, weight_bound)
    
    ## Joint denominator f2_i(A_i, S_i | Z_i, Z_Ni)
    f2_obs <- pmax(gA_obs * gS_obs, weight_bound)
    
    ########################################################
    ## Estimate ybar^DR(a, alpha)
    ##
    ## DR form:
    ##   IPW residual correction
    ##   + g-computation plug-in term
    ########################################################
    est_ybar_a_alpha <- function(a, alpha_val) {
      
      ## Policy probability for the observed S_i
      p_alpha_obs <- dbinom(S_obs, size = d_obs, prob = alpha_val)
      
      est_i <- numeric(n)
      
      for (i in seq_len(n)) {
        ## Predicted Q at the observed A_i, S_i, covariates
        new_obs <- dat[i, , drop = FALSE]
        q_obs <- predict(Qfit, newdata = new_obs, type = "response")
        
        ## IPW residual correction factor
        H_i <- as.numeric(A_obs[i] == a) * p_alpha_obs[i] / f2_obs[i]
        aug_i <- H_i * (Y_obs[i] - q_obs)
        
        ## g-computation plug-in term:
        ## average Q(a, s, Z_i, Z_Ni) over s under policy alpha
        d_i <- d_obs[i]
        s_vals <- 0:d_i
        p_s <- dbinom(s_vals, size = d_i, prob = alpha_val)
        
        newdat <- data.frame(tmpA = rep(a, length(s_vals)),
                             tmpS = s_vals)
        names(newdat) <- c(treatment_name, s_name)
        ## Hold covariates fixed at the observed values for individual i
        for (v in z_covariates) newdat[[v]] <- dat[[v]][i]
        for (v in zn_covariates) newdat[[v]] <- dat[[v]][i]
        
        q_alpha <- predict(Qfit, newdata = newdat, type = "response")
        plug_i <- sum(q_alpha * p_s)
        
        ## DR estimate contribution for individual i
        est_i[i] <- aug_i + plug_i
      }
      
      mean(est_i)
    }
    
    ########################################################
    ## Estimate ybar^DR(alpha)
    ##
    ## Same logic as above, but now average over both A_i and S_i
    ## under the joint policy for individual and neighbors.
    ########################################################
    est_ybar_alpha <- function(alpha_val) {
      
      ## Joint policy probability for the observed (A_i, S_i)
      p_joint_obs <- dbinom(S_obs, size = d_obs, prob = alpha_val) *
        ifelse(A_obs == 1, alpha_val, 1 - alpha_val)
      
      est_i <- numeric(n)
      
      for (i in seq_len(n)) {
        ## Predicted Q at the observed A_i, S_i, covariates
        new_obs <- dat[i, , drop = FALSE]
        q_obs <- predict(Qfit, newdata = new_obs, type = "response")
        
        ## IPW residual correction factor
        H_i <- p_joint_obs[i] / f2_obs[i]
        aug_i <- H_i * (Y_obs[i] - q_obs)
        
        ## Plug-in term averaged over both A_i and S_i under the policy
        d_i <- d_obs[i]
        s_vals <- 0:d_i
        p_s <- dbinom(s_vals, size = d_i, prob = alpha_val)
        
        newdat0 <- data.frame(tmpA = rep(0, length(s_vals)),
                              tmpS = s_vals)
        newdat1 <- data.frame(tmpA = rep(1, length(s_vals)),
                              tmpS = s_vals)
        names(newdat0) <- c(treatment_name, s_name)
        names(newdat1) <- c(treatment_name, s_name)
        
        ## Hold covariates fixed at the observed values for individual i
        for (v in z_covariates) {
          newdat0[[v]] <- dat[[v]][i]
          newdat1[[v]] <- dat[[v]][i]
        }
        for (v in zn_covariates) {
          newdat0[[v]] <- dat[[v]][i]
          newdat1[[v]] <- dat[[v]][i]
        }
        
        q0 <- predict(Qfit, newdata = newdat0, type = "response")
        q1 <- predict(Qfit, newdata = newdat1, type = "response")
        
        plug_i <- sum(((1 - alpha_val) * q0 + alpha_val * q1) * p_s)
        
        ## DR estimate contribution for individual i
        est_i[i] <- aug_i + plug_i
      }
      
      mean(est_i)
    }
    
    ## Compute APOs for each alpha value
    apo <- do.call(rbind, lapply(alpha, function(a0) {
      data.frame(
        alpha = a0,
        estimand = c("ybarDR(1,alpha)", "ybarDR(0,alpha)", "ybarDR(alpha)"),
        estimate = c(est_ybar_a_alpha(1, a0),
                     est_ybar_a_alpha(0, a0),
                     est_ybar_alpha(a0))
      )
    }))
    
    ## Convert APOs to causal effects
    make_effects_from_apo(
      apo = apo,
      y1_name = "ybarDR(1,alpha)",
      y0_name = "ybarDR(0,alpha)",
      ym_name = "ybarDR(alpha)"
    )
  }
  
  ## Set seed for reproducibility
  set.seed(seed)
  
  ## Point estimates on the original data
  point_eff <- fit_once_dr(data)
  
  ## Number of individuals and effect rows
  n <- nrow(data)
  n_eff <- nrow(point_eff)
  
  ## Storage for bootstrap draws
  boot_mat <- matrix(NA_real_, nrow = B, ncol = n_eff)
  
  ## ----------------------------------------
  ## Bootstrap loop:
  ## resample individuals with replacement and
  ## carry along observed summary neighbor information
  ## ----------------------------------------
  # for (b in seq_len(B)) {
  #   idx <- sample.int(n, size = n, replace = TRUE)
  #   d_b <- data[idx, , drop = FALSE]
  #   
  #   fit_b <- try(fit_once_dr(d_b), silent = TRUE)
  #   
  #   if (!inherits(fit_b, "try-error")) {
  #     boot_mat[b, ] <- fit_b$estimation
  #   }
  # }
  b <- 1
  attempt <- 0
  max_attempts <- B * 20
  
  while (b <= B && attempt < max_attempts) {
    attempt <- attempt + 1
    
    idx <- sample.int(n, size = n, replace = TRUE)
    d_b <- data[idx, , drop = FALSE]
    
    warned <- FALSE
    
    fit_b <- tryCatch(
      withCallingHandlers(
        fit_once_dr(d_b),
        warning = function(w) {
          if (grepl("fitted probabilities numerically 0 or 1 occurred",
                    conditionMessage(w))) {
            warned <<- TRUE
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) NULL
    )
    
    if (!warned && !is.null(fit_b)) {
      boot_mat[b, ] <- fit_b$estimation
      b <- b + 1
    }
  }
  
  if (b <= B) {
    warning("Only ", b - 1, " valid bootstrap samples obtained out of ", B)
  }
  
  
  ## Compute percentile CI bounds
  alpha_tail <- (1 - ci_level) / 2
  lo <- alpha_tail
  hi <- 1 - alpha_tail
  
  ## Bootstrap variances and percentile confidence limits
  var_vec <- apply(boot_mat, 2, var, na.rm = TRUE)
  lci_vec <- apply(boot_mat, 2, quantile, probs = lo, na.rm = TRUE)
  uci_vec <- apply(boot_mat, 2, quantile, probs = hi, na.rm = TRUE)
  
  ## Variance rows corresponding to each effect row
  var_rows <- point_eff
  var_rows$estimation <- var_vec
  var_rows$type[point_eff$type == "Direct"] <- "Var DE"
  var_rows$type[point_eff$type == "Indirect"] <- "Var IE"
  var_rows$type[point_eff$type == "Total"] <- "Var TE"
  var_rows$type[point_eff$type == "Overall"] <- "Var OE"
  
  ## Main rows get CI columns
  out_main <- point_eff
  out_main$lower <- lci_vec
  out_main$upper <- uci_vec
  
  ## Variance rows do not get CI columns
  out_var <- var_rows
  out_var$lower <- NA_real_
  out_var$upper <- NA_real_
  
  ## Stack estimate rows and variance rows
  out <- rbind(out_main, out_var)
  
  ## Keep requested columns
  out <- out[, c("estimation", "alpha0", "alpha1", "type", "lower", "upper")]
  
  ## Reorder rows to match the desired table layout
  out <- reorder_effect_table(out)
  
  ## Clean row names
  rownames(out) <- NULL
  
  list(
    results = out,
    point_estimates = point_eff,
    bootstrap_draws = boot_mat
  )
}
