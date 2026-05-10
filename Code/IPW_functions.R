
################################################################################
# IPW HELPER FUNCTIONS #########################################################

num_neighbors=function(net, node){
  return(length(neighbors(net, node)))
}

trt_neighbors=function(net, node){
  return(sum(subset(df, df$id %in% neighbors(net, node))$treatment))
}

avg_neighbors=function(net, node, variable){
  return(mean(subset(df, df$id %in% neighbors(net, node))[, variable]))
}


pi=function(alpha){
  return(alpha^df$na_a*(1-alpha)^(df$na-df$na_a))
}

Y_IPW=function(alpha, score){
  p=pi(alpha)
  return(c(mean(df$outcome*df$treatment*p/score), 
           mean(df$outcome*(1-df$treatment)*p/score), 
           mean(df$outcome*dbinom(df$treatment, 1, alpha)*p/score)))
}

Y_DE=function(alpha, score){
  p=pi(alpha)
  return(mean(df$outcome*df$treatment*p/score)- 
           mean(df$outcome*(1-df$treatment)*p/score))
}

Y_IE=function(alpha, score){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  return(mean(df$outcome*(1-df$treatment)*p0/score)-
           mean(df$outcome*(1-df$treatment)*p1/score))
}

Y_TE=function(alpha, score){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  return(mean(df$outcome*(df$treatment)*p0/score)-
           mean(df$outcome*(1-df$treatment)*p1/score))
}

Y_OE=function(alpha, score){
  p0=dbinom(df$treatment, 1, alpha[1])*pi(alpha[1])
  p1=dbinom(df$treatment, 1, alpha[2])*pi(alpha[2])
  return(mean(df$outcome*p0/score)-
           mean(df$outcome*p1/score))
}

Var_Y=function(a, alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p=pi(alpha)
  if (a==1){
    U_21=-colMeans((df$outcome*df$treatment*p)*grad_propensity_inv)
    B=df$outcome*df$treatment*p/score
    Y=mean(B)
  } else {
    U_21=-colMeans((df$outcome*(1-df$treatment)*p)*grad_propensity_inv)
    B=df$outcome*(1-df$treatment)*p/score
    Y=mean(B)
  }
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=df$id[df$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var_Y_margin=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p=pi(alpha)
  v_22=0
  v_21=rep(0, length(theta))
  B=df$outcome*dbinom(df$treatment, 1, alpha)*p/score
  Y=mean(B)
  for (j in 1:m){
    c=df$id[df$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*r*colSums(v_11[c,])
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(-colMeans((df$outcome*dbinom(df$treatment, 1, alpha)*p)*grad_propensity_inv), 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var_DE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p=pi(alpha)
  U_21=-colMeans((df$outcome*(2*df$treatment-1)*p)*grad_propensity_inv)
  B=df$outcome*(2*df$treatment-1)*p/score
  Y=mean(B)
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=df$id[df$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var_IE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  U_21=-colMeans(df$outcome*(1-df$treatment)*(p0-p1)*grad_propensity_inv)
  B=df$outcome*(1-df$treatment)*(p0-p1)/score
  Y=mean(B)
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=df$id[df$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var_TE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  U_21=-colMeans(df$outcome*(df$treatment*p0-(1-df$treatment)*p1)*grad_propensity_inv)
  B=df$outcome*(df$treatment*p0-(1-df$treatment)*p1)/score
  Y=mean(B)
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=df$id[df$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var_OE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p0=dbinom(df$treatment, 1, alpha[1])*pi(alpha[1])
  p1=dbinom(df$treatment, 1, alpha[2])*pi(alpha[2])
  v_22=0
  v_21=rep(0, length(theta))
  B=df$outcome*(p0-p1)/score
  Y=mean(B)
  for (j in 1:m){
    c=df$id[df$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*r*colSums(v_11[c,])
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(-colMeans(df$outcome*(p0-p1)*grad_propensity_inv), 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  return(c(Var_Y(1, alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta), 
           Var_Y(0, alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta), 
           Var_Y_margin(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta)))
}

################################################################################
# IPW 1 FUNCTIONS ##############################################################
compute_integrand <- function(theta, subdata, b, base_covariate){
  p <- 1
  theta_num <- as.numeric(theta)
  for (j in 1:nrow(subdata)){
    # Extract covariates as numeric (handle factors/characters)
    x_row <- subdata[j, base_covariate, drop = TRUE]
    x_num <- as.numeric(as.character(x_row))
    x_vec <- c(1, x_num)
    
    beta <- as.numeric(theta_num[1:length(x_vec)])
    pred <- sum(beta * x_vec)
    
    treat <- as.numeric(as.character(subdata$treatment[j]))
    prob <- plogis(pred + as.numeric(b))
    
    p <- p * dbinom(treat, size = 1, prob = prob)
  }
  # Random effect density
  p * dnorm(as.numeric(b), mean = 0, sd = as.numeric(theta_num[length(theta_num)]))
}

# propensity function
propensity <- function(theta){
  prop <- foreach(i = 1:n,
                  .export = c("base_covariate", "net0", "df", "compute_integrand"),
                  .packages = c("plyr","dplyr","igraph","numDeriv","gtools","doParallel"),
                  .combine = "c") %dopar% {
                    subdata <- subset(df, df$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
                    integrand <- function(b){
                      compute_integrand(theta, subdata, b, base_covariate)
                    }
                    as.numeric(integrate(integrand, -Inf, Inf)$value)
                  }
  return(prop)
}

# IPW_1_model function
IPW_1_model <- function(df, base_covariate, alpha){
  
  # Ensure all covariates and treatment are numeric
  df[base_covariate] <- lapply(df[base_covariate], function(x) as.numeric(as.character(x)))
  df$treatment <- as.numeric(as.character(df$treatment))
  n <- as.numeric(n)
  
  # Construct neighborhood data
  nn_data <- data.frame()
  for (i in 1:n){
    sub <- subset(df, df$id %in% neighborhood(net0, order=1, nodes = i)[[1]])
    sub$group <- i
    nn_data <- rbind(nn_data, sub)
  }
  nn_data[base_covariate] <- lapply(nn_data[base_covariate], as.numeric)
  nn_data$treatment <- as.numeric(nn_data$treatment)
  
  # Fit GLMM
  formula <- as.formula(paste("treatment ~", paste(base_covariate, collapse = "+"), "+ (1|group)"))
  M <- glmer(formula, data = nn_data, family = binomial(link= "logit"))
  
  theta <- as.numeric(c(fixef(M), as.data.frame(VarCorr(M))$sdcor))
  score <- propensity(theta)
  
  # Helper function for log_propensity
  log_propensity <- function(theta, subdata){
    integrand <- function(b){ compute_integrand(theta, subdata, b, base_covariate) }
    return(log(as.numeric(integrate(integrand, -Inf, Inf)$value)))
  }
  
  # U_11
  U_11 <- foreach(i = 1:n,
                  .export = c("base_covariate", "net0", "df", "n", "compute_integrand", "log_propensity"),
                  .packages = c("lme4","plyr","dplyr","igraph","numDeriv","gtools","doParallel"),
                  .combine = "+") %dopar% {
                    subdata <- subset(df, df$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
                    -(1/as.numeric(n)) * hessian(function(th) log_propensity(th, subdata), x = theta)
                  }
  
  # v_11
  v_11 <- foreach(i = 1:n,
                  .export = c("base_covariate", "net0", "df", "compute_integrand", "log_propensity"),
                  .packages = c("lme4","plyr","dplyr","igraph","numDeriv","gtools","doParallel"),
                  .combine = "rbind") %dopar% {
                    subdata <- subset(df, df$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
                    grad(function(th) log_propensity(th, subdata), x = theta)
                  }
  
  # grad_propensity_inv
  grad_propensity_inv <- foreach(i = 1:n,
                                 .export = c("base_covariate", "net0", "df","compute_integrand"),
                                 .packages = c("lme4","plyr","dplyr","igraph","numDeriv","gtools","doParallel"),
                                 .combine = "rbind") %dopar% {
                                   subdata <- subset(df, df$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
                                   propensity_inv <- function(th){
                                     integrand <- function(b){ compute_integrand(th, subdata, b, base_covariate) }
                                     1 / as.numeric(integrate(integrand, -Inf, Inf)$value)
                                   }
                                   grad(propensity_inv, x = theta)
                                 }
  
  # V_11
  V_11 <- foreach(j = 1:m,
                  .export = c("df"),
                  .packages = c("plyr","dplyr","igraph","numDeriv","gtools","doParallel"),
                  .combine = "+") %do% {
                    colSums(v_11[df$id[df$component==j], ]) %*% t(colSums(v_11[df$id[df$component==j], ]))
                  }
  V_11 <- V_11 * (m/n^2)
  APO=rbind(cbind(ldply(alpha, Y_IPW, score=score), alpha=alpha, 
                  type="point estimate"), 
            cbind(ldply(alpha, Var, U_11, v_11, V_11, 
                        grad_propensity_inv, score=score, theta=theta), 
                  alpha=alpha, type="variance"))
  names(APO)=c("a=1", "a=0", "margin", "alpha", "type")
  contrast=t(combn(alpha,2)) 
  
  CE=rbind(cbind(ldply(alpha, Y_DE, score=score), 
                 alpha0=alpha, alpha1=alpha, type="Direct"), 
           cbind(ldply(alpha, Var_DE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=alpha, alpha1=alpha, type="Var DE"),
           cbind(adply(contrast, 1, Y_IE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Indirect")[, -1],
           cbind(adply(contrast, 1, Var_IE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var IE")[, -1],
           cbind(adply(contrast, 1, Y_TE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Total")[, -1],
           cbind(adply(contrast, 1, Var_TE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var TE")[, -1],
           cbind(adply(contrast, 1, Y_OE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Overall")[, -1],
           cbind(adply(contrast, 1, Var_OE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var OE")[, -1]
  )
  names(CE)=c("estimation", "alpha0", "alpha1", "type")
  return(list(APO, CE))
}

################################################################################
# IPW 2 FUNCTIONS: #############################################################
propensity_2=function(M1, M2){
  p1=predict(M1, df, type = "response")
  p2=predict(M2, df, type = "response")
  return((p1^df$na_a)*((1-p1)^df$notna_a)*dbinom(df$treatment, size = 1, prob = p2))
}

IPW_2_model=function(df, M1,M2, alpha){
  
  n = as.numeric(nrow(df))
  m = as.numeric(length(unique(df$component)))
  
  score=propensity_2(M1, M2)
  theta=as.numeric(c(as.numeric(M1$coefficients),
                     as.numeric(M2$coefficients)))
  
  k1 = length(M1$coefficients)
  
  U_11 = foreach(i = 1:n, 
                 .packages = c("lme4", "plyr", "dplyr", "igraph", "numDeriv", "gtools", "doParallel", "foreach"),
                 .export = c("base_covariate", "avg_covariate"),
                 .combine = "+") %dopar% {
                   log_propensity_2_ind=function(theta){
                     
                     theta = as.numeric(theta)
                     
                     v1=as.numeric(theta[1:k1])
                     v2=as.numeric(theta[(k1+1):length(theta)])
                     
                     x1 = as.numeric(c(
                       1, 
                       df$treatment[i],
                       as.numeric(df[i, base_covariate, drop = TRUE]),
                       as.numeric(df[i, avg_covariate, drop = TRUE])
                     ))
                     
                     x2 = as.numeric(c(
                       1,
                       as.numeric(df[i, base_covariate, drop = TRUE])
                     ))
                     
                     p1=plogis(sum(v1 * x1))
                     p2=plogis(sum(v2 * x2))
                     
                     na_a = as.numeric(df$na_a[i])
                     notna_a = as.numeric(df$notna_a[i])
                     trt = as.numeric(df$treatment[i])
                     
                     prop=(p1^na_a)*
                       ((1-p1)^notna_a)*
                       dbinom(trt, size = 1, prob = p2)
                     return(log(as.numeric(prop)))
                   }
                   -(1/n)*hessian(log_propensity_2_ind, x=theta)
                 }
  v_11 = foreach(i = 1:n, 
                 .packages = c("lme4", "plyr", "dplyr", "igraph", "numDeriv", "gtools", "doParallel", "foreach"),
                 .export = c("base_covariate", "avg_covariate"),
                 .combine = "rbind") %dopar% {
                   log_propensity_2_ind=function(theta){
                     
                     theta = as.numeric(theta)
                     
                     v1=as.numeric(theta[1:k1])
                     v2=as.numeric(theta[(k1+1):length(theta)])
                     
                     x1 = as.numeric(c(
                       1, 
                       df$treatment[i],
                       as.numeric(df[i, base_covariate, drop = TRUE]),
                       as.numeric(df[i, avg_covariate, drop = TRUE])
                     ))
                     
                     x2 = as.numeric(c(
                       1, 
                       as.numeric(df[i, base_covariate, drop = TRUE])
                     ))
                     
                     p1=plogis(sum(v1 * x1))
                     p2=plogis(sum(v2 * x2))
                     
                     na_a = as.numeric(df$na_a[i])
                     notna_a = as.numeric(df$notna_a[i])
                     trt = as.numeric(df$treatment[i])
                     
                     prop=(p1^na_a)*
                       ((1-p1)^notna_a)*
                       dbinom(trt, size = 1, prob = p2)
                     return(log(as.numeric(prop)))
                   }
                   as.numeric(grad(log_propensity_2_ind, x=theta))
                 }
  grad_propensity_inv = foreach (i = 1:n, 
                                 .packages = c("lme4", "plyr", "dplyr", "igraph", "numDeriv", "gtools", "doParallel", "foreach"),
                                 .export = c("base_covariate", "avg_covariate"),
                                 .combine = "rbind") %dopar% {
                                   propensity_2_ind_inv = function(theta){
                                     
                                     theta = as.numeric(theta)
                                     
                                     v1=as.numeric(theta[1:k1])
                                     v2=as.numeric(theta[(k1+1):length(theta)])
                                     
                                     x1 = as.numeric(c(
                                       1, 
                                       df$treatment[i],
                                       as.numeric(df[i, base_covariate, drop = TRUE]),
                                       as.numeric(df[i, avg_covariate, drop = TRUE])
                                     ))
                                     
                                     x2 = as.numeric(c(
                                       1,
                                       as.numeric(df[i, base_covariate, drop = TRUE])
                                     ))
                                     
                                     p1=plogis(sum(v1 * x1))
                                     p2=plogis(sum(v2 * x2))
                                     
                                     na_a = as.numeric(df$na_a[i])
                                     notna_a = as.numeric(df$notna_a[i])
                                     trt = as.numeric(df$treatment[i])
                                     
                                     prop=(p1^na_a)*
                                       ((1-p1)^notna_a)*
                                       dbinom(trt, size = 1, prob = p2)
                                     return(1/as.numeric(prop))
                                   }
                                   as.numeric(grad(propensity_2_ind_inv, x=theta))
                                 }
  V_11 = foreach(j = 1:m, 
                 .packages = c("lme4", "plyr", "dplyr", "igraph", "numDeriv", "gtools", "doParallel", "foreach"),
                 .export = c("base_covariate", "avg_covariate"),
                 .combine = "+") %dopar% {
                   colSums(v_11[df$id[df$component==j],])%*%t(colSums(v_11[df$id[df$component==j],]))
                 }
  V_11=V_11*(m/n^2)
  APO=rbind(cbind(ldply(alpha, Y_IPW, score=score), alpha=alpha, 
                  type="point estimate"), 
            cbind(ldply(alpha, Var, U_11, v_11, V_11, 
                        grad_propensity_inv, score=score, theta=theta), 
                  alpha=alpha, type="variance"))
  names(APO)=c("a=1", "a=0", "margin", "alpha", "type")
  contrast=t(combn(alpha,2)) 
  
  CE=rbind(cbind(ldply(alpha, Y_DE, score=score), 
                 alpha0=alpha, alpha1=alpha, type="Direct"), 
           cbind(ldply(alpha, Var_DE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=alpha, alpha1=alpha, type="Var DE"),
           cbind(adply(contrast, 1, Y_IE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Indirect")[, -1],
           cbind(adply(contrast, 1, Var_IE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var IE")[, -1],
           cbind(adply(contrast, 1, Y_TE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Total")[, -1],
           cbind(adply(contrast, 1, Var_TE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var TE")[, -1],
           cbind(adply(contrast, 1, Y_OE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Overall")[, -1],
           cbind(adply(contrast, 1, Var_OE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var OE")[, -1]
  )
  names(CE)=c("estimation", "alpha0", "alpha1", "type")
  return(list(APO, CE))
}
