# Hello, world!

kappalate <- function(given_formula, data, zmodel = NULL, vce = NULL, std = NULL, which = NULL, subset = NULL, pstolerance = NULL) {

  # Ensure required packages are loaded
  if (!requireNamespace("formula.tools", quietly = TRUE)) {
    stop("Package 'formula.tools' is required. Install it using install.packages('formula.tools').")
  }
  if (!requireNamespace("geex", quietly = TRUE)) {
    stop("Package 'geex' is required. Install it using install.packages('geex').")
  }
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is required. Install it using install.packages('Formula').")
  }
  # Load libraries silently
  invisible(
    suppressMessages(
      suppressWarnings(
        try({
          library(formula.tools)
          library(geex)
          library(Formula)
        }, silent = TRUE)
      )))

  # Check whether the formula is valid and raise an error if not.
  is.formula <- function(x){
    inherits(x,"formula")
  }

  is_formula <- is.formula(given_formula)
  if (!is_formula) {
    stop("The input must be a valid IV formula of the form: 'outcome ~ exogenous variables | endogenous variable | instrumental variable'. E.g.: 'y ~ x1 + x2 | t | z'.")
  }
  length_correct <- all(length(Formula(given_formula))== c(1,3))
  if (!length_correct) {
    stop("The input must be a valid IV formula of the form: 'outcome ~ exogenous variables | endogenous variable | instrumental variable'. E.g.: 'y ~ x1 + x2 | t | z'.")
  }

  # We parse the given input arguments using the formula.tools and Formula packages.
  # The parsed IV arguments are saved as characters to access the necessary data in the df.
  tvar <- all.vars(rhs(formula(Formula(given_formula), rhs = 2)))
  zvar <- all.vars(rhs(formula(Formula(given_formula), rhs = 3)))
  if (length(tvar)!= 1 || length(zvar)  != 1){
    stop("Only one instrument z and one treatment variable t are allowed.")
    }
  # First we parse the given transformations (if any) of the outcome data
  yvar_0 <- lhs(given_formula)
  transformed_y <- eval(yvar_0, envir = data)
  data$y_col <- transformed_y
  yvar <- as.character("y_col")
  yvar_0<- deparse(yvar_0)
  if (!is.numeric(data$y_col) && !is.integer(data$y_col)) {
    stop("Error: Outcome must be of type numeric or integer.")
  }
  
  # Next for the treatment 
  tvar_0 <- rhs(formula(Formula(given_formula), rhs = 2))
  transformed_t <- eval(tvar_0, envir = data)
  data$t_col <- transformed_t
  tvar <- as.character("t_col")
  tvar_0<- deparse(tvar_0)
  if (!is.numeric(data$t_col) && !is.integer(data$t_col)) {
    stop("Error: treatment must be of type numeric or integer.")
  }
  # and instrument data
  zvar_0 <- rhs(formula(Formula(given_formula), rhs = 3))
  transformed_z <- eval(zvar_0, envir = data)
  data$z_col <- transformed_z
  zvar <- as.character("z_col")
  zvar_0<- deparse(zvar_0)
  if (!is.numeric(data$z_col) && !is.integer(data$z_col)) {
    stop("Error: Instrument must be of type numeric or integer.")
  }
  # Lastly, we parse the given transformations (if any) on the covariate data
  x_formula <- formula(Formula(given_formula), rhs = 1:1)
  x_matrix <- model.matrix(x_formula, data = data)
  x_df <- as.data.frame(x_matrix)
  x_df <- x_df[, !colnames(x_df) %in% "(Intercept)", drop = FALSE]
  xvarsips_0 <- as.character(colnames(x_df))
  names(x_df) <- paste0("x_col_", 1:ncol(x_df))
  xvarsips <- as.character(colnames(x_df))
  data <- cbind(data, x_df)
  
  # We check if only a subset of the data should be used
  # First we assert what input constitutes a valid index vector
  is_valid_index <- function(index, data) {
    if (is.null(index)) return(FALSE)
    n <- nrow(data)
    if (all(is.numeric(index))) {
      return(all(index > 0 & index <= n & index == as.integer(index)))
    } else if (is.logical(index)) {
      return(length(index) == n)
    } else {
      return(FALSE)
    }
  }

  # Then we check if a subset is given and if it is a valid index vector and locally assign the subset
  if (!is.null(subset)){
    if (is_valid_index(subset, data)){
      data <- data[subset, ]
    } else {
      stop("The specified subset is not a valid index vector. It should either be numeric e.g. subset = c(1,2,3,4) where all specified index values must be part of the dataframe or logical e.g. subset = c(TRUE, FALSE, FALSE), with length = length of the dataframe.")
    }
  }

  # Next we locally remove rows with missing data
  data <- na.omit(data[, c(yvar, xvarsips, tvar, zvar)])

  # Then we assign some standard values for options that were not set.
  if (is.null(zmodel)) {
    zmodel <- "cbps"}
  if (is.null(vce)) {
    vce <- "robust"}
  if (is.null(std)) {
    std <- "on"}
  if (is.null(which)) {
    which <- "norm"}
  if (is.null(pstolerance)) {
    pstolerance <- 1e-5}

  # Next we check some input conditions and raise errors if they are not fulfilled.
  # First for input strings etc.
  if (!zmodel %in% c("cbps", "logit", "probit")) {
    stop(sprintf("zmodel = '%s'  not allowed. Allowed values are 'cbps', 'logit', or 'probit'.", zmodel))
  }
  if (!std %in% c("on", "off")) {
    stop(sprintf("std = '%s' not allowed. Allowed values are 'on' or 'off'.", std))
  }
  if (!which %in% c("norm", "all")) {
    stop(sprintf("which = '%s'= not allowed. Allowed values are 'norm' or 'all'.", which))
  }
  if (!is.numeric(pstolerance) || length(pstolerance) != 1 || pstolerance <= 0 || pstolerance >= 1) {
    stop(sprintf("pstolerance = '%s'= not allowed. pstolerance` must be a single numeric value strictly between 0 and 1.", pstolerance))
  }

  # Then for some requirements on the treatment and instrument data:
  tvar_values <- unique(data[[tvar]])
  if (length(tvar_values) < 2) {
    stop("Treatment variable must take on at least two distinct values.")
  }
  if (length(tvar_values) > 2) {
    bintreat <- 0
  } else if (length(tvar_values) == 2) {
    if (all(tvar_values %in% c(0,1))) {
      bintreat <- 1
    } else {
      bintreat <- 0
    }
  }

  zvar_values <- unique(data[[zvar]])
  if (length(zvar_values) != 2) {
    stop("Instrument  must be binary.")
  }
  if (!all(zvar_values %in% c(0,1))) {
    stop("Instrument must only take on values zero or one")
  }

  # Next we standardize the non-binary exogenous variables if std == "on" and assign the scaled values to the data frame.
  if (std == "on"){
    for (col in xvarsips){
      unique_values <- unique(data[[col]])
      if (length(unique_values) != 2){
        data[[col]] <- scale(data[[col]])
      }
    }
  }

  # We then examine compliance with the instrument.
  dmeanz1 <- mean(as.numeric(unlist(data[data[[zvar]] == 1, tvar])), na.rm = TRUE)
  dmeanz0 <- mean(as.numeric(unlist(data[data[[zvar]] == 0, tvar])), na.rm = TRUE)
  

  if (bintreat == 1) {
    if (dmeanz0 == dmeanz1) {
      stop("Error: zero denominator - LATE not defined.")
    } else if (dmeanz1 == 1 & dmeanz0 == 0) {
      stop("Error: instrument identical to treatment.")
    } else if (dmeanz1 == 0 & dmeanz0 == 1) {
      stop("Error: instrument identical to treatment.")
    }
  }

  # We can now begin with the main estimation procedure
  # Start with estimation of the instrument propensity score when no cbps is selected
  # Note that pz <- exp(zhat) / (1 + exp(zhat)) with the z estimates in the link function is equal to using the predict function
  formula_z <- reformulate(termlabels = xvarsips, response = zvar)
  if (zmodel != "cbps") {
    if (zmodel == "logit") {
      logit_model <- glm(formula_z, data = data, family = binomial())
      bips <- coef(logit_model)
      ips <- predict(logit_model, type = "response")
    } else if (zmodel == "probit") {
      probit_model <- glm(formula_z, data = data, family = binomial(link = "probit"))
      bips <- coef(probit_model)
      ips <- predict(probit_model, type = "response")
    }

  # Now for the case with cbps
  # Determine starting values from logit
  } else {
    logit_model <- glm(formula_z, data = data, family = binomial())
    initial <- coef(logit_model)

  # Then set up the moment condition in geex following covariate balancing by (Imai & Ratkovic, 2013). Here the moment condition is E[((z - pz)/(pz(1-pz)))X] = 0, with pz <- exp(zhat) / (1 + exp(zhat)), and zhat <- X*theta.
    cbps_moment <- function(data){
      X <- model.matrix(formula_z, data = data)
      Z <- data[[zvar]]
      function(theta){
        zhat <- X %*% theta
        pz <- plogis(zhat)
        c(((Z - pz)/(pz*(1-pz)))%*%X)
      }
    }
    # The moment estimation in geex is then done via m_estimate(). If convergence fails we continue with the logistic regression estimates.
    tryCatch({
      cbps_results <- m_estimate(
      estFUN = cbps_moment,
      data   = data,
      root_control = setup_root_control(FUN = rootSolve::multiroot, start = initial)
      )
    bips <- coef(cbps_results)
    ips <- plogis(model.matrix(formula_z, data = data) %*% bips)
    }, error = function(e) {
      message(sprintf("An error occurred using M-Estimation for the covariate balancing of treatment propensities: %s Did not acheive convergence. Instead logistic regression is used to estimate treatment propensity.", e$message))
      bips <<- coef(logit_model)
      ips <<- predict(logit_model, type = "response")
      })
    }

  # We check for overlap assumptions and if any are present return a vector of indices corresponding to the observations (rows) where the assumption is violated.
  violators <- (ips < pstolerance | ips > (1 - pstolerance))
  violating_indices <- as.integer(rownames(data))[violators]
  fail <- sum(violators)
  
  if (fail > 0) {
    plural <- ifelse(fail == 1, "", "s")
    plural2 <- ifelse(fail == 1, "s", "")
    range_msg <- paste0("[",
                        format(pstolerance, scientific = FALSE),
                        ", ",
                        format(1 - pstolerance, scientific = FALSE),
                        "]")
    
    cat(paste0(fail, " observation", plural,
        " violate", plural2,
        " the overlap assumption with an instrument propensity score outside the range ",
        range_msg), ".\n")
    
    cat("Error: Instrument overlap assumption has been violated. Program will exit.\n")
    cat(paste0("Consider changing the pstolerance. Currently pstolerance = ", format(pstolerance, scientific = TRUE),  ".\n"))
    cat("To help identify the overlap violators a vector of the row indices corresponding to the observations in the passed data frame where the assumption is violated has been returned.")
    return(violating_indices)
  }  
  
  # LATE Estimation
  # Since until now we stored only the variable names we create the matrices containing the actual data points
  y_data <- data[[yvar]]
  z_data <- data[[zvar]]
  t_data <- data[[tvar]]
  # We then create the variables identifying the LATE. See Abadie's Kappa and Weighting Estimators of the Local Average Treatment Effect" by Tymon Słoczyński, Derya Uysal, & Jeffrey M. Wooldridge for the detailed theory.
  numhat <- (z_data / ips) * y_data - ((1 - z_data) / (1 - ips)) * y_data
  kappa_1 <- (z_data / ips) * t_data - ((1 - z_data) / (1 - ips)) * t_data
  kappa_0 <- (1 - t_data) * ((1 - z_data) - (1 - ips)) / (ips * (1 - ips))
  kappaw <- 1 - (t_data * (1 - z_data)) / (1 - ips) - ((1 - t_data) * z_data) / ips
  num1hat <- kappa_1 * y_data
  num0hat <- kappa_0 * y_data
  # and their sample means
  nums <- mean(numhat)
  kappa_1s <- mean(kappa_1)
  kappa_0s <- mean(kappa_0)
  kappas <- mean(kappaw)
  num1hats <- mean(num1hat)
  num0hats <- mean(num0hat)
  # The different LATE estimators a, a1, a0, a10 are then calculated as follows. See the paper for more details.
  late_a <- nums/kappas
  late_a1 <- nums/kappa_1s
  late_a0 <- nums/kappa_0s
  late_a10 <- num1hats/kappa_1s - num0hats/kappa_0s
  list_lates_a <- c(late_a,late_a1,late_a0,late_a10)
  # Next we regress y on a constant, for the cases where z==1, and z==0 individually. Then we do the same for t.
  formula_y_t <- reformulate(termlabels = c("1"), response = yvar)
  # First y(z==1).
  model_y_t_z1 <- lm(formula_y_t, data = data[data[[zvar]]==1,], weights = (1/ips[which(data[[zvar]] == 1)]))
  by1 <- coef(model_y_t_z1)
  y1hat <- predict(model_y_t_z1)
  # Then y(z==0).
  model_y_t_z0 <- lm(formula_y_t, data = data[data[[zvar]]==0,], weights = (1/(1-ips[which(data[[zvar]] == 0)])))
  by0 <- coef(model_y_t_z0)
  y0hat <- predict(model_y_t_z0)
  # t(z==1)
  if (dmeanz1 == 1){
    d1hat <- 1
  } else if (dmeanz1 == 0){
    d1hat <- 0
  } else {
    formula_t_z <- reformulate(termlabels = c("1"), response = tvar)
    model_t_z1 <- lm(formula_t_z, data = data[data[[zvar]]==1,], weights = (1/ips[which(data[[zvar]] == 1)]))
    bd1 <- coef(model_t_z1)
    d1hat <- predict(model_t_z1)
  }
  # t(z==0)
  if (dmeanz0 == 0){
    d0hat <- 0
  } else if (dmeanz0 == 1){
    d0hat <- 1
  } else {
    formula_t_z <- reformulate(termlabels = c("1"), response = tvar)
    model_t_z0 <- lm(formula_t_z, data = data[data[[zvar]]==0,], weights = (1/(1-ips[which(data[[zvar]] == 0)])))
    bd0 <- coef(model_t_z0)
    d0hat <- predict(model_t_z0)
  }

  denom1s <- mean(d1hat)
  denom0s <- mean(d0hat)
  num1s <- mean(y1hat)
  num0s <- mean(y0hat)

  num_norms <- num1s - num0s
  denom_norms <- denom1s - denom0s
  late_norm <- num_norms/denom_norms


  # Moment conditions for first stage
  if (zmodel == "logit"){
    eqips <- expression((Z - pz)%*%X)
  } else if (zmodel =="probit"){
    eqips <- expression((((Z - pz)/(pz*(1-pz)))*pd)%*%X)
  } else if (zmodel == "cbps"){
    eqips <- expression(((Z - pz)/(pz*(1-pz)))%*%X)
  }
  # Other moments
  # psi_delta
  eq_delta <- expression((Z*Y)/pz - ((1-Z)*Y)/(1-pz) - deltap)
  # psi_gamma
  eq_gamma <- expression(1 - ((1-Z)*t)/(1-pz) - Z*(1-t)/pz - gammap)
  # psi_gamma1
  eq_gamma1 <- expression((Z*t)/pz - ((1-Z)*t)/(1-pz) - gamma1)
  # psi_gamma0
  eq_gamma0 <- expression((Z*(t-1))/pz - ((1-Z)*(t-1))/(1-pz) - gamma0)
  # psi_delta1
  eq_delta1 <- expression(t*((Z-pz)/(pz*(1-pz)))*Y - delta1)
  # psi_delta0
  eq_delta0 <- expression((1-t)*(((1-Z)-(1-pz))/(pz*(1-pz)))*Y - delta0)
  # psi_mu1
  eq_mu1 <- expression((Z*(Y-mu1))/pz)
  # psi_mu0
  eq_mu0 <- expression(((1-Z)*(Y-mu0))/(1-pz))
  # psi_m1
  eq_m1 <- expression((Z*(t-m1))/pz)
  # psi_m0
  eq_m0 <- expression(((1-Z)*(t-m0))/(1-pz))

  # psi_taua
  eq_tau_a <- expression (tau_a - deltap/gammap)
  # psi_taua1
  eq_tau_a1 <- expression (tau_a1 - deltap/gamma1)
  # psi_taua0
  eq_tau_a0 <- expression (tau_a0 - deltap/gamma0)
  # psi_taua10
  eq_tau_a10 <- expression(tau_a10 - (delta1/gamma1 - delta0/gamma0))

  if (dmeanz0 != 0 & dmeanz0 != 1 & dmeanz1 != 0 & dmeanz1 != 1) {
    eq_tau_norm <- expression(tau_norm - (mu1 - mu0)/(m1 - m0))
  } else if (dmeanz0 == 0 & dmeanz1 != 0 & dmeanz1 != 1) {
    eq_tau_norm <- expression(tau_norm - (mu1 - mu0)/(m1))
  } else if (dmeanz0 != 0 & dmeanz0 != 1 & dmeanz1 == 1) {
    eq_tau_norm <- expression(tau_norm - (mu1 - mu0)/(1 - m0))
  } else if (dmeanz0 == 1 & dmeanz1 != 0 & dmeanz1 != 1) {
    eq_tau_norm <- expression(tau_norm - (mu1 - mu0)/(m1 - 1))
  } else if (dmeanz0 != 0 & dmeanz0 != 1 & dmeanz1 == 0) {
    eq_tau_norm <- expression(tau_norm - (mu1 - mu0)/(-m0))
  }

  tryCatch({
  # M-Estimation of the LATE
  # tau_a
  if (bintreat == 1){
  if (which == "all"){
  tau_a_condition <- function(data){
    X <- model.matrix(formula_z, data = data)
    Z <- data[[zvar]]
    Y <- data[[yvar]]
    t <- data[[tvar]]
    function(theta){
      zhat <- X %*% theta[1:ncol(X)]
      deltap <- theta[ncol(X)+1]
      gammap <- theta[ncol(X)+2]
      tau_a <- theta[ncol(X)+3]
      pz <- plogis(zhat)
      if (zmodel == "probit"){
        pz <- pnorm(zhat)}
      pd <- dnorm(zhat)
      c(eval(eqips), eval(eq_delta), eval(eq_gamma), eval(eq_tau_a))
    }
  }
  tau_a_estimation <- m_estimate(
    estFUN = tau_a_condition,
    data   = data,
    root_control = setup_root_control(start = c(bips,nums,kappas,late_a)))

  coef_tau_a <- coef(tau_a_estimation)
  tau_a <- coef_tau_a[length(coef_tau_a)]
  vcov_tau_a <- vcov(tau_a_estimation)
  var_tau_a <- vcov_tau_a[nrow(vcov_tau_a),ncol(vcov_tau_a)]

  # tau_a1
  tau_a1_condition <- function(data){
    X <- model.matrix(formula_z, data = data)
    Z <- data[[zvar]]
    Y <- data[[yvar]]
    t <- data[[tvar]]
    function(theta){
      zhat <- X %*% theta[1:ncol(X)]
      deltap <- theta[ncol(X)+1]
      gamma1 <- theta[ncol(X)+2]
      tau_a1 <- theta[ncol(X)+3]
      pz <- plogis(zhat)
      if (zmodel == "probit"){
        pz <- pnorm(zhat)}
      pd <- dnorm(zhat)
      c(eval(eqips), eval(eq_delta), eval(eq_gamma1), eval(eq_tau_a1))
    }
  }
  tau_a1_estimation <- m_estimate(
    estFUN = tau_a1_condition,
    data   = data,
    root_control = setup_root_control(start = c(bips,nums,kappa_1s,late_a1)))

  coef_tau_a1 <- coef(tau_a1_estimation)
  tau_a1 <- coef_tau_a1[length(coef_tau_a1)]
  vcov_tau_a1 <- vcov(tau_a1_estimation)
  var_tau_a1 <- vcov_tau_a1[nrow(vcov_tau_a1),ncol(vcov_tau_a1)]

  # tau_a0
  tau_a0_condition <- function(data){
    X <- model.matrix(formula_z, data = data)
    Z <- data[[zvar]]
    Y <- data[[yvar]]
    t <- data[[tvar]]
    function(theta){
      zhat <- X %*% theta[1:ncol(X)]
      deltap <- theta[ncol(X)+1]
      gamma0 <- theta[ncol(X)+2]
      tau_a0 <- theta[ncol(X)+3]
      pz <- plogis(zhat)
      if (zmodel == "probit"){
        pz <- pnorm(zhat)}
      pd <- dnorm(zhat)
      c(eval(eqips), eval(eq_delta), eval(eq_gamma0), eval(eq_tau_a0))
    }
  }
  tau_a0_estimation <- m_estimate(
    estFUN = tau_a0_condition,
    data   = data,
    root_control = setup_root_control(start = c(bips,nums,kappa_0s,late_a0)))

  coef_tau_a0 <- coef(tau_a0_estimation)
  tau_a0 <- coef_tau_a0[length(coef_tau_a0)]
  vcov_tau_a0 <- vcov(tau_a0_estimation)
  var_tau_a0 <- vcov_tau_a0[nrow(vcov_tau_a0),ncol(vcov_tau_a0)]
  }

  # tau_a10
  if (zmodel != "cbps"){
  tau_a10_condition <- function(data){
    X <- model.matrix(formula_z, data = data)
    Z <- data[[zvar]]
    Y <- data[[yvar]]
    t <- data[[tvar]]
    function(theta){
      zhat <- X %*% theta[1:ncol(X)]
      delta1 <- theta[ncol(X)+1]
      gamma1 <- theta[ncol(X)+2]
      delta0 <- theta[ncol(X)+3]
      gamma0 <- theta[ncol(X)+4]
      tau_a10 <- theta[ncol(X)+5]
      pz <- plogis(zhat)
      if (zmodel == "probit"){
        pz <- pnorm(zhat)}
      pd <- dnorm(zhat)
      c(eval(eqips), eval(eq_delta1),eval(eq_gamma1), eval(eq_delta0), eval(eq_gamma0), eval(eq_tau_a10))
    }
  }
  tau_a10_estimation <- m_estimate(
    estFUN = tau_a10_condition,
    data   = data,
    root_control = setup_root_control(start = c(bips, num1hats, kappa_1s, num0hats, kappa_0s, late_a10)))

  coef_tau_a10 <- coef(tau_a10_estimation)
  tau_a10 <- coef_tau_a10[length(coef_tau_a10)]
  vcov_tau_a10 <- vcov(tau_a10_estimation)
  var_tau_a10 <- vcov_tau_a10[nrow(vcov_tau_a10),ncol(vcov_tau_a10)]
  }
  }

  # tau_norm
  if (bintreat == 0 || (dmeanz0 != 0 && dmeanz0 != 1 && dmeanz1 != 0 && dmeanz1 != 1)) {
    tau_norm_condition <- function(data){
      X <- model.matrix(formula_z, data = data)
      Z <- data[[zvar]]
      Y <- data[[yvar]]
      t <- data[[tvar]]
      function(theta){
        zhat <- X %*% theta[1:ncol(X)]
        mu1 <- theta[ncol(X)+1]
        mu0 <- theta[ncol(X)+2]
        m1 <- theta[ncol(X)+3]
        m0 <- theta[ncol(X)+4]
        tau_norm <- theta[ncol(X)+5]
        pz <- plogis(zhat)
        if (zmodel == "probit"){
          pz <- pnorm(zhat)}
        pd <- dnorm(zhat)
        c(eval(eqips), eval(eq_mu1),eval(eq_mu0), eval(eq_m1), eval(eq_m0), eval(eq_tau_norm))
      }
    }
    tau_norm_estimation <- m_estimate(
      estFUN = tau_norm_condition,
      data   = data,
      root_control = setup_root_control(start = c(bips, num1s, num0s, denom1s, denom0s, late_norm)))
    }

  else if (bintreat == 1 && (dmeanz0 == 0 || dmeanz0 == 1 )){
    tau_norm_condition <- function(data){
      X <- model.matrix(formula_z, data = data)
      Z <- data[[zvar]]
      Y <- data[[yvar]]
      t <- data[[tvar]]
      function(theta){
        zhat <- X %*% theta[1:ncol(X)]
        mu1 <- theta[ncol(X)+1]
        mu0 <- theta[ncol(X)+2]
        m1 <- theta[ncol(X)+3]
        tau_norm <- theta[ncol(X)+4]
        pz <- plogis(zhat)
        if (zmodel == "probit"){
          pz <- pnorm(zhat)}
        pd <- dnorm(zhat)
        c(eval(eqips), eval(eq_mu1),eval(eq_mu0), eval(eq_m1), eval(eq_tau_norm))
      }
    }
    tau_norm_estimation <- m_estimate(
      estFUN = tau_norm_condition,
      data   = data,
      root_control = setup_root_control(start = c(bips, num1s, num0s, denom1s, late_norm)))
  }

  else if (bintreat == 1 && (dmeanz1 == 0 || dmeanz1 == 1 )){
    tau_norm_condition <- function(data){
      X <- model.matrix(formula_z, data = data)
      Z <- data[[zvar]]
      Y <- data[[yvar]]
      t <- data[[tvar]]
      function(theta){
        zhat <- X %*% theta[1:ncol(X)]
        mu1 <- theta[ncol(X)+1]
        mu0 <- theta[ncol(X)+2]
        m0 <- theta[ncol(X)+3]
        tau_norm <- theta[ncol(X)+4]
        pz <- plogis(zhat)
        if (zmodel == "probit"){
          pz <- pnorm(zhat)}
        pd <- dnorm(zhat)
        c(eval(eqips), eval(eq_mu1),eval(eq_mu0), eval(eq_m0), eval(eq_tau_norm))
      }
    }
    tau_norm_estimation <- m_estimate(
      estFUN = tau_norm_condition,
      data   = data,
      root_control = setup_root_control(start = c(bips, num1s, num0s, denom0s, late_norm)))
  }

  coef_tau_norm <- coef(tau_norm_estimation)
  tau_norm <- coef_tau_norm[length(coef_tau_norm)]
  vcov_tau_norm <- vcov(tau_norm_estimation)
  var_tau_norm <- vcov_tau_norm[nrow(vcov_tau_norm),ncol(vcov_tau_norm)]
  
  }, error = function(e) {
    stop(sprintf(
      "Error: Convergence of M-Estimation unsuccessful. Make sure that standarization of covariates is on as this makes convergence more likely in some cases. \nDetails: %s",
      e$message
    ))
  })
   
  # Results of the LATE Estimation
  if (bintreat == 0) {
    b <- matrix(tau_norm, nrow = 1, ncol = 1)
    V <- matrix(var_tau_norm, nrow = 1, ncol = 1)
    rownames(b) <- c("")
    colnames(b) <- c("tau_u")
    rownames(V) <- colnames(V) <- c("tau_u")
  }
  else if (bintreat == 1){
    if (which == "all") {
      if (zmodel != "cbps") {
        b <- matrix(c(tau_a, tau_a1, tau_a0, tau_a10, tau_norm), nrow = 1, ncol = 5)
        V <- matrix(0, nrow = 5, ncol = 5)
        diag(V) <- c(var_tau_a, var_tau_a1, var_tau_a0, var_tau_a10, var_tau_norm)
        rownames(b) <- c("")
        colnames(b) <- c("tau_a", "tau_a,1", "tau_a,0", "tau_a,10", "tau_u")
        rownames(V) <- colnames(V) <- c("tau_a", "tau_a,1", "tau_a,0", "tau_a,10", "tau_u")
        }
      else if (zmodel == "cbps") {
        b <- matrix(c(tau_a, tau_norm), nrow = 1, ncol = 2)
        V <- matrix(0, nrow = 2, ncol = 2)
        diag(V) <- c(var_tau_a, var_tau_norm)
        rownames(b) <- c("")
        colnames(b) <- c("tau_a", "tau_u")
        rownames(V) <- colnames(V) <- c("tau_a", "tau_u")
      }
    }
    else if (which == "norm") {
      if (zmodel != "cbps") {
        b <- matrix(c(tau_a10, tau_norm), nrow = 1, ncol = 2)
        V <- matrix(0, nrow = 2, ncol = 2)
        diag(V) <- c(var_tau_a10, var_tau_norm)
        rownames(b) <- c("")
        colnames(b) <- c("tau_a,10", "tau_u")
        rownames(V) <- colnames(V) <- c("tau_a,10", "tau_u")
      }
      else if (zmodel == "cbps") {
        b <- matrix(tau_norm, nrow = 1, ncol = 1)
        V <- matrix(var_tau_norm, nrow = 1, ncol = 1)
        rownames(b) <- c("")
        colnames(b) <- c("tau_u")
        rownames(V) <- colnames(V) <- c("tau_u")
      }
    }
  }


  # Return these variables
  ret <- list(
    coefficients = b,
    vcov_matrix = V,
    N = nrow(data),
    metadata = list(
      title = "LATE estimation",
      depvar = yvar_0,
      tvar = tvar_0,
      zvar = zvar_0,
      xvarsips = xvarsips_0,
      zmodel = zmodel,
      cmd = "kappalate",
      formula = given_formula
  )
  )
  # Define a new class
  class(ret) <- "kappalate"
  # return the list
  return(ret)


  }#Ends the function code
