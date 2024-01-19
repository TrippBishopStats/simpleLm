compute_alpha <- function(Y) {
  sum(Y)/length(Y)
}

compute_beta <- function(X, Y) {
  n <- length(X)
  (sum(X*Y) - 1/n*sum(X)*sum(Y))/(sum(X^2) - 1/n*sum(X)^2)
}

compute_var <- function(X,Y,alpha=NULL,beta=NULL) {

  if (is.null(alpha)) {
    alpha <- compute_alpha(Y)
  }

  if(is.null(beta)) {
    beta <- compute_beta(X,Y)
  }

  n <- length(X)

  (sum(Y^2) - 1/n*sum(Y)^2 - beta*sum(X*Y) + beta*1/n*sum(X)*sum(Y))/n
}

generate_conf_int_alpha <- function(alpha_hat, s2_hat, n, level=0.95) {
  gamma <- 1 - level
  interval <- qt(1-gamma/2, n - 2)*sqrt(s2_hat/(n-2))
  c(alpha_hat - interval, alpha_hat + interval)
}

generate_conf_int_beta <- function(beta_hat, s2_hat, X, level=0.95) {
  n <- length(X)
  gamma <- 1 - level
  X_bar <- mean(X)
  interval <- qt(1-gamma/2, n - 2)*sqrt((n*s2_hat)/((n-2)*sum((X - X_bar)^2)))
  c(beta_hat - interval, beta_hat + interval)
}

generate_conf_int_s2 <- function(s2_hat, n, level=0.95) {
  gamma <- 1 - level
  lower <- (n*s2_hat)/(qchisq(1-gamma/2, n-2))
  upper <- (n*s2_hat)/(qchisq(gamma/2, n-2))
  c(lower, upper)
}

residuals <- function(X, Y, coeffs) {
  Y - coeffs["intercept"] - coeffs["slope"]*X
}

poly_fit <- function(X, coeffs) {
  value <- coeffs[1]
  for(i in 2:length(coeffs)) {
    value <- value + X^(i-1)*coeffs[i]
  }
  value
}

compute_r2 <- function(Y, residuals) {
  RSS <- sum(residuals^2)
  Y_bar <- mean(Y)
  TSS <- sum((Y - Y_bar)^2)
  1 - RSS/TSS
}

#' @title create_lm
#' @description Computes the coefficients and residuals for a simple linear regression.
#' @author Tripp Bishop
#' @export create_lm
#' @param X A vector of the predictor variable.
#' @param Y A vector of the response variable.
#' @param level A number between 0 and 1, the level of confidence.
#' @return A list containing parameters, confidence intervals of the parameters, residuals, and the coefficient of determination.
create_lm <- function(X,Y, level = 0.95) {

  ##################################################################
  # To begin, confirm that the parameters contain reasonable data. #
  ##################################################################

  ### DEPRECATED ADD NEW ASSERTION MECHANISM

  # X %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric
  #
  # Y %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric
  #
  # assert_are_same_length(X,Y)
  #
  # level %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric %>%
  #   assert_any_are_in_range(lower = 0, upper = 1)

  #########################################################################
  # Now that the inputs have been validated, proceed with the regression. #
  #########################################################################

  beta <- compute_beta(X,Y)
  alpha <- compute_alpha(Y) - beta*mean(X)
  s2 <- compute_var(X,Y,alpha,beta)
  n <- length(X)

  conf_alpha <- generate_conf_int_alpha(alpha, s2, n, level)
  conf_beta <- generate_conf_int_beta(beta, s2, X, level)
  conf_s2 <- generate_conf_int_s2(s2, n, level)

  vars <- c(alpha, beta, s2)
  names(vars) <- c("intercept","slope","variance")
  conf_intervals <- list(`Conf alpha` = conf_alpha, `Conf beta` = conf_beta, `Conf var` = conf_s2)
  resids <- residuals(X,Y,vars)

  list(
    data = list(predictor=X, response=Y),
    coeffs = vars,
    conf_int = conf_intervals,
    residuals = resids,
    R2 = compute_r2(Y, resids)
  )
}

#' @title plot_fit
#' @description Generates a plot of model data and the best fit regression line determined by the model.
#' @author Tripp Bishop
#' @export plot_fit
#' @param model A create_lm model object.
#' @return A ggplot2 object containing a plot of the model fit.
plot_fit <- function(model) {
  tibble(predictor=model$data$predictor, response=model$data$response) %>%
    ggplot(aes(x=predictor, y=response)) +
    geom_point() +
    geom_abline(intercept = model$coeffs["intercept"], slope = model$coeffs["slope"], colour="#99CC99", size=1.25)
}

#' @title plot_predictions
#' @description Generates a plot of model predictions and optionally confidence and prediction intervals.
#' @author Tripp Bishop
#' @export plot_predictions
#' @param predictions A tibble containing model prediction data.
#' @param se The mean and prediction intervals. Default is FALSE.
#' @return A ggplot2 object containing a plot of the model predictions.
plot_predictions <- function(predictions, se = FALSE) {

  ### DEPRECATED ADD NEW ASSERTION MECHANISM
  # predictions %>% assert_is_data.frame()
  ###

  gg_obj <- predictions %>%
    ggplot(aes(x=new_values)) +
    geom_point(aes(y=predictions)) +
    geom_line(aes(y=predictions), colour="#99CC99")

  if(se) {
    gg_obj <- gg_obj +
      geom_line(aes(y=lower_mean_bounds), colour="red", size=1) +
      geom_line(aes(y=upper_mean_bounds), colour="red", size=1) +
      geom_line(aes(y=lower_pred_bounds), colour="royalblue", size=1) +
      geom_line(aes(y=upper_pred_bounds), colour="royalblue", size=1)
  }

  gg_obj
}

#' @title plot_residuals
#' @description Generates a scatter plot of the residuals generated by the specified create_lm model object.
#' @author Tripp Bishop
#' @export plot_residuals
#' @param model A create_lm model object.
#' @return A ggplot2 object containing a plot of the residual values of the create_lm model.
plot_residuals <- function(model) {
  tibble(predictor=model$data$predictor, residuals=model$residuals) %>%
    ggplot(aes(x=predictor, y=residuals)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 0, colour="#aa0099", size=1.25)
}

#' @title compute_predictions
#' @description Generates predictions for the response given a model and a vector of values of the predictor.
#' @author Tripp Bishop
#' @export compute_predictions
#' @param model A create_lm model object.
#' @param new_data The values of the predictor variable to make predictions for.
#' @param level A number between 0 and 1, the level of confidence.
#' @return A tibble containing the predictions as well as response mean and prediction intervals.
compute_predictions <- function(model, new_data, level = 0.95) {

  ##################################################################
  # To begin, confirm that the parameters contain reasonable data. #
  ##################################################################

  ### DEPRECATED ADD NEW ASSERTION MECHANISM

  # new_data %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric
  #
  # level %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric %>%
  #   assert_any_are_in_range(lower = 0, upper = 1)

  #########################################################################
  # Now that the inputs have been validated, proceed with the predictions #
  #########################################################################

  X <- model$data$predictor
  n <- length(X)
  X_bar <- mean(X)
  gamma <- 1 - level

  t <- qt(1 - gamma/2, n - 2)
  ssx <- sum((X - X_bar)**2)
  a <- sqrt((n*model$coeffs["variance"])/(n-2))
  new_diff <- (new_data - X_bar)**2
  c <- a*sqrt(1/n + new_diff/ssx)
  d <- a*sqrt(1 + 1/n + new_diff/ssx)

  mean_interval_width <- c*t
  pred_interval_width <- d*t

  preds <- model$coeffs["intercept"] + model$coeffs["slope"]*(new_data)
  lower_mean <- preds - mean_interval_width
  upper_mean <- preds + mean_interval_width
  lower_pred <- preds - pred_interval_width
  upper_pred <- preds + pred_interval_width

  tibble(
    new_values = new_data,
    predictions = preds,
    lower_mean_bounds = lower_mean,
    upper_mean_bounds = upper_mean,
    lower_pred_bounds = lower_pred,
    upper_pred_bounds = upper_pred
  )
}

#' @title polynomial_lm
#' @description Computes the coefficients and residuals for a simple polynomial linear regression.
#' @author Tripp Bishop
#' @export polynomial_lm
#' @param X A vector of the predictor variable.
#' @param Y A vector of the response variable.
#' @param order An integer greater than or equal to 2. The order of polynomial to fit to the data.
#' @return A list containing parameters, residuals, and the coefficient of determination.
polynomial_lm <- function(X, Y, order=2) {

  ##################################################################
  # To begin, confirm that the parameters contain reasonable data. #
  ##################################################################

  ### DEPRECATED ADD NEW ASSERTION MECHANISM

  # X %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric
  #
  # Y %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric
  #
  # assert_are_same_length(X,Y)
  #
  # order %>%
  #   assert_is_not_null %>%
  #   assert_is_numeric %>%
  #   is_greater_than_or_equal_to(2)

  #########################################################################
  # Now that the inputs have been validated, proceed with the regression. #
  #########################################################################

  k <- order + 1
  n <- length(X)

  # build the predictor matrix
  mat_X <- cbind(rep(1, n), X)

  for(i in 2:order) {
    X_order <- X**i
    mat_X <- cbind(mat_X, X_order)
  }

  colnames(mat_X) <- sapply(0:order, FUN = function(i) paste("x", i, sep = ''))

  mat_coeffs <- matrix(0, k, k)

  for(i in 1:k) {
    for(h in 1:k) {
      col_sum <- 0
      for(j in 1:n) {
        col_sum <-  col_sum + mat_X[j,i]*mat_X[j,h]
      }
      mat_coeffs[i,h] <- col_sum
    }
  }

  mat_soln <- numeric(k)

  for(h in 1:k) {
    col_sum <- 0
    for(j in 1:n) {
      col_sum <-  col_sum + mat_X[j,h]*Y[j]
    }
    mat_soln[h] <- col_sum
  }

  vars <- solve(mat_coeffs, mat_soln)
  Y_hat <- poly_fit(X, vars)
  resids <- Y - Y_hat

  list(
    data = list(predictor=X, response=Y),
    coeffs = vars,
    residuals = resids,
    R2 = compute_r2(Y, resids)
  )

}

#' @title poly_plot_fit
#' @description Generates a plot of model data and the best fit regression curve determined by the model.
#' @author Tripp Bishop
#' @export poly_plot_fit
#' @param model A polynomial_lm model object.
#' @return A ggplot2 object containing a plot of the model fit.
poly_plot_fit <- function(model) {

  X_bounds <- range(model$data$predictor)
  increments <- (X_bounds[2] - X_bounds[1])/100 # generate 100 points to render the best fit curve

  X_fit <- seq(from=X_bounds[1], to=X_bounds[2], by=increments)
  Y_fit <- poly_fit(X_fit, model$coeffs)
  df_fit <- tibble(x=X_fit, y=Y_fit)

  tibble(predictor=model$data$predictor, response=model$data$response) %>%
    ggplot(aes(x=predictor, y=response)) +
    geom_point() +
    geom_line(data=df_fit, aes(x=x, y=y), colour="#99CC99", size=1.25)
}

#' @title poly_plot_residuals
#' @description Generates a scatter plot of the residuals generated by the specified create_lm model object.
#' @author Tripp Bishop
#' @export poly_plot_residuals
#' @param model A polynomial_lm model object.
#' @return A ggplot2 object containing a plot of the model fit.
poly_plot_residuals <- function(model) {

  df_resids <- tibble(predictor=model$data$predictor, resids=model$residuals)

  df_resids %>%
    ggplot(aes(x=predictor, y=resids)) +
      geom_point() +
      geom_hline(yintercept = 0, colour="#aa0099", size=1)
}

