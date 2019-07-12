#' Make predictions from an 'asgl' object
#'
#' Makes predictions from a fitted 'asgl' object.
#'
#' @param object Fitted \code{asgl} object
#' @param newx Matrix of new values for \code{x} at which predictions are to
#'        be made.
#' @param s Value(s) of the penalty tuning parameter at which predictions are to
#'        be made. Default is the entire sequence used to fit the regularization
#'        path.
#' @param type Type of predicted required. For \code{"gaussian"} models,
#'        \code{"link"} and \code{"response"} give the same result: the
#'        predicted response. For \code{"binomial"} models, \code{"link"} gives
#'        the linear predictor and \code{"response"} gives the predicted
#'        probabilities.
#' @param ... Ignored.
#'
#' @return A matrix of predictions or, if \code{s} is of length 1, a vector of
#'         predictions.
#'
#' @examples
#' # linear regression
#' n <- 500; p <- 20; groupsize <- 5
#' index <- ceiling(1:p / groupsize)
#' beta <- (-2:2)
#' x <- matrix(rnorm(n * p), ncol = p, nrow = n)
#' y <- as.vector(x[,1:5] %*% beta + 0.1 * rnorm(n))
#' fit1 <- asgl(x, y, index, family = "gaussian")
#' predict(fit1, x)
#'
#' # logistic regression
#' eta <- x[, 1:5] %*% beta
#' prob <- exp(eta) / (1 + exp(eta))
#' y <- rbinom(n, 1, prob)
#' fit2 <- asgl(x, y, index, family = "binomial")
#' predict(fit2, x, s = fit2$lambda[20], type = "response")
#'
#' @export
predict.asgl <- function(object, newx, s = NULL, type = c("link", "response"),
                         ...) {

  # Validate arguments
  type <- match.arg(type)
  if (missing(newx)) {
    stop("the argument 'newx' must be provided.", call. = FALSE)
  }
  if (!is.matrix(newx)) {
    stop("the argument 'newx' must be a matrix.", call. = FALSE)
  }
  if (!is.null(s)) {
    if (!all(s %in% object$lambda)) {
      stop("the argument 's' must be NULL or a subset of the tuning parameter ",
           "sequence used to fit the model.", call. = FALSE)
    }
  }

  # Make predictions
  if (object$family == "gaussian") {
    pred <- newx %*% object$beta
  } else {
    beta <- rbind(object$intercept, object$beta)
    pred <- cbind(1, newx) %*% beta
  }

  # Subset predictions according to selected tuning parameters, if applicable
  if (!is.null(s)) {
    ind <- which(object$lambda %in% s)
    pred <- pred[, ind]
    if (length(ind) == 1) {
      pred <- as.vector(pred)
    }
  }

  # Return linear predictor / response for Gaussian model
  if (type == "link" || object$family == "gaussian") {
    return(pred)
  }

  # Return response for binomial model
  if (type == "response") {
    if (object$family == "binomial") {
      prob <- exp(pred) / (1 + exp(pred))
    }
    return(prob)
  }
}
