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
