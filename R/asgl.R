#' Fit a GLM with Adaptive Sparse Group Lasso Penalty
#'
#' Fit a generalized linear model via penalized maximum likelihood. The
#' regularization path is computed for the adaptive sparse group lasso penalty
#' along a sequence of tuning parameter values. Supports linear and logistic
#' regression.
#'
#' @param x Input matrix of dimensions \emph{n} by \emph{p}.
#' @param y Response vector of length \emph{n}. For linear models
#'        (\code{family = "gaussian"}), a numeric vector; for logistic models
#'        (\code{family = "binomial"}), a categorical factor with two levels.
#' @param index A vector of length \emph{p} defining group membership for each
#'        covariate.
#' @param family Response type. One of \code{"gaussian"} or \code{"binomial"}.
#' @param offset Optional. A vector of length \emph{n} to be included in the
#'        linear predictor
#' @param alpha Penalty mixing parameter. The penalty reduces to the adaptive
#'        group lasso if \code{alpha = 0} and to the adaptive lasso if
#'        \code{alpha = 0}.
#' @param lambda Optional. A vector of user-specified penalty tuning parameter
#'        values at which to fit model. If \code{NULL}, the sequence is chosen
#'        automatically (recommended).
#' @param lambda_min Smallest tuning parameter value, as a fraction of the
#'        largest parameter value in the sequence (the smallest value for which
#'        all coefficients are zero).
#' @param nlambda The number of tuning parameter values in the sequence.
#' @param maxit Maximum number of iterations until convergence
#' @param thresh Convergence threshold for change in model coefficient values.
#' @param gamma Backtracking parameter for use in gradient descent step
#' @param step Initial gradient descent step size.
#' @param standardize If \code{TRUE}, variables are standardized prior to model
#'        fitting.
#' @param grp_weights A vector of weights for each group of coefficients, the
#'        same length as \code{index}.
#' @param ind_weights A vector of weights for each indiviual coefficient, of
#'        length \emph{p}.
#'
#' @return An object of class "\code{asgl}"
#'
#' @examples
#' # linear regression
#' n <- 500; p <- 20; groupsize <- 5
#' index <- ceiling(1:p / groupsize)
#' beta <- (-2:2)
#' x <- matrix(rnorm(n * p), ncol = p, nrow = n)
#' y <- as.vector(x[,1:5] %*% beta + 0.1 * rnorm(n))
#' asgl(x, y, index, family = "gaussian")
#'
#' # logistic regression
#' eta <- x[, 1:5] %*% beta
#' prob <- exp(eta) / (1 + exp(eta))
#' y <- rbinom(n, 1, prob)
#' asgl(x, y, index, family = "binomial")
#'
#' # adaptive weights
#' coefs <- glm.fit(x, y, family = binomial())$coefficients
#' ind_weights <- 1 / abs(coefs)
#' grp_weights <- numeric(length(unique(index)))
#' for (i in unique(index)) {
#'   grp_weights[i] <-  1 / sqrt(sum(coefs[which(index == i)]^2))
#' }
#' asgl(x, y, index, family = "binomial",
#'      grp_weights = grp_weights, ind_weights = ind_weights)
#'
#' @useDynLib asgl, .registration = TRUE
#' @export
asgl <- function(x, y, index, family = c("gaussian", "binomial"),
                 offset = NULL, alpha = 0.95, lambda = NULL, lambda_min = 0.1,
                 nlambda = 20, maxit = 1000, thresh = 0.001, gamma = 0.8,
                 step = 1, standardize = FALSE,
                 grp_weights = NULL, ind_weights = NULL){

  # Validate input matrix
  if (!is.matrix(x)) {
    stop("the argument 'x' must be a matrix.", call. = FALSE)
  }
  dimx <- dim(x)
  if (is.null(dimx) || dimx[2] < 2) {
    stop("the argument 'x' must be a matrix with 2 or more columns.",
         call. = FALSE)
  }
  nobs <- dimx[1]
  nvar <- dimx[2]

  # Validate response vector
  if (!is.vector(y)) {
    stop("the argument 'y' must be a vector.", call. = FALSE)
  }
  leny <- length(y)
  if (leny != nobs) {
    stop(paste("the length of 'y' (", leny, ") is not equal to the number of ",
               "rows of 'x' (", nobs, ").", sep = ""), call. = FALSE)
  }

  # Validate index vector
  if (!is.vector(index)) {
    stop("the argument 'index' must be a vector.", call. = FALSE)
  }
  leni <- length(index)
  if (leni != nvar) {
    stop(paste("the length of 'index' (", leni, ") is not equal to the number ",
               "of columns of 'x' (", nvar, ").", sep = ""), call. = FALSE)
  }

  # Validate model offset
  if (is.null(offset)) {
    offset <- rep.int(0, leny)
  }
  if (!is.vector(offset)) {
    stop("the argument 'offset' must be a vector.", call. = FALSE)
  }
  leno <- length(offset)
  if (leno != nobs) {
    stop(paste("the length of 'offset' (", leno, ") is not equal to the ",
               "number of rows of 'x' (", nobs, ").", sep = ""), call. = FALSE)
  }

  # Validate adaptive weights
  lenui <- length(unique(index))
  if (is.null(grp_weights)) {
    grp_weights <- rep.int(1, lenui)
  }
  if (!is.vector(grp_weights)) {
    stop("the argument 'grp_weights' must be a vector.", call. = FALSE)
  }
  if (is.null(ind_weights)) {
    ind_weights <- rep.int(1, nvar)
  }
  if (!is.vector(ind_weights)) {
    stop("the argument 'ind_weights' must be a vector.", call. = FALSE)
  }
  lengw <- length(grp_weights)
  if (lengw != lenui) {
    stop(paste("the length of 'grp_weights' (", lengw, ") is not equal to the ",
               "number of unique elements of 'index' (", lenui, ").", sep = ""),
         call. = FALSE)
  }
  leniw <- length(ind_weights)
  if (leniw != nvar) {
    stop(paste("the length of 'ind_weights' (", leniw, ") is not equal to the ",
               "number of columns of 'x' (", nvar, ").", sep = ""),
         call. = FALSE)
  }

  # Validate model family
  family <- match.arg(family)

  # Validate additional arguments
  if (alpha < 0 || alpha > 1) {
    stop("the argument 'alpha' must be between 0 and 1.")
  }

  # ----
  # Standardize input matrix, if applicable
  X.transform <- NULL
  if (standardize) {
    x <- scale(x, center = TRUE,
               scale = apply(x, 2, function(v) sqrt(sum(v^2)/length(v))))
    X.transform <- list(X.scale = attributes(x)$`scaled:scale`,
                        X.means = attributes(x)$`scaled:center`)
  }
  # ----

  # Prepare tuning parameter sequence
  if (is.null(lambda)) {
    lambda <- get_lambda_sequence(x, y, index, family, lambda_min, nlambda,
                                  alpha, grp_weights, ind_weights)
  } else {
    nlambda <- length(lambda)
  }

  # ---- FIT THE MODEL

  ord <- order(index)
  index <- index[ord]
  x <- x[, ord]
  unOrd <- match(1:length(ord), ord)

  groups <- unique(index)
  num.groups <- length(groups)
  group.length <- as.numeric(table(index))
  range.group.ind <- c(0, cumsum(group.length))

  allbeta <- rep(0, nvar * nlambda)
  beta <- rep(0, nvar)
  range.lambda.ind <- seq(0, nvar * nlambda, nvar)
  eta <- rep(0, nobs) + offset

  if (family == "gaussian") {

    result <- .C("fit_gaussian",
              alpha = as.double(alpha),
              allbeta = as.double(allbeta),
              beta = as.double(beta),
              rangeLambdaInd = as.integer(range.lambda.ind),
              X = as.double(as.vector(x)), y = as.double(y),
              nrow = as.integer(nobs), ncol = as.integer(nvar),
              numGroup = as.integer(num.groups),
              rangeGroupInd = as.integer(range.group.ind),
              groupLen = as.integer(group.length),
              lambda = as.double(lambda),
              lambda1 = as.double(0), lambda2 = as.double(0),
              nlam = as.integer(nlambda),
              innerIter = as.integer(maxit), outerIter = as.integer(maxit),
              thresh = as.double(thresh), outerThresh = as.double(thresh),
              eta = as.double(eta), gamma = as.double(gamma),
              betaIsZero = as.integer(rep(1, num.groups)),
              step = as.double(step),
              grpWeights = as.double(grp_weights),
              indWeights = as.double(ind_weights))
  }

  if (family == "binomial") {

    intercepts <- rep(log(sum(y)) - log(nobs - sum(y)), nlambda)
    eta <- eta + intercepts[1] + offset

    result <- .C("fit_binomial",
              alpha = as.double(alpha),
              allbeta = as.double(allbeta),
              beta = as.double(beta),
              rangeLambdaInd = as.integer(range.lambda.ind),
              X = as.double(as.vector(x)), y = as.integer(y),
              nrow = as.integer(nobs), ncol = as.integer(nvar),
              numGroup = as.integer(num.groups),
              rangeGroupInd = as.integer(range.group.ind),
              groupLen = as.integer(group.length),
              lambda = as.double(lambda),
              lambda1 = as.double(0), lambda2 = as.double(0),
              nlam = as.integer(nlambda),
              innerIter = as.integer(maxit), outerIter = as.integer(maxit),
              thresh = as.double(thresh), outerThresh = as.double(thresh),
              eta = as.double(eta), gamma = as.double(gamma),
              betaIsZero = as.integer(rep(1, num.groups)),
              betaZero = as.double(intercepts[1]), step = as.double(step),
              grpWeights = as.double(grp_weights),
              indWeights = as.double(ind_weights),
              intercepts = as.double(intercepts))
  }

  result <- list(beta = result$allbeta, lambda = result$lambda,
                 intercept = result$intercepts, X.transform = X.transform)
  result$beta <- matrix(result$beta, nrow = nvar, ncol = nlambda)
  result$beta <- result$beta[unOrd, ]
  class(result) <- "asgl"
  result
}
