get_lambda_sequence <- function(x, y, index, family, lambda_min, nlambda, alpha,
                            grp_weights, ind_weights) {

  # Setup
  groups <- unique(index)
  num_groups <- length(groups)
  range_group_ind <- as.numeric(c(0, cumsum(table(index))))
  group_length <- diff(range_group_ind)

  if (family == "gaussian") {
    resp <- y
  }
  if (family == "binomial") {
    m.y <- mean(y)
    resp <- y - m.y * (1 - m.y)
  }

  lambda_max <- 0

  # Loop through groups
  for (i in 1:num_groups) {

    # Omit unpenalized groups
    w <- grp_weights[i]
    if (w == 0) next

    # For each group coefficient, find the minimum value of lambda at which
    # the coefficient is set to zero.
    ind <- groups[i]
    x_g <- x[, which(index == ind)]
    i_g <- ind_weights[index == ind]
    cor <- t(x_g) %*% resp
    ord <- order(abs(cor) / (alpha * i_g), decreasing = TRUE)
    ord_cor <- abs(cor)[ord]
    ord_i_g <- i_g[ord]
    ord_lam <- ord_cor / (alpha * ord_i_g)

    # Find the minimum value of lambda which sets all group coefficients to zero
    if (length(ord_lam) == 1) {
      max_lam <- ord_lam
    } else {
      for (j in 2:length(ord)) {
        # Compute l2-norm of soft-threshold when lambda = lam[j]
        norm <- sqrt(sum((ord_cor[1:(j - 1)] -
                            ord_cor[j] * ord_i_g[1:(j-1)] / ord_i_g[j])^2))

        if (norm >= (1 - alpha) * ord_lam[j] * w * sqrt(group_length[i])) {

          # If inequality is not satisfied, minimum lambda is somewhere in
          # interval (lam[j-1], lam[j]): find solution using quadratic formula
          our_cor <- ord_cor[1:(j-1)]
          our_i_g <- ord_i_g[1:(j-1)]
          our_range <- ord_lam[c(j, j-1)]

          A <- alpha^2 * sum(our_i_g^2) - (1 - alpha)^2 * w^2 * group_length[i]
          B <- -2 * alpha * sum(our_cor * our_i_g)
          C <- sum(our_cor^2)

          solution <- c((-B - sqrt(B^2 - 4 * A * C)) / (2 * A),
                        (-B + sqrt(B^2 - 4 * A * C)) / (2 * A))
          break
        } else {
          if (j == length(ord)) {
            our_cor <- ord_cor
            our_i_g <- ord_i_g
            our_range <- c(0, ord_lam[j])

            A <- alpha^2 * sum(our_i_g^2) -
              (1 - alpha)^2 * w^2 * group_length[i]
            B <- -2 * alpha * sum(our_cor * our_i_g)
            C <- sum(our_cor^2)

            solution <- c((-B - sqrt(B^2 - 4 * A * C)) / (2 * A),
                          (-B + sqrt(B^2 - 4 * A * C)) / (2 * A))
            break
          }
        }
      }
      max_lam <- max(solution[which(solution >= our_range[1] &
                                    solution <= our_range[2])])
    }

    # Update maximum lambda value
    if (max_lam > lambda_max) {
      lambda_max <- max_lam
    }
  }

  # Finally, obtain lambda sequence
  max_lam <- lambda_max
  min_lam <- lambda_min * max_lam
  lambdas <- exp(seq(log(max_lam), log(min_lam),
                    (log(min_lam) - log(max_lam)) / (nlambda - 1)))
  lambdas / nrow(x)
}
