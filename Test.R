SS_precompute_matrices_R <- function(theta, drift_const, h) {
  eta <- theta[1]
  alpha <- theta[6]
  beta <- theta[7]
  gamma <- theta[8]
  
  A_ss <- matrix(c(0, drift_const,
                   1, -eta),
                 nrow = 2, ncol = 2)
  I2 <- diag(2)
  A_kron <- kronecker(A_ss, I2) + kronecker(I2, A_ss)
  
  alpha_mat <- matrix(0, nrow = 4, ncol = 4)
  alpha_mat[4, 4] <- alpha
  
  block_exp <- function(TL, TR, BR) {
    n <- nrow(TL)
    P <- matrix(0, nrow = 2 * n, ncol = 2 * n)
    P[1:n, 1:n] <- TL
    P[1:n, (n + 1):(2 * n)] <- TR
    P[(n + 1):(2 * n), (n + 1):(2 * n)] <- BR
    expm(h * P)[1:n, (n + 1):(2 * n), drop = FALSE]
  }
  
  I1 <- block_exp(A_kron + alpha_mat, alpha_mat, A_kron)
  I3 <- block_exp(
    A_kron + alpha_mat,
    alpha_mat,
    kronecker(A_ss, I2) + diag(4)
  )
  I5 <- block_exp(A_kron + alpha_mat, diag(4), matrix(0, 4, 4))
  I5_G5 <- as.numeric(I5 %*% c(0, 0, 0, gamma))
  
  P4 <- matrix(0, nrow = 6, ncol = 6)
  P4[1:4, 1:4] <- A_kron + alpha_mat
  P4[4, 6] <- beta
  P4[5:6, 5:6] <- A_ss
  I4 <- expm(h * P4)[1:4, 5:6, drop = FALSE]
  
  list(
    A = A_ss,
    exp_Ah = expm(h * A_ss),
    I1 = I1,
    I3 = I3,
    I4 = I4,
    I5_G5 = I5_G5
  )
}

SS_step_omega_R <- function(Y_mid, kappa, mats) {
  x0 <- Y_mid[1] - kappa
  x1 <- Y_mid[2]
  
  x_outer <- c(x0 * x0, x0 * x1, x1 * x0, x1 * x1)
  xxs <- c(x0 * kappa, 0, x1 * kappa, 0)
  
  omega_vec <- as.numeric(
    mats$I1 %*% x_outer +
      2 * (mats$I3 %*% xxs) +
      mats$I4 %*% c(x0, x1) +
      mats$I5_G5
  )
  
  matrix(omega_vec, nrow = 2, ncol = 2, byrow = TRUE)
}

SS_f_step_R <- function(Y, h_step, theta, drift_const, kappa) {
  a <- theta[2]
  b <- theta[3]
  c <- theta[4]
  d <- theta[5]
  
  X <- Y[1]
  V <- Y[2]
  
  c(
    X,
    V + h_step * (a * X^3 + b * X^2 + c * X + d -
                    drift_const * (X - kappa))
  )
}

SS_step_R <- function(Y_old, Y_new, h, theta, drift_const,
                      kappa_pos, kappa_neg, mats) {
  kappa <- if (Y_old[1] > 0) kappa_pos else kappa_neg
  
  f_new <- SS_f_step_R(Y_new, -h / 2, theta, drift_const, kappa)
  Y_mid <- SS_f_step_R(Y_old,  h / 2, theta, drift_const, kappa)
  
  center <- c(kappa, 0)
  mu <- as.numeric(mats$exp_Ah %*% (Y_mid - center) + center)
  z <- f_new - mu
  
  omega22 <- SS_step_omega_R(Y_mid, kappa, mats)[2,2]
  #omega11 <- Omega[1, 1] # needs to be inverted
  #omega12 <- Omega[1, 2]
  #omega21 <- Omega[2, 1]
  #omega22 <- Omega[2, 2]
  
  Y_mid1 <- SS_f_step_R(Y_old,  3/2 * h / 2, theta, drift_const, kappa)
  
  omegacor22 <- SS_step_omega_R(Y_mid1, kappa, mats)[2, 2]
  #omegacor11 <- Omegacor[1, 1]
  #omegacor12 <- Omegacor[1, 2]
  #omegacor21 <- Omegacor[2, 1]
  #omegacor22 <- Omegacor[2, 2]
  
  #det_omegacor <- max(omegacor11 * omegacor22 - omegacor12 * omegacor21, 1e-8)
  quadratic <- (
    #omega22 * z[1]^2 -
    #(omega12 + omega21) * z[1] * z[2] +
    z[2]^2 / omega22 
  ) 
  
  0.5 * (2/3 * log(omegacor22) + quadratic)
}

SS_neg_log_lik_R <- function(theta, data, h) {
  X <- as.matrix(data[, c("X", "V")])
  X <- X[complete.cases(X), , drop = FALSE]
  
  if (nrow(X) < 2) {
    stop("Need at least two complete observations.")
  }
  
  data_old <- X[-nrow(X), , drop = FALSE]
  data_new <- X[-1, , drop = FALSE]
  
  a <- theta[2]
  b <- theta[3]
  c <- theta[4]
  
  m <- mean(data_old[, 1])
  m2 <- mean(data_old[, 1]^2)
  var_x <- var(data_old[, 1])
  
  drift_const <- 3 * a * m2 + 2 * b * m + c
  root <- sqrt(max((3 * a * m + b)^2 + 9 * a^2 * var_x, 0))
  kappa_pos <- (-b - root) / (3 * a)
  kappa_neg <- (-b + root) / (3 * a)
  
  mats <- SS_precompute_matrices_R(theta, drift_const, h)
  
  nlls <- vapply(
    seq_len(nrow(data_old)),
    function(i) {
      SS_step_R(
        Y_old = data_old[i, ],
        Y_new = data_new[i, ],
        h = h,
        theta = theta,
        drift_const = drift_const,
        kappa_pos = kappa_pos,
        kappa_neg = kappa_neg,
        mats = mats
      )
    },
    numeric(1)
  )
  
  sum(nlls)
}

SS_loglikelihood_R <- function(data, h, theta) {
  -SS_neg_log_lik_R(theta = theta, data = data, h = h)
}

# Same sign convention as pseudo_loglikelihood():
SS_loglikelihood_R(Ca2_final, h = 0.02, theta_test)
```