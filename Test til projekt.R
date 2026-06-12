Ca2_final

library(expm)
h <- 0.02
#theta=c(eta,a,b,c,d,alpha,beta,gamma)
theta <- c(1,2,3,4,5,6,7,8) #Data fra side 2 fra Predrags papir

A <- function(theta) {
  result <- matrix(data=c(0,theta[4],
                          1,-theta[1]), nrow = 2,ncol = 2)
  return(result)
}


A(theta)

b <- function(theta) {
  result <- rbind(-theta[5]/theta[4], 0) #Opmærksom på c\neq0
  return(result)
}
b(theta)


mu_h <- function(x, h, theta) {
  A_mat <- A(theta)
  b_vec <- b(theta)
  result <- expm(A_mat*h)%*%(x - b_vec) + b_vec
  return(result)
}


f <- function(x, h, theta) {
  result <- rbind(x[1],
                  x[2]+h*(theta[2]*x[1]^3+theta[3]*x[1]^2))
  return(result)
}

f(c(Ca2_final$X[1],Ca2_final$V[1]),h,theta)

A(theta)


#Kronecker product is sort of the same as tensor product
#M and I are defined on page 10 in Strang splitting estimator
M1 <- function(theta, A_kron) {
  alpha_mat <- matrix(0,4,4)
  alpha_mat[4,4] <- theta[6]
  # A_mat <- A(theta)
  # A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat) 
  M <- rbind(
    cbind(A_kron + alpha_mat, alpha_mat),
    cbind(matrix(0,4,4), A_kron)
  )
  return(M)
}

M2 <- function(theta, A_kron) {
  alpha_mat <- matrix(0,4,4)
  alpha_mat[4,4] <- theta[6]
  # A_mat <- A(theta)
  # A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat)
  M <- rbind(
    cbind(A_kron + alpha_mat, alpha_mat),
    cbind(matrix(0,4,4), kronecker(diag(2),(diag(2) + A_mat)))
  )
  return(M)
}

M3 <- function(theta, A_kron) {
  alpha_mat <- matrix(0,4,4)
  alpha_mat[4,4] <- theta[6]
  # A_mat <- A(theta)
  # A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat)
  M <- rbind(
    cbind(A_kron + alpha_mat, alpha_mat),
    cbind(matrix(0,4,4), kronecker((A_mat + diag(2)),diag(2)))
  )
  return(M)
}

M4 <- function(theta, A_kron) {
  alpha_mat <- matrix(0,4,4)
  alpha_mat[4,4] <- theta[6]
  # A_mat <- A(theta)
  beta_mat <- matrix(0,nrow = 4,ncol = 2)
  beta_mat[4,2] <- theta[7]
  # A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat)
  M <- rbind(
    cbind(A_kron + alpha_mat, beta_mat),
    cbind(matrix(0,nrow = 2,ncol = 4), A_mat)
  )
  return(M)
}

M5 <- function(theta, A_kron) {
  alpha_mat <- matrix(0,4,4)
  alpha_mat[4,4] <- theta[6]
  # A_mat <- A(theta)
  # A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat)
  M <- rbind(
    cbind(A_kron + alpha_mat, diag(4)),
    cbind(matrix(0,4,4), matrix(0,4,4))
  )
  return(M)
}


I1 <- function(theta,h,A_kron) {
  M1_mat <- expm(M1(theta,A_kron)*h)
  return(M1_mat[1:4,5:8])
}

I2 <- function(theta,h,A_kron) {
  M2_mat <- expm(M2(theta,A_kron)*h)
  return(M2_mat[1:4,5:8])
}

I3 <- function(theta,h,A_kron) {
  M3_mat <- expm(M3(theta,A_kron)*h)
  return(M3_mat[1:4,5:8])
}

I4 <- function(theta,h,A_kron) {
  M4_mat <- expm(M4(theta,A_kron)*h)
  return(M4_mat[1:4,5:6])
}

I5 <- function(theta,h,A_kron) {
  M5_mat <- expm(M5(theta,A_kron)*h)
  return(M5_mat[1:4,5:8])
}



I5compute <- I5(theta,h,A_kron) %*% c(0,0,0,theta[8])

Omega <- function(x,theta,h, A_kron, I5compute) { 
  B <- b(theta)
  
  result <- 
    I1(theta,h,A_kron) %*% as.vector( (x-B) %*% t(x-B) ) + 
    I2(theta,h,A_kron) %*% as.vector( (x-B) %*% t(B) ) +
    I3(theta,h,A_kron) %*% as.vector( B %*% t(x-B) ) +
    I4(theta,h,A_kron) %*% (x-B) +
    I5compute #I5(theta,h) %*% c(0,0,0,theta[8])
  
  result <- matrix(result, nrow = sqrt(length(result))) # Unvectorize
  return(result)
}


#theta = theta_test

#eta <- theta[1]
#a <- theta[2]
#alpha <- theta[6]
#beta <- theta[7]
#gamma <- theta[8]
  
#Original(should not use) 
loglik <- function(theta, h, S) {
  
  nll <- 0 #Negative loglikelihood
  #A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat)
  
  Z_tk = matrix(0, length(S[,1]), 2)
  
  for (i in 2:(length(S[,1])-2)){
    prev_data <- S[i-1,]
    curr_data <- S[i,]
    
    Z_tk[i,] <- f(curr_data, -h/2, theta) - mu_h(f(prev_data, h/2, theta), h, theta)
    omega_h <- Omega(f(prev_data, h/2, theta), theta, h)
    
    omega_34h <- Omega(f(prev_data, h/2, theta), theta, (3*h)/4) # To fix the way I have generated V, page 10 on "Strang splitting for parametric inference in second order stochast..."
    
    #omega_23h <- Omega(f(prev_data, h/2, theta), theta, (3*h)/2) # To fix the way I have generated V, page 10 on "Strang splitting for parametric inference in second order stochast..."
    
    
    logdet_omega_34h <- log(det(omega_34h)) # Eliminate bias
    
    
    nll <- nll + (4/3) * logdet_omega_34h + t(Z_tk[i,])%*%solve(omega_h)%*%Z_tk[i,]
    
    #nll <- nll + (2/3) * omega_23h[2,2] + (Z_tk[i,2]^2)/omega_h[2,2]
    
    
  }
  print(nll)
  return(nll)
  #return()
}
#theta_test

loglik(theta_test, h, cbind(Ca2_final$X, Ca2_final$V))




#theta_test

#Tried to optimize
loglik <- function(theta, h, S) {
  
  nll <- 0 #Negative loglikelihood
  
  A_mat <- A(theta)
  A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat)
  
  I5compute = I5(theta,h,A_kron) %*% c(0,0,0,theta[8])
  
  Z_tk = matrix(0, length(S[,1]), 2)
  
  for (i in 2:(length(S[,1])-2)){
    prev_data <- S[i-1,]
    curr_data <- S[i,]
    
    Z_tk[i,] <- f(curr_data, -h/2, theta) - mu_h(f(prev_data, h/2, theta), h, theta)
    omega_h <- Omega(f(prev_data, h/2, theta), theta, h, A_kron, I5compute)
    
    omega_34h <- Omega(f(prev_data, h/2, theta), theta, (3*h)/4, A_kron, I5compute) # To fix the way I have generated V, page 10 on "Strang splitting for parametric inference in second order stochast..."
    
    #omega_23h <- Omega(f(prev_data, h/2, theta), theta, (3*h)/2) # To fix the way I have generated V, page 10 on "Strang splitting for parametric inference in second order stochast..."
    
    
    logdet_omega_34h <- log(det(omega_34h)) # Eliminate bias
    
    
    nll <- nll + (4/3) * logdet_omega_34h + t(Z_tk[i,])%*%solve(omega_h)%*%Z_tk[i,]
    
    #nll <- nll + (2/3) * omega_23h[2,2] + (Z_tk[i,2]^2)/omega_h[2,2]
    
  }
  print(nll)
  return(nll)
}

loglik(theta_test, h, cbind(Ca2_final$X, Ca2_final$V))

#Another optimization
loglik <- function(theta, h, S) {
  
  nll <- 0 #Negative loglikelihood
  #A_kron <- kronecker(A_mat,diag(2)) + kronecker(diag(2),A_mat)
  
  Z_tk = matrix(0, length(S[,1]), 2)
  
  for (i in 2:(length(S[,1])-2)){
    prev_data <- S[i-1,]
    curr_data <- S[i,]
    
    Z_tk[i,] <- f(curr_data, -h/2, theta) - mu_h(f(prev_data, h/2, theta), h, theta)
    omega_h <- Omega(f(prev_data, h/2, theta), theta, h)
    
    omega_34h <- Omega(f(prev_data, h/2, theta), theta, (3*h)/4) # To fix the way I have generated V, page 10 on "Strang splitting for parametric inference in second order stochast..."
    
    #omega_23h <- Omega(f(prev_data, h/2, theta), theta, (3*h)/2) # To fix the way I have generated V, page 10 on "Strang splitting for parametric inference in second order stochast..."
    
    
    logdet_omega_34h <- log(det(omega_34h)) # Eliminate bias
    
    
    nll <- nll + (4/3) * logdet_omega_34h + t(Z_tk[i,])%*%solve(omega_h)%*%Z_tk[i,]
    
    #nll <- nll + (2/3) * omega_23h[2,2] + (Z_tk[i,2]^2)/omega_h[2,2]
    
    
  }
  print(nll)
  return(nll)
  #return()
}
#theta_test

loglik(theta_test, h, cbind(Ca2_final$X, Ca2_final$V))



theta_test
#theta_init <- c(47, -125,   23,  220,  32,  95,  -87, 1061)
#theta_init <- c(25, -172,   3,  340,  31,  86,  27, 438)
theta_init <- c(68, -98, 8, 67, -9, 64, 227, 4371)
eps <- 1e-10


result_lbfgsb <- optim(
  par = theta_init,
  fn = loglik,
  S = S,
  h = h,
  method = "L-BFGS-B",
  lower = c(0, -Inf, -Inf, -Inf, -Inf, 0, -Inf, 0),
  upper = c(Inf, 0, Inf, Inf, Inf, Inf, Inf, Inf),
  control = list(
    factr = 1e2,
    trace = 1,
    REPORT = 5,
    maxit = 500
  )
)

result_lbfgsb
print(result_lbfgsb)



theta_init <- C(47, -125,   23,  220,  32,  95,  -87, 1061)
# c(50, -200, 10, 100, 10, 30, -5, 1000)
S <- cbind(Ca2_final$X, Ca2_final$V)
h <- 0.02

result_NM <- optim(
  par = theta_init,
  fn = loglik,
  S = S,
  h = h,
  method = "Nelder-Mead",
  control = list(
    maxit = 200,
    reltol = 1e-3
  )
)

result_NM
print(result_NM)




