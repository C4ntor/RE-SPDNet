simulate_wishart_process <- function(time_grid, d = 3, v = 5, lengthscale = 1.0, scale_matrix = NULL) {
  # time_grid: vector of time points, e.g., seq(0, T, by = dt)
  # d: dimension of covariance matrices (number of variables)
  # v: degrees of freedom (number of latent Gaussian processes)
  # lengthscale: parameter for the RBF kernel
  # scale_matrix: d x d positive-definite matrix; if NULL, identity is used
  
  T_len <- length(time_grid)
  
  # Set scale matrix to identity if not provided
  if (is.null(scale_matrix)) {
    scale_matrix <- diag(d)
  }
  
  # Compute the time covariance matrix using an RBF kernel
  time_cov <- outer(time_grid, time_grid, function(t, tprime) exp(- ( (t - tprime)^2 ) / (2 * lengthscale^2)))
  
  # Preallocate an array to store the latent processes:
  # Dimensions: d (variable dimension) x v (number of latent processes) x T_len (time points)
  latent_array <- array(0, dim = c(d, v, T_len))
  
  # For each latent process (and for each dimension), simulate a realization over the time grid
  for (i in 1:v) {
    for (j in 1:d) {
      # Simulate a sample from a multivariate normal with mean zero and covariance time_cov
      latent_array[j, i, ] <- mvrnorm(n = 1, mu = rep(0, T_len), Sigma = time_cov)
    }
  }
  
  # For each time point, construct the d x v matrix X_t and compute the Wishart matrix
  Sigma <- array(0, dim = c(d, d, T_len))
  # Obtain the Cholesky factor L of the scale_matrix (ensuring positive-definiteness)
  L <- t(chol(scale_matrix))
  
  for (t in 1:T_len) {
    # Extract the d x v matrix for time t
    X_t <- latent_array[ , , t]
    # Compute the sum of outer products (this gives a Wishart-distributed matrix)
    W_t <- X_t %*% t(X_t)
    # Scale with the Cholesky factor to obtain the dynamic covariance matrix
    Sigma_t <- L %*% W_t %*% t(L)
    Sigma[,,t] <- Sigma_t
  }
  
  return(Sigma)
}

Cholmat <- function (X, tol = sqrt(.Machine$double.eps)) {
  if (!is.numeric(X)) 
    stop("argument is not numeric")
  if (!is.matrix(X)) 
    stop("argument is not a matrix")
  n <- nrow(X)
  if (ncol(X) != n) 
    stop("matrix is not square")
  if (max(abs(X - t(X))) > tol) 
    stop("matrix is not symmetric")
  D <- rep(0, n)
  L <- diag(n)
  i <- 2:n
  D[1] <- X[1, 1]
  if (abs(D[1]) < tol) 
    stop("matrix is numerically singular")
  L[i, 1] <- X[i, 1]/D[1]
  for (j in 2:(n - 1)) {
    k <- 1:(j - 1)
    D[j] <- abs(X[j, j] - sum((L[j, k]^2) * D[k]))+0.00000000001
    if (abs(D[j]) < tol) 
      stop("matrix is numerically singular")
    i <- (j + 1):n
    L[i, j] <- (X[i, j] - colSums(L[j, k] * t(L[i, k, drop = FALSE]) * 
                                    D[k]))/D[j]
  }
  k <- 1:(n - 1)
  D[n] <- abs(X[n, n] - sum((L[n, k]^2) * D[k]))
  if (abs(D[n]) < tol) 
    stop("matrix is numerically singular")
  (L %*% diag(sqrt(D)))
}
isPSD <- function(x, tol = 0.00000001) {
  if (!is.square.matrix(x)) 
    stop("argument x is not a square matrix")
  if (!isSymmetric(x)) 
    stop("argument x is not a symmetric matrix")
  if (!is.numeric(x)) 
    stop("argument x is not a numeric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow(x)
  for (i in 1:n) {
    if (abs(eigenvalues[i]) < tol) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues < 0)) {
    return(FALSE)
  }
  return(TRUE)
}
approximant <- function(A){
  B <- (A+t(A))/2
  C <- (A-t(A))/2
  BB <- tensr::polar(B)
  H <- BB$Z
  speB <- eigen(B)
  L <- speB$values
  Z <- speB$vectors
  L[L<0] <- 0
  Xf <- Z %*% diag(L) %*% t(Z)
  return(Xf)
}