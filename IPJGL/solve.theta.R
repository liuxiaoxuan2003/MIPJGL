solve.theta1 <- function(A, V1, V2, rho, n, lambda2, normtype)
{
  source('./IPJGL/Dnorm.R')
  # print('lambda2:')
  # print(lambda2)
  p <- length(V1[1, ])
  B <- 2 * lambda2 * Dnorm(V1, normtype) + 2 * lambda2 * Dnorm(V2, normtype) + 3 * rho * diag(p)

  iB <- solve(B)
  C <- A %*% iB
  results <- eigen(C, symmetric = TRUE)

  U <- results$vectors
  D <- results$values
  D <- -D + sqrt(D ^ 2 + 4 * n * diag(iB))

  theta <- 0.5 * U %*% diag(D) %*% t(U)
  return(theta)
}


solve.theta <- function(A, V, rho, n, lambda2, normtype)
{
  source('./IPJGL/Dnorm.R')
  p <- length(V[1, ])
  B <- 2 * lambda2 * Dnorm(V, normtype) + 2 * rho * diag(p)

  iB <- solve(B)
  C <- A %*% iB
  results <- eigen(C, symmetric = TRUE)

  U <- results$vectors
  D <- results$values
  D <- -D + sqrt(D ^ 2 + 4 * n * diag(iB))

  theta <- 0.5 * U %*% diag(D) %*% t(U)
  return(theta)
}
