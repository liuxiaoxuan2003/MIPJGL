theta2partial <- function(theta) {
  d <- solve(sqrt(diag(diag(theta))))
  partial <- -d %*% theta %*% d
  diag(partial)<-1
  return(partial)
}
# t = results.IPJGL$Z_arr[[i-1]]
# d <- solve(sqrt(diag(diag(t))))
# partial <- -d %*% t %*% d
# diag(partial)<-1

