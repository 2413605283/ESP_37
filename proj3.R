library(splines)


evaluate_matrices <- function(K) {
  
  knots <- seq(1, 80, length.out = K + 4)  
  
  ## outer.ok = TRUE?
  time_points <- 1:80
  X_tilde <- splineDesign(knots, time_points, outer.ok = TRUE)  
  
  
  
  
  
  
  
  d <- 1:80  
  edur <- 3.151  
  sdur <- 0.469  
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)
  
  
  X <- X_tilde %*% diag(pd)  
  
  
  S <- crossprod(diff(diag(K), diff = 2))
  
  
  return(list(X_tilde = X_tilde, pd = pd, X = X, S = S))
}


result <- evaluate_matrices(80)


print(result)
