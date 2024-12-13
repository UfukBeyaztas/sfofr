# Least squares estimation for the MSAR model
# The following function is a modified version of the LSE_estimation.R
# obtained from https://github.com/XueningZhu/MSAR_code
# Please check the following paper:
# X. Zhu, D. Huang, R. Pan, H. Wang (202), Multivariate spatial autoregressive 
# model for large scale social networks, Journal of Econometrics, 215, 591-606.

msar_fun <- function(Y, X, W, rho, B, Sige, tol = 10^{-3})
  {

  n <- nrow(W)
  p <- ncol(Y)
  vecY <- as.vector(Y)
  ww <- crossprod(W)
  q <- nrow(B)
  
  Omee <- solve(Sige)
  Del1 <- 1
  n.iter <- 1
  ind_d = 1:p^2
  ind_b = (p^2+1):(p^2+p*q)
  ind_e0 = which(lower.tri(Sige, diag = T))
  obj0 = -1
  
  while (mean(abs(Del1))>tol & n.iter<=100) {

    if (any(abs(diag(rho))>1)) {
      rho <- Diagonal(n = p, x = runif(p, -0.5, 0.5))
      B <- Matrix(runif(p*q), nrow = q, ncol = p)
    }
    vecTheta <- c(as.vector(rho), as.vector(B))
    
    grad_hessian <- grahes(rho, Y, vecY, W, ww, n, p, X, B, Sige)
    
    hessian <- grad_hessian$hessian
    Del <- solve(hessian)%*%grad_hessian$grad
    if (any(abs(diag(rho))>1))
      vecTheta <- vecTheta - Del*0.1
    else
      vecTheta <- vecTheta - Del
    rho <- Matrix(matrix(vecTheta[ind_d], nrow = p))
    B <- Matrix(matrix(vecTheta[ind_b], nrow = q))
    E <- Y - W %*% Y %*% rho - X %*% B
    Sige <- t(E) %*% E/n
    Del1 <- abs(grad_hessian$obj-obj0)
    obj0 <- grad_hessian$obj
    n.iter <- n.iter+1
  }
  
  return(list(rho = rho, B = B, Sige = Sige, theta = vecTheta, iter = n.iter))
}
