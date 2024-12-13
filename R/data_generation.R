data_generation <- function(n, nphi = 10, gpy=NULL, gpx=NULL, W.type = c("ivw","expw"),
                            rf = 0.9, sd.error=0.01, tol = 0.001, max_iter = 1000)
{

  if(!W.type %in% c("ivw","expw"))
    stop("Weight matrix type must be one of the followings: ivW and expw!")
  
  W.type <- match.arg(W.type)
  
  if(is.null(gpy))
    gpy <- seq(0, 1, length.out = 101)
  if(is.null(gpx))
    gpx <- seq(0, 1, length.out = 101)
  
  X.phi1 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1:nphi)
      phi[j, ] <- (j^-1.5) * sqrt(2) * cos(j * pi * gpx)
    return(phi)
  }
  
  X.phi2 <- function(nphi, gpx) {
    phi <- matrix(0, nphi, length(gpx))
    for (j in 1:nphi)
      phi[j, ] <- (j^-1.5) * sqrt(2) * sin(j * pi * gpx)
    return(phi)
  }
  
  rX.s <- function(nphi, gpx) {
    xsi <- rnorm(2 * nphi, 0, 1)
    X <- xsi * rbind(X.phi1(nphi, gpx), X.phi2(nphi, gpx))
    Xs <- colSums(X)
    return(Xs)
  }
  
  if(W.type == "ivw") {
    generate_W <- function(n) {
      wei <- matrix(, n, n)
      for(i in 1:n){
        for(j in 1:n){
          if(i != j){
            wei[i,j] <- 1/(1 + abs(i-j))
          }else{
            wei[i,j] <- 0
          }
        }
      }
      W <- matrix(0, n, n)
      for(i in 1:n)
        W[i,] <- wei[i,] / sum(wei[i,])
      return(W)
    }
  }else if(W.type == "expw") {
    generate_W <- function(n) {
      d <- 0.5
      wei <- matrix(0, n, n)
      
      for (i in 1:n) {
        for (j in 1:n) {
          if (i != j) {
            wei[i, j] <- exp(-d * abs(i - j))
          }
        }
      }
      
      W <- matrix(0, n, n)
      for (i in 1:n)
        W[i, ] <- wei[i, ] / sum(wei[i, ])
      return(W)
    }
  }
  
  beta <- function(s, t) 2 + s + t + 0.5 * sin(2 * pi * s * t)
  
  rho <- function(u, t) {
    rf * (1 + u * t) / (1 + abs(u - t))
  }
  
  X <- t(replicate(n, rX.s(nphi, gpx)))
  
  W <- generate_W(n)
  
  beta.ts <- outer(gpx, gpy, beta)
  Xbeta.t <- X %*% (beta.ts) * (gpx[2] - gpx[1])
  eps <- matrix(rnorm(n * length(gpy), 0, sd.error), n, length(gpy))
  G_t <- Xbeta.t + eps
  rho.ut <- outer(gpy, gpy, rho)
  
  fn <- function(f) as.matrix(W %*% f %*% rho.ut * (gpy[2] - gpy[1]))
  
  Y_0 <- G_t
  Y_1 <- G_t + fn(Y_0)
  diff <- max(abs(Y_0 - Y_1))
  Y_0 <- Y_1
  L <- 1
  
  while(diff > tol & L < max_iter) {
    Y_1 <- G_t + fn(Y_0)
    diff <- max(abs(Y_0 - Y_1))
    Y_0 <- Y_1
    L <- L + 1
  }
  
  Y <-  Y_1
  Y_true <- Y_1 - eps 
  
  return(list(Y = Y, Y_true = Y_true, X = X, W = W, rho = rho.ut, beta = beta.ts))
}
