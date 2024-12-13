# Gradient and hessian matrix
# The following function is a modified version of the LSE_grad_hessian.R
# obtained from https://github.com/XueningZhu/MSAR_code
# Please check the following paper:
# X. Zhu, D. Huang, R. Pan, H. Wang (202), Multivariate spatial autoregressive 
# model for large scale social networks, Journal of Econometrics, 215, 591-606.

grahes <- function(rho, Y, vecY, W, ww, n, p, X, B, Sige)
  {
  
  # rho     : spatial autocorrolation matrix
  # Y       : response variable
  # vecY    : vectorized response variable
  # W       : weight matrix
  # ww      : W^TW
  # n       : sample size
  # p       : number of response variables
  # X       : predictor
  # B       : parameter matrix
  # Sige    : variance-covariance matrix
  
  ifelse <- function(x, a, b) {
    if (x) {
      return(a)
    }else {
      return(b)
    }
  }
  
  q <- ncol(X)
  tW = t(W)
  
  I <- Diagonal(n*p, x = 1)
  Ip <- Diagonal(p, x = 1)
  In <- Diagonal(n, x = 1)
  Inp <- Diagonal(n*p, x = 1)
  
  IX <- kronecker(Ip, X)
  IXbeta <- as.vector(X %*% B)
  
  Omee <- solve(Sige)
  IOmee <- kronecker(Omee, In)
  
  S <- I - kronecker(t(rho), W)
  tS <- t(S)
  SYX <- S %*% vecY - IXbeta
  OmeeS <- IOmee %*% S
  tSe <- t(OmeeS)
  OmeSYX <- tSe %*% SYX
  
  Ome <- tSe %*% S
  m <- 1/diag(Ome)
  
  MY <- m*OmeSYX
  
  dww <- diag(ww)
  vMY <- as.vector(MY)
  mtSe <- m*tSe
  
  right1 <- matrix(tSe %*% SYX, ncol = p)
  right2 <- matrix(SYX, ncol = p)
  right3 <- Y
  
  left1b <- matrix(m^2*vMY, ncol = p)
  left2b <- matrix(-m*vMY, ncol = p)
  
  F_d <- matrix(0, nrow = n*p, ncol = p^2)
  Q_db2 = Matrix(0, nrow = p^2, ncol = p*q)
  
  for (j1 in 1:p) {
    for (j2 in 1:p) {
      Ij2j1 <- Matrix(0, nrow = p, ncol = p)
      Ij2j1[j2, j1] <- 1
      Ij1j2 <- t(Ij2j1)
      
      Vj1j21 <- Ij1j2 %*% Omee %*% t(rho) + rho %*% Omee %*% Ij2j1

      F_d1j <- ifelse(p==1, -m^2*as.vector((dww*right1) %*% (Vj1j21)),
                     -m^2*as.vector((dww*right1) %*% diag(diag(Vj1j21))))
      
      F_d2j <- -m*as.vector(t(W) %*% right2 %*% t(Ij1j2 %*% Omee))
      F_d3j <- -mtSe %*% as.vector(W %*% right3 %*% t(Ij2j1))
      
      
      ii <- (j2-1)*p+j1
      F_d[,ii] <- F_d1j + F_d2j + F_d3j[,1]

      Q_db1j <- t(dww*X) %*% left1b %*% (diag(Vj1j21)*Omee)-
        t(dww*t(W) %*% X) %*% left1b %*% (diag(Vj1j21)*(rho %*% Omee))
      Q_db2j <- -t(tW %*% X) %*% left2b %*% t(Omee %*% Ij2j1)
      Q_db2[ii,] <- as.vector(Q_db1j + Q_db2j)*2
    }
  }
  
  F_b <- -m*tSe %*% IX 
  grad_D <- 2*colSums(vMY*F_d)
  grad_beta <- 2*colSums(vMY*F_b)
  grad <- c(grad_D, grad_beta)
  
  left11 <- matrix(2*m^3*vMY, ncol = p)
  left12 <- matrix(-m^2*vMY, ncol = p)
  right1 <- matrix(OmeSYX, ncol = p)
  
  left2 <- matrix(-m^2*vMY, ncol = p)
  right21 <- Y
  right22 <- X %*% B
  
  left4 <- matrix(m*vMY, ncol = p)
  right4 <- Y
  
  Q_dd2 <- matrix(0, p^2, p^2)
  
  for (j1 in 1:p) {
    for (j2 in 1:p) {
      for (k1 in 1:p) {
        for (k2 in 1:p) {
          jj <- (j2-1)*p+j1
          kk <- (k2-1)*p+k1
          if (jj <= kk) {
            ii <- (kk-1)*p^2 + jj
            
            Ij2j1 <- Matrix(0, nrow = p, ncol = p)
            Ij2j1[j2, j1] <- 1
            Ij1j2 <- t(Ij2j1)
            Ik2k1 <- Matrix(0, nrow = p, ncol = p)
            Ik2k1[k2, k1] <- 1
            Ik1k2 <- t(Ik2k1)
            IIk2k1 <- kronecker(Ik2k1, In)
            
            Vj1j21 <- Ij1j2 %*% Omee %*% t(rho) + rho %*% Omee %*% Ij2j1
            Vk1k21 <- Ik1k2 %*% Omee %*% t(rho) + rho %*% Omee %*% Ik2k1
            Vjk <- Ij1j2 %*% Omee %*% Ik2k1 + Ik1k2 %*% Omee %*% Ij2j1

            right11jk <- ifelse(p==1, (dww^2*right1) %*% (Vj1j21)*diag(Vk1k21),
                               (dww^2*right1) %*% diag(diag(Vj1j21)*diag(Vk1k21)))
            right12jk <- ifelse(p==1, (dww*right1) %*% (Vjk), (dww*right1) %*% diag(diag(Vjk)))
            
            right211jk <- -(dww*tW) %*% right21 %*% t(diag(Vj1j21)*(Ik1k2 %*% Omee))
            right212jk <- -(dww*W) %*% right21 %*% t(diag(Vj1j21)*(Omee %*% Ik2k1))
            right213jk <- (dww*ww) %*% right21 %*% t(diag(Vj1j21)*Vk1k21)
            right22jk <- (dww*tW) %*% right22 %*% t(diag(Vj1j21)*t(Omee %*% Ik2k1))
            
            right211kj <- -(dww*tW) %*% right21 %*% t(diag(Vk1k21)*(Ij1j2 %*% Omee))
            right212kj <- -(dww*W) %*% right21 %*% t(diag(Vk1k21)*(Omee %*% Ij2j1))
            right213kj <- (dww*ww) %*% right21 %*% t(diag(Vk1k21)*Vj1j21)
            right22kj <- (dww*tW) %*% right22 %*% t(diag(Vk1k21)*t(Omee %*% Ij2j1))
            
            right4jk <- ww %*% right4 %*% t(Vjk)
            
            Q_dd2[jj,kk] <- sum(left11*right11jk)+sum(left12*right12jk)+
              sum(left2*right211jk) + sum(left2*right212jk) + sum(left2*right213jk)+ sum(left2*right22jk)+
              sum(left2*right211kj) + sum(left2*right212kj) + sum(left2*right213kj)+ sum(left2*right22kj)+
              sum(left4*right4jk)
          }
        }
      }
    }
  }
  
  Q_dd2[lower.tri(Q_dd2)] <- t(Q_dd2)[lower.tri(Q_dd2)]
  Q_dd <- 2*t(F_d) %*% F_d + Q_dd2*2
  Q_db <- 2*t(F_d) %*% F_b+Q_db2
  Q_bb <- 2*t(F_b) %*% F_b
  
  Q_h <- Matrix(0, nrow = p^2+p*q+p^2, ncol = p^2+p*q+p^2)
  ind_d <- 1:p^2
  ind_b <- (p^2+1):(p^2+p*q)
  
  Q_h[ind_d, ind_d] <- Q_dd
  Q_h[ind_b, ind_b] <- Q_bb
  
  Q_h[ind_d, ind_b] <- Q_db
  Q_h[ind_b, ind_d] <- t(Q_db)
  
  ind0 <- c(ind_d, ind_b)
  Q_h <- Q_h[ind0, ind0]
  
  return(list(grad = grad, hessian = Q_h, obj = mean(vMY^2)))
}
