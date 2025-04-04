getPCA <- function(data, nbasis, gp)
  {
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  bs_basis <- create.bspline.basis(rangeval = c(gp[1], gp[p]), nbasis = nbasis)
  inp_mat <- inprod(bs_basis, bs_basis)
  sinp_mat <- sqrtm(inp_mat)
  evalbase = eval.basis(gp, bs_basis)
  fdobj <- fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  
  mean_coef <- apply(t(pcaobj$coefs), 2, mean)
  sdata <- scale(t(pcaobj$coefs), scale = FALSE)
  new.data <- sdata %*% sinp_mat
  dcov <- cov(new.data)
  d.eigen <- eigen(dcov)
  var_prop <- cumsum(d.eigen$values)/sum(d.eigen$values)
  ncomp <- which(var_prop>0.95)[1]
  loads <- d.eigen$vectors[,1:ncomp]
  PCs <- solve(sinp_mat) %*% loads
  colnames(PCs) = 1:ncomp
  for(i in 1:ncomp)
    colnames(PCs)[i] = paste("PC", i, sep = "")
  PCAcoef <- fd(PCs, bs_basis)
  mean_coef <- fd(as.vector(mean_coef), bs_basis)
  pcaobj2 <- pcaobj
  pcaobj2$coefs <- t(sdata)
  PCAscore <- inprod(pcaobj2, PCAcoef)
  colnames(PCAscore) = 1:ncomp
  for(i in 1:ncomp)
    colnames(PCAscore)[i] = paste("Score", i, sep = "")
  
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, meanScore = mean_coef, ncomp = ncomp,
              bs_basis = bs_basis, evalbase = evalbase, gp = gp))
}
