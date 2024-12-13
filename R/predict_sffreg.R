predict_sffreg <- function(object, xnew, Wnew=NULL)
  {

  Wnew <- normalize_weights(Wnew)
  n <- dim(xnew)[1]
  fpcay <- object$model_info$fpcay
  fpcax <- object$model_info$fpcax
  nbasis <- object$model_info$nbasis
  model.type <- object$model_info$model.type
  scoX.test <- getPCA_test(fpcax, xnew)
  cfb <- object$model_info$cfb
  cfr <- object$model_info$cfr
  scoY <- object$model_info$fpcay$PCAscore

  if(model.type == "independent") {
    predicted.values <- ((scoX.test %*% cfb %*% t(fpcay$PCAcoef$coefs)) +
                           matrix(fpcay$meanScore$coefs, nrow = n, ncol = nbasis, byrow = T)) %*%
      t(fpcay$evalbase)
  }else if(model.type == "spatial") {
    phatvec <- solve(diag(n*dim(scoY)[2]) - t(cfr) %x% Wnew) %*%
      (diag(dim(scoY)[2]) %x% scoX.test) %*% as.vector(cfb)
    phat1 <- matrix(phatvec, ncol = dim(scoY)[2]) %*% t(fpcay$PCAcoef$coefs) +
      matrix(fpcay$meanScore$coefs, nrow = n, ncol = nbasis, byrow = T)
    predicted.values <- phat1 %*% t(fpcay$evalbase)
  }

  return(predicted.values)
}
