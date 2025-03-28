sffreg <- function(Y, X, W, nbasis=NULL, gpy=NULL, gpx=NULL,
                   model.type = c("independent","spatial"))
  {

  if(!model.type %in% c("independent","spatial"))
    stop("Model type must be one of the followings: independent and spatial!")
  if(dim(Y)[1] != dim(X)[1])
    stop("Number of observations for functional response and functional covariate
         must be equal!")
  if(dim(W)[1] != dim(W)[2])
    stop("Weight matrix W must be square matrix!")
  if(dim(W)[1] != dim(Y)[1])
    stop("Dimensions of W and Y does not match!")

  if(is.null(gpy))
    gpy <- seq(0, 1, length.out = dim(Y)[2])
  if(is.null(gpx))
    gpx <- seq(0, 1, length.out = dim(X)[2])

  if(is.null(nbasis))
    nbasis <- 4

  W <- normalize_weights(W)

  if(model.type == "independent") {
    fpcay <- getPCA(data = Y, nbasis = nbasis, gp = gpy)
    fpcax <- getPCA(data = X, nbasis = nbasis, gp = gpx)
  }else if(model.type == "spatial") {
    fpcay <- getSPCA(data = Y, nbasis = nbasis, gp = gpy, W = W)
    fpcax <- getSPCA(data = X, nbasis = nbasis, gp = gpx, W = W)
  }

  n <- dim(Y)[1]

  scoY <- fpcay$PCAscore
  scoX <- fpcax$PCAscore

  if(model.type == "independent") {
    dmodel <- lm(scoY~scoX)
    cfb <- as.matrix(dmodel$coefficients)[-1,]
    cfr <- NULL
    bhat <- fpcax$evalbase %*% (fpcax$PCAcoef$coefs %*% cfb %*%
                                  t(fpcay$PCAcoef$coefs)) %*% t(fpcay$evalbase)
    rhohat <- NULL
    yhat <- ((scoX %*% cfb %*% t(fpcay$PCAcoef$coefs)) +
               matrix(fpcay$meanScore$coefs, nrow = n, ncol = nbasis, byrow = T)) %*%
      t(fpcay$evalbase)
    residuals <- Y - yhat
  }else if(model.type == "spatial") {
    rho0 <- Matrix(0, dim(scoY)[2], dim(scoY)[2])
    B0 <- Matrix(0, dim(scoX)[2], dim(scoY)[2])
    Sige0 <- Diagonal(n =dim(scoY)[2], x = rep(1e-2, dim(scoY)[2]))

    msar_model <- msar_fun(scoY, scoX, W, rho = rho0, B = B0, Sige = Sige0)
    cfb <- msar_model$B
    cfr <- msar_model$rho

    bhat <- fpcax$evalbase %*% (fpcax$PCAcoef$coefs %*% cfb %*%
                                  t(fpcay$PCAcoef$coefs)) %*% t(fpcay$evalbase)
    bhat <- as.matrix(bhat)
    rhohat <- fpcay$evalbase %*% (fpcay$PCAcoef$coefs %*% cfr %*%
                                    t(fpcay$PCAcoef$coefs)) %*% t(fpcay$evalbase)
    rhohat <- as.matrix(rhohat)

    yhatvec <- solve(diag(n*dim(scoY)[2]) - t(cfr) %x% W) %*%
      (diag(dim(scoY)[2]) %x% scoX) %*% as.vector(cfb)
    yhat1 <- matrix(yhatvec, ncol = dim(scoY)[2]) %*% t(fpcay$PCAcoef$coefs) +
      matrix(fpcay$meanScore$coefs, nrow = n, ncol = nbasis, byrow = T)
    yhat <- yhat1 %*% t(fpcay$evalbase)
    residuals <- Y - yhat
  }

  model_info <- list()
  model_info$model.type <- model.type
  model_info$nbasis <- nbasis
  model_info$fpcay <- fpcay
  model_info$fpcax <- fpcax
  model_info$cfb <- cfb
  model_info$cfr <- cfr

  return(list(fitted.values = yhat, bhat = bhat, rhohat = rhohat,
              residuals = residuals, model_info = model_info))

}
