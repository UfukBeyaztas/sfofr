\name{predict_sffreg}
\alias{predict_sffreg}
\title{Prediction for Spatial Function-on-Function Regression Models}
\description{Generates predictions for new functional predictors using a fitted Spatial Function-on-Function Regression (SFoFR) model. The function supports both independent and spatial regression types and handles spatial weight matrices for test samples when spatial regression is used.}
\usage{
predict_sffreg(object, xnew, Wnew = NULL)
}
\arguments{
  \item{object}{An \code{sffreg} object containing a fitted Spatial Function-on-Function Regression model.}
  \item{xnew}{A matrix of new functional predictors for the test sample, where rows correspond to observations and columns represent function evaluations at discrete grid points.}
  \item{Wnew}{An optional spatial weight matrix for the test sample. Required if \code{model.type = "spatial"}.}
}
\details{The function uses the fitted \code{sffreg} object to predict the functional response for new functional predictors based on the selected model type (\code{"independent"} or \code{"spatial"}).

\itemize{
  \item \code{"independent"}: Predicts the functional response by projecting the new functional predictor onto the eigenfunctions of the functional principal components (FPC) used in the original model.
  \item \code{"spatial"}: Incorporates spatial autocorrelation by utilizing the spatial weight matrix (\code{Wnew}) in conjunction with the spatial functional principal components (SFPC). The predictions are computed by solving a linear system involving the Kronecker product of \code{Wnew} and the spatial coefficients.
}

The function ensures that the predictions align with the dimensionality and structure of the original functional response.
}
\value{A matrix of predicted functional responses, where rows correspond to observations in the test sample and columns represent predicted function values at the evaluation grid points.}
\seealso{\code{\link{sffreg}}, \code{\link{data_generation}}}
\examples{
## Simulate data and fit SFoFR model
n <- 500
data <- data_generation(n, W.type = "ivw", rf = 0.8)
fit <- sffreg(data$Y, data$X, data$W, model.type = "spatial")

## Generate new functional predictors
data_test <- data_generation(n, W.type = "ivw", rf = 0.8)

## Predict using spatial model
Wnew <- data_test$W  # Use the same weight matrix for simplicity
predictions_spatial <- predict_sffreg(fit, data_test$X, Wnew = Wnew)
}
