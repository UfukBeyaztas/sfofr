\name{sffreg}
\alias{sffreg}
\title{Spatial Function-on-Function Regression Model}
\description{Fits a Spatial Function-on-Function Regression (SFoFR) model, where a functional response is regressed on a functional predictor with spatial autocorrelation captured through a spatial weight matrix. The model leverages functional principal component (FPC) and spatial functional principal component (SFPC) decompositions for dimension reduction.}
\usage{
sffreg(Y, X, W, nbasis = NULL, gpy = NULL, gpx = NULL,
model.type = c("independent", "spatial"))
}
\arguments{
  \item{Y}{A matrix of functional responses, where rows correspond to observations and columns represent function evaluations at discrete grid points.}
  \item{X}{A matrix of functional predictors, structured similarly to \code{Y}.}
  \item{W}{A spatial weight matrix encoding the spatial relationships among observations.}
  \item{nbasis}{Number of basis functions for B-spline basis expansion. If \code{NULL}, a default value is used.}
  \item{gpy}{Grid points for the functional response. Defaults to 101 equally spaced points over \eqn{[0,1]}.}
  \item{gpx}{Grid points for the functional predictor. Defaults to 101 equally spaced points over \eqn{[0,1]}.}
  \item{model.type}{Specifies the regression model type. Options are:
    \itemize{
      \item \code{"independent"}: Assumes no spatial autocorrelation.
      \item \code{"spatial"}: Includes spatial autocorrelation using SFPC-based techniques.
    }
  }
}
\details{The SFoFR model is formulated as:
\deqn{
Y_i(t) = \sum_{i^{\prime}=1}^n w_{i i^{\prime}} \int_0^1 Y_{i^{\prime}}(u) \rho(u,t) du + \int_0^1 X_i(s) \beta(s,t) ds + \epsilon_i(t),
}
where:
\itemize{
  \item \eqn{w_{ii'}} are elements of the spatial weight matrix \eqn{\bm{W}}.
  \item \eqn{\rho(u,t)} is the spatial autocorrelation function.
  \item \eqn{\beta(s,t)} is the bivariate regression coefficient function.
  \item \eqn{\epsilon_i(t)} is the random functional error term.
}

Using FPC and SFPC decompositions, the infinite-dimensional SFoFR model is reduced to a multivariate spatial autoregressive (MSAR) model:
\deqn{
\widetilde{Y} = (\bm{\rho}^\top \otimes \bm{W}) \widetilde{Y} + \bm{X}^{\top*} \widetilde{\beta} + \widetilde{e},
}
where \eqn{\widetilde{Y}}, \eqn{\widetilde{\beta}}, and \eqn{\widetilde{e}} are vectorized representations of the response, regression coefficients, and residuals, respectively. The parameters \eqn{\bm{\rho}} and \eqn{\bm{\beta}} are estimated via a least-squares objective function, ensuring identifiability under specific regularity conditions.

The estimates of the spatial autocorrelation and regression coefficient functions are:
\deqn{
\widehat{\rho}(u,t) = \phi^\top(u) \widehat{\bm{\rho}} \phi(t), \quad \widehat{\beta}(s,t) = \psi^\top(s) \widehat{\bm{\beta}} \phi(t),
}
where \eqn{\phi(t)} and \eqn{\psi(s)} are SFPC and FPC eigenfunctions, respectively.
}
\value{A list with the following components:
  \item{fitted.values}{Predicted functional response values.}
  \item{bhat}{Estimated regression coefficient matrix.}
  \item{rhohat}{Estimated spatial autocorrelation matrix.}
  \item{residuals}{Residuals of the fitted model.}
  \item{model_info}{A list containing model-specific details, including the type of model, basis functions, and principal components.}
}
\references{
Khoo, T. H., Pathmanathan, D. and Dabo-Niang, S. (2023), ‘Spatial autocorrelation of global stock exchanges using functional areal spatial principal component analysis’, Mathematics 11(2), 674.
Zhu, X., Cai, Z. and Ma, Y. (2022), ‘Network functional varying coefficient model’, Journal of the American Statistical Association: Theory and Methods 117(540), 2074–2085.
}
\author{
Ufuk Beyaztas, Han Lin Shang, Gizel Bakicierler, Abhijit Mandal, Roger S. Zoh, and Carmen D. Tekwe
}
\seealso{
\code{\link{data_generation}}, \code{\link{getSPCA}}, \code{\link{predict_sffreg}}
}
\examples{
## Generate data and fit SFoFR model
n <- 500
data <- data_generation(n, W.type = "expw", rf = 0.8)
fit <- sffreg(data$Y, data$X, data$W, model.type = "spatial")
}
