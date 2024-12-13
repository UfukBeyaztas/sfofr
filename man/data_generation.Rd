\name{data_generation}
\alias{data_generation}
\title{Generate Spatially Correlated Functional Data}
\description{Simulates spatially correlated functional data for use in spatial function-on-function regression models. The function generates functional predictors, responses, spatial weight matrices, and the underlying regression and spatial autocorrelation structures, based on specified formulas and parameters.}
\usage{
data_generation(n, nphi = 10, gpy = NULL, gpx = NULL, W.type = c("ivw", "expw"),
                rf = 0.9, sd.error = 0.01, tol = 0.001, max_iter = 1000)
}
\arguments{
  \item{n}{Number of samples to generate.}
  \item{nphi}{Number of basis functions used to construct the functional predictor.}
  \item{gpy}{Grid points for the functional response. Defaults to 101 equidistant points over the interval [0,1] if not provided.}
  \item{gpx}{Grid points for the functional predictor. Defaults to 101 equidistant points over the interval [0,1] if not provided.}
  \item{W.type}{Type of spatial weight matrix. Options are "ivw" (inverse distance weighting) and "expw" (exponential decay weighting).}
  \item{rf}{Spatial dependence parameter \eqn{\alpha}, which regulates the strength of spatial effects.}
  \item{sd.error}{Standard deviation of the error term added to the functional response.}
  \item{tol}{Tolerance for the Neumann series approximation used to compute the spatial effects.}
  \item{max_iter}{Maximum number of iterations for the Neumann series approximation.}
}
\details{The function generates data based on the spatial function-on-function regression model: \deqn{
Y_i(t) = \sum_{i^{\prime}=1}^n w_{i i^{\prime}} \int_0^1 Y_{i^{\prime}}(u) \rho(u,t) du + \int_0^1 X_i(s) \beta(s,t) ds + \epsilon_i(t),
}where:
\itemize{
  \item \eqn{w_{ii'}} are elements of a spatial weight matrix \eqn{\bm{W}} that encodes spatial relationships between locations.
  \item \eqn{\rho(u,t)} is the spatial autocorrelation function defined as \eqn{\rho(u,t) = \alpha \frac{1 + u t}{1 + |u - t|}}, where \eqn{\alpha} controls the spatial dependence strength.
  \item \eqn{\beta(s,t)} is the bivariate regression coefficient function, specified as \eqn{\beta(s,t) = 2 + s + t + 0.5 \sin(2 \pi s t)}.
  \item \eqn{\epsilon_i(t)} is a random error term.
}

The functional predictor \eqn{X_i(s)} is constructed as: \deqn{
X_i(s) = \sum_{k=1}^{10} \frac{1}{k^{3/2}} \left( \nu_{i1,k} \sqrt{2} \cos(k \pi s) + \nu_{i2,k} \sqrt{2} \sin(k \pi s) \right),
}
where \eqn{\nu_{i1,k}} and \eqn{\nu_{i2,k}} are independent standard normal random variables.

Two types of spatial weight matrices (\eqn{\bm{W}}) are supported:
\itemize{
  \item \emph{Inverse distance weighting (ivw)}: \eqn{w_{ii'} = \frac{1}{1 + |i - i'|}}, normalized to ensure each row sums to 1.
  \item \emph{Exponential decay weighting (expw)}: \eqn{w_{ii'} = \exp(-d |i - i'|)}, with a decay parameter \eqn{d = 0.5}, also normalized row-wise.
}
The functional response is computed iteratively using the Neumann series approximation until convergence, ensuring accurate spatial effects representation.
}
\value{
Returns a list containing:
  \item{Y}{Simulated functional response matrix \eqn{Y_i(t)}.}
  \item{Y_true}{Functional response matrix without the error term.}
  \item{X}{Simulated functional predictor matrix \eqn{X_i(s)}.}
  \item{W}{Spatial weight matrix.}
  \item{rho}{Spatial autocorrelation matrix \eqn{\rho(u,t)}.}
  \item{beta}{True regression coefficient matrix \eqn{\beta(s,t)}.}
}
\references{
Hoshino, T. (2024), Functional spatial autoregressive models, Technical report, arXiv. URL: \url{https://arxiv.org/pdf/2402.14763}.
}
\author{
Ufuk Beyaztas, Han Lin Shang, Gizel Bakicierler, Abhijit Mandal, Roger S. Zoh, and Carmen D. Tekwe
}
\note{
Ensure that the inputs, particularly \code{W.type}, are specified correctly to avoid errors. The spatial dependence strength (\code{rf}) must satisfy the contraction condition \eqn{\Vert \rho \Vert_{\infty} < 1 / \Vert \bm{W} \Vert_{\infty}} to guarantee convergence.
}

\seealso{
\code{\link{getSPCA}}, \code{\link{sffreg}}, \code{\link{predict_sffreg}}
}
\examples{
## Example usage of the data_generation function
n <- 100
data <- data_generation(n, W.type = "ivw", rf = 0.8)
str(data)
}
