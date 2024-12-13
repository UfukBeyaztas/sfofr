\name{getSPCA}
\alias{getSPCA}
\title{Spatial Principal Component Analysis for Functional Data}
\description{Performs Spatial Principal Component Analysis (SFPC) on functional data with a spatial structure. The method computes spatial functional principal components (SFPCs) and their corresponding scores, capturing both variability and spatial dependence in the data.}
\usage{
getSPCA(data, nbasis, gp, W)
}
\arguments{
  \item{data}{A matrix of functional data, where rows correspond to observations and columns represent functional values evaluated at discrete grid points.}
  \item{nbasis}{Number of basis functions to be used for the B-spline basis expansion.}
  \item{gp}{Grid points for the functional data, representing the domain over which the functions are observed.}
  \item{W}{Spatial weight matrix, encoding the spatial relationships among the observations.}
}
\details{This function implements Spatial Principal Component Analysis (SFPC), which extends the classical FPCA to account for spatial dependencies among functional data. The methodology is grounded in the eigen-decomposition of the covariance function, augmented with a spatial weight matrix to capture spatial autocorrelation.

The covariance function of a centered functional variable \eqn{x(t)} is decomposed using the Karhunen-Loève theorem:
\deqn{
G(s,t) = \sum_{k=1}^{\infty} \delta_k \eta_k(s) \eta_k(t),
}
where \eqn{\eta_k(t)} are the eigenfunctions and \eqn{\delta_k} are the eigenvalues. In SFPC, these eigenfunctions are computed with additional constraints that maximize spatial autocorrelation as measured by Moran's I statistic:
\deqn{
\mathrm{Var}(\widetilde{\bm{\kappa}}) I(\widetilde{\bm{\kappa}}) = n^{-1} \widetilde{\bm{\chi}}^\top \bm{D}^\top \bm{W} \bm{D} \widetilde{\bm{\chi}},
}
where \eqn{\bm{W}} is the spatial weight matrix and \eqn{\bm{D}} is the matrix of basis expansion coefficients.

The eigenfunctions obtained from SFPC represent spatially smoothed principal components, and their corresponding scores summarize the joint variation and spatial structure of the functional data.
}
\value{
A list with the following components:
  \item{PCAcoef}{ A \code{fd} object representing the estimated spatial principal components (eigenfunctions).}
  \item{PCAscore}{  A matrix of SFPC scores for each observation.}
  \item{meanScore}{  A \code{fd} object representing the mean function of the functional data.}
  \item{ncomp}{  Number of components retained, based on the proportion of variance explained.}
  \item{bs_basis}{  B-spline basis object used for the functional representation.}
  \item{evalbase}{  Evaluated basis matrix for grid points \code{gp}.}
  \item{gp}{  The grid points over which the functions are evaluated.}
  \item{W}{  The spatial weight matrix used for SFPC analysis.}
}
\references{
Khoo, T. H., Pathmanathan, D. and Dabo-Niang, S. (2023), ‘Spatial autocorrelation of global stock exchanges using functional areal spatial principal component analysis’, Mathematics 11(2), 674.
}
\author{
Ufuk Beyaztas, Han Lin Shang, Gizel Bakicierler, Abhijit Mandal, Roger S. Zoh, and Carmen D. Tekwe
}
\seealso{
\code{\link{data_generation}}, \code{\link{sffreg}}, \code{\link{predict_sffreg}}
}
\examples{
n <- 100
data <- data_generation(n, W.type = "ivw", rf = 0.8)
sfpc <- getSPCA(data$Y, nbasis = 10, gp = seq(0, 1, length.out = 101), W = data$W)
}
