% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tensor_clustering.R
\name{sim_dTBM}
\alias{sim_dTBM}
\title{Simulation of degree-corrected tensor block models}
\usage{
sim_dTBM(
  seed = NA,
  imat = FALSE,
  asymm = FALSE,
  p,
  r,
  core_control = c("random", "control"),
  delta = NULL,
  s_min = NULL,
  s_max = NULL,
  dist = c("normal", "binary"),
  sigma = 1,
  theta_dist = c("abs_normal", "pareto", "non"),
  alpha = NULL,
  beta = NULL
)
}
\arguments{
\item{seed}{number, random seed for generating data}

\item{imat}{logic variable, if "TRUE", generate matrix data; if "FALSE", generate order-3 tensor data}

\item{asymm}{logic variable, if "TRUE", clustering assignment differs in different modes; if "FALSE", all the modes share the same clustering assignment}

\item{p}{vector, dimension of the tensor/matrix observation}

\item{r}{vector, cluster number on each mode}

\item{core_control}{character, the way to control the generation of core tensor/matrix; see "details"}

\item{delta}{number, Frobenius norm of the slices in core tensor if \code{core_control = "control"}}

\item{s_min}{number, value of off-diagonal elements in original core tensor/matrix if \code{core_control = "control"}}

\item{s_max}{number, value of diagonal elements in original core tensor/matrix if \code{core_control = "control"}}

\item{dist}{character, distribution of tensor/matrix observation; see "details"}

\item{sigma}{number, standard deviation of Gaussian noise if \code{dist = "normal"}}

\item{theta_dist}{character, distribution of degree heterogeneity; see "details"}

\item{alpha}{number, shape parameter in pareto distribution if \code{theta_dist = "pareto"}}

\item{beta}{number, scale parameter in pareto distribution if \code{theta_dist = "pareto"}}
}
\value{
a list containing the following:

\code{Y} {array ( if \code{imat = FALSE} )/matrix ( if \code{imat = TRUE} ), simulated tensor/matrix observations with dimension \code{p}  }

\code{X} {array ( if \code{imat = FALSE} )/matrix ( if \code{imat = TRUE} ), mean tensor/matrix of the observation, i.e., the expectation of \code{Y}}

\code{S} {array ( if \code{imat = FALSE} )/matrix ( if \code{imat = TRUE} ), core tensor/matrix recording the block effects with dimension \code{r}}

\code{theta} {a list of vectors, degree heterogeneity on each mode}

\code{z} {a list of vectors, clustering assignment on each mode}
}
\description{
Generate order-3 tensor/matrix observations with degree heterogeneity under degree-corrected tensor block models.
}
\details{
The general tensor observation is generated as

\code{Y = S x1 Theta1 M1 x2 Theta2 M2 x3 Theta3 M3 + E,}

where \code{S} is the core tensor, \code{Thetak} is a diagonal matrix with elements in the \code{k}-th vector of \code{theta},
\code{Mk} is the membership matrix based on the clustering assignment in the \code{k}-th vector of \code{z} with \code{r[k]} clusters,
\code{E} is the mean-zero noise tensor, and \code{xk} refers to the matrix-by-tensor product on the \code{k}-th mode, for \code{k = 1,2,3}.

If \code{imat = TRUE}, \code{Y,S,E} degenerate to matrix and \code{Y = Theta1 M1 S M2^T Theta2^T + E}.

If \code{asymm = FALSE}, \code{Thetak = Theta} and \code{Mk = M} for all \code{k = 1,2,3}.

\code{core_control} specifies the way to generate \code{S}:

If \code{core_control = "control"}, first generate \code{S} as a diagonal tensor/matrix with diagonal elements \code{s_max} and off-diagonal elements \code{s_min};
then scale the original core such that Frobenius norm of the slices equal to \code{delta}, i.e, \code{delta = sqrt(sum(S[1,,]^2))} or \code{delta = sqrt(sum(S[1,]^2))};
ignore the scaling if \code{delta = NULL}; option \code{"control"} is only applicable for symmetric case \code{asymm = FALSE}.

If \code{core_control = "random"}, generate \code{S} with random entries following uniform distribution U(0,1).

\code{dist} specifies the distribution of \code{E}: \code{"normal"} for Gaussian and \code{"binary"} for Bernoulli distribution; \code{sigma} specifices the standard deviation if \code{dist = "normal"}.

\code{theta_dist} firstly specifies the distribution of \code{theta}: \code{"non"} for constant 1, \code{"abs_normal"} for absoulte normal distribution, \code{"pareto"} for pareto distribution; \code{alpha, beta} specify the shape and scale parameter if \code{theta_dist = "pareto"};
then scale \code{theta} to have mean equal to one in each cluster.
}
\examples{

test_data = sim_dTBM(seed = 1, imat = FALSE, asymm = FALSE, p = c(50,50,50), r = c(3,3,3),
                    core_control = "control", s_min = 0.05, s_max = 1,
                    dist = "normal", sigma = 0.5,
                    theta_dist = "pareto", alpha = 4, beta = 3/4)
}
