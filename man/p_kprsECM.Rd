% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ECM_functions.r
\name{p_kprsECM}
\alias{p_kprsECM}
\title{p_kprsECM computes the first derivative as a function of s using an analytical formula}
\usage{
p_kprsECM(s, N1, N2, sigma_eigv)
}
\arguments{
\item{s}{The value needed for the first derivative to equal log(lambdaECM)}

\item{N1}{Number of rows in matrix 1}

\item{N2}{Number of rows in matrix 2}

\item{sigma_eigv}{A vector of eigenvalues of (Sigma11)*(inv-Sigma22)}
}
\description{
p_kprsECM computes the first derivative as a function of s using an analytical formula
}
