% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects_and_functions.r
\name{p_lmglm}
\alias{p_lmglm}
\title{p_lmglm, computation of log(mgf) of the log(Wilks (Lambda) Statistic)}
\usage{
p_lmglm(s_hat, n, m, omega)
}
\arguments{
\item{s_hat, }{The value needed for the first derivative to equal log(Wilks), found by function p_find_s}

\item{n}{The error Df of the one-way MANOVA analysis considered}

\item{m}{The hypothesis Df of the one-way MANOVA analysis considered}

\item{omega}{The a vector of eigenvalues of the Wilks Non-Centrality Parameter 
corresponding to one independent variable.}
}
\description{
p_lmglm, computation of log(mgf) of the log(Wilks (Lambda) Statistic)
}
