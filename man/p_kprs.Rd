% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects_and_functions.r
\name{p_kprs}
\alias{p_kprs}
\title{p_kprs computes the first derivative as a function of s using an analytical formula}
\usage{
p_kprs(s, n, m, omega)
}
\arguments{
\item{s}{The value needed for the first derivative to equal log(Wilks)}

\item{n}{The error Df of the one-way MANOVA analysis considered}

\item{m}{The hypothesis Df of the one-way MANOVA analysis considered}

\item{omega}{A vector of eigenvalues of the Wilks Non-Centrality Parameter corresponding
to one independent variable.}
}
\description{
p_kprs computes the first derivative as a function of s using an analytical formula
}
