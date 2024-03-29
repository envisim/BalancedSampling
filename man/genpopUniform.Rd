% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/populations.R
\name{genpopUniform}
\alias{genpopUniform}
\alias{genpopPoisson}
\title{Generate populations}
\usage{
genpopUniform(size, dims = 2L, from = 0, to = 1)

genpopPoisson(
  parents,
  children,
  dims = 2L,
  from = 0,
  to = 1,
  distribution = function(n) rnorm(n, 0, 0.02),
  mirror = TRUE
)
}
\arguments{
\item{size}{The size of the population}

\item{dims}{The number of auxiliary variables}

\item{from}{A number or a vector of size \code{dims} with the minimum values}

\item{to}{A number or a vector of size \code{dims} with the maximum values}

\item{parents}{The number of parent locations}

\item{children}{A number or a vector of size \code{parents} with the mean number of
children to be spawned.}

\item{distribution}{A function taking a number as a variable, returning the
offset from the parent location.}

\item{mirror}{If \code{TRUE}, the population is mirrored to be inside \code{from} and \code{to}.}
}
\description{
Generate uniform and poisson cluster process populations

If \code{from} and \code{to} is used with \code{genpopPoisson} together with \code{mirror}, the
population will be bounded within these values.
For the \code{genpopUniform}, these numbers represent the minimum and maximum
values of the uniform distribution.
}
\section{Functions}{
\itemize{
\item \code{genpopPoisson()}: Poisson cluster process

}}
\examples{
\dontrun{
set.seed(12345);
x = genpopUniform(120, 2L);
N = nrow(x);
n = 60;
prob = rep(n / N, N);
s = lpm2(prob, x);
b = sb(prob, x, s);
}

\dontrun{
set.seed(12345);
x = genpopPoisson(70, 50, 2L);
N = nrow(x);
n = 60;
prob = rep(n / N, N);
s = lpm2(prob, x);
b = sb(prob, x, s);
}

}
