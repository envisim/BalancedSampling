% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cube.R
\name{cube}
\alias{cube}
\alias{cubestratified}
\title{The Cube method}
\usage{
cube(prob, x, eps = 1e-12)

cubestratified(prob, x, integerStrata, eps = 1e-12)
}
\arguments{
\item{prob}{A vector of length N with inclusion probabilities.}

\item{x}{An N by q matrix of balancing auxiliary variables.}

\item{eps}{A small value used to determine when an updated probability is
close enough to 0.0 or 1.0.}

\item{integerStrata}{An integer vector of length N with stratum numbers.}
}
\value{
A vector of selected indices in 1,2,...,N.
}
\description{
Selects balanced samples with prescribed inclusion probabilities
from a finite population using the fast flight Cube Method.
}
\details{
If \code{prob} sum to an integer n, and \code{prob} is included as the first
balancing variable, a fixed sized sample (n) will be produced.
\subsection{Stratified cube}{

For \code{cubestratified}, \code{prob} is automatically inserted as a balancing variable.

The stratified version uses the fast flight Cube method and pooling of
landing phases.
}
}
\section{Functions}{
\itemize{
\item \code{cubestratified()}: 

}}
\examples{
\dontrun{
set.seed(12345);
N = 1000;
n = 100;
prob = rep(n/N, N);
x = matrix(runif(N * 2), ncol = 2);
s = cube(prob, x);
plot(x[, 1], x[, 2]);
points(x[s, 1], x[s, 2], pch = 19);

set.seed(12345);
N = 1000;
n = 100;
prob = rep(n/N, N);
x = matrix(runif(N * 2), ncol = 2);
strata = c(rep(1L, 100), rep(2L, 200), rep(3L, 300), rep(4L, 400));
s = cubestratified(prob, x, strata);
plot(x[, 1], x[, 2]);
points(x[s, 1], x[s, 2], pch = 19);

set.seed(12345);
prob = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
N = length(prob);
x = matrix(runif(N * 2), ncol = 2);
ep = rep(0L, N);
r = 10000L;
for (i in seq_len(r)) {
  s = cube(prob, cbind(prob, x));
  ep[s] = ep[s] + 1L;
}
print(ep / r);
}

}
\references{
Deville, J. C. and Tillé, Y. (2004).
Efficient balanced sampling: the cube method.
Biometrika, 91(4), 893-912.

Chauvet, G. and Tillé, Y. (2006).
A fast algorithm for balanced sampling.
Computational Statistics, 21(1), 53-62.

Chauvet, G. (2009).
Stratified balanced sampling.
Survey Methodology, 35, 115-119.
}
\seealso{
Other sampling: 
\code{\link{hlpm2}()},
\code{\link{lcube}()},
\code{\link{lpm}()},
\code{\link{scps}()}
}
\concept{sampling}
