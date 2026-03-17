\name{print.rbfmvar}
\alias{print.rbfmvar}
\title{Print Method for rbfmvar Objects}
\usage{
\method{print}{rbfmvar}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"rbfmvar"}.}

\item{...}{Further arguments passed to or from other methods (unused).}
}
\value{
Invisibly returns \code{x}.
}
\description{
Prints a formatted summary of RBFM-VAR estimation results.
}
\examples{
set.seed(1)
n  <- 60
Y  <- cbind(cumsum(rnorm(n)), cumsum(rnorm(n)))
res <- rbfmvar(Y, lags = 2)
print(res)
}
