\name{sbpval}
\alias{sbpval}
\alias{sbpval.summary.scaleboot}
\alias{sbpval.summary.scalebootv}
\title{Extract P-values}
\description{
  \code{sbpval} extracts p-values from \code{"summary.scaleboot"} or
    \code{"summary.scalebootv"} objects.
}
\usage{
sbpval(x, ...)

\method{sbpval}{summary.scaleboot}(x,sd=FALSE,best=FALSE,...)

\method{sbpval}{summary.scalebootv}(x,...)

}
\arguments{
  \item{x}{an object used to select a method.}
  \item{sd}{logical. Should standard errors be returned as well?}
  \item{best}{logical. Should only the best model be returned?}
  \item{...}{further arguments passed to or from other methods.}  
}
\details{
  This method is used only to extract previously calculated p-values from the
  summary object.
}
\value{
  The \code{sbpval} method for the class \code{"summary.scaleboot"} returns a
  matrix of p-values. Each row is a vector of
  \eqn{p_k} for \eqn{k} as specified in the \code{summary} method.
  If \code{sd=TRUE}, then it returns a list with two components:
  \code{estimate} for the p-values and \code{sd} for its standard errors.
  If \code{best=TRUE}, then only the result for the best fit model is returned.

  The \code{sbpval} method for the class \code{"summary.scalebootv"} returns a
  list vector with each component obtained by applying \code{sbpval} to
  each \code{"scaleboot"} object.
}
\author{Hidetoshi Shimodaira}
\seealso{\code{\link{summary.scaleboot}}.}
\examples{
data(mam15)
a <- mam15.relltest[["t4"]] # an object of class "scaleboot"
b <- summary(a) # calculate p-values
b # print the p-values
sbpval(b) # extract a matrix p-values
sbpval(b,sd=TRUE) # with sd
sbpval(b,best=TRUE) # extract a vector of p-values
sbpval(b,best=TRUE,sd=TRUE) # with sd
}
\keyword{ models }
\keyword{ nonparametric }
