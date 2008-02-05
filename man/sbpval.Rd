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

\method{sbpval}{summary.scaleboot}(x,sd=FALSE,select=c("average","best","all"),...)

\method{sbpval}{summary.scalebootv}(x,...)

}
\arguments{
  \item{x}{an object used to select a method.}
  \item{sd}{logical. Should standard errors be returned as well?}
  \item{select}{character. If "average" or "best", only the p-values of
    corresponding models are returned. If "all", then p-values of all
    the models are returned.}
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
sbpval(b) # extract a vector of p-values which are averaged by Akaike weights.
sbpval(b,sd=TRUE) # with sd
sbpval(b,select="all") # extract a matrix of p-values
sbpval(b,select="all",sd=TRUE) # with sd
}
\keyword{ models }
\keyword{ nonparametric }
