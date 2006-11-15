\name{summary.scaleboot}
\alias{summary.scaleboot}
\alias{summary.scalebootv}
\alias{print.summary.scaleboot}
\alias{print.summary.scalebootv}
\title{P-value Calculation for Multiscale Bootstrap}
\description{
  \code{summary} method for class \code{"scaleboot"} and \code{"scalebootv"}.
}
\usage{

\method{summary}{scaleboot}(object,models=names(object$fi),k=1:3,s=1,sp=-1,...)

\method{summary}{scalebootv}(object,models=attr(object,"models"),k=1:3,...)

}
\arguments{
  \item{object}{an object used to select a method. }
  \item{models}{character vector of model names. If numeric,
    \code{names(object$fi)[models]} is used for each \code{"scaleboot"} object.}
  \item{k}{numeric vector of \eqn{k} for calculating p-values.}
  \item{s}{\eqn{\sigma_0^2}}
  \item{sp}{\eqn{\sigma_p^2}}
  \item{...}{further arguments passed to and from other methods.}  
}

\details{ For each model, a class of approximately unbiased p-values,
indexed by \eqn{k=1,2,...}, is calculaed.  The p-values are named 
\code{k.1}, \code{k.2}, ..., where \eqn{k=1} (\code{k.1}) corresponds to
the ordinary bootstrap probability, and \eqn{k=2} (\code{k.2})
corresponds to the third-order accurate p-value of Shimodaira (2002). As the
\eqn{k} value increases, the bias of testing decreases, although the
p-value becomes less stable numerically and the monotonicity of rejection
regions becomes worse. Typically, \eqn{k=3} provides a reasonable
compromise. The \code{\link{sbpval}} method is available to extract p-values from
the \code{"summary.scaleboot"} object.

The p-value is defined as
  \deqn{ p_k = 1 - \Phi\left( \sum_{j=0}^{k-1}
\frac{(\sigma_p^2-\sigma_0^2)^j}{j!}  \frac{d^j \psi(x|\beta)}{d
x^j}\Bigr|_{\sigma_0^2} \right),} where \eqn{\psi(\sigma^2|\beta)} is the
model specification function, \eqn{\sigma_0^2} is the evaluation point
for the Taylor series, and \eqn{\sigma_p^2} is an additional
parameter. Typically, we do not change the default values
\eqn{\sigma_0^2=1} and \eqn{\sigma_p^2=-1}.

The p-values are justified only for good fitting models. By default,
the model which minimizes the AIC value is selected. We can modify the
AIC value by using the \code{sbaic} function. We also diagnose the
fitting by using the \code{plot} method.
}
\value{
  \code{summary.scaleboot} returns
  an object of the class \code{"summary.scaleboot"}, which is inherited
  from the class \code{"scaleboot"}. It is a list containing all the components of class
  \code{"scaleboot"} and the following components:
  \item{pv}{matrix of p-values of size \code{length(models)} *
    \code{length(k)} with elements \eqn{p_k}.}
  \item{pe}{matrix of standard errors of p-values.}
  \item{best}{ a list consisting of components \code{model} for the best
  fitting model name, \code{aic} for its AIC value, \code{pv} for a
  vector of p-values, and \code{pe} for a vector of standard errors.}
  \item{parex}{a list of components \code{k}, \code{s}, and \code{sp}.}
}
\author{Hidetoshi Shimodaira}
\seealso{\code{\link{sbfit}}, \code{\link{sbpsi}}, \code{\link{sbpval}},
   \code{\link{sbaic}}.}
\examples{
data(mam15)
## For a single hypothesis
a <- mam15.relltest[["t4"]] # an object of class "scaleboot"
summary(a) # calculate and print p-values (k=1:3)
summary(a,k=1:4) # up to "k.4" p-value.

## For multiple hypotheses
b <- mam15.relltest[1:15] # an object of class "scalebootv"
summary(b) # calculate and print p-values (k=1:3)
summary(b,k=1:4) # up to "k.4" p-value.

}
\keyword{ models }
\keyword{ nonparametric }
