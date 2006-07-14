\name{sbpsi}
\alias{sbpsi}
\alias{sbpsi.poly}
\alias{sbpsi.sing}
\alias{sbmodelnames}
\title{Model Specification Functions}
\description{
  \code{sbpsi.poly} and \code{sbpsi.sing} are \eqn{\psi} functions to
  specify a polynomial model and a singular model, respectively.
}
\usage{

sbpsi.poly(beta,s=1,k=1,sp=-1,aux=NULL,check=FALSE)

sbpsi.sing(beta,s=1,k=1,sp=-1,aux=NULL,check=FALSE)

sbmodelnames(m=3,poly=m,sing=m)
}
\arguments{
  \item{beta}{numeric vector of parameters;
    \eqn{\beta_0}=\code{beta[1]}, \eqn{\beta_1}=\code{beta[2]},...
    \eqn{\beta_{m-1}}=\code{beta[m]}, where \eqn{m} is the number of
    parameters.
  }
  \item{s}{\eqn{\sigma_0^2}.}
  \item{k}{numeric to specify the order of derivatives.}
  \item{sp}{\eqn{\sigma_p^2}.}
  \item{aux}{auxiliary parameter. Currently, not used.}
  \item{check}{logical for boundary check.}
  \item{m}{numeric to specify the maximum number of parameters.}
  \item{poly}{maximum number of parameters in poly models.}
  \item{sing}{maximum number of parameters in sing models.}
}
\details{
  For \eqn{k=1}, \code{sbpsi} functions return their \eqn{psi} function
  values at \eqn{\sigma^2=\sigma_0^2}.  Currently, two types of
  \code{sbpsi} functions are
  implemented. \code{sbpsi.poly} defines the polynomial model;
  \deqn{\psi(\sigma^2 | \beta) =
    \sum_{j=0}^{m-1} \beta_j \sigma^{2j}.}
  \code{sbpsi.sing} defines the singular model;
  \deqn{\psi(\sigma^2 | \beta) = \beta_0 +
    \sum_{j=1}^{m-2} \beta_j \sigma^{2j} (1+a)/(1+a \sigma),}
  where \eqn{a=\beta_{m-1}/(1-\beta_{m-1})} so that \eqn{0\le a \le
    \infty} for \eqn{0\le\beta_{m-1}\le1}.

  For \eqn{k>1}, \code{sbpsi} functions return values extrapolated at
  \eqn{\sigma^2=\sigma_p^2} using derivatives up to order \eqn{k-1}
  evaluated at \eqn{\sigma^2=\sigma_0^2};
  \deqn{\sum_{j=0}^{k-1} \frac{(\sigma_p^2-\sigma_0^2)^j}{j!}
    \frac{d^j \psi(x|\beta)}{d x^j}\Bigr|_{\sigma_0^2},}
  which reduces to \eqn{\psi(\sigma_0^2|\beta)} for \eqn{k=1}.
}
\value{
  \code{sbpsi.poly} and \code{sbpsi.sing} are examples of a sbpsi
  function; users can develop their own sbpsi functions for better
  model fitting by preparing \code{sbpsi.xxx} and \code{sbini.xxx}
  functions for model \code{xxx}.
  If check=FALSE, a sbpsi function returns
  the \eqn{\psi} function value or the extrapolation value.
  If check=TRUE, a sbpsi function returns NULL when all
  the elements of beta are included in the their valid
  intervals. Otherwise, a sbpsi function returns a list with components
  \code{beta} for the parameter value being modified to be on a boundary
  of the interval and \code{mask} of a logical vector indicating which
  elements are not on the boundary.

  \code{sbmodelnames} returns a character vector of model names.
}
\author{Hidetoshi Shimodaira}
\seealso{\code{\link{sbfit}}.}
\keyword{ models }