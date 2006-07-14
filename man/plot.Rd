\name{plot.scaleboot}
\alias{plot.scaleboot}
\alias{lines.scaleboot}
\alias{sblegend}
\alias{plot.scalebootv}
\title{Plot Diagnostics for Multiscale Bootstrap}
\description{
  \code{plot} method for class \code{"scaleboot"}.
}
\usage{

\method{plot}{scaleboot}(x, models = names(x$fi), xval = c("inverse","sigma","square"),
     yval = c("zvalue", "pvalue", "psi"), xlab = NULL,
     ylab = NULL,log.xy = "", xlim = NULL, ylim = NULL,
     add = F, length.x = 300, col =1:6, lty = 1:5, lwd = par("lwd"),
     pch = 1, cex = 1, pt.col = col[1],pt.lwd = lwd[1], ...)

\method{plot}{scalebootv}(x,models=attr(x,"models"),...)

\method{lines}{scaleboot}(x,z,models=names(x$fi), length.x=z$length.x,
      col=z$col,lty=z$lty,lwd=z$lwd,... )

sblegend(x="topright",y=NULL,z,inset=0.1,...)

}
\arguments{
  \item{x}{an object used to select a method. 
    For \code{sblegend},
    \code{x} is a numeric or character such as "lefttop" or "righttop", which
    is passed to \code{legend}.
  }
  \item{models}{character vector of model names. If numeric,
    \code{names(object$fi)[models]} is used for each \code{"scaleboot"}
    object.}
  \item{xval}{specifies x-axis. "inverse" for \eqn{1/\sigma} or
    \code{1/sqrt(sa[i])}, "sigma" for \eqn{\sigma} or \code{sqrt(sa[i])}.}
  \item{yval}{specifies y-axis. "zvalue" for
    \eqn{\psi(\sigma^2|\beta)/\sigma} or \code{qnorm(1-bp[i])},
    "pvalue" for \eqn{1-\Phi(\psi(\sigma^2|\beta)/\sigma)} or \code{bp[i]},
    "psi" for \eqn{\psi(\sigma^2|\beta)} or \code{sqrt(sa[i])*qnorm(1-bp[i])}.}
  \item{xlab}{label for x-axis.}
  \item{ylab}{label for y-axis.}
  \item{log.xy}{character to specify log-scale. "", "x", "y", or "xy".}
  \item{xlim}{range for x-axis.}
  \item{ylim}{range for y-axis.}
  \item{add}{logical for adding another plot.}
  \item{length.x}{the number of segments to draw curves.}
  \item{col}{color for model curves.}
  \item{lty}{lty for model curves.}
  \item{lwd}{lwd for model curves.}
  \item{pch}{pch for bp points.}
  \item{cex}{cex for bp points.}
  \item{pt.col}{col for bp points.}
  \item{pt.lwd}{lwd for bp points.}
  \item{...}{further arguments passed to or from other methods.}  

  \item{z}{output from previous \code{plot.scaleboot}.}
  \item{y}{numeric passed to \code{legend}.}
  \item{inset}{inset distance from the margins, which is passed to
    \code{legend}.}
}
\details{
  \code{plot} method plots bootstrap probabilities and call \code{lines}
  method, which draws fitted curves for models.
}
\author{Hidetoshi Shimodaira}
\seealso{\code{\link{sbfit}},\code{\link{legend}}.}
\examples{
data(mam15)
## a single plot
a <- mam15.relltest[["t4"]] # an object of class "scaleboot"
sblegend("topleft",z=plot(a)) # x=1/sigma, y=zvalue
sblegend("topright",z=plot(a,xval="sigma",log="x")) # x=log(sigma)
sblegend("bottomright",z=plot(a,xval="sigma",
               yval="pval",log="xy")) # x=log(sigma), y=log(pvalue)
## multiple plots
b <- mam15.relltest[1:15] # an object of class "scalebootv"
plot(b) # x=1/sigma, y=zvalue
plot(b,xval="sigma",log="x") # x=log(sigma)

}
\keyword{ hplot }
\keyword{ models }
