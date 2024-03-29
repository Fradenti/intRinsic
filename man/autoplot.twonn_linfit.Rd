% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoplot_twonn.R
\name{autoplot.twonn_linfit}
\alias{autoplot.twonn_linfit}
\title{Plot the output of the \code{TWO-NN} model estimated via least squares}
\usage{
\method{autoplot}{twonn_linfit}(object, title = "TWO-NN Linear Fit", ...)
}
\arguments{
\item{object}{object of class \code{twonn_linfit}, the output of the
\code{twonn} function when \code{method = "linfit"}.}

\item{title}{string used as title of the plot.}

\item{...}{other arguments passed to specific methods.}
}
\value{
a \code{\link[ggplot2]{ggplot2}} object displaying the goodness of
the linear fit of the TWO-NN model.
}
\description{
Use this method without the \code{.twonn_linfit} suffix.
The function returns the representation of the linear
regression that is fitted with the \code{linfit} method.
}
\seealso{
\code{\link{twonn}}

Other autoplot methods: 
\code{\link{autoplot.Hidalgo}()},
\code{\link{autoplot.gride_bayes}()},
\code{\link{autoplot.twonn_bayes}()},
\code{\link{autoplot.twonn_mle}()}
}
\concept{autoplot methods}
