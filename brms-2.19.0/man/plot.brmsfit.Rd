% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.brmsfit}
\alias{plot.brmsfit}
\title{Trace and Density Plots for MCMC Draws}
\usage{
\method{plot}{brmsfit}(
  x,
  pars = NA,
  combo = c("dens", "trace"),
  N = 5,
  variable = NULL,
  regex = FALSE,
  fixed = FALSE,
  theme = NULL,
  plot = TRUE,
  ask = TRUE,
  newpage = TRUE,
  ...
)
}
\arguments{
\item{x}{An object of class \code{brmsfit}.}

\item{pars}{Deprecated alias of \code{variable}.
Names of the parameters to plot, as given by a
character vector or a regular expression.}

\item{combo}{A character vector with at least two elements.
Each element of \code{combo} corresponds to a column in the resulting
graphic and should be the name of one of the available
\code{\link[bayesplot:MCMC-overview]{MCMC}} functions
(omitting the \code{mcmc_} prefix).}

\item{N}{The number of parameters plotted per page.}

\item{variable}{Names of the variables (parameters) to plot, as given by a
character vector or a regular expression (if \code{regex = TRUE}). By
default, a hopefully not too large selection of variables is plotted.}

\item{regex}{Logical; Indicates whether \code{variable} should
be treated as regular expressions. Defaults to \code{FALSE}.}

\item{fixed}{(Deprecated) Indicates whether parameter names
should be matched exactly (\code{TRUE}) or treated as
regular expressions (\code{FALSE}). Default is \code{FALSE}
and only works with argument \code{pars}.}

\item{theme}{A \code{\link[ggplot2:theme]{theme}} object
modifying the appearance of the plots.
For some basic themes see \code{\link[ggplot2:ggtheme]{ggtheme}}
and \code{\link[bayesplot:theme_default]{theme_default}}.}

\item{plot}{Logical; indicates if plots should be
plotted directly in the active graphic device.
Defaults to \code{TRUE}.}

\item{ask}{Logical; indicates if the user is prompted
before a new page is plotted.
Only used if \code{plot} is \code{TRUE}.}

\item{newpage}{Logical; indicates if the first set of plots
should be plotted to a new page.
Only used if \code{plot} is \code{TRUE}.}

\item{...}{Further arguments passed to
\code{\link[bayesplot:MCMC-combos]{mcmc_combo}}.}
}
\value{
An invisible list of
  \code{\link[gtable:gtable]{gtable}} objects.
}
\description{
Trace and Density Plots for MCMC Draws
}
\examples{
\dontrun{
fit <- brm(count ~ zAge + zBase * Trt
           + (1|patient) + (1|visit),
           data = epilepsy, family = "poisson")
plot(fit)
## plot population-level effects only
plot(fit, variable = "^b_", regex = TRUE)
}

}
