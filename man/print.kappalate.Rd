\name{print.kappalate}
\alias{print.kappalate}
\title{Print Method for kappalate Objects}

\description{
Prints detailed results for objects of class \code{kappalate}, including coefficient estimates, standard errors, z-values, p-values, and 95% confidence intervals in a Stata-like format.
}

\usage{
\method{print}{kappalate}(x, ...)
}

\arguments{
  \item{x}{An object of class \code{kappalate}, returned by the \code{kappalate} function.}
  \item{...}{Further arguments (not used).}
}

\examples{
# SIPP data, as in the Stata help file (requires the haven package)
if (requireNamespace("haven", quietly = TRUE)) {
  data <- haven::read_dta(system.file("extdata", "sipp.dta", package = "kappalate"))
  data <- as.data.frame(data)
  data <- data[!(is.na(data$kwage) | is.na(data$educ) | data$rsncode == 999), ]
  data$lwage <- log(data$kwage)

  fit <- kappalate(lwage ~ age_5 | nvstat | rsncode, data = data, zmodel = "cbps", which = "all")
  print(fit)
}
}
