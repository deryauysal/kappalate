\name{kappalate}
\alias{kappalate}
\title{Estimating the local average treatment effect (LATE) using Abadie's kappa approach and other weighting estimators}
\description{
\code{kappalate} estimates the local average treatment effect (LATE) using methods based on inverse probability weighting, most of which can be motivated by Abadie's (2003) kappa theorem. See Słoczyński, Uysal, and Wooldridge (2024) for a detailed treatment of the underlying theoretical results.}
\usage{
kappalate(given_formula, data, zmodel = "cbps", std = "on", which = "norm", subset = NULL)
}
\arguments{
  \item{given_formula}{
    A valid IV formula of the form: 'depvar ~ indepvars | treatment | instrument'.
    E.g., 'y ~ x1 + x2 | t | z'.
    Allows standard formula transformations directly e.g., 'log(y+1) ~ x1 + x2 | t | z'. 
    The components of this formula are as follows:
    \code{depvar}: The outcome variable.
    \code{indepvars}: A number of control variables.
    \code{treatment}: The treatment variable.
    \code{instrument}: The instrumental variable, which must be binary (0, 1). It is assumed to be valid conditional on \code{treatment} and needs to satisfy the monotonicity assumption.
    }
  \item{data}{Dataframe containing the variables used in the given_formula.}
  \item{zmodel}{Specifies the approach to estimate the instrument propensity score. Options are \code{"logit"}, \code{"probit"}, or \code{"cbps"} (default is \code{"cbps"}).}
  \item{std}{Determines whether nonbinary covariates are standardized before estimation. Options are \code{"on"} or \code{"off"}. Default is \code{"on"} which is stronly encouraged since convergence will often not be achieved otherwise.}
  \item{which}{Specifies whether to calculate all estimates or only normalized estimates. Options are \code{"all"} or \code{"norm"} (default is \code{"norm"}).} 
  \item{subset}{Specify only a subset of the dataframe to be used in the estimation. It should either be a numeric vector e.g. subset = c(1,2,3,4) where all specified index values must be part of the dataframe or  a logical vector e.g. subset = c(TRUE, FALSE, FALSE), with length = length of the dataframe.")}
  \item{pstolerance}{Specifies the tolerance for the overlap assumptions check. Must be a single numeric value strictly between 0 and 1. Default is \code{1e-5}. If the overlap assumption is violated the program will return vector of the row indices corresponding to the observations in the passed data frame where the assumption is violated.}
}

\details{
\code{kappalate} displays up to five estimates of the LATE, depending on the options chosen.
The naming convention for the alternative estimates follows Słoczyński, Uysal, and Wooldridge (2024), who recommend normalized estimators (option \code{which = "norm"}), and especially the estimator in Uysal (2011), referred to as \code{tau_u}, with covariate balancing propensity scores (option \code{zmodel = "cbps"}).
We also recommend standardizing nonbinary covariates before estimation (option \code{std = "on"}), as this may improve the performance of R's optimization process without affecting the estimation results.

\code{kappalate} allows for multivalued treatments, even though this case is not explicitly discussed by Słoczyński, Uysal, and Wooldridge (2024).
For multivalued treatments, \code{tau_u} can estimate the average causal response (ACR), a variant of the LATE, as per Frölich (2007).
}

\value{
Returns a list of class \code{kappalate} containing:
\item{coefficients}{A vector of LATE estimates.}
\item{vcov_matrix}{Variance-covariance matrix of the estimates.}
\item{N}{Number of observations used in the estimation.}
\item{metadata}{Contains information on the used estimation specifications.}
}

\examples{
# use the SIPP data
install.packages("haven")  # Install if you haven't
library(haven)

data <- read_dta("~/Documents/GitHub/kappalate/stata/sipp.dta")  # Replace with your actual file path

# Drop rows where kwage is missing (NA) or educ is missing (NA) or rsncode == 999
data <- data[!(is.na(data$kwage) | is.na(data$educ) | data$rsncode == 999), ]

# Create a new variable lwage as the natural log of kwage
data$lwage <- log(data$kwage)

kappalate_cbps <- kappalate(lwage ~ age_5 | nvstat | rsncode  , data = data, zmodel = "cbps", std = "on",  which = "all")
kappalate_logit <- kappalate(lwage ~ age_5 | nvstat | rsncode  , data = data, zmodel = "logit", std = "on",  which = "all")

print(kappalate_cbps)
print(kappalate_logit)
}

\references{
  Abadie, Alberto (2003). "Semiparametric Instrumental Variable Estimation of Treatment Response Models." \emph{Journal of Econometrics}, 113(2), 231–263.

  Frölich, Markus (2007). "Nonparametric IV Estimation of Local Average Treatment Effects with Covariates." \emph{Journal of Econometrics}, 139(1), 35–75.

  Heiler, Phillip (2022). "Efficient Covariate Balancing for the Local Average Treatment Effect." \emph{Journal of Business & Economic Statistics}, 40(4), 1569–1582.

  Słoczyński, Tymon, S. Derya Uysal, and Jeffrey M. Wooldridge (2024). "Abadie's Kappa and Weighting Estimators of the Local Average Treatment Effect." \emph{Journal of Business & Economic Statistics}, forthcoming.
  Available at \url{https://doi.org/10.1080/07350015.2024.2332763}.

  Uysal, S. Derya (2011). "Three Essays on Doubly Robust Estimation Methods." PhD dissertation, University of Konstanz.
}
