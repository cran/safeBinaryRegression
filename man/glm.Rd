\name{glm}
\Rdversion{1.1}
\alias{glm}

\title{
Fitting Generalized Linear Models
}

\description{
This function overloads the \code{\link[stats]{glm}} function so that a check for the existence of the maximum likelihood estimate is computed before fitting a \sQuote{glm} with a binary response.
}

\usage{
glm(formula, family = gaussian, data, weights, subset, na.action, start = NULL, etastart, mustart, offset, control = glm.control(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ..., separation = c("find", "test"))
}

\arguments{
The arguments are identical to the arguments of the \code{\link[stats]{glm}} function provided in the \sQuote{stats} package with the exception of

\item{separation}{either \dQuote{find} or \dQuote{test}.  Both options prevent the model from being fit to binary data when the maximum likelihood estimate does not exist.  Additionally, when \code{separation = "find"}, the terms separating the sample points are identified when the maximum likelihood estimate is found not to exist.}

The following arguments are passed to the \code{\link[stats]{glm}} function:

\item{formula}{see \code{\link[stats]{glm}}}
\item{family}{see \code{\link[stats]{glm}}}
\item{data}{see \code{\link[stats]{glm}}}
\item{weights}{see \code{\link[stats]{glm}}}
\item{subset}{see \code{\link[stats]{glm}}}
\item{na.action}{see \code{\link[stats]{glm}}}
\item{start}{see \code{\link[stats]{glm}}}
\item{etastart}{see \code{\link[stats]{glm}}}
\item{mustart}{see \code{\link[stats]{glm}}}
\item{offset}{see \code{\link[stats]{glm}}}
\item{control}{see \code{\link[stats]{glm}}}
\item{model}{see \code{\link[stats]{glm}}}
\item{method}{see \code{\link[stats]{glm}}}
\item{x}{see \code{\link[stats]{glm}}}
\item{y}{see \code{\link[stats]{glm}}}
\item{contrasts}{see \code{\link[stats]{glm}}}
\item{...}{see \code{\link[stats]{glm}}}
}

\details{
This function checks for the existence of the maximum likelihood estimate before the \sQuote{glm} function is used to fit binary regression models by solving the linear program proposed in Konis (2007).
}

\value{
See the return value for the \code{\link[stats]{glm}} function.
}

\references{
Kjell Konis (2007). Linear programming algorithms for detecting separated data in binary logistic regression models. DPhil, University of Oxford \url{http://ora.ouls.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a}
}

\author{
Kjell Konis \email{kjell.konis@epfl.ch}
}

\seealso{
\code{\link[stats]{glm}}.
}

\examples{
## A set of 4 completely separated sample points ##
x <- c(-2, -1, 1, 2)
y <- c(0, 0, 1, 1)

\dontrun{glm(y ~ x, family = binomial)}


## A set of 4 quasicompletely separated sample points ##
x <- c(-2, 0, 0, 2)
y <- c(0, 0, 1, 1)

\dontrun{glm(y ~ x, family = binomial)}
}

\keyword{regression}
