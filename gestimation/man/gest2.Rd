\name{gest2}
\alias{gest2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gest2(y.mod, data, blip.mod, treat.mod, ecfact.mod, keep = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.mod}{
%%     ~~Describe \code{y.mod} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{blip.mod}{
%%     ~~Describe \code{blip.mod} here~~
}
  \item{treat.mod}{
%%     ~~Describe \code{treat.mod} here~~
}
  \item{ecfact.mod}{
%%     ~~Describe \code{ecfact.mod} here~~
}
  \item{keep}{
%%     ~~Describe \code{keep} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y.mod, data, blip.mod, treat.mod, ecfact.mod, keep = TRUE) 
{
    y <- model.response(model.frame(y.mod, data))
    n <- length(y)
    k <- length(blip.mod)
    if (length(treat.mod) != k || length(ecfact.mod) != k) {
        stop("Lengths of blip.mod, treat.mod and ecfact.mod must all equal number of time intervals.")
    }
    kappa <- rep(0, n)
    psi.hat <- list()
    if (keep) {
        blip <- list()
        regret.sum <- list()
        residuals <- list()
        fitted <- list()
    }
    for (j in k:1) {
        z <- model.matrix(blip.mod[[j]], data)
        a <- model.response(model.frame(treat.mod[[j]], data))
        s <- a * z
        E.s <- fitted(glm(treat.mod[[j]], data = data, family = binomial)) * 
            z
        w <- model.matrix(ecfact.mod[[j]], data)
        w.qr <- qr(w)
        temp1 <- t(s - E.s) \%*\% qr.resid(w.qr, y + kappa)
        temp2 <- t(s - E.s) \%*\% qr.resid(w.qr, s)
        psi.hat[[j]] <- drop(qr.solve(temp2, temp1))
        temp3 <- z \%*\% psi.hat[[j]]
        gamma.hat <- a * temp3
        if (keep) {
            blip[[j]] <- gamma.hat
            regret.sum[[j]] <- kappa
            residuals[[j]] <- qr.resid(w.qr, y - gamma.hat + 
                kappa)
            fitted[[j]] <- qr.fitted(w.qr, y - gamma.hat + kappa) - 
                kappa + gamma.hat
        }
        kappa <- kappa + pmax(temp3, 0) - gamma.hat
    }
    obj <- list()
    class(obj) <- "gest"
    obj$call <- match.call(expand.dots = TRUE)
    obj$samplesize <- n
    obj$intervals <- k
    obj$psi.hat <- psi.hat
    if (keep) {
        obj$blip <- blip
        obj$regret.sum <- regret.sum
        obj$residuals <- residuals
        obj$fitted <- fitted
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
