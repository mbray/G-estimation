\name{gest.boot}
\alias{gest.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gest.boot(data, gest.fn, B = 200, verbose = FALSE, store = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{gest.fn}{
%%     ~~Describe \code{gest.fn} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{store}{
%%     ~~Describe \code{store} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (data, gest.fn, B = 200, verbose = FALSE, store = FALSE, 
    ...) 
{
    psi.boot <- list()
    for (i in 1:B) {
        data.boot <- data[sample(1:nrow(data), replace = TRUE), 
            ]
        psi.boot[[i]] <- gest.fn(data = data.boot, ...)$psi.hat
        if (verbose) 
            cat(i, "..")
    }
    if (verbose) 
        cat("\n")
    psi.boot <- do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), 
        psi.boot)
    psi.covmat <- lapply(psi.boot, var)
    obj <- gest.fn(data = data, ...)
    obj$psi.covmat <- psi.covmat
    if (store) {
        obj$psi.boot <- psi.boot
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
