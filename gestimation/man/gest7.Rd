\name{gest7}
\alias{gest7}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gest7(y.mod, data, blip.mod, treat.mod, psi.initial = NULL, verbose = FALSE, ...)
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
  \item{psi.initial}{
%%     ~~Describe \code{psi.initial} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (y.mod, data, blip.mod, treat.mod, psi.initial = NULL, 
    verbose = FALSE, ...) 
{
    psi <- gest6(y.mod, data, blip.mod, treat.mod, psi.initial, 
        ...)$psi.hat
    if (verbose) 
        cat("Initial: ", unlist(psi), "\n")
    psi.new <- gest6(y.mod, data, blip.mod, treat.mod, psi, ...)$psi.hat
    iter <- 1
    while (any(abs((unlist(psi.new) - unlist(psi))/unlist(psi)) > 
        1e-06)) {
        if (verbose) 
            cat("Iteration: ", unlist(psi.new), "\n")
        psi <- psi.new
        psi.new <- gest6(y.mod, data, blip.mod, treat.mod, psi, 
            ...)$psi.hat
        if (iter >= 50) 
            break
        iter <- iter + 1
    }
    gest6(y.mod, data, blip.mod, treat.mod, psi, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
