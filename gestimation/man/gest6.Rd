\name{gest6}
\alias{gest6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gest6(y.mod, data, blip.mod, treat.mod, psi.initial = NULL, ...)
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
    ...) 
{
    if (is.null(psi.initial)) {
        ecfact.mod <- list(~l0, ~l0 + a1 + l1)
        gest2(y.mod, data, blip.mod, treat.mod, ecfact.mod, ...)
    }
    else {
        psi10 <- psi.initial[[1]][1]
        psi11 <- psi.initial[[1]][2]
        psi20 <- psi.initial[[2]][1]
        psi21 <- psi.initial[[2]][2]
        ecfact.mod <- list(~l0 * I(psi10 + psi11 * l0 > 0), ~l0 * 
            I(a1 == 0 & psi10 + psi11 * l0 > 0) + l0 * I(a1 == 
            1 & psi10 + psi11 * l0 < 0) + l1 * I(psi20 + psi21 * 
            l1 < 0))
        gest2(y.mod, data, blip.mod, treat.mod, ecfact.mod, ...)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
