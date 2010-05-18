\name{simulate.dataset}
\alias{simulate.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulate.dataset(n = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (n = 500) 
{
    l0 <- rnorm(n, 450, 150)
    l1 <- rnorm(n, 1.25 * l0, 60)
    p1 <- expit(2 - 0.006 * l0)
    p2 <- expit(0.8 - 0.004 * l1)
    a1 <- rbinom(n, 1, p1)
    a2 <- pmax(a1, rbinom(n, 1, p2))
    psi10 <- 250
    psi11 <- -1
    psi20 <- 720
    psi21 <- -2
    gamma1 <- a1 * (psi10 + psi11 * l0)
    gamma2 <- a2 * (psi20 + psi21 * l1)
    mu1 <- pmax(psi10 + psi11 * l0, 0) - gamma1
    mu2 <- pmax(psi20 + psi21 * l1, 0) - gamma2
    y <- rnorm(n, 400 + 1.6 * l0, 60) - mu1 - mu2
    df1 <- data.frame(l0, a1, l1, a2, y)
    return(df1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
