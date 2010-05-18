\name{myplot}
\alias{myplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
myplot(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, ...) 
{
    require(lattice)
    mypanel <- function(x, y, ...) {
        if (is.factor(x)) {
            panel.bwplot(x, y, ..., horizontal = FALSE)
            xx <- 1:nlevels(x)
            yy <- tapply(y, x, median)
            panel.segments(xx[1:(nlevels(x) - 1)], yy[1:(nlevels(x) - 
                1)], xx[2:(nlevels(x))], yy[2:(nlevels(x))])
        }
        else {
            panel.xyplot(x, y, col = "gray60", ...)
            panel.loess(x, y, col = "blue", ...)
        }
        panel.abline(h = 0, col = "red")
    }
    xyplot(x, panel = mypanel, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
