\name{obj}
\alias{obj}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(obj)}
\format{
  The format is:
List of 8
 $ call      : language gest2(y.mod = y.mod, data = df1, blip.mod = blip.mod, treat.mod = treat.mod,      ecfact.mod = ecfact.mod)
 $ samplesize: int 500
 $ intervals : int 2
 $ psi.hat   :List of 2
  ..$ : Named num [1:2] 231.935 -0.973
  .. ..- attr(*, "names")= chr [1:2] "(Intercept)" "l0"
  ..$ : Named num [1:2] 704.74 -1.96
  .. ..- attr(*, "names")= chr [1:2] "(Intercept)" "l1"
 $ blip      :List of 2
  ..$ : num [1:500, 1] 0 0 113 -259 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
  ..$ : num [1:500, 1] 0 0 300 -613 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
 $ regret.sum:List of 2
  ..$ : num [1:500, 1] 0 0 0 613 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
  ..$ : num [1:500] 0 0 0 0 0 0 0 0 0 0 ...
 $ residuals :List of 2
  ..$ : num [1:500, 1] 34.28 21.39 -86.54 4.17 -49.96 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
  ..$ : num [1:500, 1] 47.6 54.1 -81.7 -94.1 -56.3 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
 $ fitted    :List of 2
  ..$ : num [1:500, 1] 844 1081 670 332 1474 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
  ..$ : num [1:500, 1] 830 1048 665 430 1481 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:500] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
 - attr(*, "class")= chr "gest"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(obj)
## maybe str(obj) ; plot(obj) ...
}
\keyword{datasets}
