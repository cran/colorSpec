\encoding{UTF-8}
\name{lms2000}
\alias{lms2000.1nm}
\alias{lms2000}
\docType{data}
\title{Cone Fundamentals - 2-degree (2000)}

\description{
\tabular{ll}{
\code{lms2000.1nm} \tab the Stockman and Sharpe (2000) 2° cone fundamentals from 390 to 830 nm, at 1nm intervals \cr
}
}

\format{
A \bold{colorSpec} object organized as a matrix with 3 columns:
  \tabular{ll}{
    {\code{long}}   \tab {the long wavelength responsivity} \cr
    {\code{medium}} \tab {the medium wavelength responsivity} \cr
    {\code{short}}  \tab {the short wavelength responsivity} \cr
  }
}


\source{
http://www.cvrl.org/cones.htm
}

\references{
Stockman, A., Sharpe, L. T., & Fach, C. C. (1999). 
The spectral sensitivity of the human short-wavelength cones. 
Vision Research, 39, 2901-2927.

Stockman, A., & Sharpe, L. T. (2000). 
Spectral sensitivities of the middle- and long-wavelength sensitive cones 
derived from measurements in observers of known genotype. 
Vision Research, 40, 1711-1737.

}

\seealso{
\code{\link{lms1971}} 
}

\examples{
summary(lms2000.1nm)
white.point = product( D65.1nm, lms2000.1nm, wave='auto' )
}
\keyword{datasets}
\keyword{eyes}
