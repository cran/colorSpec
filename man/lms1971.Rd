\encoding{UTF-8}
\name{lms1971}
\alias{lms1971.5nm}
\alias{lms1971}
\docType{data}

\title{Cone Fundamentals - 2-degree (1971)}
\description{
\tabular{ll}{
\code{lms1971.5nm} \tab the Vos & Walraven (1971) 2° cone fundamentals from 380 to 780 nm, at 5nm intervals \cr
}}

\format{
A \bold{colorSpec} object organized as a matrix with 3 columns:
  \tabular{ll}{
    {\code{long}}   \tab {the long wavelength responsivity} \cr
    {\code{medium}} \tab {the medium wavelength responsivity} \cr
    {\code{short}}  \tab {the short wavelength responsivity} \cr
  }
}

\source{
\url{http://www.cvrl.org/database/text/cones/vw.htm}
}

\references{
Vos, J. J. & Walraven, P. L.
On the derivation of the foveal receptor primaries. 
\bold{Vision Research}.
11 (1971) pp. 799-818.

}

\seealso{
\link{lms2000}
}

\examples{
summary(lms1971.5nm)
white.point = product( D65.1nm, lms1971.5nm, wave='auto' )
}
\keyword{datasets}
\keyword{eyes}
