\name{HigherPasserines}
\alias{HigherPasserines}
\docType{data}

\title{
Cone Fundamentals for the Higher Passerines
}

\description{
\tabular{ll}{
\code{HigherPasserines} \tab Tetrachromatic Cone Fundamentals of Higher Passerine Birds \cr
}}

\format{
A \bold{colorSpec} object organized as a matrix with the 4 spectra:
  \tabular{ll}{
    {\code{UV}} \tab {the UV wavelength responsivity} \cr
    {\code{Short}} \tab {the short wavelength responsivity} \cr
    {\code{Medium}} \tab {the medium wavelength responsivity} \cr 
    {\code{Long}} \tab {the long wavelength responsivity} \cr
  }
The {\code{wavelength}} is from 300 to 700 nm, at 1nm intervals.
}

\source{
\url{https://onlinelibrary.wiley.com/doi/10.1111/j.1095-8312.2005.00540.x/suppinfo}
}

\references{
Endler & Mielke.
Comparing entire colour patterns as birds see them.
Biological Journal of the Linnean Society.
Volume 86, Issue 4, pages 405-431, December 2005.
Original Name of File: BIJ_540_Endler_Mielke_OnlineAppendix.txt.
}

\seealso{
\code{\link{lms2000}} 
}

\examples{
summary(HigherPasserines)
}

\keyword{datasets}
\keyword{eyes}