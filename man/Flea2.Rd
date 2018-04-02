% File Flea2.Rd

\name{Flea2.RGB}
\alias{Flea2.RGB}
\alias{Flea2}
\docType{data}

\title{
Flea2 Camera FL2-14S3C from Point Grey
}
\description{
\tabular{ll}{
\code{Flea2.RGB} \tab an RGB responder to light, from 360 to 800 nm at 10 nm intervals \cr
}}

\format{
A \bold{colorSpec} object with \code{quantity} equal to \code{'energy->electrical'} and 3 spectra:
\code{Red}, \code{Green}, and \code{Blue}.
\cr
}

\details{
This data is read from the file \bold{Flea2-spectral.txt} which was digitized from the
plot in \bold{Flea2-spectral.png}.
\cr
}

\source{
\url{http://www.ptgrey.com/support/downloads/10117}
\cr
}

\seealso{
\code{\link{quantity}},
vignette \bold{blueflame}
}

\examples{
#  Make a scanner from a tungsten source and a Flea2 camera
Flea2.scanner = product( A.1nm, "VARMATERIAL", Flea2.RGB, wavelength=420:680 )
Flea2.scanner = calibrate( Flea2.scanner )
}

\keyword{datasets}
\keyword{cameras}
