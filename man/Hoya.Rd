% File Hoya.Rd

\name{Hoya}
\alias{Hoya}
\docType{data}

\title{
standard Hoya filters
}

\description{
\tabular{ll}{
\code{Hoya} \tab 4 standard Hoya filters; from 300 to 750 nm at 10nm intervals.\cr
}
}

\format{
A \bold{colorSpec} object with \code{\link{quantity}} equal to \code{'transmittance'} and 4 spectra:

\tabular{ll}{
\code{R-60}   \tab   long-pass red filter with cutoff about 600nm\cr
\code{G-533}  \tab   band-pass green filter with peak about 533nm\cr
\code{B-440}  \tab   band-pass blue filter with peak about 440nm\cr
\code{LB-120} \tab   Light-balancing Blue filter with mired shift equal to -120 \cr
}
}

\source{
http://www.hoyaoptics.com/
}

\seealso{
\code{\link{quantity}}
}

\examples{
#   compute response of ACES scanner to the Hoya filters
product( Hoya, scanner.ACES, wave='auto' )
}

\keyword{datasets}
\keyword{materials}
