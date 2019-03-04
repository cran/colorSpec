% File computeCCT.Rd
\encoding{UTF-8}

\name{computeSSI}
\title{Compute the Spectrum Similarity Index of light spectra}
\alias{computeSSI}
\alias{computeSSI.colorSpec}


\description{
Compute the Spectrum Similarity Index (SSI), an index between 0 and 100,
of a \bold{colorSpec} object with \code{type} equal to \code{'light'}.
It compares a test spectrum with a reference spectrum (an ideal).
The value 100 means a perfect match to the reference, 
and a smaller value mean a poorer match (similar to CRI).
Only values in the interval [375,675] nm are used; for details see \cite{Holm}.
}

\usage{
\S3method{computeSSI}{colorSpec}( x, reference=NULL, digits=0, isotherms='mccamy', locus='robertson' )
}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type} equal to \code{'light'}, and M test spectra }
  
\item{reference}{a \bold{colorSpec} object with \code{type} equal to \code{'light'},
and either 1 or M reference spectra.
\code{reference} can also be \code{NULL} (the default),
which means to generate each reference spectrum from the corresponding test spectrum.}

\item{digits}{the number of digits after the decimal point in the returned vector.
According to \cite{Holm} the output should be rounded to the nearest integer,
which corresponds to \code{digits=0}.
To return full precision, set \code{digits=Inf}.}

\item{isotherms}{this is only used when \code{reference=NULL}.
It is passed to \code{\link{computeCCT}()}
in order to compute the CCT of each test spectrum.}

\item{locus}{this is only used when \code{reference=NULL}.
It is passed to \code{\link{computeCCT}()}
in order to compute the CCT of each test spectrum.}

}



\details{
If \code{reference} contains a single spectrum, 
then each test spectrum is compared to that one reference spectrum.
If \code{reference} contains M spectra,
then the i'th test spectrum is compared to the i'th reference spectrum.

If \code{reference=NULL} then for each test spectrum
the CCT is computed and used to compute a reference spectrum with the same CCT.
It is either a Planckian (black-body) or daylight illuminant, see \cite{Holm} for details.
The test spectrum and auto-computed reference spectrum are then compared.
}

\value{
\code{computeSSI()} returns a numeric vector of length M,
where M is the number of spectra in \code{x}.
The vector's \code{names} is set from \code{specnames(x)}
and a compact code for the corresponding reference spectrum.

If the \code{type} of \code{x} is not \code{'light'},
or \code{reference} is not valid,
then the function returns NULL.
}




\references{
J. Holm and T. Maier and P. Debevec and C. LeGendre and J. Pines and J. Erland and G. Joblove and S. Dyer and B. Sloan and J. di Gennaro and D. Sherlock.
\emph{A Cinematographic Spectral Similarity Index}.
SMPTE 2016 Annual Technical Conference and Exhibition.
pp. 1-36.
(2016).
}

\seealso{
\code{\link{type}()}, 
\code{\link{computeCCT}()},
\code{\link{planckSpectra}()},    
\code{\link{daylightSpectra}()}, 
\code{\link{specnames}()}
}

\examples{
computeSSI( planckSpectra( 1000*(2:6) )  )
##  P2000_SSI[2027K] P3000_SSI[3057K] P4000_SSI[D4063] P5000_SSI[D5061] P6000_SSI[D6063] 
##                99               98               93               92               92 
}


\keyword{light}

