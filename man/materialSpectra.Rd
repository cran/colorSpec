% File materialSpectra.Rd

\name{materialSpectra}
\alias{lensAbsorbance}
\alias{neutralMaterial}
\alias{materialSpectra}

\title{compute standard material spectra}

\description{
Compute hypothetical neutral gray material reflectance/transmittance,
and absorbance of the human lens, as a function of age.
}

\usage{
neutralMaterial( gray=1, wavelength=380:780 )
lensAbsorbance( age=32, wavelength=400:700 )
}

\arguments{
\item{gray}{a vector of gray levels, in the interval [0,1].  
        gray=1 represents the \emph{perfect reflecting diffuser}. }
        
\item{age}{a vector of ages in years; all ages must be >= 20. }    
        
\item{wavelength}{a vector of wavelengths for the returned object}  
}
    
\value{
\code{neutralMaterial} returns a
\bold{colorSpec} object with \code{\link{quantity}} equal to \code{'reflectance'}.
The reflectance of each spectrum is constant.
The number of spectra in the object is the number of gray levels = \code{length(gray)}.\cr

\code{lensAbsorbance} returns a
\bold{colorSpec} object with \code{\link{quantity}} equal to \code{'absorbance'}.
The absorbance model for the human lens is taken from \cite{Pokorny}.
The number of spectra in the object is the number of ages = \code{length(age)}.\cr

}

\references{
Pokorny, Joel, Vivianne C. Smith, and Margaret Lutze.
Aging of the Human Lens.
Applied Optics.
Vol. 26, No. 8.
15 April 1987.
Table I.  Page 1439.

}

\seealso{
\link{lightSpectra}, \code{\link{quantity}}
}

\examples{
#   make a perfect reflecting diffuser (PRD)
prd = neutralMaterial( 1 )

#   make a perfect transmitting filter (PTF)
ptf = prd
quantity(ptf) = 'transmittance'

#   compare transmittance at 3 ages: 20, 32, and 80 years
plot( linearize(lensAbsorbance( c(20,32,80) )), col='black', lty=1:3 )
}

\keyword{materials}
