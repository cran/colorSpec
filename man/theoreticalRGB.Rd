% File theoreticalRGB.Rd

\name{theoreticalRGB}
\alias{theoreticalRGB}
\alias{BT.709.RGB}
\alias{Adobe.RGB}
\alias{ACES.RGB}
\docType{data}

\title{
Theoretical RGB Cameras - BT.709.RGB, Adobe.RGB, and ACES.RGB
}

\description{
These are 3 built-in \bold{colorSpec} objects, 
with \code{\link{quantity}} equal to \code{'energy->electrical'}.

\describe{
\item{\code{BT.709.RGB}}{a theoretical RGB responder to light.  
The 3 responsivity spectra are constructed so that the RGBs from this theoretical camera, when displayed on an sRGB display after correct EOTF adjustment, would emit light with the same XYZs as the captured scene (up to a constant multiple). 
All three responsivities have negative lobes.}

\item{\code{Adobe.RGB}}{a theoretical RGB responder to light.  
The 3 responsivity spectra are constructed so that the RGBs from this theoretical camera, when displayed on an Adobe RGB display after correct EOTF adjustment, would emit light with the same XYZs as the captured scene (up to a constant multiple).
All three responsivities have negative lobes.}

\item{\code{ACES.RGB}}{a theoretical RGB responder to light.  
Unlike the two above cameras, the responsivities are non-negative and so this camera
could be built, in principle.
These are the ACES RICD (Reference Input Capture Device)  spectral sensitivities.}
}
}

\format{
All are \bold{colorSpec} objects with \code{\link{quantity}} equal to \code{'energy->electrical'} and 3 spectra:
\code{r}, \code{g}, and \code{b}.
The wavelengths are 360 to 830 nm at 1 nm intervals.
}


\details{
All responsitivity spectra are linear combinations of the spectra in \code{\link{xyz1931.1nm}}.
These 3 theoretical cameras satisfy the \emph{Maxwell-Ives criterion} by construction.
\cr
For \code{BT.709.RGB} and \code{Adobe.RGB},
the responsivities are scaled so the response to \code{\link{D65.1nm}} is RGB=(1,1,1).
These responsivities have negative lobes.
\cr
The BT.709 primaries and white point are the same as those of sRGB (though the EOTF functions are different).  \cr
Adobe RGB and sRGB share the same Red, Blue, and White chromaticities, and only differ in the Green.
This implies that for both cameras the Green output is 0 at Red and Blue, and 1 at White.
This in turn implies that the Green output is identical for both cameras for all input spectra,
and so the Green responsivity spectra are identical for both cameras.
\cr
For \code{ACES.RGB} the responsivities are area normalized
as in Annex C of S-2008-001.
They are scaled so that the response to Illuminant E is RGB=(1,1,1).
For an example of white-balancing, as in Annex B,
see the examples below.
}


\references{
Poynton, Charles.
\bold{Digital Video and HD - Algorithms and Interfaces}.
Morgan Kaufmann.
Second Edition. 2012.
Figure 26.5 on page 302.

\bold{Academy Color Encoding Specification (ACES)}.
S-2008-001.
2011.
Annex B, pp. 23-25.
Annex C, pp. 26-33.
}


\seealso{
\code{\link{quantity}()},
\code{\link{D65.1nm}},
\code{\link{xyz1931.1nm}},
\code{\link{ptransform}()},
\code{\link{calibrate}()},
vignette \href{../doc/blueflame.pdf}{\bold{Blue Flame and Green Comet}}
}

\examples{
#######    BT.709.RGB is created using the following recipe  ########
P = matrix( c(0.64,0.33,NA,  0.3,0.6,NA, 0.15,0.06,NA ), 3, 3, byrow=TRUE )
rownames(P) = c('R','G','B')    
BT.709.RGB  = ptransform( xyz1931.1nm, P, D65.1nm )  
quantity(BT.709.RGB) = "energy->electrical"

#######    Adobe.RGB recipe is the same, except for the matrix P  ########
P = matrix( c(0.64,0.33,NA,  0.21,0.71,NA, 0.15,0.06,NA ), 3, 3, byrow=TRUE )
rownames(P) = c('R','G','B')    
Adobe.RGB  = ptransform( xyz1931.1nm, P, D65.1nm )  
quantity(Adobe.RGB) = "energy->electrical"

#######  white-balancing ACES.RGB for CIE Standard Illuminant D60 ########
# in a scene illuminated by daylight illuminant D6003,
# and with a perfect-reflecting-diffuser in that scene,
# object 'camera1' would have response RGB=(1,1,1) for that diffuser.
D6003 = daylightSpectra( 6000*1.4388/1.4380, wavelength=wavelength(ACES.RGB) )
camera1 = calibrate( ACES.RGB, D6003, 1, method='scaling' )
}


\keyword{datasets}
