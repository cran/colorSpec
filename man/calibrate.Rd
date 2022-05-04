% File calibrate.Rd

\name{calibrate}
\alias{calibrate}
\alias{calibrate.colorSpec}

\title{make a linear modification to a colorSpec responder}

\description{
make a linear modification to a \bold{colorSpec} responder with M spectra,
so a specific stimulus (a single spectrum) creates a specific response (an M-vector).
It is generalized form of \emph{white balance}.\cr
The options are complicated, but in all cases the returned object is
\code{\link{multiply}(x,gmat)} where \code{gmat} is an internally calculated MxM matrix - called the \emph{gain matrix}.
Stated another way, the spectra in the output are linear combinations of spectra in
the input \code{x}.\cr
In case of ERROR, a message is logged and the original \code{x} is returned.
}

\usage{
\S3method{calibrate}{colorSpec}( x, stimulus=NULL, response=NULL, method=NULL )

}

\arguments{
  \item{x}{a \bold{colorSpec} responder with M spectra.
            The \code{type} must be \code{'responsivity.light'} or \code{'responsivity.material'}.}

  \item{stimulus}{a \bold{colorSpec} object with a single spectrum, with \code{type} either
  \code{'light'} or \code{'material'} to match \code{x}.
  The wavelength sequence of \code{stimulus} must be equal to that of \code{x}.\cr
  If \code{stimulus} is \code{NULL}, then an appropriate default is chosen, see \bold{Details}.
  }

  \item{response}{an M-vector, or a scalar which is then replicated to length M.
  Normally all entries are not \code{NA}, but it is OK to have exactly one that is not \code{NA}.
  In this special case, a single scaling factor is computed from that non-\code{NA} coordinate,
  and then applied to \emph{all} M coordinates; the \code{method} must be \code{'scaling'}.
  This is useful for the recommended method for calibration in ASTM E308-01 section 7.1.2.
  The same type of scaling is also recommended method in CIE 15: Technical Report section 7.1.
  In this case \code{response=c(NA,100,NA)} so the special coordinate is the luminance Y.
  See the \bold{Examples} below and the vignettes
  \href{../doc/gallery.pdf}{\bold{Viewing Object Colors in a Gallery}}
  and  
  \href{../doc/lens-aging.pdf}{\bold{The Effect of the Aging Human Lens on Color Vision}}.\cr
  All entries in \code{response}, that are not \code{NA}, must be positive.\cr  
  If \code{response} is \code{NULL}, then an appropriate default \emph{may be} chosen, see \bold{Details}.
  }

  \item{method}{an MxM \emph{adaptation matrix}.
  \code{method} can also be \code{'scaling'} and it is then set to the MxM identity matrix,
  which scales each responsivity spectrum in \code{x} independently.\cr
  If M=3, \code{method} can also be \code{'Bradford'}, \code{'Von Kries'}, \code{'MCAT02'}, or \code{'Bianco+Schettini'},
  and it is then set to the popular corresponding \emph{chromatic adaptation matrix}.
  For these special matrices, the spectra in \code{x} are \bold{not} scaled independently;
  there is "cross-talk".\cr
  If \code{method} is \code{NULL}, then an appropriate default is chosen, see \bold{Details}.
  }
}

\value{
   a \bold{colorSpec} object equal to \code{\link{multiply}(x,gmat)}
   where \code{gmat} is an internally calculated MxM matrix.
   The \code{\link{quantity}()} and \code{\link{wavelength}()} are preserved.\cr
   Note that \code{gmat} is not the same as the the MxM \emph{adaptation matrix}.
   To inspect \code{gmat} execute \code{\link{summary}()} on the returned object.
   If \code{method} is \code{'scaling'} then \code{gmat} is diagonal and the
   diagonal entries are the M gain factors needed to achieve the calibration.\cr
   Useful data is attached as attribute \code{"calibrate"}.
}


\details{
  If \code{stimulus} is \code{NULL}, it is set to
  \code{\link{illuminantE}()} or \code{\link{neutralMaterial}()} to match \code{x}.

  If \code{response} is \code{NULL} and the response of \code{x} is \code{electrical} or \code{action},
  then \code{response} is set to an M-vector of all 1s.
  If \code{response} is \code{NULL} and the response of \code{x} is \code{neural},
  then this is an ERROR and the user is prompted to supply a specific \code{response}.

  If \code{method} is \code{NULL}, its assignment is complicated.\cr
  If M=3 and the response of \code{x} is \code{neural}, 
  and the \code{specnames} of \code{x} partially match \code{c('x','y','z')} (case-insensitive),
  and none of the components of \code{response} are \code{NA},
  then the neural response is assumed to be human,   and the \code{method} is set to \code{'Bradford'}.\cr
  Otherwise \code{method} is set to \code{'scaling'}.
}

\note{
Chromatic adaptation transforms, such as \code{'Bradford'},
do not belong in the realm of spectra,
for this is not really a spectral calculation.
For more about this subject see the explanation in \cite{Digital Color Management},
Chapter 15 - Myths and Misconceptions.
These sophisticated adaptation transforms are provided in \code{calibrate()} because it is possible and convenient.
}

\references{
ASTM E308-01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
2001.

CIE 15: Technical Report: Colorimetry, 3rd edition.
CIE 15:2004.

Edward J. Giorgianni and Thomas E. Madden.
\bold{Digital Color Management: Encoding Solutions.} 2nd Edition
John Wiley. 2009.
Chapter 15 - Myths and Misconceptions.
}


\seealso{
\code{\link{is.regular}()},
\code{\link{multiply}()},
\code{\link{quantity}()},
\code{\link{wavelength}()},
\code{\link{colorSpec}},
\code{\link{summary}()},
\code{\link{illuminantE}()},
\code{\link{neutralMaterial}()},
\code{\link{product}()}
}

\examples{
wave = 380:780

# make an art gallery illuminated by illuminant A, and with tristimulus XYZ as output
gallery = product( A.1nm, 'artwork', xyz1931.1nm, wave=wave )

#  calibrate simplistically,
#  so the perfect reflecting diffuser has the standard XYZ coordinates for Illuminant A
#  using the convention that Y=100 (instead of Y=1)
A = 100 * spacesXYZ::standardXYZ('A')
A
##         X   Y      Z
##  A 109.85 100 35.585


gallery.cal1 = calibrate( gallery, response=A, method='scaling' )

#  calibrate following the ASTM and CIE recommendation
gallery.cal2 = calibrate( gallery, response=c(NA,100,NA), method='scaling' )

#   make the Perfect Reflecting Diffuser for testing
prd = neutralMaterial( 1, wave=wave ) ; specnames(prd) = 'PRD'

#   compare responses to the PRD for gallery.cal1 and gallery.cal2
white.1 = product( prd, gallery.cal1 )
white.2 = product( prd, gallery.cal2 )
white.1 ; white.2 ; white.1 - white.2 

##           X   Y      Z
##  PRD 109.85 100 35.585
##             X   Y        Z
##  PRD 109.8488 100 35.58151
##                X             Y           Z
##  PRD 0.001210456 -2.842171e-14 0.003489313


# make an RGB flatbead scanner from illuminant F11 and a Flea2 camera
scanner = product( subset(Fs.5nm,'F11'), 'paper', Flea2.RGB, wave='auto')
# adjust RGB gain factors (white balance) so the perfect reflecting diffuser yields RGB=(1,1,1)
scanner = calibrate( scanner )

# same flatbead scanner, but this time with some "white headroom"
scanner = product( subset(Fs.5nm,'F11'), 'paper', Flea2.RGB, wave='auto' )
scanner = calibrate( scanner, response=0.95 )
scanner
}

\keyword{colorSpec}
