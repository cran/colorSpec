% File calibrate.Rd

\name{calibrate}
\alias{calibrate}
\alias{calibrate.colorSpec}

\title{make a linear modification to a colorSpec responder}

\description{
make a linear modification to a \bold{colorSpec} responder with M spectra
so a specific single spectrum (the stimulus) creates a given specific response.
It is generalized \emph{white balance}.\cr
The options are complicated, but in all cases the returned object is
\code{multiply(x,mat)} where  \code{mat} is an internally calculated MxM matrix.
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
  
  \item{stimulus}{a \bold{colorSpec} object with a single spectrum, with \code{type}
  \code{'light'} or \code{'material'} to match \code{x}.  
  The wavelength sequence of \code{stimulus} must be equal to that of \code{x}.\cr
  If \code{stimulus} is \code{NULL}, then an appropriate default is chosen, see \bold{Details}.
  
  }
  
  \item{response}{an M-vector, or a scalar which is replicated to length M.
  All entries in \code{response} must be positive.\cr
  If \code{response} is \code{NULL}, then an appropriate default \emph{may be} chosen, see \bold{Details}.
  }
  
  \item{method}{an MxM \emph{adaption matrix}. 
  \code{method} can also be \code{'scaling'} and it is then set to the MxM identity matrix,
  which scales each responsivity spectrum in \code{x} independently. \cr
  If M=3, \code{method} can also be \code{'Bradford'}, \code{'Von Kries'}, or \code{'MCAT02'},
  and it is then set to the popular corresponding \emph{chromatic adaption matrix}.
  For these special matrices, the spectra in \code{x} are \bold{not} scaled independently;
  there is "cross-talk".\cr
  If \code{method} is \code{NULL}, then an appropriate default is chosen, see \bold{Details}.
  }
}
    
\value{
   a \bold{colorSpec} object equal to \code{multiply(x,mat)} 
   where \code{mat} is an internally calculated MxM matrix.  
   The \code{\link{quantity}} and \code{\link{wavelength}} are preserved.\cr
   Note that \code{mat} is not the same as the the  MxM \emph{adaption matrix}.
   To inspect \code{mat} execute \code{\link{summary}} on the returned object.
   If \code{method} is \code{'scaling'} then \code{mat} is diagonal and the
   diagonal entries are the M gain factors needed to achieve the calibration.\cr
   Useful data is attached as attribute \code{"calibrate"}.
}


\details{
  If \code{stimulus} is \code{NULL}, it is set to
   \code{\link{illuminantE}()} or \code{\link{neutralMaterial}()} to match \code{x}.
  
  If \code{response} is \code{NULL} and the response of \code{x} is \code{electrical} or \code{action},
  then \code{response} is set to an M-vector of all 1s.
  If \code{response} is \code{NULL} and the response of \code{x} is \code{neural},
  then this is an ERROR and the user is prompted to supply a specific \code{response}.\cr
 
  If \code{method} is \code{NULL} and M=3 and the response of \code{x} is \code{neural},
  then the neural response is assumed to be human,
  and the \code{method} is set to  the 3x3 Bradford matrix.\cr
  Otherwise \code{method} is set to  the MxM identity matrix,
  which scales each responsivity spectrum in \code{x} independently.
  In cameras this is usally called \emph{white balance},
  and so \code{calibrate()} can be considered a generalization of white balance.
}

\note{
Chromatic adaption transforms, such as \code{'Bradford'},
do not belong in the realm of spectra;
this is not really a spectral calculation.
For more about this subject see the explanation in \cite{Digital Color Management}.
This adaption option is provided in \code{calibrate} because it is possible and convenient.
}

\references{
Edward J. Giorgianni and Thomas E. Madden.
\bold{Digital Color Management: Encoding Solutions.} 2nd Edition
John Wiley. 2009.
Chapter 15 - Myths and Misconceptions.
}


\seealso{
\code{\link{is.regular}},
\code{\link{quantity}},
\code{\link{wavelength}},
\code{\link{colorSpec}},
\code{\link{summary}},
\code{\link{illuminantE}},
\code{\link{neutralMaterial}},
\code{\link{product}}
}

\examples{
# make an art gallery illuminated by illuminant A, and with tristimulus XYZ as output
gallery = product( A.1nm, 'artwork', xyz1931.1nm, wave='auto')      
# chromatically adapt the output XYZs to D50 white point, using Bradford matrix
gallery.D50 = calibrate( gallery, response=officialXYZ('D50') )

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
