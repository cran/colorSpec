% File bandSpectra.Rd

\name{bandSpectra}
\alias{bandMaterial}
\alias{bandRepresentation}
\alias{bandRepresentation.colorSpec}

\title{Compute Band-based Material Spectra, and Bands for Existing Material Spectra}

\description{
A band-based material spectrum is a superimposition of bandpass filters, and (optionally) a bandstop filter.
The 2 functions in this topic convert a vector of numbers between 0 and 1 to a band representation,
and back again.
}

\usage{
bandMaterial( lambda, wavelength=380:780 )

\S3method{bandRepresentation}{colorSpec}( x )
}

\arguments{
\item{lambda}{a numeric Mx2 matrix with wavelength pairs in the rows,
or a vector that can be converted to such a matrix, by row.
The two wavelengths in a row (the \emph{transition wavelengths}) define either a bandpass or bandstop filter,
and all the rows are superimposed to define the transmittance spectrum of the final material.
If the 2 wavelengths are denoted by \eqn{\lambda_1} and \eqn{\lambda_2},
and \eqn{\lambda_1 < \lambda_2} then the filter is a bandpass filter.
If the 2 wavelengths are swapped, then the spectrum is "flipped"
and is a bandstop filter, and the band "wraps around" from long wavelengths to short.
There can be at most 1 bandstop filter in the matrix, otherwise it is an error.
The bands must be pairwise disjoint, otherwise it is an error.
To get a material with transmittance identically 0, set \code{lambda} to a 0x2 matrix.
To get a material with transmittance identically 1, 
set \code{lambda} to a 1x2 matrix with
\eqn{\lambda_1=\beta_0} and \eqn{\lambda_2=\beta_N}, where \eqn{N} is the number of wavelengths.
See vignette \href{../doc/convexity.pdf}{\bold{Convexity and Transitions}}
for the definition of \eqn{\beta_0} and  \eqn{\beta_N}
and other mathematical details.
\cr
\code{lambda} can also be a list of such matrices, which are processed separately, see \bold{Value}.
}
    
\item{wavelength}{a vector of wavelengths for the returned object}  

\item{x}{a \bold{colorSpec} object with \code{type} equal to \code{'material'}}
}
    

\value{

\code{bandMaterial()} returns a
\bold{colorSpec} object with \code{\link{quantity}} equal to \code{'transmitance'}.
If \code{lambda} is a matrix, then the object has 1 spectrum.
If \code{lambda} is a list of matrices with length N, then the object has N spectra.

\code{bandRepresentation()}
returns a list of matrices with 2 columns.
There is a matrix in the list for each spectrum in \code{x}.
}

\details{
\code{bandRepresentation()} is a right-inverse of \code{bandMaterial()}, 
see \bold{Examples} and the test script \code{test-bands.R}. 
For more mathematical details, see the vignette
\href{../doc/convexity.pdf}{\bold{Convexity and Transitions}}.
}



\seealso{
\code{\link{rectangularMaterial}()},
vignette \href{../doc/convexity.pdf}{\bold{Convexity and Transitions}}
}

\examples{
#  make a vector superimposing a bandpass and a bandstop filter, and of the proper length 401
vec = c( rep(1,100), 0.5, rep(0,40), .25, rep(1,50), 0.9, rep(0,100), 0.4, rep(1,107) )

#	 convert that vector to a colorSpec object, with a single spectrum
spec = colorSpec( vec, wavelength=380:780, quantity='transmittance', specnames='sample' )

#	 extract and print the 2 bands
lambda = bandRepresentation( spec ) ;  print(lambda)

##  $sample
##      lambda1 lambda2
##  BS   673.10   480.0
##  BP1  521.25   572.4

#  convert the 2 bands  (the transition wavelengths) back to a vector of length 401
#  and compare with the original vector
delta = vec - coredata( bandMaterial(lambda) )

range(delta)
##  [1] -9.092727e-14  2.275957e-14
}

\keyword{materials}
