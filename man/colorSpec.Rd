% File colorSpec.Rd

\name{colorSpec}
\title{constructing and testing colorSpec Objects}
\alias{colorSpec}
\alias{is.colorSpec}

\description{
  The function \code{colorSpec} is used to construct \bold{colorSpec} objects.

  \code{is.colorSpec} tests whether an object is a valid \bold{colorSpec} object.
}
\usage{
colorSpec( data, wavelength, quantity='auto', organization='auto' )

is.colorSpec(x)
}

\arguments{
  \item{data}{a vector or matrix of the spectrum values.
    In case \code{data} is a vector, there is a single spectrum and the number
    of points in that spectrum is the length of the vector.
    In case \code{data} is a matrix, the spectra are stored in the columns,
    so the number of points in each spectrum is the number of rows.
    It is OK for the matrix to have only 0 or 1 column.
    The column names (if any) are taken as the spectrum names.
    If no names are given, or if there are duplicate names,
    then \code{'S1', 'S2', ...} are used.
    Names can also be assigned after construction too; see \code{\link{specnames}}.
    Row names are ignored.
    Compare this function with \link{ts}.}
    
  \item{wavelength}{a numeric vector of wavelengths for all the spectra.
    The length of this vector must be equal to \code{\link[base]{NROW}(data)}.
    The sequence must be increasing.}
    
  \item{quantity}{a character string giving the \code{quantity} of all spectra;
    see \code{\link{quantity}} for a list of possible values.  
    In case of \code{'auto'}, a guess is made from the column names.
    This guess can be overridden later.}
    
  \item{organization}{a character string giving the desired organization 
    of the returned \bold{colorSpec} object.
    In case of \code{'auto'}, the organization is \code{'vector'} or \code{'matrix'}
    depending on \code{data}.
    Other possible organizations are \code{'df.col'} or \code{'df.row'}.
    The organization can be changed later, see \code{\link{organization}}
    for discussion of all 4 possible organizations.}
    
  \item{x}{an \R object to test for validity.}
    }
    
\value{
  \code{colorSpec} returns a \bold{colorSpec} object, or \code{NULL} in case of ERROR.
  
  \code{is.colorSpec} returns \code{TRUE} or \code{FALSE}.
           If \code{FALSE} it logs helpful reasons that \code{x} is invalid.
           
}    


\details{
  A \bold{colorSpec} object is either a \code{vector}, \code{matrix}, or \code{data.frame}.
  It is of S3 class \code{'colorSpec'} with these extra attributes:
  
  \describe{
  \item{\code{wavelength}}{a numeric vector of wavelengths for all the spectra.  If the organization of the object is \code{df.col}, then this is absent.}

  \item{\code{quantity} }{a character string that gives the physical quantity of all spectra, see \code{\link{quantity}} for a list of possible values.}

  \item{\code{metadata}}{a list for user-defined data.  The names \code{'path'}, \code{'header'} and \code{'date'} are already reserved; see \code{\link{metadata}}.}

  \item{\code{step.wl}}{step between adjacent wavelengths in nm. This is assigned only when the wavelengths are regular; see \code{\link{is.regular}}. }

  \item{\code{specname} }{only assigned when the \code{organization} is \code{'vector'}, in which case it is equal to the single character string name of the single spectrum.  See \code{\link{specnames}}.}
  
  \item{\code{sequence}}{only assigned when the object was returned from \code{\link{product}}. It is a list of the \bold{colorSpec} terms in this product.}
  
 \item{\code{calibration}}{only assigned when the object was returned from \code{\link{calibrate}}.}
    
}
}


\seealso{
  \code{\link{wavelength}},
  \code{\link{quantity}},
  \code{\link{metadata}},
  \code{\link{step.wl}},
  \code{\link{specnames}},
  \code{\link{is.regular}},
  \code{\link{coredata}}  
}

\examples{
#  make a synthetic Gaussian bandpass filter

center = 600
wave   = 400:700
trans  = exp( -(wave-center)^2 / 20^2 )

filter.bp   = colorSpec( trans, wave, 'transmittance' )

organization( filter.bp )  # returns: [1] "vector"

specnames( filter.bp ) = "myfilter"

# and now plot it
plot( filter.bp )
}

\keyword{colorSpec}
