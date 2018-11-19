% File colorSpec.Rd

\name{colorSpec}
\title{constructing and testing colorSpec Objects}
\alias{colorSpec}
\alias{is.colorSpec}
\alias{as.colorSpec}
\alias{as.colorSpec.default}

\description{
The function \code{colorSpec()} is the constructor for \bold{colorSpec} objects.

\code{is.colorSpec()} tests whether an object is a valid \bold{colorSpec} object.\cr
\code{as.colorSpec()} converts other variables to a \bold{colorSpec} object,
and is designed to be overridden by other packages.
}



\usage{
colorSpec( data, wavelength, quantity='auto', organization='auto', specnames=NULL )

is.colorSpec(x)

\S3method{as.colorSpec}{default}( \ldots )
}

\arguments{
  \item{data}{a vector or matrix of the spectrum values.
    In case \code{data} is a vector, there is a single spectrum and the number
    of points in that spectrum is the length of the vector.
    In case \code{data} is a matrix, the spectra are stored in the columns,
    so the number of points in each spectrum is the number of rows in \code{data},
    and the number of spectra is the number of columns in \code{data}.
    It is OK for the matrix to have only 0 or 1 column.}

  \item{wavelength}{a numeric vector of wavelengths for all the spectra, in nm.
    The length of this vector must be equal to \code{\link[base]{NROW}(data)}.
    The sequence must be increasing.
    The \code{wavelength} vector can be changed after construction.}
    
  \item{quantity}{a character string giving the \code{quantity} of all spectra
    in \code{data}; see \code{\link{quantity}} for a list of possible values.  
    In case of \code{'auto'}, a guess is made from the \code{specnames}.
    The \code{quantity} can be changed after construction.}
    
  \item{organization}{a character string giving the desired organization 
    of the returned \bold{colorSpec} object.
    In case of \code{'auto'}, the organization is \code{'vector'} or \code{'matrix'}
    depending on \code{data}.
    Other possible organizations are \code{'df.col'} or \code{'df.row'};
    see \code{\link{organization}} for discussion of all 4 possible organizations.
    The organization can be changed after construction.
    }
    
  \item{specnames}{a character vector with length equal to the number of spectra,
  and with no duplicates.
  If \code{specnames} is \code{NULL} and \code{data} is a vector,
  then \code{specnames} is set to \code{deparse(substitute(data))}.
  If \code{specnames} is \code{NULL} and \code{data} is a matrix,
  then \code{specnames} is set to \code{colnames(data)}.
  If \code{specnames} is \emph{still} not a character vector with the right length,
  or if there are duplicate names, then \code{specnames} is set to
  \code{'S1', 'S2', ...} with a warning message.
  The \code{specnames} vector can be changed after construction.
  }  
  
  \item{x}{an \R object to test for being a valid \bold{colorSpec} object.}
  
  \item{...}{arguments for use in other packages.}  
}
    
\value{
\code{colorSpec()} returns a \bold{colorSpec} object, or \code{NULL} in case of ERROR.
Compare this function with \code{stats::\link{ts}()}.
  
\code{is.colorSpec()} returns \code{TRUE} or \code{FALSE}.
It does more than check the \code{class},
and also checks \code{wavelength}, \code{quantity}, and \code{organization}.
If \code{FALSE}, it logs (at \code{loglevel='DEBUG'})
the reason that \code{x} is invalid.
           
\code{as.colorSpec.default()} issues an ERROR message and returns \code{NULL}
}    


\details{
  Under the hood, a \bold{colorSpec} object is either a \code{vector}, \code{matrix}, or \code{data.frame}.
  It is of S3 class \code{'colorSpec'} with these extra attributes:
  
  \describe{
  \item{\code{wavelength}}{a numeric vector of wavelengths for all the spectra.  If the organization of the object is \code{'df.col'}, then this is absent.}

  \item{\code{quantity} }{a character string that gives the physical quantity of all spectra, see \code{\link{quantity}()} for a list of possible values.}

  \item{\code{metadata}}{a named list for user-defined data.  The names \code{'path'}, \code{'header'} and \code{'date'} are already reserved; see \code{\link{metadata}()}.}

  \item{\code{step.wl}}{step between adjacent wavelengths in nm. This is assigned only when the wavelengths are regular; see \code{\link{is.regular}()}. }

  \item{\code{specname} }{only assigned when the \code{organization} is \code{'vector'}, in which case it is equal to the single character string name of the single spectrum.
  Note that the word \code{specname} is singular. Also see \code{\link{specnames}()}.}
}
And there are a few more attributes that exist only in special cases;
see the \bold{colorSpec User Guide}.
}


\seealso{
  \code{\link{wavelength}},
  \code{\link{quantity}},
  \code{\link{organization}},  
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

filter.bp   = colorSpec( trans, wave, quantity='transmittance', specnames='myfilter' )

organization( filter.bp )  # returns:  "vector"

# and now plot it
plot( filter.bp )
}

\keyword{colorSpec}
