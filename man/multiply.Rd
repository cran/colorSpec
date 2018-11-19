% File multiply.Rd

\name{multiply}
\title{multiply a colorSpec object by scalar, vector, or matrix}
\alias{multiply}
\alias{normalize}
\alias{multiply.colorSpec}
\alias{normalize.colorSpec}

\description{
multiply spectra by coefficients and return modified object
}
\usage{
\S3method{multiply}{colorSpec}( x, s )

\S3method{normalize}{colorSpec}( x, norm='L1' )
}

\arguments{
  \item{x}{a \bold{colorSpec} object with M spectra}
  
  \item{s}{a scalar, an M-vector, or an MxP matrix. 
           In the case of a matrix, assigning \code{colnames(s)} is recommended;
           see \bold{Details}.}  
  
  \item{norm}{one of \code{'L1'}, \code{'L2'}, or \code{'Linf'}, specifying one of the standard vector norms \eqn{L^1, L^2, or L^{inf}}. \cr
  \code{norm} can also be a numeric wavelength (e.g. 560 nm), and then the spectrum is scaled to have value 1 at this wavelength.  Of course, this is not a true vector norm. }
}
    
\value{
  \code{multiply} returns a \bold{colorSpec} object with the matrix of spectra
  of \code{x} multiplied by \code{s}.

  \code{normalize} returns a \bold{colorSpec} object with each spectrum
  of \code{x} scaled to have given norm equal to 1.
  
  In both functions, the \code{\link{quantity}} and \code{\link{wavelength}} are preserved.  
}


\details{
  For \code{multiply()}:\cr
  If \code{s} is an  MxP  matrix, say \bold{S}, and one thinks of the spectra as
  organized in an  NxM  matrix \bold{X}, 
  then the new spectra are defined by the matrix \bold{XS}, which is  NxP.
  If the P column names of \code{s} are set, then they are copied to the 
  \code{specnames} of the output.
  Otherwise, default spectrum names are assigned as in   
  \code{\link{colorSpec}()}, with a warning. \cr
  If \code{s} is an M-vector, then \bold{S}=\code{diag(s)} is computed and used in the previous sentence.
  This has the effect of multiplying spectrum i by \code{s}[i].\cr
  If \code{s} is a scalar then every spectrum is multiplied by \code{s}.\cr
  The multiplication may produce negative entries, but no check is made for this.\cr  
  WARNING: An M-vector and an Mx1 matrix may yield quite different results.

  For \code{normalize()}:\cr
  \code{normalize()} calls \code{multiply()} with \code{s} = an M-vector.
  If the norm of a spectrum is 0, then it is left unchanged.   
}

\note{
If \code{x} is organized as a matrix, and \code{s} is a scalar,
the one can use the simpler and equivalent \code{s*x}.
}


\seealso{
  \code{\link{product}()},
  \code{\link{quantity}()},
  \code{\link{wavelength}()},
  \code{\link{specnames}()},  
  \code{\link{colorSpec}()}
}

\keyword{colorSpec}
