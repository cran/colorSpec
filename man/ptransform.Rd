% File ptransform.Rd

\name{ptransform}
\alias{ptransform}
\alias{ptransform.colorSpec}

\title{make a linear transformation to a colorSpec responder}

\description{
apply a linear transformation to a \bold{colorSpec} responder with M spectra,
so that multiples of M given \emph{primary} vectors
are transformed to the standard basis of \eqn{R^M}.
And a given \emph{white} vector is transformed to the M-vector of all 1s.\cr
The returned object is always
\code{multiply(x,A)} where  \code{A} is an internally calculated MxM matrix.
The name \code{ptransform} is short for \emph{projective transformation}.\cr
In case of ERROR, a message is logged and NULL returned.
}

\usage{
\S3method{ptransform}{colorSpec}( x, primary, white, digits=Inf )

}

\arguments{
  \item{x}{a \bold{colorSpec} responder with M spectra.  
    \code{type(x)} must be \code{'responsivity.light'} or \code{'responsivity.material'}.
  }
  
  \item{primary}{an MxM matrix whose rows define the M primary vectors
    in the response space of \code{x}.
    It is OK for each row to have a single value that is NA.
    In this case the NA value is changed so that the sum of the row is 1.
    This is done because typically the rows represent chromaticities
    in the response space of \code{x}.
    After this change, the rows of \code{primary} must form a basis of \eqn{R^M}.\cr
    \code{rownames(primary)} must be defined; when M=3 they are typically
    \code{c('R','G','B')}.
  }
  
  \item{white}{an M-vector in the response space of \code{x},
  that is typically the ideal white point of a color display.
  When \code{white} is expressed in the basis defined by \code{primary},
  the coordinates must all be non-zero.\cr
  \code{white} can also be a \bold{colorSpec} object with a single spectrum
  suitable as stimulus for \code{x};
  in this case the vector \code{white} is set to
  \code{ product( white, x, wavelength='auto' )}.
  }
  
  \item{digits}{if a positive integer,
  then \code{white} is rounded to this number of decimal digits,
  but in a non-standard way, see \bold{Details}.  
  This is typically done so the internally calculated MxM matrix \code{A} agrees
  with that from a color standard, see \bold{Examples}.
  }
}
    
\value{
   a \bold{colorSpec} object equal to \code{multiply(x,A)} 
   where \code{A} is an internally calculated MxM matrix.  
   The \code{\link{quantity}} and \code{\link{wavelength}} are preserved.
   The \code{specnames} of the returned object are set to \code{tolower( rownames(primary) )}.\cr
   The user may want to change the \code{\link{quantity}} of the returned object;
   see \bold{Examples}.
}


\details{
The formal mathematical requirements for \code{primary} and \code{white} are:
\itemize{ 
  \item The rows of \code{primary} must form a basis of \eqn{R^M}.
  Equivalently, if \eqn{P} denotes the matrix \code{primary},
  then \eqn{P} is invertible.
  \item The coordinates of \code{white} in this basis are all non-zero.
  Equivalently, if \eqn{x} is the solution of \eqn{ x P = white},
  then every component of x is non-zero.
  }
Assuming both of these are true, then there is a unique matrix \eqn{A} so that
\itemize{ 
  \item \eqn{A} transforms a multiple of the \eqn{i}'th row of \eqn{P} to the \eqn{i}'th
  standard basis vector of \eqn{R^M}.
  \item \eqn{A} transforms \code{white} to the M-vector of all 1s.
}
This statement is essentially the fundamental theorem of (analytic) projective geometry;
see \emph{Bumcrot} page 87, and \emph{Semple} page 398.
The rows of \eqn{P} plus \eqn{white} define a \emph{projective frame};
the former are the \emph{fundamental points} and the latter is the \emph{unit point}.

If \code{digits} is a positive integer,
the chromaticity of \code{white} is computed by dividing \code{white} by \code{sum(white)}.
The latter must be non-zero, or else it is an ERROR.
This chromaticity is rounded to \code{digits} decimal digits, while preserving the sum of 1.
This \emph{rounded chromaticity} is non-zero, and defines a line through 0.
The vector \code{white} is projected onto this line to get the new and rounded \code{white},
with the rounded chromaticity.
See \bold{Examples}.
}



\references{
Bumcrot, Robert J.
\bold{Modern Projective Geometry}.
Holt, Rinehart, and Winston. 1969.

IEC 61966-2-1:1999.
Multimedia systems and equipment - Colour measurement and management.
Part 2-1: Colour management - Default RGB colour space - sRGB.
\url{https://webstore.iec.ch/publication/6169}

Semple, J. G. and G. T. Kneebone.
\bold{Algebraic Projective Geometry}.
Oxford. 1952.
}




\seealso{
\code{\link{quantity}},
\code{\link{wavelength}},
\code{\link{colorSpec}},
\code{\link{multiply}},
\code{\link{product}}
}

\examples{
############ Example for sRGB   ###########

# assign the standard sRGB primaries
P = matrix( c(0.64,0.33,NA,  0.3,0.6,NA, 0.15,0.06,NA ), 3, 3, byrow=TRUE )
rownames(P) = c('R','G','B')
P
#   [,1] [,2] [,3]
# R 0.64 0.33   NA
# G 0.30 0.60   NA
# B 0.15 0.06   NA

white = product( D65.1nm, xyz1931.1nm, wave='auto' )
white
#           X        Y        Z
# D65 100.437 105.6708 115.0574

white/sum(white)
#             X         Y         Z
# D65 0.3127269 0.3290232 0.3582499    

# But the sRGB standard D65 is xy=(0.3127,0.3290)
# so when the next line is executed,
# the calculated 3x3 matrix will *NOT* agree with the sRGB standard
y = ptransform( xyz1931.1nm, P, white, digits=Inf )

product( D65.1nm, y, wave='auto' )
#     R G B
# D65 1 1 1      # this is exactly what we want, but the internal 3x3 matrix is a little off

# now repeat, but this time round the white chromaticity to
# xy=(0.3127,0.3290) in order to get the matrix right
y = ptransform( xyz1931.1nm, P, white, digits=4 )

rgb = product( D65.1nm, y, wave='auto' )
rgb
#            R        G        B
# D65 1.000238 1.000053 0.999835   # off in the 4'th digit  (WARN: this is linear RGB)

255 * rgb
#            R        G        B
# D65 255.0607 255.0134 254.9579   # good enough for 8-bit RGB

65535 * rgb
#            R        G        B
# D65 65550.59 65538.44 65524.18   # NOT good enough for 16-bit RGB  

# So for 16-bit RGB, one can get the white RGB right, or the 3x3 matrix right, but not both !


############ Example for ProPhoto RGB   ###########

# assign the standard ProPhoto RGB primaries
P = matrix( c(0.7347,0.2653,NA,  0.1596,0.8404,NA, 0.0366,0.0001,NA ), 3, 3, byrow=TRUE )
rownames(P) = c('R','G','B')
P
#     [,1]   [,2] [,3]
# R 0.7347 0.2653   NA
# G 0.1596 0.8404   NA
# B 0.0366 0.0001   NA

white = product( D50.5nm, xyz1931.5nm, wave='auto' )
white
#            X       Y        Z
# D50 101.2815 105.042 86.67237

white / sum(white)
#             X         Y         Z
# D50 0.3456755 0.3585101 0.2958144  

# but the ProPhoto RGB standard is xy=(0.3457,0.3585);  proceed anyway
y = ptransform( xyz1931.5nm, P, white, digits=Inf )

product( D50.5nm, y, wave='auto' )
#     R G B
# D50 1 1 1     #  this is exactly what we want, but the internal 3x3 matrix is a little off

# the following line is an equivalent way to compute y.
# pass D50.5nm directly as the 'white' argument
y = ptransform( xyz1931.5nm, P, D50.5nm )
}


\keyword{colorSpec}
