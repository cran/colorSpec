% File emulate.Rd

\name{emulate}
\alias{emulate}
\alias{emulate.colorSpec}

\title{modify a colorSpec responder to emulate (approximate) another responder}

\description{
The two possible modifications are:
\itemize{
\item pre-multiplication by a transmitting filter
\item post-multiplication by a matrix
}
Both of these are optional.
If neither of these modifications is enabled, the original \code{x} is returned.\cr
}

\usage{
\S3method{emulate}{colorSpec}( x, y, filter=FALSE, matrix=TRUE )

}

\arguments{
  \item{x}{a \bold{colorSpec} responder with M spectra, to be modified. 
            The \code{type} must be \code{'responsivity.light'} or \code{'responsivity.material'}}
  
  \item{y}{a \bold{colorSpec} responder with N spectra, to be emulated by a modified \code{x}.
  It must have the same \code{type} and \code{wavelength} vector as \code{x}
  }
  
  \item{filter}{enable filter pre-multiplication.}
  
  \item{matrix}{enable matrix post-multiplication.
  If \code{matrix=TRUE} then the computed matrix \code{A} is MxN. }
}
    
\value{
a \bold{colorSpec} object close to \code{y}, as in \bold{Details}.
The \code{quantity} is the same as \code{y}.
The \code{specnames()} are the same as those of \code{y},
except that \code{".em"} is appended to each one.
The function attaches attribute \code{"emulate"},
whose value is a list containing \code{filter} and/or \code{A} as appropriate.
}


\details{
If \code{filter=FALSE} and \code{matrix=TRUE} then
the returned value is \code{multiply(x,A)},
where the matrix \code{A} is computed to minimize the difference with
\code{y}, in the least squares sense (Frobenius matrix norm).
The function \code{MASS::ginv()} is used here.

If \code{filter=TRUE} and \code{matrix=FALSE} then
the returned value is \code{product(filter,x)},
where the object \code{filter} is computed to minimize the difference with
\code{y}, in the least squares sense (Frobenius matrix norm).
This calculation is fairly straightforward,
but requires that the responsivity of \code{x} does not vanish at any wavelength.
It also requires that \code{M=N}.
The computed filter may be unrealistic, i.e. the transmittance may be > 1.
If this happens a \code{WARN} message is issued.

If \code{filter=TRUE} and \code{matrix=TRUE} then
the returned value is \code{product(filter,multiply(x,A))},
where \code{(filter,A)} are chosen with the above minimization criterion.
If \code{N=1} then we must have \code{M=1} as well;
the calculation is trivial and the emulation is exact.
If \code{N} \eqn{\ge} \code{2},
the calculation is iterative - solving alternatively for 
\code{filter} and \code{A} until convergence.
The function \code{MASS::ginv()} is used on each iteration.
This is a bilinear optimization.
If convergence fails, it is an error and the function returns \code{NULL}.
If convergence succeeds, there is 1 degree of freedom in the \code{(filter,A)} pair.
If one is scaled by a positive constant, the other can be scaled by the inverse,
and the returned object is the same.
The filter is scaled so the maximum transmittance is 1.

If \code{filter=FALSE} and \code{matrix=FALSE} then
the original \code{x} is returned, with a \code{WARN} message.
}



\seealso{
\code{\link{wavelength}},
\code{\link{type}},
\code{\link{quantity}},
\code{\link{multiply}},
\code{\link{product}},
\code{\link{specnames}}
}

\section{Examples}{
see the vignette \bold{Emulation of one Camera by another Camera}
}

\keyword{colorSpec}
