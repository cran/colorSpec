% File canonicalOptimalColors.Rd
\encoding{UTF-8}

\name{canonicalOptimalColors}
\title{compute the Canonical Optimal Colors}
\alias{canonicalOptimalColors}
\alias{canonicalOptimalColors.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} 
equal to \code{'responsivity.material'}.
The set of all possible material reflectance functions (or transmittance functions) 
is convex, closed, and bounded (in any reasonable function space),
and this implies that the set of all possible output responses 
from \code{x} is also convex, closed, and bounded.
The latter set is called the \emph{object-color solid} or \emph{Rösch Farbkörper} for \code{x}.
A color on the boundary of the \emph{object-color solid} is called an \emph{optimal color} for \code{x}.
The corresponding transmittance spectrum is called an  \emph{optimal spectrum} for \code{x}.
The special points \bold{W} (the response to the perfect reflecting diffuser)
and \bold{0} (the response to the perfect absorbing diffuser) are optimal.

Currently the function only works if the number of spectra in 
\code{x} is 3 (e.g. RGB or XYZ).
In this case the \emph{object-color solid} is a zonohedron whose boundary
is the union of parallelograms, which may be coplanar.
These parallelograms are indexed by distinct pairs of the wavelengths of \code{x};
if \code{x} has N wavelengths, then there are N*(N-1) parallelograms.
The center of each parallelogram is called a \emph{canonical optimal color}.
Interestingly, the special points \bold{W} and \bold{0} are not canonical.
}


\usage{
\S3method{canonicalOptimalColors}{colorSpec}( x, lambda, spectral=FALSE )

}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type}
equal to \code{'responsivity.material'} and 3 spectra}

\item{lambda}{a numeric Mx2 matrix whose rows contain distinct pairs of wavelengths of \code{x},
               or a numeric vector that can be converted to such a matrix, by row.
               If any entry in \code{lambda} is not a wavelength of \code{x},
               it is an error. }
          
\item{spectral}{if \code{TRUE}, the function returns a \bold{colorSpec} object
with the optimal spectra, see \bold{Value}.
}

}
    
    
    
\value{
If argument \code{spectral=FALSE},
\code{canonicalOptimalColors()} returns a \code{data.frame} with a row for each row in \code{lambda}.
The columns in the output are:

\item{lambda}{the given matrix argument \code{lambda}} 

\item{optimal}{
the computed optimal colors - an Mx3 matrix
}
  
\item{transitions}{the number of transitions in the optimal spectrum, this is a positive even number}

If \code{rownames(lambda)} is not \code{NULL}, they are copied to the row names of the output.

If argument \code{spectral=TRUE},
it returns a \bold{colorSpec} object with quantity 
\code{'transmittance'}.
This object contains the optimal spectra,
and the above-mentioned \code{data.frame} can then be obtained by applying
\code{\link{extradata}()} to the returned object.

In case of global error, the function returns \code{NULL}.
}


\details{
The 3 responsivities are regarded not as continuous functions, 
but as step functions.
This implies that the color solid is a zonohedron.
In the preprocessing phase the zonohedral representation is calculated.
The faces of the zonohedron are either parallelograms,
or \emph{compound faces} that can be partitioned into parallelograms.
The centers of all these parallelograms are the canonical optimal colors.
\cr
The optimal spectra take value 1/2 at the 2 given wavelengths, and 0 or 1 elsewhere.
If the 2 wavelengths are \eqn{\lambda_1} and \eqn{\lambda_2},
and \eqn{\lambda_1 < \lambda_2} then the spectrum is approximately a bandpass filter.
If the 2 wavelengths are swapped, then the spectrum is "flipped"
and is approximately a bandstop filter.
}



\seealso{
\code{\link{probeOptimalColors}()},
\code{\link{bandRepresentation}()},
\code{\link{scanner.ACES}},
\code{\link{extradata}()},
\code{\link{type}},
vignette \bold{convexity},
}

\references{
Centore, Paul.
\emph{A zonohedral approach to optimal colours}.
\bold{Color Research & Application}.
Vol. 38.
No. 2.
pp. 110-119.
April 2013.

Logvinenko, A. D.
An object-color space. 
\bold{Journal of Vision}.
9(11):5, 1-23, (2009).
\url{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
doi:10.1167/9.11.5.

Schrödinger, E. (1920). 
Theorie der Pigmente von grösster Leuchtkraft. 
\bold{Annalen der Physik}.
62, 603-622.

West, G. and M. H. Brill.
Conditions under which Schrödinger object colors are optimal.
\bold{Journal of the Optical Society of America}.
73. pp. 1223-1225. 1983.
}


\examples{
wave    = seq(400,700,by=5)
D50.eye = product( D50.5nm, 'material', xyz1931.1nm, wavelength=wave )
canonicalOptimalColors( D50.eye, c(500,600, 550,560, 580,585) )
##    lambda.1 lambda.2   optimal.x   optimal.y   optimal.z transitions
##  1      500      600 47.02281830 80.07281030  4.33181530           2
##  2      550      560  5.18490614 10.09045773  0.06121505           2
##  3      580      585 26.91247649 21.49031008  0.03457904           6
}


\keyword{colorSpec}
