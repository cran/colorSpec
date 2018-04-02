% File probeOptimalColors.Rd
\encoding{UTF-8}

\name{probeOptimalColors}
\title{compute optimal colors by ray tracing}
\alias{probeOptimalColors}
\alias{probeOptimalColors.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} 
equal to \code{'responsivity.material'}.
The set of all possible material reflectance functions (or transmittance functions) 
is convex, closed, and bounded (in any reasonable function space),
and this implies that the set of all possible output responses 
from \code{x} is also convex, closed, and bounded.
The latter set is called the \emph{object-color solid} or \emph{Rösch Farbkörper} for \code{x}.
A color on the boundary of the \emph{object-color solid} is called an \emph{optimal color}.
The special points \bold{W} (the response to the perfect reflecting diffuser)
and \bold{0} are on the boundary of this set.
The interior of the line segment of neutrals joining \bold{0} to \bold{W} is in the interior of the
\emph{object-color solid}.
It is natural to parameterize this segment from 0 to 1 (from \bold{0} to \bold{W}).

A ray \eqn{r} that is based at a point on the interior of the neutral line segment
must intersect the boundary of the \emph{object-color solid} in a unique optimal color.
The purpose of the function \code{probeOptimalColors()} is to compute that intersection point.

The function only works as stated if:
\enumerate{
\item the number of spectra in \code{x} is 3 (e.g. RGB or XYZ)
\item the chromaticity diagram of \code{x} is convex and well-ordered (injective with no reversals)
}
The 1st condition makes the situation simple enough to deal with.
The 2nd condition implies that a reflectance function is optimal iff it
takes the values 0 or 1, and has 0, 1, or 2 transitions;
see \cite{Logvinenko} or \cite{West} for the proof of this.
The proof in \cite{Schrödinger} is flawed because he ignored condition 2.
This 2-transition condition also simplifies the situation.
Fortunately, the CIE chromaticity diagrams (both 1931 and 1964) are convex
and injective.
For counter-examples see the \bold{References}.
If we define a \emph{Schrödinger object color} to be a color
defined by a reflectance function with  0, 1, or 2 transitions,
then it would be accurate to state that
\code{probeOptimalColors()} computes \emph{Schrödinger object colors}.
}


\usage{
\S3method{probeOptimalColors}{colorSpec}( x, gray, direction, tol=1.e-6, aux=FALSE )

}

\arguments{
  \item{x}{a \bold{colorSpec} object with \code{type}
  equal to \code{'responsivity.material'} and 3 spectra}
  
  \item{gray}{vector of numbers in the open interval (0,1) that define neutral grays on the
              line segment from black to white; this neutral gray point is the basepoint of a probe ray}
  
  \item{direction}{ a vector of 3 numbers, or a matrix with 3 columns, that define 1 or more directions
              for the probe rays.  A vector of 3 numbers is changed to a matrix with 1 row.}
              
  \item{tol}{error tolerance for the intersection of probe and object-color boundary}  
               
  \item{aux}{a logical that specifies whether to return extra performance data; see \bold{Details}}               
}
    
\value{
\code{probeOptimalColors()} returns a \code{data.frame} with a row for each traced ray.
There are \code{length(gray)} * \code{ncol(direction)} rays.
The columns in the output are:
  \item{gray}{the graylevel defining the \eqn{basepoint} of the ray. \eqn{basepoint = gray*W} }
  \item{direction}{the \eqn{direction} of the ray}  
  \item{s}{computed scalar so that \eqn{basepoint + s*direction} is optimal}
  \item{optimal}{the optimal color on the boundary; \eqn{optimal = basepoint + s*direction} }
  \item{lambda}{lambda.1 and lambda.2 at the 2 transitions, in nm. lambda.1 < lambda.2 => bandpass,
                and lambda.1 > lambda.2 => bandstop.}
  \item{dol}{\code{delta} and \code{omega} - the Logvinenko parameters \eqn{(\delta,\omega)} for optimal colors, plus \code{lambda} in nm. 
  \eqn{\omega} is the reparameterization of \eqn{\lambda} ;  see \bold{References}}
                
And if aux is \code{TRUE}, these auxiliary columns related to performance:   
  \item{time_grid}{time to find initial estimate point on boundary, in seconds}
  \item{iters}{number of interations of Newton's Method to find the ray intersection}
  \item{btracks}{total # of backtracks in "damped" Newton's method}
  \item{time_newt}{time spent in Newton iterations, in seconds}
  \item{error}{root-finding error, in coordinates of the optimal color. Always less than argument \code{tol}}
  
  If an individual ray could not be traced (see \bold{Known Issues}), 
  the row contains \code{NA} in appropriate columns.\cr
  In case of global error, the function returns \code{NULL}.
}    


\details{
Each gray level and each direction defines a ray.
So the total number of rays traced is \code{length(gray)} * \code{ncol(direction)}.
The intersection problem is reduced to a 2-dimension root finding problem
which is solved using Newton's Method.
The initial estimate is found by precomputing a fine quadrilateral mesh
over the optimal colors.
\cr
The responsivity functions may be negative,
but may not all simultaneously vanish at any wavelength.
In that case the mapping from the \eqn{\omega} and \eqn{\delta} sphere 
to the output response space is not injective.
}

\section{Known Issues}{
The optimal color boundary is not differentiable at \bold{0} and \bold{W}.
There may be numerical iteration failures near these 2 points.
}

\seealso{
  \code{\link{type}},
  vignette \bold{optimals}
}

\references{
Logvinenko, A. D.
An object-color space. 
\bold{Journal of Vision}.
9(11):5, 1-23, (2009).
\url{http://journalofvision.org/9/11/5/}.
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


\keyword{colorSpec}
