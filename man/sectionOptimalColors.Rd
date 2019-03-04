% File sectionOptimalColors.Rd
\encoding{UTF-8}

\name{sectionOptimalColors}
\title{compute sections of an optimal color surface by hyperplanes}
\alias{sectionOptimalColors}
\alias{sectionOptimalColors.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} 
equal to \code{'responsivity.material'}.
The set of all possible material reflectance functions (or transmittance functions) 
is convex, closed, and bounded (in fact they form a cube),
and this implies that the set of all possible output responses 
from \code{x} is also convex, closed, and bounded.
The latter set is called the \emph{object-color solid} or \emph{Rösch Farbkörper} for \code{x}.
If the dimension of the response of \code{x} is 2, 
this solid is a convex polygon
that is centrally symmetric - a \emph{zonogon}.
If the dimension of the response of \code{x} is 3 (e.g. RGB or XYZ), 
this solid is a special type of centrally symmetric convex polyhedron
called a \emph{zonohedron}, see \cite{Centore}.
This function only supports dimensions 2 and 3.
Denote this object-color solid by \bold{Z}.

A color on the boundary of \bold{Z} is called an \emph{optimal color}.
Consider the intersection of a hyperplane with the boundary of \bold{Z}.
Let the equation of the hyperplane be given by:
\deqn{ <v,normal> = \beta }
where \eqn{normal} is orthogonal to the hyperplane,
and \eqn{\beta} is the plane constant, and \eqn{v} is a variable vector.
The purpose of the function \code{sectionOptimalColors()} 
is to compute the intersection set.

In dimension 2 this hyperplane is a line, and the intersection is
generically 2 points, and 1 point if the line only intersects the boundary
(we ignore the special case when the intersection is an edge of the polygon).

In dimension 3 this hyperplane is a 2D plane, and the intersection is
generically a polygon, and 1 point if the line only intersects the boundary
(we ignore the special case when the intersection is a face of the zonohedron).

Of course, the intersection can also be empty.
}


\usage{
\S3method{sectionOptimalColors}{colorSpec}( x, normal, beta )
}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type}
equal to \code{'responsivity.material'} and M spectra, where M=2 or 3.}

\item{normal}{a nonzero vector of dimension M, that is the normal to a hyperplane}

\item{beta}{a vector of numbers of positive length.
The number \code{beta[k]} defines the plane \code{<v,normal> = beta[k]}. }.
}
    
    
    
\value{
The function returns a list with an item for each value in vector \code{beta}.
Each item in the output is a list with these items:

\item{beta}{the value of the plane constant \eqn{\beta}} 
\item{section}{an NxM matrix, where N is the number of points in the section,
and M is the dimension of \code{normal}.
If the intersection is empty, then N=0.}

In case of global error, the function returns \code{NULL}.
}


\details{
Consider first the case that the dimension of \code{x} is 3,
so that \bold{Z} is a zonohedron.
In the preprocessing phase the zonohedral representation is calculated.
The faces of \bold{Z} are either parallelograms,
or \emph{compound faces} that are partitioned into parallelograms.
The centers of all these parallelograms are computed,
along with their extent in direction \eqn{normal}.
For a given plane \eqn{<v,normal>=\beta},
the parallelograms that intersect the plane are extracted.
The boundary of each parallelogram intersects the plane in 2 points (in general)
and one of those points is computed.
The set of all these points is then sorted into proper order around the boundary.
\cr
In the case that the dimension of \code{x} is 2,
so that \bold{Z} is a zonogon,
the parallelograms are replaced by line segments (edges),
and the processing is much easier.
}

\section{WARNING}{
The preprocessing calculation of the zonohedron dominates the total time.
And this time goes up rapidly with the number of wavelengths.
We recommend using a wavelength step of 5nm, as in the \bold{Examples}.
For best results, batch a lot of \code{beta}s into a single function call
and then process the output.
\cr
Moreover, the preprocessing time is dominated by the partitioning
of the compound faces into parallelograms.
This is made worse by an \code{x} whose spectral responses have little overlap,
as in \code{\link{scanner.ACES}}.
In these cases, try a larger step size, and then reduce.
Optimizing these compound faces is a possible topic for the future.
}

\seealso{
vignette \bold{optimals},
\code{\link{probeOptimalColors}()}
}

\references{
Centore, Paul.
\emph{A Zonohedral Approach to Optimal Colours}.
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
}


\examples{
wave = seq(420,680,by=5)
Flea2.scanner = product( A.1nm, "material", Flea2.RGB, wavelength=wave )
seclist = sectionOptimalColors( Flea2.scanner, normal=c(0,1,0), beta=10 )
length( seclist[[1]]$section )
seclist[[1]]$section[ 1:5, ]
## [1] 207   # the polygon has 207 vertices, and the first 5 are:
##            Red Green      Blue
##  [1,] 109.2756    10 3.5391342
##  [2,] 109.5729    10 2.5403628
##  [3,] 109.8078    10 1.7020526
##  [4,] 109.9942    10 1.0111585
##  [5,] 110.1428    10 0.4513051
}


\keyword{colorSpec}
