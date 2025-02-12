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
is convex, closed, and bounded,
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
Let the equation of a hyperplane be given by:
\deqn{ <v,normal> = \beta }
where \eqn{normal} is orthogonal to the hyperplane,
and \eqn{\beta} is the plane constant, and \eqn{v} is a variable.
The purpose of the function \code{sectionOptimalColors()} 
is to compute the intersection of the hyperplane and \eqn{\partial}\bold{Z}.

In dimension 3 this hyperplane is a 2D plane,
and the intersection is generically a convex polygon.
If the plane is a supporting plane, then the intersection is a face
(a vertex, an edge, or a facet) of \bold{Z};
in this case the function computes the center of the face.

In dimension 2 this hyperplane is a line, and the intersection is generically 2 points.
If the line is supporing line, then the intersection is a vertex or an edge of \bold{Z};
in the case of an edge the function computes the center of the edge.

Of course, the intersection can also be empty.

The function is essentially a wrapper around 
\code{\link[zonohedra:section.zonohedron]{zonohedra::section.zonohedron()}}
and
\code{\link[zonohedra:section.zonogon]{zonohedra::section.zonogon()}}.
}


\usage{
\S3method{sectionOptimalColors}{colorSpec}( x, normal, beta )
}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type}
equal to \code{'responsivity.material'} and M spectra, where M=2 or 3.}

\item{normal}{a nonzero vector of dimension M, that is the normal to the family of parallel hyperplanes}

\item{beta}{a vector of numbers of positive length.
The number \code{beta[k]} defines the k'th plane \code{<v,normal> = beta[k]}. }.
}
    
    
    
\value{
The function returns a list with an item for each value in vector \code{beta}.
Each item in the output is a list with these items:

\item{beta}{the value of the plane constant \eqn{\beta}} 
\item{section}{an NxM matrix, where N is the number of points in the section,
and M is the dimension of \code{normal}.
The points of the section are the rows of the matrix.
If the intersection is empty, then N=0.}

In case of global error, the function returns \code{NULL}.
}


\details{
In the case that the dimension of \code{x} is 3, then \bold{Z} is a zonohedron.
For processing details see:
\code{\link[zonohedra:section.zonohedron]{zonohedra::section.zonohedron()}}.

In the case that the dimension of \code{x} is 2, then \bold{Z} is a zonogon.
For processing details see:
\code{\link[zonohedra:section.zonogon]{zonohedra::section.zonogon()}}.
}



\seealso{
vignette \href{../doc/optimals.pdf}{\bold{Plotting Chromaticity Loci of Optimal Colors}},
\code{\link{probeOptimalColors}()},
\code{\link[zonohedra:section.zonohedron]{zonohedra::section.zonohedron()}},
\code{\link[zonohedra:section.zonogon]{zonohedra::section.zonogon()}}
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
9(11):5, 1-23, (2009).\cr
\code{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
doi:10.1167/9.11.5.
}


\examples{
wave = seq(420,680,by=5)
Flea2.scanner = product( A.1nm, "material", Flea2.RGB, wavelength=wave )
seclist = sectionOptimalColors( Flea2.scanner, normal=c(0,1,0), beta=10 )
nrow( seclist[[1]]$section )
##  [1] 89 


seclist[[1]]$section[ 1:5, ]
##  the polygon has 89 vertices, and the first 5 are:
##           Red Green      Blue
##  233 109.2756    10 3.5391342
##  185 109.5729    10 2.5403628
##  136 109.8078    10 1.7020526
##  86  109.9942    10 1.0111585
##  35  110.1428    10 0.4513051
}


\keyword{colorSpec}
