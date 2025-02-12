\name{inside}

\alias{insideOptimalColors}
\alias{insideOptimalColors.colorSpec}


\title{
test whether points are inside the surface of \emph{optimal colors}
}

\description{
Test points for being inside the surface of optimal colors,
which is the boundary of the \emph{object color solid}.
It is essentially a wrapper around
\code{\link[zonohedra:inside]{zonohedra::inside()}}.
}

\usage{
\S3method{insideOptimalColors}{colorSpec}( x, p )
}

\arguments{
\item{x}{
a \bold{colorSpec} object - with 1, 2, or 3 spectra
}

\item{p}{
an NxM numeric matrix, where M is the number of spectra.
The points to be tested are in the rows.
\code{p} can also be a numeric vector that can be converted
to such a matrix, by row.
}
}


\value{
\code{insideOptimalColors()} returns a \code{data.frame} with N rows and these columns:

\item{p}{the given point}

\item{distance}{the signed distance from the point to the boundary surface of optimal colors.
When \code{distance} < 0, the point is in the interior, and the distance is the
true Euclidean distance.
For boundary points, \code{distance} is 0.  When \code{distance} > 0, it may be larger
than the true Euclidean distance, so \code{distance} is really a pseudo-distance.}

\item{inside}{whether the point is inside the zonotope. \code{inside = distance<=0}}


If the row names of \code{p} are unique,
they are copied to the row names of the output.
\cr
In case of error, the function returns \code{NULL}.
}

\seealso{
\code{\link[zonohedra:inside]{zonohedra::inside()}}
}


\examples{
wave = seq(400,800,by=4)

D50.eye = product( D50.5nm, 'varmat', xyz1931.1nm, wave=wave )

white = colSums( 4 * as.matrix(D50.eye) )

insideOptimalColors( D50.eye,  c(30,50,70,  0,0,0, white, white+1 )  )

##          p.1       p.2       p.3 inside  distance
##  1  30.00000  50.00000  70.00000   TRUE -7.528005
##  2   0.00000   0.00000   0.00000   TRUE  0.000000
##  3 101.24793 105.04021  86.52988   TRUE  0.000000
##  4 102.24793 106.04021  87.52988  FALSE  1.259912
}



