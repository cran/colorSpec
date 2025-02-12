% File sectionSchrodingerColors.Rd
\encoding{UTF-8}

\name{sectionSchrodingerColors}
\title{compute sections of a Schrodinger color surface by hyperplanes}
\alias{sectionSchrodingerColors}
\alias{sectionSchrodingerColors.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} equal to \code{'responsivity.material'}.
The set of all possible material reflectance functions (or transmittance functions)
that take the value 0 or 1, and with 2 or 0 transitions is called the \emph{2-transition spectrum space}.
When there are 2 transitions, there are 2 types of spectra: \emph{bandpass} and \emph{bandstop}.
When there are 0 transitions, the spectrum is either identically 0 or identically 1.
When \code{x} is applied to this space, the corresponding surface in response space
is called the \emph{2-transition surface}.  
The special points \bold{0} and \bold{W} (the response to the perfect reflecting diffuser)
are on this surface.
The surface is symmetrical about the neutral gray midpoint \bold{G}=\bold{W/2}.
Following \cite{West and Brill}, colors on the surface are called \emph{Schrödinger colors}.
Denote the surface by \eqn{S_2}.
For details see:
\code{\link[zonohedra:raytrace2trans]{zonohedra::raytrace2trans()}}.
This function only supports \code{x} with 3 channels.

Let the equation of a hyperplane be given by:
\deqn{ <v,normal> = \beta }
where \eqn{normal} is orthogonal to the hyperplane,
and \eqn{\beta} is the plane constant, and \eqn{v} is a variable.
The purpose of the function \code{sectionSchrodingerColors()} 
is to compute the intersection of the hyperplane and \eqn{S_2}.
The intersection is generically a single polygon, which is not necessarily convex.

In pathological cases, the intersection can be many polygons, but this function
only returns one of them, with a warning that there are more.
Of course, the intersection can also be empty.

The function is essentially a wrapper around 
\code{\link[zonohedra:section.zonohedron]{zonohedra::section2trans()}}.
}


\usage{
\S3method{sectionSchrodingerColors}{colorSpec}( x, normal, beta )
}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type}
equal to \code{'responsivity.material'} and 3 spectra.}

\item{normal}{a nonzero vector of dimension 3, that is the normal to the family of parallel hyperplanes}

\item{beta}{a vector of numbers of positive length.
The number \code{beta[k]} defines the k'th plane \code{<v,normal> = beta[k]}. }.
}
    
    
    
\value{
The function returns a list with an item for each value in vector \code{beta}.
Each item in the output is a list with these items:

\item{beta}{the value of the plane constant \eqn{\beta}} 
\item{section}{an Nx3 matrix, where N is the number of points in the section.
The points of the section are the rows of the matrix.
If the intersection is empty, then N=0.}

In case of global error, the function returns \code{NULL}.
}


\seealso{
vignette \href{../doc/optimals.pdf}{\bold{Plotting Chromaticity Loci of Optimal and Schrodinger Colors}},
\code{\link{sectionOptimalColors}()},
\code{\link[zonohedra:section.zonohedron]{zonohedra::section2trans()}}
}

\references{
West, G. and M. H. Brill.
Conditions under which Schrödinger object colors are optimal.
\bold{Journal of the Optical Society of America}.
73. pp. 1223-1225. 1983.
}


\examples{
wave = seq(420,680,by=5)
Flea2.scanner = product( A.1nm, "material", Flea2.RGB, wavelength=wave )
seclist = sectionSchrodingerColors( Flea2.scanner, normal=c(0,1,0), beta=10 )
nrow( seclist[[1]]$section )
##  [1] 106


seclist[[1]]$section[ 1:5, ]
##  the polygon has 106 vertices, and the first 5 are:
##              Red Green     Blue
##  [1,]  0.4054689    10 28.09329
##  [2,]  5.2833637    10 27.76564
##  [3,] 10.9896471    10 27.45328
##  [4,] 17.2285339    10 27.17679
##  [5,] 23.7280652    10 26.94359
}


\keyword{colorSpec}
