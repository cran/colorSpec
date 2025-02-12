% File plotOptimals.Rd
\encoding{UTF-8}

\name{plotOptimals}
\title{Plot Optimal Colors}
\alias{plotOptimals3D}
\alias{plotOptimals3D.colorSpec}
\alias{plotOptimals2D}
\alias{plotOptimals2D.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} 
equal to \code{'responsivity.material'} and 3 responsivity spectra.
The function \code{plotOptimals3D()} 
makes a plot of the \emph{object-color solid} for \code{x}.
This solid is a \emph{zonohedron} in 3D.
The 3D drawing package \pkg{rgl} is required.
\cr
Consider a \bold{colorSpec} object \code{x} with \code{type} 
equal to \code{'responsivity.material'} and 2 responsivity spectra.
The function \code{plotOptimals2D()} 
makes a plot of the \emph{object-color solid} for \code{x}.
This solid is a \emph{zonogon} in 2D.
The 3D drawing package \pkg{rgl} is \emph{not} required.
\cr
The set of all possible material reflectance functions (or transmittance functions) 
is convex, closed, and bounded (in any reasonable function space),
and this implies that the set of all possible output responses 
from \code{x} is also convex, closed, and bounded.
The latter set is called the \emph{object-color solid},
or \emph{Rösch Farbkörper}, for \code{x}.
A color on the boundary of the \emph{object-color solid} is called an \emph{optimal color}.

These functions are essentially wrappers around
\code{\link[zonohedra:plot.zonogon]{zonohedra::plot.zonogon()}} and
\code{\link[zonohedra:plot.zonohedron]{zonohedra::plot.zonohedron()}}. 
}


\usage{
\S3method{plotOptimals3D}{colorSpec}( x, ... )

\S3method{plotOptimals2D}{colorSpec}( x, ... )

}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type} equal to
\code{'responsivity.material'} and 2 or 3 spectra, as appropriate.
}

\item{...}{more arguments passed to
\code{\link[zonohedra:plot.zonogon]{zonohedra::plot.zonogon()}} or
\code{\link[zonohedra:plot.zonohedron]{zonohedra::plot.zonohedron()}} 
as appropriate.
For \code{plotOptimals3D()}, examples are \code{type} and \code{both}.
For \code{plotOptimals2D()}, examples are \code{orientation} and \code{elabels}.
}

}

    
\value{
The functions return \code{TRUE} or \code{FALSE}.
}    


\section{Note}{
If all responsivity functions of \code{x} are non-negative,
the \emph{object-color solid} of \code{x} is inside the box.
If the responsivity functions of \code{x} have negative lobes,
the \emph{object-color solid} of \code{x} extends outside the box.
Indeed, the box may actually be \emph{inside} the optimals.
}


\seealso{
\code{\link{type}()},
\code{\link{probeOptimalColors}()},
\code{\link{sectionOptimalColors}()},  
\code{\link[zonohedra:plot.zonogon]{zonohedra::plot.zonogon()}},
\code{\link[zonohedra:plot.zonohedron]{zonohedra::plot.zonohedron()}}, 
vignette \href{../doc/optimals.pdf}{\bold{Plotting Chromaticity Loci of Optimal and Schrodinger Colors}}
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
Journal of Vision.
9(11):5, 1-23, (2009).\cr
\code{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
doi:10.1167/9.11.5.

West, G. and M. H. Brill.
Conditions under which Schrödinger object colors are optimal.
\bold{Journal of the Optical Society of America}.
73. pp. 1223-1225. 1983.
}


\examples{\donttest{
human = product( D50.5nm, 'slot', xyz1931.5nm, wave=seq(400,770,by=5) )
plotOptimals3D( human )

plotOptimals2D( subset(human,2:3) )     # y and z only

scanner = product( D50.5nm, 'slot', BT.709.RGB, wave=seq(400,770,by=5) )
plotOptimals3D( scanner )
}}

\keyword{colorSpec}
