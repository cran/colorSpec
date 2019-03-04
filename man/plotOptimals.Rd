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
For more discussion see \code{\link{sectionOptimalColors}()}.

}


\usage{
\S3method{plotOptimals3D}{colorSpec}( x, size=50, type='w', both=TRUE )

\S3method{plotOptimals2D}{colorSpec}( x )

}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type} equal to
\code{'responsivity.material'} and 2 or 3 spectra, as appropriate.
}

\item{size}{an integer giving the number of wavelengths at which to resample \code{x}.
To skip resampling, set \code{size=NA}.
}

\item{type}{\code{type='w'} for a wireframe plot of the parallelogram faces.
\code{type='p'} for a point plot with points at the centers of the parallelograms.
}

\item{both}{the color solid is symmetric about its center, so only half of it must
be computed.
If \code{both=TRUE} it plots one half in black and the other half in red.
If \code{both=FALSE} it only plots one half in black.
}
}

    
\value{
The functions return \code{TRUE} or \code{FALSE}.
}    


\section{Details for 3D}{
If \code{n} is the number of wavelengths,
the number of parallelogram faces of the zonohedron is \code{n*(n-1)}.
The time to compute these faces increase with \code{n} even faster,
so that is why the default \code{size=50} is a fairly small number.
It was chosen to be a reasonable compromise between detail and performance.
\cr
In addition to the wireframe or points,
it draws the box with opposite vertices at the "poles" \bold{0} and \bold{W}
and the diagonal segment of neutral grays that connects \bold{0} and \bold{W}.
}

\section{Details for 2D}{
If \code{n} is the number of wavelengths,
the number of edges in the zonogon is \code{2*n}.
Computing these edges is fast and visualization is easy,
so there are no plotting options at this time.
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
vignette \bold{optimals}
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
9(11):5, 1-23, (2009).
\url{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
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
