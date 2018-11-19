% File plotOptimals.Rd
\encoding{UTF-8}

\name{plotOptimals}
\title{Plot Optimal Colors}
\alias{plotOptimals3D}
\alias{plotOptimals3D.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} 
equal to \code{'responsivity.material'} and 3 responsivity spectra.
The function \code{plotOptimals3D()} 
makes a wireframe plot of the \emph{object-color solid} for \code{x}.
The 3D drawing package \pkg{rgl} is required.
\cr
The set of all possible material reflectance functions (or transmittance functions) 
is convex, closed, and bounded (in any reasonable function space),
and this implies that the set of all possible output responses 
from \code{x} is also convex, closed, and bounded.
The latter set is called the \emph{object-color solid},
or \emph{Rösch Farbkörper}, for \code{x}.
A color on the boundary of the \emph{object-color solid} is called an \emph{optimal color}.
The special points \bold{W} (the response to the perfect reflecting diffuser)
and \bold{0} are on the boundary of this set.
The interior of the line segment of neutrals joining \bold{0} to \bold{W} is in the interior of the set.
For more discussion see \code{\link{probeOptimalColors}()}.

}


\usage{
\S3method{plotOptimals3D}{colorSpec}( x, size=c(33,33) )

}

\arguments{
  \item{x}{a \bold{colorSpec} object with \code{type} equal to
  \code{'responsivity.material'} and 3 spectra }
  
  \item{size}{an integer 2-vector with the number of meridians and parallels to plot}
               
}
    
\value{
\code{TRUE} or \code{FALSE}
}    


\details{
The boundary of the \emph{object-color solid} for \code{x} has
a natural parameterization by \eqn{\omega} and \eqn{\delta},
which are analogous to longitude and latitude for the sphere.
See \cite{Logvinenko} for more details.
These 2 parameters define reflectance spectra with 2 or fewer transitions
between 0 and 1.
By default, the function draws a wireframe with 33  meridians
and 33 parallels.
\cr
In addition it draws the box with opposite vertices at the "poles" \bold{0} and \bold{W}
and the diagonal segment of neutral grays that connects \bold{0} and \bold{W}.
It draws a small ball at the midpoint; 
the \emph{Rösch Farbkörper} is symmetrical about this midpoint.

}

\section{Note 1}{
As both \bold{References} point out, the plotted set of 2-transition colors is only optimal
under certain conditions on \code{x}, 
which are fortunately true for the human responsivity functions \code{xyz1931.1nm}.
But they are certainly not true for \emph{all} \code{x},
so the plotted surface seen here might really be sub-optimal.
\cr
For general \code{x} it may require reflectance
functions with 3,4,5,... transitions to define the optimals.

}


\section{Note 2}{
If all responsivity functions of \code{x} are non-negative,
the \emph{object-color solid} of \code{x} is inside the box.
If the responsivity functions of \code{x} have negative lobes,
the \emph{object-color solid} of \code{x} extends outside the box.
Indeed, the box may actually be \emph{inside} the optimals.
\cr
The responsivity functions cannot all simultaneously vanish at any wavelength.
In that case the mapping from the \eqn{\omega} and \eqn{\delta} sphere 
to the output response space is not injective.
}


\section{Future Work}{
For an object \code{x} with 2 spectra,
it would not take much work to write a new function
\code{plotOptimals2D()} that plots the 1-transition colors for \code{x}.
These are the short-pass and long-pass colors,
also known as \emph{edge colors} or \emph{Kantenfarben}.\cr
With a lot more work it would be possible to
plot the true optimals for \code{x}, with any number of transitions between 0 and 1.

}



\seealso{
\code{\link{type}},
\code{\link{probeOptimalColors}}  
vignette \bold{optimals}
}



\references{
Logvinenko, A. D.
An object-color space. 
Journal of Vision.
9(11):5, 1-23, (2009).
\url{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
doi:10.1167/9.11.5.

West, G. and M. H. Brill.
Conditions under which Schrödinger object colors are optimal.
Journal of the Optical Society of America.
73. pp. 1223-1225. 1983.
}


\examples{\donttest{
# requires package rgl
library( rgl )
human = product( D50.5nm, 'slot', xyz1931.5nm, wave=400:770 )
plotOptimals3D( human )
scanner = product( D50.5nm, 'slot', BT.709.RGB, wave=400:770 )
plotOptimals3D( scanner )

}}

\keyword{colorSpec}
