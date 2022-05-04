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

Currently the function only works if the number of spectra in 
\code{x} is 3 (e.g. RGB or XYZ).

Before \pkg{colorSpec} v 0.8-1 this function used a 2D root-finding
method that could only find
optimal colors whose spectra contain 0, 1, or 2 transitions.
But starting with v0.8-1, we have switched to zonohedral representation
of the object-color solid, which makes it possible to discover more than 2 transitions.
The inspiration for this change is the article by \cite{Centore}.
To \emph{inspect} these computed spectra, the argument \code{spectral} must be set
to \code{TRUE}.
}


\usage{
\S3method{probeOptimalColors}{colorSpec}( x, gray, direction, aux=FALSE, spectral=FALSE, tol=1.e-6 )

}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type}
equal to \code{'responsivity.material'} and 3 spectra}

\item{gray}{vector of numbers in the open interval (0,1) that define neutral grays on the
          line segment from black to white; this neutral gray point is the basepoint of a probe ray}

\item{direction}{a numeric Nx3 matrix with directions of the probe rays in the rows,
               or a numeric vector that can be converted to such a matrix, by row.}
          
   
\item{aux}{a logical that specifies whether to return extra performance and diagnostic data; see \bold{Details}}

\item{spectral}{if \code{TRUE}, the function returns a \bold{colorSpec} object
with the optimal spectra, see \bold{Value}.
}

\item{tol}{error tolerance for the intersection of probe and object-color boundary}
}
    
    
    
\value{
If argument \code{spectral=FALSE},
\code{probeOptimalColors()} returns a \code{data.frame} with a row for each traced ray.
There are \code{length(gray)} * \code{nrow(direction)} rays.
The columns in the output are:

\item{gray}{the graylevel defining the \eqn{basepoint} of the ray. 
\eqn{basepoint = gray*W} 
}

\item{direction}{the \eqn{direction} of the ray} 

\item{s}{computed scalar so that \eqn{basepoint + s*direction} is optimal}

\item{optimal}{
the optimal color on the boundary; \eqn{optimal = basepoint + s*direction} 
}
  
\item{lambda}{
lambda.1 and lambda.2 at the 2 transitions, in nm. 
lambda.1 < lambda.2 => bandpass,
and lambda.1 > lambda.2 => bandstop.
It will happen that the optimal spectrum has more than 2 transitions;
in this case both lambdas are set to \code{NA}.
}

\item{dol}{
\code{delta} and \code{omega} - the Logvinenko parameters \eqn{(\delta,\omega)} for optimal colors, plus \code{lambda} (\eqn{\lambda}) in nm. 
\eqn{\omega} is the reparameterization of \eqn{\lambda} ;  see \cite{Logvinenko}.
If there are more than 2 transistions, these are set to \code{NA}.
}
            
If aux is \code{TRUE}, these auxiliary columns related to performance and diagnostics are added:   

\item{timetrace}{time to trace the ray, in seconds}

\item{parallelograms}{
# of parallelograms in the (possibly compound) face. 1 means just a single parallelogram.
}

\item{tested}{
# of parallelograms actually tested for ray intersection. This only has meaning for compound faces.
}

\item{alpha}{the 2 coordinates of the intersection point inside the parallelogram}


If argument \code{spectral=TRUE},
\code{probeOptimalColors()} returns a \bold{colorSpec} object with quantity 
\code{'reflectance'}.
This object contains the optimal spectra, and can be used to inspect the spectra
with more than 2 transitions, which will happen.
The above-mentioned \code{data.frame} can then be obtained by applying
\code{\link{extradata}()} to the returned object.

If an individual ray could not be traced (which should be rare), 
the row contains \code{NA} in appropriate columns.\cr
In case of global error, the function returns \code{NULL}.
}    


\details{
Each gray level and each direction defines a ray.
So the total number of rays traced is \code{length(gray)} * \code{nrow(direction)}.
The 3 responsivities are regarded not as continuous functions, 
but as step functions.
This implies that the color solid is a zonohedron.
In the preprocessing phase the zonohedral representation is calculated.
The faces of the zonohedron are either parallelograms,
or \emph{compound faces} that can be partitioned into parallelograms.
The centers of all these parallelograms are computed, along with their normals
and plane constants.
\cr
This representation of the color solid is very strict regarding the
2-transition assumption.
During use, one can count on there being some spectra with more than two transitions.
Forcing the best 2-transition spectrum is a possible topic for the future.
}

\section{WARNING}{
The preprocessing calculation of the zonohedron dominates the total time.
And this time goes up rapidly with the number of wavelengths.
We recommend using a wavelength step of 5nm, as in the \bold{Examples}.
For best results, batch a lot of rays into a single function call
and then process the output.
\cr
Moreover, the preprocessing time is dominated by the partitioning
of the compound faces into parallelograms.
This is made worse by spectral responses with little overlap,
as in \code{\link{scanner.ACES}}.
In these cases, try a larger step size, and then reduce.
Optimizing these compound faces is a possible topic for the future.
}

\seealso{
\code{\link{type}},
vignette \href{../doc/optimals.pdf}{\bold{Plotting Chromaticity Loci of Optimal Colors}},
\code{\link{scanner.ACES}},
\code{\link{extradata}()}
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
9(11):5, 1-23, (2009).\cr
\code{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
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
probeOptimalColors( D50.eye, c(0.2,0.5,0.9), c(1,2,1, -1,-2,-1) )

##    gray direction.1 direction.2 direction.3         s  optimal.1  optimal.2
##  1  0.2           1           2           1 32.306207  52.533143  85.612065
##  2  0.2          -1          -2          -1  8.608798  11.618138   3.782055
##  3  0.5           1           2           1 20.993144  71.560483  94.485416
##  4  0.5          -1          -2          -1 20.993144  29.574196  10.512842
##  5  0.9           1           2           1  4.333700  95.354911 103.165832
##  6  0.9          -1          -2          -1 35.621938  55.399273  23.254556

##     optimal.3 lambda.1 lambda.2    dol.delta    dol.omega   dol.lambda
##  1  49.616046 451.8013 598.9589   0.63409966   0.48287469 536.97618091
##  2   8.701041 636.3031 429.4659   0.08458527   0.99624955 674.30015903
##  3  64.267740 441.9105 615.0822   0.78101041   0.49048222 538.73234859
##  4  22.281453 615.0822 441.9105   0.21898959   0.99048222 662.20606601
##  5  82.227974 422.9191 648.7404   0.95800430   0.49825407 540.49590064
##  6  42.272337 593.2415 455.2425   0.42035428   0.97962398 650.57382749


# create a 0-1 spectrum with 2 transitions
rectspec = rectangularMaterial( lambda=c(579.8697,613.7544), alpha=1, wave=wave )

# compute the corresponding color XYZ
XYZ = product( rectspec, D50.eye )
XYZ
##                             X        Y          Z
##  BP_[579.87,613.754] 33.42026 21.96895 0.02979764

# trace a ray from middle gray through XYZ
white.XYZ   = product( neutralMaterial(1,wave=wave), D50.eye )
direction   = XYZ - white.XYZ/2

res = probeOptimalColors( D50.eye, 0.5, direction, aux=FALSE )
res$s         
##  1.00004   the ray has gone past the original color to the boundary

res$optimal
##              X        Y          Z
##  [1,] 33.41958 21.96774 0.02808178

res$lambda    
##  NA NA     because there are more than 2 transitions in the true optimal

# since s=1.00004 > 1,
# XYZ is actually in the interior of the color solid, and not on the boundary.
# The boundary is a little-bit further along the ray,
# and the corresponding spectrum has more than 2 transitions.
}


\keyword{colorSpec}
