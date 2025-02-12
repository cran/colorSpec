% File computeADL.Rd
\encoding{UTF-8}

\name{computeADL}
\title{compute ADL coordinates by ray tracing}
\alias{computeADL}
\alias{computeADL.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} equal to \code{responsivity.material}.
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
For details see:
\code{\link[zonohedra:raytrace2trans]{zonohedra::raytrace2trans()}}.

This surface is a subset of the \emph{object-color solid} or \emph{Rösch Farbkörper} for \code{x}.
Points on the boundary of of the solid are called \emph{optimal colors},
see \code{\link{probeOptimalColors}()}.
In most cases, a point on the 2-transition surface (a Schrödinger color) is an optimal color,
but this is not always true, see the vignette
\href{../doc/convexity.pdf}{\bold{Convexity and Transitions}}.

Now consider a color response \bold{R} not equal to \bold{G}.
There is a ray based at \bold{G} and passing through \bold{R}
that intersects the 2-transition surface
at a \emph{Schrödinger color} \bold{B}
with Logvinenko coordinates \eqn{(\delta,\omega)}.
If these 2 coordinates are combined with \eqn{\alpha}, where
\bold{R} = \bold{G} + \eqn{\alpha}(\bold{B}-\bold{G}), 
it yields the \emph{Logvinenko coordinates}
\eqn{(\alpha,\delta,\omega)} of \bold{R}.
These coordinates are also denoted by ADL; see \bold{References}.
A color is Schrödinger iff \eqn{\alpha{=}1}.

The coordinates of \bold{0} are \eqn{(\alpha,\delta,\omega)}=(1,0,NA).
The coordinates of \bold{W} are \eqn{(\alpha,\delta,\omega)}=(1,1,NA).
The coordinates of \bold{G} are \eqn{(\alpha,\delta,\omega)}=(0,NA,NA).
}


\usage{
\S3method{computeADL}{colorSpec}( x, response )
}

\arguments{
  \item{x}{a \bold{colorSpec} object with \code{type} equal to \code{responsivity.material} and 3 spectra }
  
  \item{response}{a numeric Nx3 matrix with responses in the rows, or a numeric vector
  that can be converted to such a matrix, by row.}
}
    
\value{
\code{computeADL()} returns a \code{data.frame} with a row for each response.
The columns in the data frame are:
\item{response}{the input response vector}
\item{lambda}{lambda.1 and lambda.2 at the 2 transitions, in nm.
lambda.1 < lambda.2  => bandpass,
and lambda.1 > lambda.2 => bandstop.}
\item{ADL}{the computed ADL coordinates of the response vector}  
\item{omega}{the reparameterized \eqn{\lambda} in the interval [0,1];  see \cite{Logvinenko}}

If an individual ray could not be traced,
the row contains \code{NA} in appropriate columns.
\cr
In case of global error, the function returns \code{NULL}.
}    


\details{
For each response, a ray is computed and the ray tracing is
done by \code{\link{probeOptimalColors}()}.

}


\section{WARNING}{
Since this function is really a simple wrapper around 
\code{\link[zonohedra:raytrace2trans]{zonohedra::raytrace2trans()}}.
Please see the warnings about star-shaped regions there.
}


\seealso{
\code{\link{type}()},
\code{\link{probeOptimalColors}()},
\code{\link[zonohedra:raytrace2trans]{zonohedra::raytrace2trans()}},
vignette \href{../doc/convexity.pdf}{\bold{Convexity and Transitions}}
}

\examples{
D50.eye = product( D50.5nm, 'varmat', xyz1931.1nm, wave=seq(360,830,by=5) )
computeADL( D50.eye, c(30,50,70) )

##    response.X response.Y response.Z lambda.1 lambda.2 ADL.alpha ADL.delta  ADL.lambda
##  1         30         50         70 427.2011 555.5262 0.7371480 0.5384106 473.3594244

##     omega
## 0.3008816
 
## since alpha < 1, XYZ=c(30,50,70) is *inside* the Schrodinger color surface of D50.eye
}


\references{
Logvinenko, A. D.
An object-color space. 
\bold{Journal of Vision}.
9(11):5, 1-23, (2009).\cr
\code{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
doi:10.1167/9.11.5.

Godau, Christoph and Brian Funt.
XYZ to ADL: Calculating Logvinenko's Object Color Coordinates.
Proceedings Eighteenth IS&T Color Imaging Conference.
San Antonio. Nov 2009.

West, G. and M. H. Brill.
Conditions under which Schrödinger object colors are optimal.
\bold{Journal of the Optical Society of America}.
73. pp. 1223-1225. 1983.
}


\keyword{colorSpec}
