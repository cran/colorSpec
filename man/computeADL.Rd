% File computeADL.Rd
\encoding{UTF-8}

\name{computeADL}
\title{compute ADL coordinates by ray tracing}
\alias{computeADL}
\alias{computeADL.colorSpec}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} equal to \code{responsivity.material}.
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
The solid is symmetrical about the neutral gray midpoint \bold{G}=\bold{W/2}.

Now suppose that \code{x} has 3 spectra (3 responses)
and consider a color response \bold{R} not equal to \bold{G}.
There is a ray based at \bold{G} and passing through \bold{R}
that intersects the boundary of the 
\emph{object-color solid} at an \emph{optimal color} \bold{B} on the boundary
with Logvinenko coordinates \eqn{(\delta,\omega)}.
If these 2 coordinates are combined with \eqn{\alpha}, where
\bold{R} = \eqn{(1-\alpha)}\bold{G} + \eqn{\alpha}\bold{B}, 
it yields the \emph{Logvinenko coordinates}
\eqn{(\alpha,\delta,\omega)} of \bold{R}.
These coordinates are also denoted by ADL; see \bold{References}.
A response is in the \emph{object-color solid} iff \eqn{\alpha \le 1}.
A response is  \emph{optimal} iff \eqn{\alpha=1}.

The coordinates of \bold{0} are \eqn{(\alpha,\delta,\omega)}=(1,0,0).
The coordinates of \bold{W} are \eqn{(\alpha,\delta,\omega)}=(1,1,0).
The coordinates of \bold{G} are undefined.
}


\usage{
\S3method{computeADL}{colorSpec}( x, response )

}

\arguments{
  \item{x}{a \bold{colorSpec} object with \code{type} equal to \code{responsivity.material} and 3 spectra }
  
  \item{response}{ a vector of 3 numbers, or a matrix with 3 columns, that define 1 or more responses.
A 3-vector is re-shaped to a 1x3 matrix.}              
}
    
\value{
\code{computeADL()} returns a \code{data.frame} with a row for each response.
The columns in the data frame are:
  \item{response}{the input response vector}
  \item{ADL}{the computed ADL coordinates of the response vector}  
  \item{omega}{the reparameterized \eqn{\lambda} in the interval [0,1];  see \bold{References}}
  \item{lambda}{lambda.1 and lambda.2 at the 2 transitions, in nm.
  lambda.1 < lambda.2  => bandpass,
  and lambda.1 > lambda.2 => bandstop.}

  If an individual ray could not be traced, the row contains \code{NA} in appropriate columns.\cr
  In case of global error, the function returns \code{NULL}.
}    


\details{
For each response, a ray is computed and the ray tracing is
done by \code{\link{probeOptimalColors}()}.

}


\section{Known Issues}{
The optimal color boundary is not differentiable at \bold{0} and \bold{W}.
There may be numerical iteration failures if the response
is near the neutral axis.

}


\seealso{
\code{\link{type}},
\code{\link{probeOptimalColors}},
vignette \bold{optimals}
}

\examples{
D50.eye = product( D50.5nm, 'varmat', xyz1931.1nm, wave='auto' )
computeADL( D50.eye, c(30,50,70) )
# response.x response.y response.z   ADL.alpha   ADL.delta  ADL.lambda     omega lambda.1 lambda.2
#         30         50         70   0.7364348   0.5384243 473.3909184 0.3008561 427.1431 555.5176
#since alpha < 1, this response is *inside* the object-color solid
}


\references{
Logvinenko, A. D.
An object-color space. 
\bold{Journal of Vision}.
9(11):5, 1-23, (2009).
\url{http://journalofvision.org/9/11/5/}.
doi:10.1167/9.11.5.

Godau, Christoph and Brian Funt.
XYZ to ADL: Calculating Logvinenko's Object Color Coordinates.
Proceedings Eighteenth IS&T Color Imaging Conference.
San Antonio. Nov 2009.

}


\keyword{colorSpec}
