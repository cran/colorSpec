\name{insideSchrodingerColors}

\alias{insideSchrodingerColors}
\alias{insideSchrodingerColors.colorSpec}


\title{
test whether points are inside the surface of Schrödinger colors
}

\description{
Consider a \bold{colorSpec} object \code{x} with \code{type} equal to \code{'responsivity.material'}.
The set of all possible material reflectance functions (or transmittance functions)
that take the value 0 or 1, and with 2 or 0 transitions is called the \emph{2-transition spectrum space}.
When there are 2 transitions, there are 2 types of spectra: \emph{bandpass} and \emph{bandstop}.
When there are 0 transitions, the spectrum is either identically 0 or identically 1.
When \code{x} is applied to this space, the corresponding surface in response space
is called the \emph{2-transition surface}.
The surface is closed.
The special points \bold{0} and \bold{W} (the response to the perfect reflecting diffuser)
are on this surface.
The surface is symmetrical about the neutral gray midpoint \bold{G}=\bold{W/2}.
Following \cite{West and Brill}, colors on the surface are called \emph{Schrödinger colors}.
Denote the surface by \eqn{S_2}.

This function computes whether given points are inside \eqn{S_2}.
It only supports \code{x} with 3 channels.
It is essentially a wrapper around
\code{\link[zonohedra:inside2trans]{zonohedra::inside2trans()}}.
}

\usage{
\S3method{insideSchrodingerColors}{colorSpec}( x, p )
}

\arguments{
\item{x}{
a \bold{zonohedron} object
}

\item{p}{
an Nx3 numeric matrix. The points to be tested are in the rows.
\code{p} can also be a numeric vector that can be converted
to such a matrix, by row.
}
}


\details{
If the surface has no self-intersections,
the the definition of whether a point \code{p} is
"inside" is fairly straightforward:
it is where the linking number of
\code{p} and the surface is non-zero.
In fact, if it is non-zero then it must be +1 or -1.
The \emph{linking number} is analogous the \emph{winding number} in 2D.

Unfortunately, there is currently no test for whether the
surface \emph{has} self-intersections,
For a bad surface with self-intersections,
the linking number might be any integer.
Since there is no such test,
we simply use the same non-zero linking number rule always.

The computed \code{linkingnumber} is returned so that the user
can apply the non-zero rule, or the even-odd rule,
as appropriate for their situation.
These 2 rules are analogous to the two winding number rules
used for polygons in computer graphics,
see \bold{Point in polygon}.

The case where a point is \emph{on} the surface
(i.e. the \code{distance} to the surface is 0) is problematic.
The \code{linkingnumber} is then undefined,
and we currently set \code{inside} to be undefined as well.
Thus \code{inside} should be interpreted as \emph{strictly inside}.
However, in some situations, the user may want to consider
\code{inside} to be \code{TRUE} in this problematic case.
Or the user may want to consider points that are within
a very small epsilon of the surface,
where roundoff might have occurred, to have \code{inside=FALSE}
or \code{inside=NA}.
So the both the computed \code{linkingnumber} and \code{distance} are
returned so the user can use them
to make their own definition of what "inside" means.
}


\value{
\code{insideSchrodingerColors()} returns a \code{data.frame}
with N rows and these columns:

\item{p}{the given point}

\item{distance}{the distance from the point to the surface.
This is the true Euclidean distance,
and not a signed distance.
If the point is on the surface, the \code{distance} should be 0
up to numerical precision.
}

\item{linkingnumber}{the linking number of the point and the surface.
If the point is \emph{on} the surface (\code{distance==0}),
the (mathematical) linking number is undefined,
and the computed \code{linkingnumber} is \code{NA} (integer).
}

\item{inside}{whether the point is inside the surface; a logical.
This is currently set to \code{linkingnumber != 0}.
If the \code{linkingnumber} is \code{NA} (integer),
then \code{inside} is \code{NA} (logical).
}  

\item{timecalc}{the time to do the calculations, in seconds}

If the row names of \code{p} are unique,
they are copied to the row names of the output.
\cr
In case of error, the function returns \code{NULL}.
}



\references{
\bold{Point in polygon --- Wikipedia, The Free Encyclopedia}.
\url{https://en.wikipedia.org/w/index.php?title=Point_in_polygon&oldid=1139808558}.
2023.

\bold{Spivak, Michael}.
\emph{A Comprehensive Introduction to Differential Geometry}.
Volume 1.
3rd edition.
Publish or Perish.
1999.

\bold{West, G. and M. H. Brill.}
Conditions under which Schrödinger object colors are optimal.
\bold{Journal of the Optical Society of America}.
73. pp. 1223-1225. 1983.
}

\seealso{
\code{\link[zonohedra:inside2trans]{zonohedra::inside2trans()}}
}


\examples{
D50.eye = product( D50.5nm, 'varmat', xyz1931.1nm, wave=seq(360,830,by=5) )

insideSchrodingerColors( D50.eye,  c(30,50,70,  0,0,0)  )

##    p.1 p.2 p.3 distance linkingnumber inside     timecalc
##  1  30  50  70 7.516577             1   TRUE 4.858999e-04
##  2   0   0   0 0.000000            NA     NA 9.399839e-06
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{adaptation}% use one of  RShowDoc("KEYWORDS")

