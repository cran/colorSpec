% File interpolate.Rd

\name{interpolate}
\title{interpolate spectra}
\alias{interpolate}
\alias{interpolate.colorSpec}

\description{
  interpolate along a 1-parameter path of spectra
}
\usage{
\S3method{interpolate}{colorSpec}( x, p, pout, pname=deparse(substitute(p)) )

}

\arguments{
\item{x}{a \bold{colorSpec} object, typically with multiple spectra}

\item{p}{a numeric vector with \code{length(p)==numSpectra(x)}. The value \code{p[i]} is associated with the \code{i}'th spectrum in \code{x}.}

\item{pout}{a numeric vector of parameter values at which interpolation of the spectra in \code{x} take place}

\item{pname}{the name of the parameter \code{p}}
}


\value{
\code{interpolate(x)} returns a \bold{colorSpec} object \code{y}
with a spectrum for each value in \code{pout}.
The organization of \code{y} is \code{'df.row'},
and \code{extradata(y)} has a single column which is a copy of \code{pout}.
The name of the column is \code{pname}.
The names in \code{specnames(y)} are \code{<pname>=<pout>}.
Other properties of \code{y},
e.g. \code{\link{wavelength}}, \code{\link{quantity}}, ...,
are the same as \code{x}.\cr
In case of ERROR, the function returns \code{NULL}.
}    


\details{
Each spectrum in \code{x} can be thought of as a point in a high-dimensional space,
and each point has a real-valued parameter associated with it.
The function performs natural spline interpolation on these points,
one coordinate at a time.
For each wavelength value it calls \code{\link{spline}}
with \code{method='natural'}.
}


\seealso{
\code{\link{organization}},   
\code{\link{wavelength}},
\code{\link{extradata}},
\code{\link[stats]{spline}}
}

\examples{
path = system.file( "extdata/stains/PhenolRed-Fig7.txt", package="colorSpec" )
wave = 350:650
phenolred = readSpectra( path, wavelength=wave )
pH = as.numeric( sub( '[^0-9]+([0-9]+)$', '\\\\1', specnames(phenolred) ) )
pHvec = seq(min(pH),max(pH),by=0.05)
phenolinterp = interpolate( phenolred, pH, pHvec )
}

\keyword{colorSpec}
