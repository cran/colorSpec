% File resample.Rd

\name{resample}
\title{resample a colorSpec Object to new wavelengths}
\alias{resample}
\alias{resample.colorSpec}

\description{
  interpolate or smooth to new wavelengths.  Simple extrapolation is also performed.
}
\usage{
\S3method{resample}{colorSpec}( x, wavelength, method='auto', span=0.02 )

}

\arguments{
\item{x}{a \bold{colorSpec} object}

\item{wavelength}{vector of new wavelengths, in nanometers}

\item{method}{interpolation methods available are 
    \code{'sprague'}, \code{'spline'}, and \code{'linear'}.
    Also available is \code{'auto'} which means to use \code{'sprague'}
    if \code{x} is regular, and \code{'spline'} otherwise.
    This is the CIE recommendation.\cr
    Smoothing methods available are \code{'loess'}.
    See \bold{Details}.
    }

\item{span}{smoothing argument passed to \code{\link[stats]{loess}}
            before interpolation, and not used by other methods.
           The default value \code{span=0.02} is suitable for
           \code{.scope} spectra but may be too small in many other cases.}  
}

\value{
\code{resample(x)} returns a \bold{colorSpec} object 
with the new \code{wavelength}.
Other properties, e.g. \code{\link{organization}}, \code{\link{quantity}}, ...,
are preserved.\cr
In case of ERROR, the function returns \code{NULL}.
           
}    


\details{
For method \code{'sprague'} the quintic polynomial in \cite{De Kerf} is used.
Six weights are applied to nearby data values: 3 in front and 3 behind.
The \code{'sprague'} method is only supported when \code{x} is regular.

For method \code{'spline'} the function \code{\link[stats]{spline}}
is called with \code{method='natural'}.
Four weights are applied to nearby data values: 2 in front and 2 behind.
The \code{'spline'} method is supported even when \code{x} is irregular.

For method \code{'linear'} the function \code{\link[stats]{approx}} is called.

For method \code{'loess'} the function \code{\link[stats]{loess}} is called
with the given \code{span} parameter.
Smoothing is most useful for noisy data, e.g. raw data from a spectrometer.
I have found that \code{span=0.02} works well for \code{.scope} files, 
but this may be too small in other cases,
when it triggers an error in \code{loess()}. 

For the non-linear methods,
undershoot can sometimes create negative values, and some effort is made to clip
these to 0, but only when appropriate.
}


\references{
De Kerf, Joseph L. F.
\bold{The interpolation method of Sprague-Karup}.
\emph{Journal of Computational and Applied Mathematics}.
volume I, no 2, 1975.
equation (S).
}


\seealso{
\code{\link{organization}},   
\code{\link{quantity}},
\code{\link{wavelength}},
\code{\link{is.regular}},
\code{\link[stats]{spline}},
\code{\link[stats]{approx}},
\code{\link[stats]{loess}}
}

\examples{
path = system.file( "extdata/sources/pos1-20x.scope", package='colorSpec' )
y = readSpectra( path )
# plot noisy data in gray
plot( y, col='gray' )
# plot smoothed plot in black on top of the noisy one to check quality
plot( resample( y, 200:880, meth='loess', span=0.02 ), col='black', add=TRUE )
}

\keyword{colorSpec}
