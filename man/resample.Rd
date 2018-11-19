% File resample.Rd

\name{resample}
\title{resample a colorSpec Object to new wavelengths}
\alias{resample}
\alias{resample.colorSpec}

\description{
  interpolate or smooth to new wavelengths.  
  Simple extrapolation and clamping is also performed.
}
\usage{
\S3method{resample}{colorSpec}( x, wavelength, method='auto', span=0.02, extrapolation='const', clamp='auto' )

}

\arguments{
\item{x}{a \bold{colorSpec} object}

\item{wavelength}{vector of new wavelengths, in nanometers}

\item{method}{interpolation methods available are 
    \code{'sprague'}, \code{'spline'}, and \code{'linear'}.
    Also available is \code{'auto'} which means to use \code{'sprague'}
    if \code{x} is regular, and \code{'spline'} otherwise.
    % This is the CIE recommendation.\cr
    An available smoothing method is \code{'loess'}.
    See \bold{Details}.
    }

\item{span}{smoothing argument passed to \code{\link[stats]{loess}()}
            during interpolation, and not used by other methods.
           The default value \code{span=0.02} is suitable for
           \code{.scope} spectra but may be too small in many other cases.}  
           
\item{extrapolation}{extrapolation methods available are
    \code{'const'} and \code{'linear'}.
    These can be abbreviated to the initial letter.
    Also available is a numeric value, which is used for simple padding.
    See \bold{Details}.\cr
    Also available is a vector or list of length 2 that combines 2 of the above.
    The first item is used on the low side (shorter wavelengths),
    and the second item is used on the high side (longer wavelengths).
    }

\item{clamp}{clamp methods available are
    \code{'auto'}, \code{TRUE}, and \code{FALSE}.
    Also available is a numeric vector of length 2,
    which defines the clamping interval.    
    See \bold{Details}.
    }
}

\value{
\code{resample(x)} returns a \bold{colorSpec} object 
with the new \code{wavelength}.
Other properties, e.g. \code{\link{organization}}, \code{\link{quantity}}, ...,
are preserved.\cr
In case of ERROR, the function returns \code{NULL}.
           
}    


\details{
If method is \code{'sprague'}, the quintic polynomial in \cite{De Kerf} is used.
Six weights are applied to nearby data values: 3 on each side.
The \code{'sprague'} method is only supported when \code{x} is regular.

If method is \code{'spline'}, the function \code{stats::\link{spline}()}
is called with \code{method='natural'}.
% Four weights are applied to nearby data values: 2 on each side.
The \code{'spline'} method is supported even when \code{x} is irregular.

If method is \code{'linear'}, the function \code{stats::\link{approx}()} is called.
Two weights are applied to nearby data values: 1 on each side.
The \code{'linear'} method is supported even when \code{x} is irregular.

If method is \code{'loess'}, the function \code{stats::\link[stats]{loess}()} is called
with the given \code{span} parameter.
Smoothing is most useful for noisy data, e.g. raw data from a spectrometer.
I have found that \code{span=0.02} works well for Ocean Optics \code{.scope} files, 
but this may be too small in other cases, which triggers an error in \code{stats::loess()}. 
The \code{'loess'} method is supported even when \code{x} is irregular.

If extrapolation is \code{'const'},
the extreme values at each end are simply extended.
If extrapolation is \code{'linear'},
the line defined by the 2 extreme values at each end is used for extrapolation.
If the ultimate and penultimate wavelengths are equal,
then this line is undefined and the function reverts to \code{'const'}.

If clamp is \code{'auto'},
output values are clamped to the physically realizable interval
appropriate for \code{x}.
This is the interval [0,1]  when
\code{quantity(x)} is \code{'reflectance'} or \code{'transmittance'},
and the interval [0,\eqn{\infty})  otherwise.
Exception: If an input spectrum in \code{x} violates a limit,
then clamping the output spectrum to this limit is \emph{NOT} enforced.
This happens most frequenty for theoretical (or matrixed) cameras,
such as \code{\link{BT.709.RGB}}.

If clamp is \code{TRUE}, the result is the same as \code{'auto'},
but with no exceptions.
If clamp is \code{FALSE}, then no clamping is done.

If clamp is a numerical interval, then clamping is done to that interval,
with no exceptions.
The two standard intervals mentioned above can be expressed in \bold{R}
as \code{c(0,1)} and \code{c(0,Inf)} respectively.
}


\references{
De Kerf, Joseph L. F.
\bold{The interpolation method of Sprague-Karup}.
\emph{Journal of Computational and Applied Mathematics}.
volume I, no 2, 1975.
equation (S).
}


\seealso{
\code{\link{organization}()},   
\code{\link{quantity}()},
\code{\link{wavelength}()},
\code{\link{is.regular}()},
\code{\link{theoreticalRGB}},
\code{\link{spline}()},
\code{\link{approx}()},
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
