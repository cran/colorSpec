% File plot.Rd

\name{plot}
\alias{plot.colorSpec}

\title{plot spectra}

\description{
plot the spectra in a \bold{colorSpec} object as lines or points
}

\usage{
\S3method{plot}{colorSpec}( x, color=NULL, subset=NULL, main=TRUE, legend=TRUE, CCT=FALSE, add=FALSE, ... )
}

\arguments{
\item{x}{a colorSpec object}

\item{color}{
If \code{color=NULL} then colors are computed from the spectra themselves.
If \code{type(x)} is \code{'material'} the color is computed using
illuminant \code{D65.1nm} and responder \code{BT.709.RGB} with no further normalization.
Otherwise the spectrum color is faked by changing its \code{quantity} to \code{'energy'}
and taking the \code{\link{product}} with \code{BT.709.RGB}.
The resulting RGBs are normalized to have a maximum of 1.
This RGB normalization is done \emph{before} processing the \code{subset} argument.\cr
If \code{color='auto'} then a suitable set of colors is generated using
\code{\link{colorRamp}()}.\cr
Otherwise \code{color} is passed on to \code{\link{lines.default}()}
as the \code{col} argument,  e.g. \code{col='black'}.
}

\item{subset}{specifies a subset of \code{x} to plot;
see \code{\link{subset}()} for acceptable arguments.}

\item{main}{If \code{main=TRUE} then a main title is generated from the file \code{'path'} in
the \code{metadata} list, or from \code{deparse(substitute(x))}.
If \code{main=FALSE} then no main title is displayed.
And if \code{main} is a string then that string is used as the main title.
}

\item{legend}{If \code{legend=TRUE} then a pretty legend using \code{\link{specnames}(x)}
is placed in the \code{'topright'} corner of the plot.  
If \code{legend} is a string it is interpreted as naming a corner of the 
plot and passed as such to the function \code{\link[graphics]{legend}}.
If \code{legend=FALSE} then no legend is drawn.
}

\item{CCT}{
If \code{CCT=TRUE} and the \code{type} of \code{x} is \code{'light'} then
the CCT of each spectrum is added to the legend; see \code{\link{computeCCT}()}.
The package \bold{spacesXYZ} must be installed.
}

\item{add}{
If \code{add=TRUE} then the lines are added to an existing plot,
and these arguments are ignored:
\code{main}, \code{ylab}, \code{xlim}, \code{ylim},
and \code{log}; see \bold{Details}.
}

\item{...}{other graphical parameters, see \bold{Details}}

}
\value{
\code{TRUE} or \code{FALSE}
}

\details{
Commonly used graphical parameters are:

\describe{

\item{\code{type}}{passed to \code{\link{lines.default}()}, 
with default \code{type='l'}.
Other valid values are \code{'p'} (points), \code{'b'}, \code{'c'}, \code{'o'}, 
\code{'h'}, \code{'S'}, \code{'s'}, and \code{'n'},
see \code{\link[graphics]{plot}()} for their meanings.
\cr
An additional \code{type='step'} is available.
This option draws each spectrum as a \emph{step function},
similar to \code{'S'} and \code{'s'},
except that the jumps are \emph{between} the wavelengths
(with appropriate extensions at min and max wavelengths).
The function \code{\link{segments}()} is used for the drawing.
For \code{type='step'}, \code{lwd} and \code{lty} should be vectors of length 1 or 2.
If the length of \code{lwd} is 1, then horizontal segments are draw with that width,
but vertical segments are not drawn.
If the length of \code{lwd} is 2, 
then vertical segments are draw with width \code{lwd[2]}.
If the length of \code{lty} is 2, then the styles are applied to the
horizontal and vertical segments in that order.
If the length of \code{lty} is 1, then that style is applied to both
horizontal and vertical segments.
For examples of this plotting option, see the vignette 
\href{../doc/convexity.pdf}{\bold{Convexity and Transitions}}.}

\item{\code{lwd, lty}}{passed to \code{\link{lines.default}()},
except when \code{type='step'} when they are passed to \code{\link{segments}()}.
In the former case these can be vectors, and components are passed
sequentially to each spectrum, similar to \code{\link{matplot}()}.
In the latter case, see the description in \code{type}.
The default value for both is 1.}

\item{\code{pch}}{passed to \code{\link{lines.default}()},
but it only has meaning when \code{type='p'}, \code{'b'}, or \code{'o'}.
This can be a vector, and components are passed sequentially to each spectrum.}

\item{\code{ylab}}{If \code{ylab} is a string 
then it is passed on to \code{\link{plot.default}()},
otherwise suitable default string is generated.}

\item{\code{xlim, ylim}}{If \code{xlim} and \code{ylim} are 2-vectors, they are passed
to \code{\link{plot.default}}. 
If one of the components is \code{NA} then a suitable default is supplied.}

\item{\code{log}}{passed on to \code{\link{plot.default}()}. 
Care must be taken for \code{y} because
many spectra are 0 at some wavelengths, and even negative. Use \code{ylim} in such cases.}
}
}

\examples{
plot( 100 * BT.709.RGB )
plot( xyz1931.1nm, add=TRUE, lty=2, legend=FALSE )
}

\seealso{
\code{\link{computeCCT}()},
\code{\link{subset}()},
\code{\link[graphics]{lines}()},
\code{\link[graphics]{segments}()},
\code{\link[graphics]{plot}()},
\code{\link[graphics]{matplot}()},
\code{\link[grDevices]{colorRamp}()}
}

\keyword{colorSpec}

