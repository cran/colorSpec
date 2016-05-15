% File plot.Rd

\name{plot}
\alias{plot.colorSpec}

\title{plot spectra}

\description{
plot the spectra in a \bold{colorSpec} object as lines
}

\usage{
\S3method{plot}{colorSpec}( x, color=NULL, subset=NULL, main=TRUE, legend=TRUE, CCT=FALSE, add=FALSE, ... )
}

\arguments{
\item{x}{a colorSpec object}

\item{color}{
If \code{color=NULL} then colors are computed from the spectra themselves.
If the \code{type} of \code{x} is \code{'material'} the color is computed using
illuminant \code{D65.1nm} and responder \code{BT.709.RGB} with no further normalization.
Otherwise the spectrum color is faked by changing its \code{quantity} to \code{'power'}
and taking the \code{\link{product}} with \code{BT.709.RGB}.
The resulting RGBs are normalized to have a maximum of 1.
This RGB normalization is done \emph{before} processing the \code{subset} argument.\cr
If \code{color='auto'} then a suitable set of colors is generated using
\code{\link[grDevices]{colorRampPalette}}.\cr
Otherwise \code{color} is passed along to \code{\link{lines}} as the \code{col} argument, 
e.g. \code{col='black'}.
}

\item{subset}{specifies a subset of \code{x} to plot; see \code{\link{subset}} for acceptable arguments.}

\item{main}{If \code{main=TRUE} then a main title is generated from the file \code{'path'} in
the \code{metadata} list, or from \code{deparse(substitute(x))}.
If \code{main=FALSE} then no main title is displayed.
And if \code{main} is a string then that string is used as the main title.
}

\item{legend}{If \code{legend=TRUE} then a pretty legend using \code{\link{specnames}}
is placed in the \code{'topright'} corner of the plot.  
If \code{legend} is a string it is interpreted as naming a corner of the 
plot and passed as such to the function \code{\link[graphics]{legend}}.
If \code{legend=FALSE} then no legend is drawn.
}

\item{CCT}{
If \code{CCT=TRUE} and the \code{type} of \code{x} is \code{'light'} then
the CCT of each spectrum is added to the legend; see \code{\link{computeCCT}}.
}

\item{add}{
If \code{add=TRUE} then lines are added to an existing plot,
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
\item{\code{lty, lwd}}{passed on to \code{\link{lines}}, with no checks}

\item{\code{pch}}{If \code{pch} is an integer, it is passed on to \code{\link{points}}
which adds points on top of the lines}

\item{\code{ylab}}{If \code{ylab} is a string 
then it is passed on to \code{\link{plot.default}},
otherwise suitable default string is generated.}

\item{\code{xlim, ylim}}{If \code{xlim} and \code{ylim} are 2-vectors, they are passed
to \code{\link{plot.default}}. 
If one of the components is \code{NA} then a suitable default is supplied.}

\item{\code{log}}{passed on to \code{\link{plot.default}}. 
Care must be taken for \code{y} because
many spectra are 0 at some wavelengths, and even negative. Use \code{ylim} in such cases.}
}
}

\examples{
print( xyz1931.1nm )

xyz1931.1nm     # same thing, just calls print()
}

\seealso{
\code{\link[base]{print}},
\code{\link[base]{summary}}  
in base.

}

\keyword{colorSpec}

