% File plotPatchesRGB.Rd

\name{plotPatchesRGB}
\title{Plot Patches defined by Linear RGB}
\alias{plotPatchesRGB}


\description{
RGB patches are a very common way of comparing color renderings.
This function draws rectangular patches,
and can also draw triangles formed by omitting one vertex from the rectangle.
}

\usage{
plotPatchesRGB( obj, normalize=FALSE, gamma='sRGB', background='gray50',
                labels=TRUE, shape='full', add=FALSE )
}

\arguments{

\item{obj}{
an Mx3 matrix of linear RGBs for M patches, with assigned \code{rownames}.
Or a \code{data.frame} containing a unique matrix column whose name
starts with the string \code{'RGB'}.
If the \code{data.frame} also has columns \code{LEFT,TOP,WIDTH,HEIGHT} then these
are used to place the patches.
The Y coordinate increases going down the page.}

\item{normalize}{If \code{TRUE} then all RGBs are scaled so their maximum is 1}

\item{gamma}{
either the string \code{'sRGB'} or a positive number giving the gamma function of the display;
see \code{\link{DisplayRGBfromLinearRGB}()} for details.}

\item{background}{
the color for the background behind all the patches.
If it is a character string, it is passed directly to \code{par()} as parameter \code{bg}.
If it is a numeric vector of length 3, it is is treated as linear RGB and passed
to \code{DisplayRGBfromLinearRGB()} with the given \code{gamma} (just like the matrix \code{obj}).
If it is a numeric vector of length 1, it is interpreted as linear graylevel
and replicated to length 3.}
  
\item{labels}{if \code{TRUE} then the \code{rownames} are drawn
using \code{\link[graphics]{text}()}.}

\item{shape}{
If \code{shape='full'} (the default) then the full rectangle is drawn.
If \code{'half'} then the rectangle is shrunk to 1/2 size, and with the same center.
If \code{shape} is one of \code{'left'}, \code{'right'}, 
\code{'bottom'}, \code{'top'} then only a 
half-rectangle is drawn, but keeping the specified side.
If \code{shape} is one of \code{'topleft'}, \code{'topright'}, 
\code{'bottomleft'}, \code{'bottomright'},
then only a triangular half of the rectangle is drawn, 
but keeping the specified vertex. }
               
\item{add}{if \code{TRUE} then the patches are added to an existing plot.}
}

\value{ \code{TRUE} if successful, and \code{FALSE} otherwise}

\details{
If \code{obj} is a matrix, or a data.frame without columns \code{LEFT,TOP,WIDTH,HEIGHT},
then the patches are drawn vertically stacked and abutting from top to bottom.
See the \bold{blueflame} vignette for examples.\cr
And see the \bold{gallery} vignette for examples of a data.frame \emph{with} those columns.
}

\note{
The gamma function of the display, aka electro-optical conversion function or EOCF,
maps the interval [0,1] to itself in a non-linear fashion (unless \code{gamma=1}).
}

\seealso{
\code{\link{DisplayRGBfromLinearRGB}}
}

\keyword{RGB}

