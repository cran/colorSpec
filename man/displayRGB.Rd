% File displayRGB.Rd

\name{DisplayRGB}
\title{Compute Display RGB from Linear RGB}
\alias{DisplayRGBfromLinearRGB}


\description{
All RGB displays have a non-linear "gamma function" of some sort.
This function converts from linear RGB to an RGB appropriate for 
the gamma function of the display.
}

\usage{
DisplayRGBfromLinearRGB( RGB, gamma='sRGB' )

}

\arguments{
  \item{RGB}{linear RGB values organized as a vector or matrix of any size; 
                  all 3 channels are treated the same way so size does not matter}
  \item{gamma}{either the string \code{'sRGB'} or a positive number giving the gamma of the display.}}

\value{
The function first clamps the input \code{RGB} to the interval [0,1].
If \code{gamma='sRGB'} (not case-sensitive) it then maps [0,1] to [0,1] using the special
piecewise-defined sRGB function, see \cite{Wikipedia}.
In case \code{gamma} is a positive number, the function raises all values to the power \code{1/gamma}.
The dimensions and names of the input are copied to the output. \cr
In case of error, the function returns the clamped input values.
}


\source{
Wikipedia.
\bold{sRGB}.
http://en.wikipedia.org/wiki/sRGB
}

\seealso{
  \code{\link{RGBfromXYZ}}
}

\examples{
DisplayRGBfromLinearRGB( c(0.2, 0.5) )
# [1] 0.4845292 0.7353570     #  this is display sRGB, in [0,1]

DisplayRGBfromLinearRGB( c(-0.1, 0.2, 0.5, 1), 2.2 )
# [1] 0.0000000 0.4811565 0.7297401 1.0000000    #  gamma=2.2

x = seq( 0, 1, len=101)
plot( x, DisplayRGBfromLinearRGB(x), type='l' )
}

\keyword{RGB}

