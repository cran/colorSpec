% File chop.Rd

\name{chop}
\title{chop spectra into low and high parts}
\alias{chop}
\alias{chop.colorSpec}

\description{
  chop all spectra in a \bold{colorSpec} object into low and high parts at a blending interval
}
\usage{
\S3method{chop}{colorSpec}( x, interval, adj=0.5 )

}

\arguments{
  \item{x}{a \bold{colorSpec} object}
  
  \item{interval}{a numeric vector with length 2 giving the endpoints of the interval, in nm}
  
  \item{adj}{a number in [0,1] defining weights of low and high parts over the interval}
    
}
    
\value{
  \code{chop(x)} returns a \bold{colorSpec} object with twice the number of spectra in \code{x}
  and with \code{\link{organization}} equal to \code{'matrix'}.
  The names of the new spectra are formed by appending ".lo" and ".hi" to the original spectrum names.
  
}    


\details{
  For each spectrum, the low and high parts sum to the original spectrum.
  The low part vanishes on the right of the interval, and the high part vanishes on the left.

}


\seealso{
  \code{\link{organization}},  
}

\examples{
# chop blue butane flame into diatomic carbon and hydrocarbon parts
path = system.file( "extdata/sources/blueflame.txt", package="colorSpec" )
blueflame = readSpectra( path, seq(375,650,0.5) )
plot( chop( blueflame, interval=c(432,435), adj=0.8 ) )

# chop 'white' LED into blue and yellow parts
path = system.file( "extdata/sources/Gepe-G-2001-LED.sp", package="colorSpec" )
LED = readSpectra( path )
plot( chop( LED, c(470,495) ) )
}

\keyword{colorSpec}
