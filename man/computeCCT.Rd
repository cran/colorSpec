% File computeCCT.Rd
\encoding{UTF-8}

\name{computeCCT}
\title{Compute Correlated Color Temperature (CCT) of Light Spectra}
\alias{computeCCT}
\alias{computeCCT.colorSpec}
\alias{CCTfromXYZ}
\alias{CCTfromuv}

\description{
  Compute the CCT, in K,
  of a \bold{colorSpec} object with \code{type} equal to \code{'light'}
}
\usage{
\S3method{computeCCT}{colorSpec}( x, method='robertson', strict=FALSE, c2=1.4388e7 )

CCTfromXYZ( XYZ, method='robertson', strict=FALSE, c2=1.4388e7 )

CCTfromuv( uv, method='robertson', strict=FALSE, c2=1.4388e7 )
}

\arguments{
\item{x}{an \bold{colorSpec} \R object with \code{type} equal to \code{'light'}, and M spectra }
  
\item{XYZ}{a 3-vector with XYZ tristimulus values (CIE 1931)}

\item{uv}{a 2-vector with uv chromaticity values (CIE UCS 1960)}

\item{method}{For \code{'robertson'} a lookup table with linear interpolation is used.
For \code{'lm'} the closest point on the Planckian locus is calculated directly
using the Levenberg-Marquardt method. See \bold{Details}.
Partial matching of the value is implemented.
}

\item{strict}{The CIE considers the CCT of a chromaticity to be meaningful only
if its distance to the Planckian locus is less than 0.05 [in CIE UCS 1960].
If \code{strict=FALSE} then this condition is ignored.
Otherwise, if the the condition is not satisfied,
the functions returns \code{NA_real_}.
}

\item{c2}{For points on the Planckian locus, the function \code{\link{planckSpectra}()}
is called, and this argument is passed to it.
It is only used when \code{method='lm'}.
See \code{\link{planckSpectra}} for the definition,
and the special value \code{c2='calc'}.
}
}

\details{
In \code{computeCCT()}, for each spectrum,
XYZ is computed using \code{\link{xyz1931.1nm}},
and the result passed to \code{CCTfromXYZ()},
and then to \code{CCTfromuv()}, along with the other three arguments. 
If the \code{quantity} of \code{x} is \code{'photons'} (actinometric) each spectrum is
converted to \code{'energy'} (radiometric) on the fly.\cr  
In \code{CCTfromuv()}, Robertson's Method uses a lookup table with 31 entries
and can compute CCTs in the interval [1666.7K,Inf], 
see \bold{References} and \bold{Source}.
If \code{method='lm'} the result of Robertson's Method is used as an initial point
for Levenberg-Marquardt, which then "polishes" it.
The package \bold{minpack.lm} is required here.
The difference between the 2 methods is small.
To enable diagnostic information about the polishing,
execute \code{cs.options(loglevel='INFO')} before calling one of these functions.
}

\value{
\code{computeCCT()} returns a numeric vector of length M,
where M is the number of spectra in \code{x}.
The vector's \code{names} is set to \code{specnames(x)}.\cr
If the \code{type} of \code{x} is not \code{'light'},
then a warning is issued and all values are \code{NA_real_}.\cr
\code{CCTfromXYZ()} and \code{CCTfromuv()} return a single number,
but in case of ERROR, they return \code{NA_real_}.
}

\note{
The lookup table on page 228 in Wyszecki & Stiles (see \bold{References})
contains an error at 325 mired,
which was corrected by Bruce Lindbloom (see \bold{Source}).

}

\source{
\url{http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_T.html}

}

\references{
Robertson, A. R.
Computation of correlated color temperature and distribution temperature.
Journal of the Optical Society of America.
58. pp. 1528-1535 (1968).

Wyszecki, Günther and W. S. Stiles.
\bold{Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.}
John Wiley & Sons, 1982.
Table 1(3.11). pp. 227-228.
}

\seealso{
  \code{\link{type}}, 
  \code{\link{quantity}},
  \code{\link{xyz1931}},
  \code{\link{planckSpectra}},    
  \code{\link{specnames}}
}

\examples{
computeCCT( D65.1nm )               # returns 6502.068
computeCCT( D65.1nm, meth='lm' )    # returns 6502.711
computeCCT( A.1nm )                 # returns 2855.656
computeCCT( A.1nm, meth='lm' )      # returns 2855.583
CCTfromXYZ( c(1,1,1) )              # returns 5454.028

moon = readSpectra( system.file( "extdata/sources/moonlight.txt", package='colorSpec' ) )
computeCCT( moon )                  # returns 4482.371
}


\keyword{light}

