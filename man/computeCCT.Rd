% File computeCCT.Rd
\encoding{UTF-8}

\name{computeCCT}
\title{Compute Correlated Color Temperature (CCT) of Light Spectra}
\alias{computeCCT}
\alias{computeCCT.colorSpec}


\description{
  Compute the CCT, in K,
  of a \bold{colorSpec} object with \code{type} equal to \code{'light'}
}
\usage{
\S3method{computeCCT}{colorSpec}( x, isotherms='robertson', locus='robertson', strict=FALSE )
}

\arguments{
\item{x}{an \bold{colorSpec} \R object with \code{type} equal to \code{'light'}, and M spectra }

\item{isotherms}{A character vector whose elements match one
of the available isotherm families:
\code{'robertson'}, \code{'mccamy'}, and \code{'native'}.
Matching is partial and case-insensitive.
When more than one family is given, a matrix is returned, see \bold{Value}.
When \code{isotherms='native'} the isotherms are defined implicitly
as lines perpendicular to the locus, see \bold{Details} in 
\code{\link[spacesXYZ:CCT]{spacesXYZ::CCTfromXYZ}()}.
The character \code{NA} (\code{NA_character_}) is taken as
a synonym for \code{'native'}.
}

\item{locus}{
valid values are \code{'robertson'} and \code{'precision'}, see above.
Matching is partial and case-insensitive.
}

\item{strict}{The CIE considers the CCT of a chromaticity \code{uv} to be meaningful only
if the distance from \code{uv} to the Planckian locus is less than or equal to 0.05 
[in CIE UCS 1960].
If \code{strict=FALSE}, then this condition is ignored.
Otherwise, the distance is computed along the corresponding isotherm,
and if it exceeds 0.05 the returned CCT is set to \code{NA}.}

}


\details{
In \code{computeCCT()}, for each spectrum,
XYZ is computed using \code{\link{xyz1931.1nm}},
and the result passed to \code{spacesXYZ::CCTfromXYZ()}.
If the \code{quantity} of \code{x} is \code{'photons'} (actinometric) each spectrum is
converted to \code{'energy'} (radiometric) on the fly.\cr 
}

\value{
\code{computeCCT()} returns a numeric vector of length M,
where M is the number of spectra in \code{x}.
The vector's \code{names} is set to \code{specnames(x)}.\cr
If the \code{type} of \code{x} is not \code{'light'},
then a warning is issued and all values are \code{NA_real_}.\cr
}


\references{
McCamy, C. S.
\emph{Correlated color temperature as an explicit function of chromaticity coordinates}.
Color Research & Application.
Volume 17.
Issue 2.
pages 142-144.
April 1992.


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
\code{\link{type}()}, 
\code{\link{quantity}()},
\code{\link{xyz1931}},
\code{\link{planckSpectra}()},    
\code{\link{specnames}()},
\code{\link[spacesXYZ:CCT]{spacesXYZ::CCTfromXYZ}()}
}

\examples{
computeCCT( D65.1nm )                       # returns 6502.068
computeCCT( D65.1nm, isotherms='native' )   # returns 6503.323
computeCCT( A.1nm )                         # returns 2855.656
computeCCT( A.1nm, isotherms='native' )     # returns 2855.662
computeCCT( A.1nm, isotherms='mccamy' )     # returns 2857.188

moon = readSpectra( system.file( "extdata/sources/moonlight.txt", package='colorSpec' ) )
computeCCT( moon )                  # returns 4482.371
}


\keyword{light}

