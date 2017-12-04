% File computeCCT.Rd
\encoding{UTF-8}

\name{computeCCT}
\title{Compute Correlated Color Temperature (CCT) of Light Spectra}
\alias{computeCCT}
\alias{computeCCT.colorSpec}
\alias{CCTfromXYZ}

\description{
  Compute the CCT in Kelvin degrees of a \bold{colorSpec} object with \code{type} equal to \code{'light'}
}
\usage{
\S3method{computeCCT}{colorSpec}( x )

CCTfromXYZ( XYZ )

}
\arguments{
  \item{x}{an \bold{colorSpec} \R object with \code{type} equal to \code{'light'}, and M spectra }
  
  \item{XYZ}{a  3-vector with XYZ tristimulus values }
    
}

\details{
  In \code{computeCCT}, for each spectrum, XYZ is computed using \code{\link{xyz1931.1nm}},
  and the result passed to \code{CCTfromXYZ}. 
  If the \code{quantity} of \code{x} is \code{'photons'} (actinometric) each spectrum is
  converted to \code{'power'} (radiometric) on the fly.\cr  
  In \code{CCTfromXYZ}, the CCT is computed using Robertson's Method,
  which can compute CCTs in the interval [1666.7K,Inf], see \bold{References}.\cr

}

\value{
  \code{computeCCT} returns a numeric vector of length M, where M is the number of spectra in \code{x}.
  The vector's \code{names} is set to \code{specnames(x)}.\cr
  If the \code{type} of \code{x} is not \code{'light'}, then a warning is issued and all values are \code{NA}.\cr
  \code{CCTfromXYZ} returns a single number, but if the chromaticity of XYZ is too far from the daylight locus, then it returns \code{NA}.
  
}

\note{
The table in the \bold{References} contains an error at 325 mired,
which was corrected in \bold{Source}.

}

\source{
http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_T.html

}

\references{
Günther Wyszecki and W. S. Stiles.
Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.
John Wiley & Sons, 1982.
Table 1(3.11). pp. 227-228.
}

\seealso{
  \code{\link{type}}, 
  \code{\link{quantity}},
  \code{\link{xyz1931}},  
  \code{\link{specnames}}
}

\examples{
computeCCT( D65.1nm )       # returns 6502.068
computeCCT( A.1nm )         # returns 2855.656
CCTfromXYZ( c(1,1,1) )      # returns 5454.028

moon = readSpectra( system.file( "extdata/sources/moonlight.txt", package='colorSpec' ) )
computeCCT( moon )          # returns 4482.371
}


\keyword{light}

