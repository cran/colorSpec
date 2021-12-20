% File computeCRI.Rd
\encoding{UTF-8}

\name{computeCRI}
\title{Compute Color Rendering Index (CRI) of Light Spectra}
\alias{computeCRI}
\alias{computeCRI.colorSpec}

\description{
Compute the CIE 1974 color rendering index (CRI) of a light spectrum, 
called the \emph{the test illuminant}.\cr
From the given spectrum a \emph{reference illuminant} is selected with the same CCT 
(Correlated Color Temperature).
A selected set of 8 color samples is rendered in XYZ (1931) with both illuminants
and 8 color differences are computed in a special \bold{CIEUVW} color space.
For each color difference a CRI is computed, where 100 is a perfect color match.
The final CRI is the average of these 8 CRI values.
}

\usage{
\S3method{computeCRI}{colorSpec}( x, adapt=TRUE, attach=FALSE, tol=5.4e-3 )
}

\arguments{
  \item{x}{an \bold{colorSpec} \R object with \code{type} equal to \code{'light'}, and exactly 1 spectrum }
  
  \item{adapt}{if \code{TRUE}, then a special chromatic adaption is performed, see \bold{Details} }
    
  \item{attach}{if \code{TRUE}, then a large list of intermediate calculations
  is attached to the returned number, as attribute \code{data}. 
  This attached list includes data for all special 14 color samples, although
  the last 6 do not affect the returned CRI.}
        
  \item{tol}{for the CRI to be meaningful the chromaticities of the test and reference illuminants must be sufficiently close in the CIE
  1960 uniform chromaticity space.
  If the tolerance is exceeded, the function returns \code{NA}.
  The default \code{tol=5.4e-3} is the one recommended by the CIE, 
  but the argument allows the user to override it.}
}

\details{
The CCT of \code{x} is computed by \code{computeCCT()} with default options.
\cr
If \code{adapt} is \code{TRUE} the 8 test \code{uv} points
are chromatically adapted from the test illuminant to the reference illuminant
using a special von Kries type transformation; see \cite{Oleari} and
\cite{Wikipedia}.
The test UVW values are computed relative to the reference illuminant.
\cr
If \code{adapt} is \code{FALSE} the 8 test \code{uv} points
are \emph{not} chromatically adapted,
and the test UVW values are computed relative to the test illuminant.
}

\value{
\code{computeCRI()} returns a single number \eqn{\le} 100.
In case of ERROR it returns \code{NA}. 
If \code{attach} is \code{TRUE} a large list
of intermediate calculations is attached to the returned number.
}

\source{
The test color reflectance spectra are taken from:\cr
http://www.lrc.rpi.edu/programs/nlpip/lightinganswers/lightsources/scripts/NLPIP_LightSourceColor_Script.m

}

\references{
Oleari, Claudio, Gabriele Simone.
\bold{Standard Colorimetry: Definitions, Algorithms and Software}.
John Wiley. 2016.
pp. 465-470.


Günther Wyszecki and W. S. Stiles.
\bold{Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition}.
John Wiley & Sons, 1982.
Table 1(3.11). p. 828.

Wikipedia.
\bold{Color rendering index}.
\url{https://en.wikipedia.org/wiki/Color_rendering_index}

Hunt, R. W. G. and M. R. Pointer.
\bold{Measuring Colour}.  4th edition.
John Wiley & Sons.
2011.
Appendix 7.
}

\seealso{
\code{\link{type}}, 
\code{\link{xyz1931}},  
\code{\link{computeCCT}}
}

\examples{
computeCRI( subset(Fs.5nm,'F2') )       # returns 64.15195
computeCRI( subset(Fs.5nm,'F4') )       # returns 51.36348
}


\keyword{light}
