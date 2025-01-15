% File computeCRI.Rd
\encoding{UTF-8}

\name{computeCRI}
\title{Compute Color Rendering Index (CRI) of Light Spectra}
\alias{computeCRI}
\alias{computeCRI.colorSpec}
\alias{computeCRIdata}
\alias{computeCRIdata.colorSpec}

\description{
Compute the CIE 1974 color rendering index (CRI) of a light spectrum, 
called the \emph{the test illuminant}.\cr
From the given spectrum a \emph{reference illuminant} is selected with the same CCT 
(Correlated Color Temperature).
A selected set of 8 color samples is rendered in XYZ (1931) with both illuminants
and 8 color differences are computed in CIE 1964 UVW color space.
For each color sample a CRI is computed, where 100 is a perfect color match.
The final CRI is the average of these 8 CRI values.
}

\usage{
\S3method{computeCRI}{colorSpec}( x, CCT=NULL, adapt=TRUE, tol=5.4e-3, attach=FALSE )
\S3method{computeCRIdata}{colorSpec}( x, CCT=NULL, adapt=TRUE, tol=5.4e-3 )
}

\arguments{
  \item{x}{a non-empty \bold{colorSpec} object with \code{type} equal to \code{'light'}.
  The spectra in \code{x} are the \emph{test illuminants}.
  For \code{computeCRIdata()} there must be exactly 1 spectrum.
  For \code{computeCRI()} there can be multiple spectra.}
  
  \item{CCT}{the Correlated Color Temperature of the reference illuminant.
  If \code{CCT=NULL} (the default) then \code{CCT} is computed from the test illuminant.
  The user can override this with a target CCT, e.g. the advertised CCT for a particular light bulb.
  For \code{computeCRIdata()} there must be exactly one \code{CCT}.
  For \code{computeCRI()} the length of the vector must be either 1 or the number of spectra
  in \code{x}.}
  
  
  \item{adapt}{if \code{TRUE}, then a special chromatic adaption is performed, see \bold{Details} }

  \item{tol}{for the CRI to be meaningful the chromaticities of the test and reference illuminants must be sufficiently close in the CIE
  1960 uniform chromaticity space.
  If the tolerance is exceeded, the CRI is set to \code{NA_real_}.
  The default \code{tol=5.4e-3} is the one recommended by the CIE, 
  but the argument allows the user to override it.
  To ignore this test, set \code{tol=Inf}.}
  
  \item{attach}{if \code{TRUE} and there is exactly 1 spectrum in \code{x},
  then the list of intermediate calculations returned by \code{computeCRIdata()}
  is attached to the returned number, as attribute \code{data}.
  This attached list includes data for all 14 color samples,
  though only the first 8 are used to compute the CRI. If there is not exactly 1 spectrum in \code{x},
  then \code{attach} is ignored.}  
}

\details{
If not \code{NULL}, the CCT of \code{x} is computed by \code{computeCCT()} with default options.
\cr
When computing XYZs, the wavelengths of \code{x} and the color matching functions
of \code{xyz1931.1nm} are used.
\cr
If \code{adapt} is \code{TRUE} the 8 color sample \code{uv} points
are chromatically adapted from the test illuminant to the reference illuminant
using a special von Kries type transformation; see \cite{Oleari} and
\cite{Wikipedia}.
The color sample UVW values are computed with the reference illuminant.
\cr
If \code{adapt} is \code{FALSE} the 8 color sample \code{uv} points
are \emph{not} chromatically adapted,
and the color sample UVW values are computed with the test illuminant.
}

\value{
\code{computeCRI()} returns a vector of CRI values with length equal to the
number of spectra in \code{x}.  All values are \eqn{\le} 100.
In case of ERROR the CRI value is \code{NA_real_}. 
If \code{attach} is \code{TRUE} and \code{x} has exactly one spectrum,
a large list of intermediate calculations is attached to the returned number.

\code{computeCRIdata()} returns a list of intermediate calculations,
for the single spectrum in \code{x}.
The items in the list are:
\tabular{ll}{
\code{CCT} \tab the Correlated Color Temperature of the reference illuminant. \cr
\code{illum.ref} \tab the reference illuminant, which is a \bold{colorSpec}  object. \cr
\code{table1} \tab a data frame with 2 rows, with CIE data for the test and reference illuminants (see below). \cr
\code{Delta_uv} \tab  the distance between the illuminants in the CIE 1960 uniform chromaticity space. \cr
\code{table2} \tab a data frame with 14 rows, with CIE XYZ data of the 14 color samples. \cr
\code{table3} \tab a data frame with 14 rows, with CIE uv data of the 14 color samples. If argument \code{adapt} is \code{FALSE}, then \code{table3=NULL}. \cr
\code{table4} \tab a data frame with 14 rows, with UVW data of the 14 color samples. \cr
}

The columns of \code{table1} are:
\tabular{ll}{
\code{XYZ} \tab the CIE XYZ of the 2 illuminants. \cr
\code{xy} \tab the CIE xy of the 2 illuminants. \cr
\code{uv} \tab the CIE 1960 uv of the 2 illuminants. \code{Delta_uv} is the distance between these 2 points. \cr
}
The first row of \code{table1} is the given test illuminant, and the second row is the reference illuminant.
The initial letter of the rowname of the reference indicates the type:
the letter \code{P} means a Planckian illuminant, and the letter \code{D} means a Daylight illuminant.

The columns of \code{table2} are:
\tabular{ll}{
\code{referen} \tab the CIE XYZ of the color sample, as illuminated by the reference illuminant. \cr
\code{test} \tab the CIE XYZ of the color sample, as illuminated by the test illuminant. \cr
}

The columns of \code{table3} are:
\tabular{ll}{
\code{before} \tab the CIE 1960 uv of the color sample, as illuminated by the test illuminant. \cr
\code{after} \tab the \code{before} values, after adaptation to the reference illuminant. \cr
\code{difference} \tab \code{after} - \code{before} \cr
}

The columns of \code{table4} are:
\tabular{ll}{
\code{referen} \tab the UVW of the color sample, as illuminated by the reference illuminant. \cr
\code{test} \tab  the UVW of the color sample, as illuminated by the test illuminant. \cr
\code{DeltaE} \tab the distance between the test and reference UVWs \cr
\code{CRI} \tab the CRI of the color sample, which is a function of \code{DeltaE}.  \cr
}
The final CRI is the average of the \code{CRI} of the first 8 samples in \code{table4},
and the remaining samples are ignored.
}


\source{
The reflectance spectra of the 14 color samples are taken from:\cr
http://www.lrc.rpi.edu/programs/nlpip/lightinganswers/lightsources/scripts/NLPIP_LightSourceColor_Script.m

The wavelength vector is 360nm to 830nm with 5nm step.
The same data over 380nm to 780nm is in Appendix 7 of \cite{Hunt and Pointer}.
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
\code{\link{type}()}, 
\code{\link{xyz1931}},  
\code{\link{computeCCT}()}
}

\examples{
computeCRI( Fs.5nm )
##       F1       F2       F3       F4  F5  F6       F7       F8       F9      F10      F11      F12 
## 75.82257 64.15195 56.68144 51.36348  NA  NA 90.18452 95.50431 90.29347 81.03585 82.83362 83.06373 

computeCRI( Fs.5nm, adapt=FALSE )
##       F1       F2       F3       F4  F5  F6       F7       F8       F9      F10      F11      F12 
## 77.73867 65.38567 57.20553 50.65979  NA  NA 90.18551 95.96459 90.27063 82.86106 82.86306 83.10613 

computeCRI( subset(Fs.5nm,'F2') )
##       F2 
## 64.15195 

computeCRI( subset(Fs.5nm,'F2'), CCT=4200 )
##       F2
## 63.96502

computeCRIdata( subset(Fs.5nm,'F2') )   # returns a very large list, ending with CRI = 64.15195 
}


\keyword{light}
