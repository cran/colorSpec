% File materialSpectra.Rd

\name{materialSpectra}
\alias{materialSpectra}
\alias{neutralMaterial}
\alias{rectangularMaterial}
\alias{lensAbsorbance}

\title{compute standard material spectra}

\description{
Compute neutral gray material constant reflectance/transmittance,
and rectangular spectra.
Also compute absorbance of the human lens, as a function of age.
}

\usage{
neutralMaterial( gray=1, wavelength=380:780 )
rectangularMaterial( lambda, alpha=1, wavelength=380:780 )

lensAbsorbance( age=32, wavelength=400:700 )
}

\arguments{
\item{gray}{a numeric N-vector of gray levels, in the interval [0,1].  
    \code{gray=1} represents the \emph{Perfect Reflecting Diffuser}. }
    
\item{lambda}{a numeric Nx2 matrix with wavelength pairs in the rows,
    or a vector that can be converted to such a matrix, by row.
    The two wavelengths are the two \emph{transition wavelengths} of the returned spectrum,
    see \bold{Details}.}
    
\item{alpha}{a numeric N-vector of \emph{chromatic amplitudes} in the interval [-1,1].
    N must be equal to \code{nrow(lambda)}.
    \code{alpha} can also be a single number, 
    which is then replicated to length \code{nrow(lambda)}.
    The \emph{chromatic amplitude} is defined by \cite{Logvinenko}
    and controls the size of both transitions, see \bold{Details}.}
              
\item{age}{a numeric N-vector of ages in years; all ages must be \eqn{\ge} 20. }    
        
\item{wavelength}{a vector of wavelengths for the returned object}  
}
    
\value{
\code{neutralMaterial()} returns a
\bold{colorSpec} object with \code{\link{quantity}} equal to \code{'reflectance'}.
The reflectance of each spectrum is constant and taken from \code{gray}.
There are N spectra in the object - one for each gray level.

\code{rectangularMaterial()} returns a
\bold{colorSpec} object with \code{\link{quantity}} equal to \code{'transmitance'}.
The transmitance of each spectrum is a step function with 0, 1 or 2 transitions (jumps)
defined by the corresponding row in \code{lambda}.
If \code{rownames(lambda)} is not \code{NULL}, they are copied to \code{specnames}
of the output.
Otherwise the \code{specnames} are computed from the shape of the spectrum
using these acronyms:
LP (long-pass), SP (short-pass), BP (band-pass), BS (band-stop),
and N (neutral, in case \code{alpha==0}).

\code{lensAbsorbance()} returns a
\bold{colorSpec} object with \code{\link{quantity}} equal to \code{'absorbance'}.
The absorbance model for the human lens is taken from \cite{Pokorny}.
There are N spectra in the object - one for each age (N=\code{length(age)}).
}

\details{
A \emph{rectangular spectrum}, or \emph{rectangular metamer}, is easiest to define
when \eqn{\alpha=1} and \eqn{\lambda_1 < \lambda_2}.
In this case it is a band-pass filter with transmittance=1 for 
\eqn{\lambda \in [\lambda_1 , \lambda_2]} and transmittance=0 otherwise.
To create a long-pass filter, just set \eqn{\lambda_2} to \code{Inf},
or any large wavelength outside the spectrum range;
and similarly for a short-pass filter.\cr
When \eqn{0<\alpha<1} the spectrum is a weighted mixture of this band-pass filter
with a perfect neutral gray filter with transmittance=0.5 at all \eqn{\lambda},
using \eqn{\alpha} and \eqn{1-\alpha} as the two weights.
The minimum transmittance is \eqn{(1-\alpha)/2} and the maximum is \eqn{(1+\alpha)/2},
and their difference, the \emph{chromatic amplitude}, is \eqn{\alpha}.
It is still a band-pass filter.\cr
If \eqn{\alpha=0} the spectrum is a perfect neutral with transmittance=0.5.
\cr
To "flip" the spectrum to its complement (change band-pass to band-stop, etc.), 
change \eqn{\alpha} to a negative number, or swap \eqn{\lambda_1} and \eqn{\lambda_2}.
If \eqn{\lambda_1==\lambda_2} then the spectrum is undefined and a warning is issued
(unless \eqn{\alpha=0}).
}

\section{Logvinenko}{
It is clear that there are 3 degrees-of-freedom in the
spectra returned by \code{rectangularMaterial()}.
\emph{Logvinenko} shows that these spectra in fact form a 3D ball,
which he calls the \emph{rectangle color atlas}.
He also shows that if a material responder satisfies  the 2-transition condition,
then these spectra uniquely generate \emph{all} 
colors in the corresponding object color solid.
For more on this, see the vignette \href{../doc/inversion.pdf}{\bold{Estimating a Spectrum from its Response}}.
}


\section{Ostwald}{
Every spectrum returned by \code{rectangularMaterial()} is an Ostwald ideal spectrum.
In Ostwald's terminology, 
the \emph{color content} = \emph{chromatic amplitude} = \eqn{\alpha}.
And the \emph{black content} = \emph{white content} = \eqn{(1-\alpha)/2}.
Note that the sum of these 3 contents is 1.
However, Ostwald allows \emph{black content} and \emph{white content} to be unequal,
as long as the sum of the 3 contents is 1, and all are non-negative.
Thus there is one extra degree-of-freedom for Ostwald's ideal spectra,
for a total of 4 degrees-of-freedom.
If an additional argument (or arguments) were added to \code{rectangularMaterial()},
then it could return all Ostwald ideal spectra.
}

\references{
Foss, Carl E.  and Dorothy Nickerson and Walter C. Granville.
Analysis of the Ostwald Color System.
\emph{J. Opt. Soc. Am.}.
vol. 34.
no. 7.
pp. 361-381.
July, 1944.


Logvinenko, A. D.
An object-color space. 
\bold{Journal of Vision}.
9(11):5, 1-23, (2009).\cr
\code{https://jov.arvojournals.org/article.aspx?articleid=2203976}.
doi:10.1167/9.11.5.

Pokorny, Joel, Vivianne C. Smith, and Margaret Lutze.
Aging of the Human Lens.
\bold{Applied Optics}.
Vol. 26, No. 8.
15 April 1987.
Table I.  Page 1439.

}

\seealso{
\link{lightSpectra},
\code{\link{quantity}()},
\code{\link{specnames}()},
\code{\link{computeADL}()},
vignette \href{../doc/inversion.pdf}{\bold{Estimating a Spectrum from its Response}}
}

\examples{
#   make a perfect reflecting diffuser (PRD)
prd = neutralMaterial( 1 )

#   make a perfect transmitting filter (PTF)
ptf = prd
quantity(ptf) = 'transmittance'

#   make a band-stop filter (for interval [500,550])
#   with 1% transmittance in the band, and 99% outside the band
bs = rectangularMaterial( c(500,550), -0.98, 400:700 )
bs = rectangularMaterial( c(550,500),  0.98, 400:700 )  # equivalent to previous line

#   compare transmittance at 3 ages: 20, 32, and 80 years
plot( linearize(lensAbsorbance( c(20,32,80) )), col='black', lty=1:3 )
}

\keyword{materials}
