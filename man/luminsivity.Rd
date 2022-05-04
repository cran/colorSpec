\encoding{UTF-8}
\name{luminsivity}
\alias{luminsivity.1nm}
\alias{luminsivity}
\docType{data}
\title{Luminous Efficiency Functions (photopic and scotopic)}

\description{
\tabular{ll}{
\code{luminsivity.1nm} \tab    Four luminous efficiency functions, from 360 to 830 nm, at 1nm step
}
}

\format{
A \bold{colorSpec} object, with quantity \code{'energy->neural'}, and with 4 spectra:

\describe{
\item{\code{photopic1924}}{The luminous efficiency function adopted by the CIE in 1924,
and defining the \emph{standard photopic observer}.
It is only to be used when light levels are high enough
that the sensitivity of the eye is mediated by cones, and not rods.
It is the same as the \emph{y-bar} function in \code{xyz1931.1nm}.
It is used to define the \emph{candela} in the \emph{International System} (SI)
and is the only one of these functions to appear in the SI.
It was downloaded from
\url{http://www.cvrl.org/database/data/lum/vl1924e_1.csv}
where it is defined from 360 to 830 nm.
}

\item{\code{scotopic1951}}{The luminous efficiency function adopted by the CIE in 1951,
and defining the \emph{standard scotopic observer}.
It is only to be used when light levels are low enough to exclude the activation of cones.
It has no effective role in colorimetry.
It was downloaded from
\url{http://www.cvrl.org/database/data/lum/scvle_1.csv}
where it is defined from 380 to 780 nm.
It has been padded with 0s to 360 to 830 nm.
}

\item{\code{photopic1978}}{The luminous efficiency function for photopic vision,
with adjustments in the blue region by Judd (1951) and Vos (1978).
It was published by the CIE in 1988.
It was downloaded from
\url{http://www.cvrl.org/database/data/lum/vme_1.csv}
where it is defined from 380 to 780 nm.
It has been padded with 0s to 360 to 830 nm.
}

\item{\code{photopic2008}}{The CIE (2008) \emph{physiologically-relevant} luminous efficiency function for photopic vision,
by Stockman, Jagle, Pirzer, & Sharpe.
It was downloaded from
\url{http://www.cvrl.org/database/data/lum/linCIE2008v2e_1.csv}
where it is defined from 390 to 830 nm.
It has been padded with 0s to 360 to 830 nm.
}
}
}


\note{
Luminsivity is a self-coined \emph{portmanteau word}:
\code{luminsivity} = \code{luminous} * \code{responsivity}.
The word is unrelated to \emph{emissivity}.
The term \emph{luminous responsivity} is not common, but appears on page 15 of \emph{Grum}.
The term \emph{luminous efficiency function} is standard, but too long.
The term \emph{luminosity function} is common, but \emph{luminosity} is ambiguous and also appears in astronomy and scattering theory.\cr
The object \code{luminsivity.1nm} is used by the function \code{\link{photometric}()}.
}


\source{
Colour & Vision Research Laboratory.
Institute of Opthalmology.
University College London. UK.
\url{http://www.cvrl.org/}
}

\references{
Grum, Franc and Richard J. Becherer.
\bold{Radiometry}.
Optical Radiation Measurements, Volume 1.
Academic Press.
1979.

Stockman, A., Jagle, H., Pirzer, M., & Sharpe, L. T. (2008).
The dependence of luminous efficiency on chromatic adaptation.
Journal of Vision, 8, 16:1, 1-26.
}

\seealso{
\code{\link{xyz1931.1nm}},
\code{\link{photometric}} 
}

\examples{
summary(luminsivity.1nm)
product( D65.1nm, luminsivity.1nm, wave='auto' )
}
\keyword{datasets}
\keyword{eyes}
