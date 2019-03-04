% File officialXYZ.Rd

\name{standardRGB}
\title{Convert from XYZ to some standard RGB spaces}
\alias{RGBfromXYZ}


\description{
To display an XYZ value, it typically must be converted to a standard RGB space.
This is the function to do it.

}

\usage{
RGBfromXYZ( XYZ, space )
}

\arguments{
  \item{XYZ}{a 3-vector, or a matrix with 3 columns with XYZs in the rows}
  \item{space}{the name of the RGB space - either \code{'sRGB'} or \code{'Adobe RGB'}. 
                    The match is case-insensitive, and spaces in the string are ignored.}
}

\details{
The input \code{XYZ} is multiplied by the appropriate 3x3 conversion matrix
(for \bold{sRGB} or \bold{Adobe RGB}).
These matrices are taken from \cite{Lindbloom} and not from the corresponding \cite{Wikipedia} articles; 
for the reason why see \bold{Note}.

}

\note{
An RGB space is normally defined by the xy chromaticities of the 3 primaries and the white point.
We follow \cite{Lindbloom} in using the 'official' XYZ of the white point from ASTM E308.
Using this XYZ of the white point makes the color space a little more consistent
with other areas of color.\cr
For example, from IEC 61966-2-1 we have
D65 xyY=(0.3127,0.3290,1) ->  XYZ=(0.9504559,1,1.0890578).
But from ASTM E308, D65 XYZ=(0.95047,1,1.08883), which is a little different.

}

\value{
An Mx3 matrix where M is the number of rows in \code{XYZ}, or M=1 if XYZ is a 3-vector.
Each row of the returned matrix is filled with linear RGB in the appropriate RGB space.
Values outside the unit cube are not clamped.
To compute non-linear display RGB, see \code{\link{DisplayRGBfromLinearRGB}()}. \cr
In case of error the function returns \code{NULL}.
  
}


\section{WARNING}{
This function is deprecated.
New software should use \code{spacesRGB::RGBfromXYZ()} instead.
The new function returns "signal RGB" instead of linear RGB.
}


\source{
IEC 61966-2-1:1999.
Multimedia systems and equipment - Colour measurement and management.
Part 2-1: Colour management - Default RGB colour space - sRGB.
\url{http://webstore.iec.ch/publication/6169}

Lindbloom, Bruce.
RGB/XYZ Matrices.
\url{http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html}

Wikipedia.
\bold{sRGB}.
\url{http://en.wikipedia.org/wiki/SRGB}

Wikipedia.
\bold{Adobe RGB}.
\url{http://en.wikipedia.org/wiki/Adobe_RGB_color_space}

}

\seealso{
  \code{\link{D65}},
  \code{\link{officialXYZ}},
  \code{\link{DisplayRGBfromLinearRGB}}
}

\examples{
RGBfromXYZ( officialXYZ('D65'), 'sRGB' )
#      R G B
# [1,] 1 1 1    # not really 1s, but difference < 1.e-7

RGBfromXYZ( c(.3127,0.3290,0.3583)/0.3290, 'sRGB' )
#              R        G       B
# [1,] 0.9998409 1.000023 1.00024    difference > 1.e-5
}

\keyword{RGB}

