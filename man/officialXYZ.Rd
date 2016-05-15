% File officialXYZ.Rd
\encoding{UTF-8}

\name{officialXYZ}
\title{Query the Official XYZ values for Standard Illuminants}
\alias{officialXYZ}


\description{
In careful calcuations, it is sometime helpful to have the 'official' values of XYZ,
i.e. with the right number of decimal places.

}

\usage{
officialXYZ( name )
}

\arguments{
  \item{name}{a subvector of
    \code{c('A','B','C','D50','D50.ICC','D55','D65','D75',}
    \code{'E','F2','F7','F11')}, which are the names of some standard illuminants}
}

\details{
All XYZ values are taken from the ASTM publication in \bold{References},
except \code{B} which is taken from \cite{Wyszecki & Stiles}
and \code{D50.ICC} which is taken from ICC publications.
The latter is different than that of ASTM.

}

\value{
  An Mx3 matrix where M is the length of \code{name}.  
  Each row filled with the official XYZ, but if the illuminant name is not recognized the row is all \code{NA}s.
  The matrix \code{rownames} are set to \code{name}, and \code{colnames} to \code{c('X','Y','Z')}.
}

\note{
The input names are case-sensitive.
The output XYZ is normalized so that \code{Y=1}.

}


\references{
ASTM E 308 - 01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
(2001).

Günther Wyszecki and W. S. Stiles.
Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition.
John Wiley & Sons, 1982.
Table I(3.3.8) p. 769.
}

\seealso{
  \code{\link{ABC}}, 
  \code{\link{D50}},  
  \code{\link{D65}},
  \code{\link{Fluorescents}},
  \code{\link{illuminantE}}
}

\examples{
officialXYZ( c('A','D50','D50.ICC','D65') ) 
#                 X Y         Z
# A       1.0985000 1 0.3558500
# D50     0.9642200 1 0.8252100
# D50.ICC 0.9642029 1 0.8249054
# D65     0.9504700 1 1.0888300
}

\keyword{eyes}

