% File colorSpec-package.Rd
\docType{package}
\name{colorSpec-package}
\alias{colorSpec-package}
\title{Package colorSpec - Color Calculations with Emphasis on Spectral Data}

\description{
Package \pkg{colorSpec} is for working with color spectral data in \R.

}

\details{
\bold{Features}:
\enumerate{
\item a clear classification of the commonly seen spectra into 4 types - depending on the vector space to which they belong
\item  flexible organization for the spectra in memory, using an S3 class - \bold{colorSpec}
\item  a product algebra for the \bold{colorSpec} objects
\item  uniform handling of biological eyes, electronic cameras, and general action spectra
\item  a few advanced calculations, such as computing optimal colors  (aka Macadam Limits)
\item  built-in essential tables, such as the CIE illuminants and color matching functions
\item  a package logging system with log levels taken from the popular \bold{Log4J}
\item  support for reading a few spectrum file types, including CGATS
\item  bonus files containing some other interesting spectra
\item  minimal dependencies on other \R packages

}


\bold{Non-features}:
\enumerate{
\item there is no support for many common 3D color spaces, such as CIELAB, HSL, etc..  
For these spaces see packages \pkg{colorspace}  and  \pkg{colorscience}.
\item  there are few non-linear operations 
\item  there is little support for scientific units; for these see package \pkg{colorscience}
\item photons are parameterized by wavelength in nm; 
other wavelength units and alternative parameterizations,
such as wavenumber and electronvolt, are not available
}

Regarding the non-linear operations in 2, 
the only such operations are conversion of linear
RGB to display RGB, conversion of absorbance to transmittance,
and the reparameterized wavelength in \code{\link{computeADL}}.
The electronic camera model is purely linear
with no dark current offset or other deviations.

Many ideas are taken from packages
\pkg{hyperSpec}, \pkg{hsdar}, \pkg{pavo}, and \pkg{zoo}.
}


\author{
Glenn Davis <gdavis@gluonics.com>
}

\seealso{

\code{\link{colorSpec}} for the S3 class provided by this package.

\bold{colorSpec-guide.pdf} in the doc directory.  
Go to the main index page and look for 
\bold{User guides, package vignettes and other documentation.}
}

\keyword{package}

