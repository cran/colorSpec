% File print.Rd

\name{print}
\alias{summary.colorSpec}
\alias{print.colorSpec}

\title{Convert colorSpec object to readable text}

\description{
display a \bold{colorSpec} object as readable text.
Output goes to \code{stdout()}.
}

\usage{
\S3method{print}{colorSpec}( x, ...)

\S3method{summary}{colorSpec}( object, long=TRUE, ... )

}


\arguments{
\item{x}{a colorSpec object}

\item{object}{a colorSpec object}

\item{long}{logical indicating whether to print metadata, calibration, and product terms}

\item{...}{further arguments ignored}

}
\value{
Both functions return (invisibly) the character vector that was just printed.
}

\details{
\code{print()} prints a summary of the wavelength vector, and names of all spectra.
For each spectrum it prints the range of values, \code{LambdaMax}, and
\code{extradata} if any.\cr
If \code{long} is \code{FALSE}, \code{summary()} prints the same as \code{print()},
and if \code{long} is \code{TRUE} it also prints metadata and data about product terms (if any).
The function \code{print()} simply calls \code{summary()} with \code{long=FALSE}.

}

\examples{
print( xyz1931.1nm )

xyz1931.1nm     # same thing, just calls print()
}

\seealso{
\code{\link{extradata}},
\code{\link[base]{print}},
\code{\link[base]{summary}},
\code{\link{stdout}}
}

\keyword{colorSpec}

