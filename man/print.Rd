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

\item{long}{logical indicating whether to print \code{metadata}, \code{calibrate},
\code{sequence}, and \code{emulate} data}

\item{...}{further arguments ignored}

}
\value{
Both functions return (invisibly) the character vector that was just printed
to \code{stdout()}.
}

\details{
If \code{long=FALSE},
\code{summary()} prints a summary of the wavelength vector, and names of all spectra.
For each spectrum it prints the range of values, \code{LambdaMax}, and
\code{extradata} if any.
If \code{long=TRUE} it also prints data listed above (if any).\cr
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

