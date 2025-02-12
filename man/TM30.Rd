% File TM30.Rd

\name{ANSI/IES TM-30}

\alias{computeTM30}
\alias{computeTM30.colorSpec}
\alias{plot.TM30}
\alias{print.TM30}

\title{Compute ANSI/IES TM-30 data, and Create Color Rendition Reports}

\description{
Compute TM-30 data from a \bold{colorSpec} object with the spectrum of a light source
(\code{type()} is \code{'light'}).
And plot full reports from this data.
The full TM-30 standard is in the \bold{References}.

These functions were guided by the
\bold{IES Spectral Calculator},
and checked against the downloadable .XLSX spreadsheets and .PDF reports available there.
}

\usage{
\S3method{computeTM30}{colorSpec}( x, reference=NULL )

\S3method{plot}{TM30}( x, omi=c(0.75,0.2,1.25,0.25), ... )

\S3method{print}{TM30}( x, ... )
}


\arguments{
\item{x}{
For \code{computeTM30()}, a colorSpec object with a single spectrum whose \code{type} is \code{'light'}.
\cr
For \code{plot()} and \code{print()}, an object returned by \code{computeTM30()},
with class \code{'TM30'}.
}

\item{reference}{
a \pkg{colorSpec} object used as the reference light source for \code{x}.
If it is \code{NULL}, this reference is computed from the CCT of \code{x},
according to the instructions in \bold{ANSI/IES TM-30-20}.
}

\item{omi}{a vector of 4 numbers, defining the outer margins of the plot, in inches;
for details see \code{\link[graphics]{par}()}.
For the TM-30 report, text is written in the top and bottom outer margins
using \code{\link[graphics]{title}()}, so these margins may require some fine tuning.
}

\item{...}{further arguments ignored}
}

\details{
\code{computeTM30()} follows the recipes in \bold{ANSI/IES TM-30-20}.

\code{plot()} reproduces the full reports at the \bold{IES Spectral Calculator}.

\code{print} prints a short 9-line summary of the important TM-30-20 measures.

}


\value{
\code{computeTM30()} returns a list of data with class \code{'TM30'}.
The data is designed to be a reflection of the data in the .XLSX spreadsheets created
by the \bold{IES Spectral Calculator}.

\code{plot()} and \code{print()} return \code{TRUE} or \code{FALSE}, invisibly.
}

\references{
\bold{ANSI/IES TM-30-20}.
\bold{Technical Memorandum: IES Method for Evaluating Light Source Color Rendition}.
\url{https://store.ies.org/product/technical-memorandum-ies-method-for-evaluating-light-source-color-rendition/}

\bold{IES Spectral Calculator}.
\url{https://www.ies.org/standards/standards-toolbox/tm-30-spectral-calculator/}
}

\examples{
\dontrun{
F1 = subset( Fs.5nm, 1 )

pdf( "report.pdf", width=8.5, height=11 )
plot( computeTM30(F1) )
dev.off()

svg( "report.svg", width=8.5, height=11 )
plot( computeTM30(F1) )
dev.off()
}
}

\seealso{
\code{\link[graphics]{par}()},
\code{\link[graphics]{title}()}
}

\keyword{colorSpec}

