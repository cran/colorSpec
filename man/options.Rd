\name{cs.options}
\alias{cs.options}
\alias{cs.getOptions}
\alias{cs.setOptions}
\title{Functions to set and retrieve colorSpec package options}

\description{
Allow the user to set and examine a variety of options for \bold{colorSpec}.

}


\usage{
cs.options( ... )
}

\arguments{
\item{...}{
named arguments are set; 
unnamed arguments are retrieved, and must be single character strings.
See example.}
}

\value{
returns a list with all options mentioned.\cr
If no arguments are given, then it returns a list of all options.
}

\section{Options used in logging}{

\tabular{lll}{
\bold{Name}         \tab \bold{Default Value}                   \tab \bold{Description} \cr
\code{loglevel}     \tab \code{'WARN'}                          \tab controls the amount of logging information produced. \cr
. \tab . \tab . \cr
\code{logformat}    \tab \code{'\%t \%l \%n::\%f(). \%m'}  \tab the format of each log line. See \code{\link{logging}} page.\cr
. \tab . \tab . \cr
\code{stoponerror}  \tab \code{TRUE}                            \tab stop when a logging event has level ERROR; a logical. \cr
}
When setting \code{loglevel} an initial letter is sufficient.
}


\examples{
#   set 2 options and retrieve the 3rd
cs.options( loglevel="DEBUG", stoponerror=FALSE, "logformat" )
}


\seealso{
\link{logging}
}


\keyword{options}

