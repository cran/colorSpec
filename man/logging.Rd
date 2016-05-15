\name{logging}
\alias{logging}

\title{Logging in colorSpec package}

\description{
There is some flexibility in the \bold{colorSpec} logging level and format.
There is no flexibility in the output, it all goes to the console using
\code{cat} and \code{print}, but see \code{\link[base]{sink}}.

}

\details{
\code{loglevel}\cr
The levels are: "FATAL", "ERROR", "WARN", "INFO", "DEBUG", and "TRACE" - the usual ones from \bold{Log4J}.
The initial level is "WARN".\cr
A log event with level "ERROR" stops execution if the option \code{stoponerror} is \code{TRUE}.\cr
A "FATAL" log event (e.g. internal error), always stops execution.

\code{logformat}\cr
The format is given by a string with standard \bold{Log4J} \emph{conversion specifications}:

\tabular{lll}{
\code{\%t}  \tab    the date/time of the logging event. \tab \code{\%t} can be followed by standard \code{strftime} specs in braces; see example.\cr

\code{\%l}  \tab    the level of the logging event      \tab  \cr

\code{\%n}  \tab    namespace where event occurred      \tab  \cr

\code{\%f}  \tab    function where event occurred       \tab  \cr

\code{\%m}  \tab    the message itself                  \tab  \cr

}
}



\examples{
cs.options( logformat="\%t{\%H:\%M:\%OS3} \%l \%n::\%f(). \%m" )
}

\seealso{
\link{cs.options}, \link[base]{sink}
}


\keyword{logging}

