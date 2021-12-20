\name{logging}
\alias{logging}

\title{Logging in colorSpec package}

\description{
There is some flexibility in the \bold{colorSpec} logging level and format.
Logging output goes to \code{stderr()}, just like the message stream;
but see \code{\link[base]{sink}()} (and the pitfalls of using it).\cr

}

\section{Logging Options}{

\describe{
\item{\code{colorSpec.loglevel}}{The levels are:
\code{"FATAL", "ERROR", "WARN", "INFO", "DEBUG",} and \code{"TRACE"} - the usual ones from \bold{Log4j}.
The initial level is \code{"WARN"}.
A \code{"FATAL"} event usually means an internal package error.
When setting \code{colorSpec.loglevel} an initial letter is sufficient.}


\item{\code{colorSpec.logformat}}{The format is given by a string with standard \bold{Log4j} \emph{conversion specifications}:
\tabular{lll}{
\code{\%t}  \tab    the date/time of the logging event. \tab \code{\%t} can be followed by standard \code{strftime} specs in braces; see example.\cr

\code{\%l}  \tab    the level of the logging event      \tab  \cr

\code{\%n}  \tab    namespace where event occurred      \tab  \cr

\code{\%f}  \tab    function where event occurred       \tab  \cr

\code{\%m}  \tab    the message itself                  \tab  \cr
}
}


\item{\code{colorSpec.stoponerror}}{If the this option is \code{TRUE} (the default),
a log event with level \code{"ERROR"} stops execution;
otherwise, execution keeps going.
For interactive use, \code{TRUE} is probably better.
For long batch jobs, \code{FALSE} might be appropriate,
since then a single error may not force a complete repeat.\cr
A \code{"FATAL"} event always stops execution.}
}
}


\references{
Wikipedia.
\bold{Log4j}.
\url{https://en.wikipedia.org/wiki/Log4j}
}

\seealso{
\code{\link{options}}, 
\code{\link{cs.options}}, 
\code{\link{sink}}, 
\code{\link{stderr}}
}


\examples{
options( colorSpec.logformat="\%t{\%H:\%M:\%OS3} \%l \%n::\%f(). \%m", colorSpec.stoponerror=TRUE )

# or equivalently
cs.options( logformat="\%t{\%H:\%M:\%OS3} \%l \%n::\%f(). \%m", stop=TRUE )
}



\keyword{logging}

