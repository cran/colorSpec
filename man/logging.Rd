\name{logging}
\alias{logging}

\title{Logging in colorSpec package}

\description{
Logging is done using the \pkg{logger} package.
Logging output goes to \code{stderr()}, just like the message stream;
but see \code{\link[base]{sink}()} (and the pitfalls of using it).\cr

}

\section{Logging Options}{

\describe{

\item{\code{colorSpec.stoponerror}}{If the this option is \code{TRUE} (the default),
a log event with level \code{ERROR} stops execution;
otherwise, execution keeps going.
For interactive use, \code{TRUE} is probably better.
For long batch jobs, \code{FALSE} might be appropriate,
since then a single error may not force a complete repeat.\cr
A \code{FATAL} event always stops execution.}
}
For examples on changing this option, see \code{\link{cs.options}}.
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
options( colorSpec.stoponerror=TRUE )

# or equivalently
cs.options( stop=TRUE )
}



\keyword{logging}

