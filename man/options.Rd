\name{cs.options}
\alias{cs.options}
\title{Functions to set and retrieve colorSpec package options}

\description{
\bold{colorSpec} has one option.
The option is stored in the \bold{R} global list and is:

\code{colorSpec.stoponerror}

For details on what it does see \code{\link{logging}}.

It can be set using the built-in function \code{base::\link{options}()}.
When \R starts up,
an option can be set using a call to \code{base::options()} in the file \bold{Rprofile.site}.
If \bold{colorSpec} is later loaded, the value of the option will not be changed.
If an option has not been assigned, then it is created with a default value.

The little helper function \code{cs.options()} makes setting the options a little easier
in a few ways:
\itemize{  
\item it automatically prepends the string \code{'colorSpec.'}
\item partial matching of the option name is enabled
\item an error is issued when the option value has the wrong type
}
}


\usage{
cs.options( ... )
}

\arguments{
\item{...}{
named arguments are set; 
unnamed arguments are ignored with a warning.
See \bold{Examples}.}
}

\value{
returns a list with all the \bold{colorSpec} options.
}


\seealso{
\code{\link{logging}},
\code{\link{options}}
}



\examples{
cs.options( stop=TRUE )             # no problems, no warning
cs.options( stop='TRUE' )           # warns that value has the wrong type
cs.options( stop=FALSE, "DEBUG" )   # warns that the 2nd argument has no name
}


\keyword{options}

