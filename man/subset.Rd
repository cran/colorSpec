% File subset.Rd

\name{subset}
\title{extract a subset of a colorSpec Object}
\alias{subset.colorSpec}
\alias{subset}

\description{
  extract a subset of the spectra in a \bold{colorSpec} object.\cr  
  The subset can be specified by indexes, by a logical vector, or by a regular expression matching the
  \code{\link{specnames}}
}

\usage{
\S3method{subset}{colorSpec}( x, subset, ... )

}

\arguments{
  \item{x}{a \bold{colorSpec} object}
  
  \item{subset}{an integer vector, a logical vector, or a regular expression}
  
  \item{...}{additional arguments ignored}
    
}
    
\value{
  \code{subset(x)} returns a \bold{colorSpec} object with the same \code{\link{organization}} as \code{x}.
Exception: if the organization of \code{x} is \code{'vector'} and the subset is empty,
then the returned object is a matrix with 0 columns.
}    


\details{
If \code{subset} is an integer vector, each integer must be between 1 and M,
where M the number of spectra in \code{x}.
No duplicates are allowed.
The number of spectra returned is equal to \code{length(subset)}.
It is OK for the \code{length} to be 0, in which case the function returns the empty subset.

If \code{subset} is a logical vector, its length must be equal to M.
The number of spectra returned is equal to the number of \code{TRUE}s in \code{subset}.

If \code{subset} is a regular expression, the number of spectra returned is equal to
the number of \code{specnames(x)} matched by the expression. 
}

\note{
\code{subset} can also be used for re-ordering the spectra; 
just set argument \code{subset} to the desired permutation vector.
}


\seealso{
  \code{\link{organization}}
}

\examples{
tritanope = subset( lms2000.1nm, 1:2 )  # keep long and medium cone fundamentals, but drop the short

sml2000.1nm = subset( lms2000.1nm, 3:1 ) #   reorder from short to long
}

\keyword{colorSpec}
