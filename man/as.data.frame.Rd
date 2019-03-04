\name{as.data.frame}
\alias{as.data.frame.colorSpec}

\title{Convert a colorSpec Object to a data.frame}
\description{
convert a \bold{colorSpec} object to a \code{data.frame}
}

\usage{
\S3method{as.data.frame}{colorSpec}( x, row.names=NULL, optional=FALSE, organization='auto', ... )
}

\arguments{
  \item{x}{a \bold{colorSpec} object}
  \item{organization}{The organization of the returned \code{data.frame},
  which can be \code{'row'}, \code{'col'}, or \code{'auto'}.
  If \code{'auto'}, then \code{'row'} or \code{'col'} is selected automatically,
  see \bold{Details}}
  \item{row.names}{ignored}
  \item{optional}{ignored}
  \item{...}{extra arguments ignored}  
}

\details{
If \code{organization} is \code{'auto'},
and the organization of \code{x} is \code{'df.row'},
then \code{organization} is set to \code{'row'}
and the returned \code{data.frame} has the spectra in the rows.
Otherwise the returned \code{data.frame} has the spectra in the columns.
}

\value{
If the returned \code{data.frame} has the spectra in the rows,
then the spectra are in a matrix in the last column (with name \code{spectra}),
and any existing \code{extradata} are also returned in the other columns.
The wavelengths are only present in character form,
as the \code{colnames} of the matrix.\cr
If the returned \code{data.frame} has the spectra in the columns,
then the wavelengths are in the first column, and the spectra are in the other columns.
}


\seealso{
\code{\link{as.matrix}},
\code{\link{extradata}}
}


\keyword{colorSpec}
