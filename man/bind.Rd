% File bind.Rd

\name{bind}
\title{Combine colorSpec Objects}
\alias{bind}
\alias{bind.colorSpec}

\description{
  Take a sequence of \bold{colorSpec} objects and combine their spectra
}
\usage{
\S3method{bind}{colorSpec}( \ldots )

}

\arguments{
  \item{...}{\bold{colorSpec} objects with the same \code{\link{wavelength}} and \code{\link{quantity}}, and with distinct \code{\link{specnames}} (no duplicates)}
    
}
    
\value{
\code{bind()} returns a \bold{colorSpec} object, or \code{NULL} in case of ERROR.  If the bind is
successful, the number of spectra in the output object is the sum of the number of spectra
in the input objects.
           
}    


\details{
The \code{\link{organization}} of the returned object is the most complex of those in the inputs,
where the order of complexity is:

\code{'matrix'}  <  \code{'df.col'}  <  \code{'df.row'}

If the selected \code{organization} is \code{'df.row'} the \code{extradata} is combined in a way that preserves all the columns.
Missing data is filled with \code{NA}s, analogous to
\code{plyr::rbind.fill()}.

The \code{\link{metadata}} of the returned object is copied from the first object in the input list.
  
}


\seealso{
\code{\link{wavelength}},
\code{\link{quantity}},
\code{\link{specnames}},    
\code{\link{organization}},  
\code{\link{extradata}},    
\code{\link{metadata}},
\code{\link[plyr]{rbind.fill}}  
}

\examples{
Rosco = readSpectra( system.file( 'extdata/objects/Rosco.txt', package='colorSpec' ) )
Rosco = resample( Rosco, wave=wavelength(Hoya) )
numSpectra(Hoya)        # prints 4
numSpectra(Rosco)       # prints 42

filters = bind( Hoya, Rosco )
numSpectra(filters)     # prints 46

colnames( extradata(Hoya) )
## [1] "SAMPLE_NAME"  "FresnelReflectance"  "Thickness"

colnames( extradata(Rosco) )
## [1] "Model"  "SampleID"  "SAMPLE_NAME"  "Substrate"  "RefractiveIndex"  "Thickness"

##  The columns in common are "SAMPLE_NAME" and "Thickness"


colnames( extradata(filters) )
## [1] "FresnelReflectance" "Model" "RefractiveIndex" "SAMPLE_NAME"
## [5] "SampleID" "Substrate" "Thickness"
##
## "SAMPLE_NAME" and "Thickness" are combined in the usual way
## The other columns are present, and missing data is filled with NAs
}

\keyword{colorSpec}
