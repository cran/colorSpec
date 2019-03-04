% File readSpectra.Rd

\name{readSpectra}
\title{read colorSpec objects from files}
\alias{readSpectra}
\alias{readAllSpectra}
\alias{readSpectraXYY}
\alias{readSpectraSpreadsheet}
\alias{readSpectrumScope}
% \alias{readSpectraExcel}
\alias{readSpectraControl}
\alias{readSpectraCGATS}

\description{
These functions read \bold{colorSpec} objects from files.
In case of \code{ERROR}, they return \code{NULL}.
There are 5 different file formats supported; see \bold{Details}.
}

\usage{
readSpectra( pathvec, ... )

readSpectraXYY( path )
readSpectraSpreadsheet( path )
readSpectrumScope( path )
% readSpectraExcel( path, worksheet=NULL )      # requires package xlsx
readSpectraCGATS( path )
readSpectraControl( path )

}


\arguments{
  \item{pathvec}{a character vector to (possibly) multiple files.
    The file extension and a few lines from each file are read
    and a guess is made regarding the file format.}

  \item{...}{optional arguments passed on to \code{\link{resample}()}.
    The most important is \code{wavelength}.
    If these are missing then \code{\link{resample}()} is not called.}

  \item{path}{a path to a single file with the corresponding format:
    \code{XYY}, \code{Spreadsheet}, \code{Scope}, 
    \code{CGATS}, or \code{Control}.  See Details.
    If the function cannot recognize the format, it returns NULL.}

}

    
\value{
\code{readSpectra()} returns a single \bold{colorSpec} object or \code{NULL} in case of ERROR.
If there are multiple files in \code{pathvec} and they cannot be combined using \code{bind()} because their wavelengths are different, it is an ERROR.
To avoid this ERROR, the \code{wavelength} argument can be used for
resampling to common wavelengths.
If there are multiple files, the \code{\link{organization}} of the returned object is
\code{'df.row'} and the first column is the \code{path}
from which the spectrum was read.
  
The functions \code{readSpectraXYY()}, \code{readSpectraSpreadsheet()}, 
and \code{readSpectraScope()}, return a single \bold{colorSpec} object,
or \code{NULL} in case of ERROR.

The functions \code{readSpectraCGATS()} and \code{readSpectraControl()}
are more complicated.
These 2 file formats can contain multiple spectra with different \code{wavelength} sequences,
so both functions return a \emph{list} of \bold{colorSpec} objects,
even when that list has length 1.
If no spectral objects are found, they return \code{NULL}.

If \code{readSpectra()} calls \code{readSpectraCGATS()} or \code{readSpectraControl()}
and receives a list of \bold{colorSpec} objects,
\code{readSpectra()} attempts to \code{\link{bind}()} them into a single object.
If they all have the same wavelength vector, then the \code{bind()} succeeds and the
single \bold{colorSpec} object is returned.
Otherwise the \code{bind()} fails, and it is an ERROR.
To avoid this error \code{readSpectra()} can be called with a
\code{wavelength} argument. 
The multiple spectra are resampled using \code{\link{resample}()} and
\emph{then} combined using \code{bind()},
which makes it much more convenient to read such files.
}



\details{
\code{readSpectra()} reads the first few lines of the file in order
to determine the format,
and then calls the corresponding format-specific function.
If \code{readSpectra()} cannot determine the format, it returns NULL.
The 5 file formats are:
  
\code{XYY} \cr         
There is a column header line matching \code{^(wave|wv?l)} (not case sensitive)
followed by the the names of the spectra. 
All lines above this one are taken to be \code{metadata}.
The separarator on this header line can be space, tab, or comma; 
the line is examined and the separator found there is used in the lines with data below. 
The organization of the returned object is \code{'df.col'}.  
This is probably the most common file format;
see the sample file \code{ciexyz31_1.csv}.\cr

\code{Spreadsheet}  \cr
There is a line matching \code{"^(ID|SAMPLE|Time)"}. This line and lines below must be tab-separated.  Fields matching \code{'^[A-Z]+([0-9.]+)nm$'} are taken to be spectral data and other fields are taken to be \code{extradata}.  All lines above this one are taken to be \code{metadata}.
The organization of the returned object is \code{'df.row'}.
This is a good format for automated acquisition, using a spectrometer, of many spectra.
See the sample file \code{N130501.txt} from \bold{Wolf Faust}.\cr
  
\code{Scope}  \cr
This is a file format used by \bold{Ocean Optics} spectrometer software.  
There is a line \cr
\code{>>>>>Begin Processed Spectral Data<<<<<}\cr
followed by wavelength and energy separated by a tab.  
There is only 1 spectrum per file.  
The \code{organization} of the returned object is \code{'vector'}.
See the sample file \code{pos1-20x.scope}. \cr

\code{CGATS}  \cr
This is a complex format that is best understood by looking at some samples, such as\cr
\code{extdata/objects/Rosco.txt};
see also the \bold{References}.
The function \code{\link{readCGATS}()} is first called to get all the tables,
and then for each table the column names are examined.
There are 2 conventions for presenting the spectral data:
\itemize{  
\item In the standard convention the fields \code{SPECTRAL_DEC} or \code{SPECTRAL_PCT}
have the spectral values.
The former is the true value, and the latter is the true value x 100.
Each value column is preceded a corresponding wavelength column,
which has the field name \code{SPECTRAL_NM}.
Note that these field names are highly duplicated.
In principle, this convention allows each record in a CGATS table to have
a different wavelength vector.
However, this complication is rejected by \code{readSpectraCGATS()},
which treats it as an ERROR.
\item In the non-standard convention the field names that match the pattern\cr
\code{"^(nm|SPEC_|SPECTRAL_)[_A-Z]*([0-9.]+)$"}
are considered to be spectral value data, and other fields are considered \code{extradata}.
The wavelength is the numerical value of the 2nd parenthesized expression \code{([0-9.]+)}
in nanometers.
Note that every record in this CGATS table has the same wavelength vector.
Although this convention is non-standard, it appears in files from many companies,
including X-Rite.
}
If a \code{data.frame} has spectral data, it is converted to a \bold{colorSpec} object
and placed in the returned list.
The organization of the resulting \bold{colorSpec} object is \code{'df.row'}.
If the \code{data.frame} of \code{extradata} contains a column 
\code{SAMPLE_NAME}, \code{SAMPLE_ID}, \code{SampleID}, or \code{Name}, 
(examined in that order), then that column is taken to be the \code{specnames}
of the object.
If a table has no spectral data, then it is ignored.
If the CGATS file has no tables with spectral data, then it is an ERROR
and the function returns \code{NULL}.
  
\code{Control} \cr
This is a personal format used for digitizing images of plots from manufacturer
datasheets and academic papers. It is structured like a \code{.INI} file.  
There is a \code{[Control]}
section establishing a simple linear map from pixels to the wavelength and spectrum quantities.
Only 3 points are really necessary.
It is OK for there to be a little rotation of the plot axes relative to the image.
This is followed by a section for each spectrum, in XY pixel units only.
Conversion to wavelength and spectral quantities is done on-the-fly after they are read.
Extrapolation can be a problem, especially when the value is near 0.
To force constant extrapolation (see \code{\link{resample}()}),
repeat the control point (knot) at the endpoint.
See the sample file \code{Lumencor-SpectraX.txt}.
The organization of the returned objects is \code{'vector'}.\cr
}

\note{
During import, each read function tries to guess the \code{quantity} 
from spectrum names or other cues.  
For example the first line in \bold{N130501.txt} is \code{IT8.7/1},
which indicates that the quantity is \code{'transmittance'}
(a reflective target is denoted by \code{IT8.7/2}).  
If a confident guess cannot be made, it makes a wild guess and issues a warning.  
If the quantity is incorrect, one can assign the correct value after import.  
Alternatively one can add a line to the header part of the file with the 
keyword \code{'quantity'} followed by some white-space and the correct value.  
It is OK to put the value in quotes.
See example files under folder \bold{extdata}.
}

\references{
CGATS.17 Text File Format.
\url{http://www.colorwiki.com/wiki/CGATS.17_Text_File_Format}.

ANSI/CGATS.17. 
Graphic technology - Exchange format for colour and process control data using XML or ASCII text.
\url{https://webstore.ansi.org/}
2009.

ISO/28178. Graphic technology - Exchange format for colour and process control data using XML or ASCII text.
\url{https://www.iso.org/standard/44527.html}.
2009.
\cr
}


\seealso{
\code{\link{wavelength}},
\code{\link{quantity}},
\code{\link{metadata}},
\code{\link{resample}},
\code{\link{bind}},  
\code{\link{readCGATS}}  
}

\examples{
#   read file with header declaring the quantity to be "photons->neural"
bird = readSpectra( system.file( "extdata/eyes/BirdEyes.txt", package='colorSpec' ) )
quantity(bird)   # [1] "photons->neural"
}

\keyword{colorSpec}
