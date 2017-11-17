% File readCGATS.Rd

\name{readCGATS}
\title{read tables from files in ANSI/CGATS.17 format}
\alias{readCGATS}

\description{
This text format supports a preamble followed by N tables, where N \eqn{\ge} 1.
Each table can have a separate header.
A table may or may not contain spectral data, see \bold{Note}.
The function converts each table to a data.frame with attributes; see Details.
}

\usage{
readCGATS( path, collapsesingle=FALSE )

}


\arguments{
  \item{path}{the path name of a single file, in CGATS format}

  \item{collapsesingle}{If \code{path} has only one table (N=1) and \code{collapsesingle} is \code{TRUE},
  then return the single data.frame (instead of a list with 1 data.frame).
  If \code{path} has multiple tables (N \eqn{\ge} 2), then \code{collapsesingle} is ignored.
  }
}

    
\value{
  \code{readCGATS} returns a list of data.frames - 
  one data.frame for each table found in \code{path}.
  The list and each individual data.frame have attributes, see \bold{Details}.

  If \code{path} has only a single table (the majority of files have only 1) 
  and \code{collapsesingle} is \code{TRUE}, then the attributes of the list
  are copied to those of the data.frame, and the data.frame is then returned.
  }



\details{
  The returned list is given attributes: \code{"path"},  \code{"preamble"}, 
  and (if present) \code{"date"}, \code{"created"}, \code{"originator"}, and \code{"file_descriptor"}.
  The attribute values are all character vectors.
  The value of attribute \code{"path"} is the argument \code{path}.
  The length of \code{"preamble"} is (typically) greater than 1, and the others have length 1.
  Each line of the preamble is a keyword-value pair.
  The keyword \code{ORIGINATOR} is converted to attribute \code{"originator"}.
  The keyword \code{FILE_DESCRIPTOR} is converted to attribute \code{"file_descriptor"}.
  The keyword \code{CREATED} is converted to attributes  \code{"created"} and \code{"date"}.
  
  Each data.frame in the list is assigned attributes: \code{"header"}, 
  and (if present) \code{"descriptor"}.
  The length of \code{"header"} is (typically) greater than 1, and \code{"descriptor"} has length 1.
  Each line of the table header is a keyword-value pair.
  The keywords \code{DESCRIPTOR} and  \code{TABLE_DESCRIPTOR}  
  are converted to attribute \code{"descriptor"}.
  
  The CGATS standard allows \emph{both} spaces and tabs as whitespace.
  This function is stricter; when splitting lines,
  it allows both characters, but only one at a time.
  The line after \code{BEGIN_DATA_FORMAT} is examined.
  If a tab is present then the tab is used for splitting both the
  \code{BEGIN_DATA_FORMAT} line and the \code{BEGIN_DATA} block.
  Othewise the space character is used.
}

\note{
No attempt is made to recognize those tables that contain spectral data.
For conversion of spectral data to \bold{colorSpec} objects, see \code{\link{readSpectraCGATS}}.
}

\references{

ANSI/CGATS.17. 
Graphic technology - Exchange format for colour and process control data using XML or ASCII text.
https://webstore.ansi.org/
2009.

ISO/28178. Graphic technology - Exchange format for colour and process control data using XML or ASCII text.
https://www.iso.org/standard/44527.html.
2009.

CGATS.17 Text File Format.   http://www.colorwiki.com/wiki/CGATS.17_Text_File_Format.
}


\seealso{
\code{\link{readSpectraCGATS}}
}

\examples{
#   read file with 2 tables of non-spectral data
A70 = readSpectra( system.file( "tests/A70.ti3", package='colorSpec' ) )
length(A70)         # [1] 2   # the file has 2 tables
ncol( A70[[1]] )    # [1] 7   # the 1st table has 7 columns
ncol( A70[[2]] )    # [1] 4   # the 2nd table has 4 columns
}

\keyword{colorSpec}
