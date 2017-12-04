% File readCGATS.Rd

\name{readCGATS}
\title{read tables from files in ANSI/CGATS.17 format}
\alias{readCGATS}

\description{
The CGATS text format supports a preamble followed by N tables, where N \eqn{\ge} 1.
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
  The name of the table is lost.
  
  If there is an error in any table, then the function returns \code{NULL}.
  }



\details{
  The returned list is given attributes: \code{"path"},  \code{"preamble"}, 
  and (if present) \code{"date"}, \code{"created"}, \code{"originator"}, and \code{"file_descriptor"}.
  The attribute values are all character vectors.
  The value of attribute \code{"path"} is the argument \code{path},
  and the other values are extracted from \code{"preamble"}.
  The length of \code{"preamble"} is (typically) greater than 1, and the others have length 1.
  Each line of the preamble is a keyword-value pair.
  The keyword \code{ORIGINATOR} is converted to attribute \code{"originator"}.
  The keyword \code{FILE_DESCRIPTOR} is converted to attribute \code{"file_descriptor"}.
  The keyword \code{CREATED} is converted to attributes  \code{"created"} and \code{"date"}.
  The list is also given \code{names}.
  If the keyword \code{TABLE_NAME} is present in the table header, then its value is used.
  Otherwise the names are \code{"TABLE_1"}, \code{"TABLE_2"}, ...
  
  Each data.frame in the list is assigned attributes: \code{"header"}, 
  and (if present) \code{"descriptor"}.
  The length of \code{"header"} is (typically) greater than 1, and \code{"descriptor"} has length 1.
  Each line of the table header is a keyword-value pair.
  The keywords \code{DESCRIPTOR} and  \code{TABLE_DESCRIPTOR}  
  are converted to attribute \code{"descriptor"}.
  
  For the lines between \code{BEGIN_DATA} and \code{END_DATA},
  two conventions for separating the values are supported:
  \itemize{  
  \item  In the standard convention, fields are separated by contiguous spaces or tabs, and character strings (which may have embedded spaces or even tabs) 
  are enclosed by double-quotes.
  This is is the convention in the CGATS standard.
  The function \code{\link{scan}()} is used here.
  \item In the non-standard convention, fields are separated by a \emph{single} tab,
  and character strings (which may have embedded spaces but not tabs) are \emph{not} enclosed by double-quotes.
  This convention is often easier to work with in spreadsheet software.
  The function \code{\link{strsplit}()} is used here.
  }
  The function \code{readCGATS()} selects the separation convention 
  by examining the line after \code{BEGIN_DATA_FORMAT}.
  If this line is split by a single tab and the number of fields
  matches that given on the \code{NUMBER_OF_FIELDS} line,
  then the non-standard convention is selected;
  otherwise, the standard convention is selected.
}

\note{
In the \code{BEGIN_DATA_FORMAT} line(s),
field names may not be quoted and may not have embedded spaces.
\cr
The CGATS standard allows duplicated field names,
and \code{readCGATS()} returns them as they appear, with no attempt to append numbers
in order to make them unique.
Examples of field names which may be duplicated are:
\code{SPECTRAL_NM}, \code{SPECTRAL_DEC}, and \code{SPECTRAL_PCT};
for more on these see \code{\link{readSpectraCGATS}()}.
\cr
No attempt is made to recognize those tables that contain spectral data.
For conversion of spectral data to \bold{colorSpec} objects, see \code{\link{readSpectraCGATS}()}.
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
\code{\link{readSpectraCGATS}},
\code{\link{scan}},
\code{\link{strsplit}},
\code{\link{names}}
}

\examples{
#   read file with 2 tables of non-spectral data
A70 = readSpectra( system.file( "tests/A70.ti3", package='colorSpec' ) )
length(A70)         # [1] 2   # the file has 2 tables
ncol( A70[[1]] )    # [1] 7   # the 1st table has 7 columns
ncol( A70[[2]] )    # [1] 4   # the 2nd table has 4 columns
}

\keyword{colorSpec}
