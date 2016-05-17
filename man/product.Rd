% File product.Rd
\encoding{UTF-8}

\newcommand{\times}{\if{latex}{\out{$\times$}}\ifelse{html}{\out{&times;}}}
\newcommand{\bullet}{\if{latex}{\out{$\bullet$}}\ifelse{html}{\out{&bull;}}}
\newcommand{\mapsto}{\if{latex}{\out{$\mapsto$}}\ifelse{html}{\out{&#x21a6;}}}

\name{product}
\title{Compute the product of colorSpec objects}
\alias{product}
\alias{product.colorSpec}

\description{
Take a sequence of \bold{colorSpec} objects and compute their product.
The product is associative.
Only certain types of sequences are allowed, see \bold{Details}.
}

\usage{
\S3method{product}{colorSpec}( \ldots )

}

\arguments{
  \item{...}{Unnamed arguments are \bold{colorSpec} objects,
             and possibly a single character string, see \bold{Details}. 
             Named arguments (if any) are passed to \code{\link{resample}}.\cr
             The most important named argument is \code{wavelength}.
             The default \code{wavelength='identical'} means that all the \bold{colorSpec} objects must have the same wavelength sequence; if they do not it is an ERROR.\cr
  \code{wavelength} can be a new wavelength sequence, and all the objects are then
   \code{\link{resample}}d at these new wavelengths.\cr
  \code{wavelength} can also be \code{'auto'} or \code{NULL} which means to compute a suitable wavelength
  sequence from those of the objects, see \bold{Details}.
It is OK to abbreviate the string \code{wavelength} (e.g. to \code{wave}); see \bold{Examples}.\cr
Other possible named arguments are \code{method} and \code{span}.
}
}

    
\value{
\code{product()} returns a \bold{colorSpec} object or a matrix, see \bold{Details}.
In case of a \bold{colorSpec} object, the \code{\link{organization}} is
\code{'matrix'} or \code{'vector'}; any \code{\link{extradata}} is lost.
However, all terms in the product are saved in \code{attr(*,'sequence')}.
One can use \code{\link{str}} to inspect this attribute. \cr
All terms are converted to \code{\link{radiometric}} on-the-fly and the returned
\bold{colorSpec} object is also radiometric.\cr
In case of ERROR it returns \code{NULL}.
           
}    


\details{
To explain the allowable product sequences  
it is helpful to introduce some simple notation for the objects:

 \tabular{lll}{
 \bold{notation} \tab \bold{colorSpec} \code{type} \tab \bold{description of the object} \cr
 \eqn{L} \tab \code{light} \tab a light source \cr
 \eqn{M} \tab \code{material} \tab a material \cr
 \eqn{R_L} \tab \code{responsivity.light} \tab a light responder (aka detector)\cr
 \eqn{R_M} \tab \code{responsivity.material} \tab a material responder\cr 
 }
 
It is also helpful to define a sequence of positive integers
to be \emph{conformable} iff it has at most one value greater than 1. 
For example,
a sequence of all 1s is conformable.  A sequence of all \eqn{q}'s is conformable.
The sequences \code{c(1,3)} and \code{c(1,1,4,1,1,4,1)} are conformable,
but \code{c(1,1,4,1,3,4,1)} is not.
 
There are 6 types of sequences for which the product is defined:
\cr
 
1. \if{html}{\out{&nbsp;&nbsp;}} \eqn{ M_1 * M_2 * ... * M_m }  \mapsto{} \eqn{M'} \cr  
The product of \eqn{m} materials is another material.  
Think of a stack of \eqn{m} transmitting filters effectively forming a new filter.
If we think of each object as a matrix (with the spectra in the columns),
then the product is element-by-element using \R's \code{*} - the Hadamard product.
The numbers of spectra in the terms must be conformable.
If some objects have 1 spectrum and all the others have \eqn{q},
then the column-vector spectrums are repeated \eqn{q} times to form a
matrix with \eqn{q} columns.
If the numbers of spectra are not conformable, 
it is an ERROR and the function returns \code{NULL}.
\cr
As an example, suppose \eqn{M_1} has 1 spectrum and \eqn{M_2} has \eqn{q} spectra,
and \eqn{m=2}.
Then the product is a material with \eqn{q} spectra.
Think of an IR-blocking filter followed by the RGB filters in a 3-CCD camera.
\cr \cr
 
2. \if{html}{\out{&nbsp;&nbsp;}}  \eqn{ L * M_1 * M_2 * ... * M_m }  \mapsto{} \eqn{L'} \cr  
The product of a light source followed by \eqn{m} materials is a light source.  
Think of a light source
followed by a stack of \eqn{m} transmitting filters, effectively forming a new light source.
The numbers of spectra in the terms must be conformable as in sequence 1,
and the matrices are multiplied element by element.\cr
As an example, suppose \eqn{L} has 1 spectrum and \eqn{M_1} has \eqn{q} spectra,
and \eqn{m=1}.
Then the product is a light source with \eqn{q} spectra.
Think of a light source followed by a filter wheel with \eqn{q} filters.
\cr \cr
 
3. \if{html}{\out{&nbsp;&nbsp;}}  \eqn{ M_1 * M_2 * ... * M_m * R_L }  \mapsto{} \eqn{R_L'} \cr  
The product of \eqn{m} materials followed by a light responder, is a light responder.
Think of a stack of \eqn{m} transmitting filters in front of a camera, effectively forming a new camera.
The numbers of spectra in the terms must be conformable as in sequence 1,
and the matrices are multiplied element by element.\cr
As an example, suppose \eqn{R_L} has 1 spectrum and \eqn{M_1} has \eqn{q} spectra,
and \eqn{m=1}.
Then the product is a light responder with \eqn{q} spectra.
Think of a 3-CCD camera in which all 3 CCDs have exactly the same responsivity
and so can be modeled with a single object \eqn{R_L}.
\cr \cr
 
  
4. \if{html}{\out{&nbsp;&nbsp;}}  \eqn{ L * M_1 * M_2 * ... * M_m * R_L }  \mapsto{} \eqn{matrix} \cr  
The product of a light source, followed by \eqn{m} materials, 
followed by a light responder, is a matrix!
The numbers of spectra in the terms must splittable into
a conformable left part (\eqn{L'} from sequence 2.) 
and a conformable right part (\eqn{R_L'} from sequence 3.).
There is a row for each spectrum in \eqn{L'},
and a column for each spectrum in \eqn{R_L'}.
Suppose the element-by-element product of the left part is 
\eqn{n}\times{}\eqn{p} 
and the element-by-element product of the right part is 
and \eqn{n}\times{}\eqn{q},
where \eqn{n} is the number of wavelengths.
Then the output matrix is the usual matrix product \code{\%*\%}
of the transpose of the left part times and right part,
which is \eqn{p}\times{}\eqn{q}.
\cr
As an example, think of a light source followed by a 
reflective color target with 24 patches
followed by an RGB camera.
The sequence of spectra is \code{c(1,24,3)}
which is splittable into \code{c(1,24)} and \code{c(3)}.
The product matrix is 24\times{}3.
See the \bold{gallery} vignette for a worked-out example.\cr
Note that is OK for there to be no materials in this product;
it is OK if \eqn{m=0}.
See the \bold{blueflame} vignette for a worked-out example.
\cr \cr

5. \if{html}{\out{&nbsp;&nbsp;}} \eqn{L * M_1 * ... *} \bullet{} \eqn{* ... * M_m * R_L }  \mapsto{} \eqn{R_M'} \cr  
This is the strangest product.
The bullet symbol \bullet{} means that a variable material is inserted at that slot
in the sequence (or light path).
For each material spectrum inserted there is a response from \eqn{R_L}.
Therefore the product of this sequence is a material responder \eqn{R_M}.
Think of a light source \eqn{L} going through 
a transparent object  \bullet{}  on a flatbed scanner and into a camera \eqn{R_L}.
For more about the mathematics of this product, 
see the \bold{colorSpec-guide.pdf} in the doc directory.
These material responder spectra are the same as the
\emph{effective spectral responsivities} in \cite{Digital Color Management}.
The numbers of spectra in the terms must be conformable as in sequence 1.
\cr
In the function \code{product()} the location of the \bullet{} is marked
by any character string whatsoever - it's up to the user who might choose
something that describes the typical material.
For example one might choose:\cr
\code{scanner = product( A.1nm, 'photo', Flea2.RGB, wave='auto') }\cr
to model a scanner that is most commonly used to scan photographs.
Other possible strings could be \code{'artwork'}, \code{'varmaterial'},
or even \code{'slot'}.
See the \bold{gallery} vignette for a worked-out example.
\cr \cr
  
6. \if{html}{\out{&nbsp;&nbsp;}}  \eqn{ M_1 * M_2 * ... * M_m * R_M } \mapsto{} \eqn{matrix} \cr  
The product of \eqn{m} materials followed by a material responder, is a matrix !
The sequence of numbers of spectra must be splittable into left and right
parts as in sequence 4, and the product matrix is formed the same way.
One reason for computing this matrix in 2 steps is that one can
\code{\link{calibrate}} the material responder separately in a customizable way.
See the \bold{gallery} vignette for a worked-out example with a flatbed scanner.
\cr \cr  

Note that sequences 4. and 6. are the only ones that 
use the usual matrix product \code{\%*\%}.
  
  
The argument \code{wavelength} can also be \code{'auto'} or \code{NULL}.
In this case the intersection of all the wavelength ranges of the objects is computed.
If the intersection is empty, it is an ERROR and the function returns \code{NULL}.
The wavelength step \code{\link{step.wl}} 
is taken to be the smallest over all the object wavelength sequences.
If the minimum \code{\link{step.wl}} is less than 1 nanometer,
it is rounded off to the nearest power of 2 (e.g 1, 0.5, 0.25, ...).
\cr
Although the function signature shows the \pkg{colorSpec} objects, 
followed by \code{wavelength}, followed by more arguments,
actually they can come in any order.
The unnamed arguments are taken to be \pkg{colorSpec} objects
and the named arguments are taken to be arguments for \code{resample}.
}

\references{
Edward J. Giorgianni and Thomas E. Madden.
\bold{Digital Color Management: Encoding Solutions.} 2nd Edition
John Wiley. 2009.
Figure 10.11a. page 141.

Wikipedia.
\bold{Hadamard product (matrices)}.
http://en.wikipedia.org/wiki/Hadamard_product_\%28matrices\%29

}


\seealso{
  \code{\link{wavelength}},
  \code{\link{type}},
  \code{\link{resample}},  
  \code{\link{calibrate}},
  \code{\link{radiometric}},  
  \code{\link{step.wl}}
}

\examples{
#  sequence 1.
path = system.file( "extdata/filters/Midwest-SP700-2014.txt", package='colorSpec' )
blocker.IR = readSpectra( path )
product( blocker.IR, Hoya, wave='auto' )


#  sequence 2.
product( subset(solar.irradiance,1), atmosphere2003, blocker.IR, Hoya, wave='auto' )


#  sequence 3.
plumbicon = readSpectra( system.file( "extdata/cameras/plumbicon30mm.txt", package='colorSpec' ) )
product( blocker.IR, subset(Hoya,1:3), plumbicon, wave='auto' )


#  sequence 4.
product( D65.1nm, Flea2.RGB )  # a 1x3 matrix, no materials
product( D65.1nm, neutralMaterial(0.01), Flea2.RGB, wave='auto' ) # a 1x3 matrix, 1 material
path = system.file( "extdata/sources/Lumencor-SpectraX.txt", package='colorSpec' )
lumencor = readSpectra( path, wave=340:660 )
product( lumencor, Flea2.RGB, wave='auto' )   # a 7x3 matrix, no materials


#   sequence 5.
#   make an RGB scanner
bluebalancer = subset(Hoya,'LB')
# combine tungsten light source A.1nm with blue light-balance filter
# use the string 'artwork' to mark the variable material slot
scanner = product( A.1nm, bluebalancer, 'artwork', Flea2.RGB, wave='auto' )


#  sequence 6.
scanner = calibrate( scanner )
target = readSpectra( system.file( "extdata/targets/N130501.txt", package='colorSpec') )
product( target, scanner, wave='auto' )  #  a 288x3 matrix

}

\keyword{colorSpec}
