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
Only certain types of sequences are allowed.
The return value can be a new \bold{colorSpec} object or a matrix;
see \bold{Details}.
}

\usage{
\S3method{product}{colorSpec}( \ldots )

}

\arguments{
  \item{...}{unnamed arguments are \bold{colorSpec} objects,
             and possibly a single character string, see \bold{Details}. 
             Possible named arguments are:

\describe{
\item{\code{wavelength}}{The default \code{wavelength='identical'} means that all the \bold{colorSpec} objects must have the same wavelength sequence; if they do not it is an ERROR. \code{wavelength} can be a new wavelength sequence, and all the objects are then
\code{\link{resample}}d at these new wavelengths. \code{wavelength} can also be \code{'auto'} or \code{NULL} which means to compute a suitable wavelength sequence from those of the objects, see \bold{Details}. It is OK to abbreviate the string \code{wavelength}
(e.g. to \code{wave}); see \bold{Examples}.
It is OK for the wavelength sequence to be irregular;
when the return value is a matrix
the integration weights the spectrum values appropriately.\cr
}

\item{ \code{method}, \code{span}, \code{extrapolation}, \code{clamp} }{passed to \code{resample()} with no checking or changes\cr
}

\item{ \code{integration}}{only applies when the return type is matrix.
The default option is \code{'rectangular'}, which means to weight the spectrum
value equally at all wavelengths; this is the ASTM E308-01 recommendation.
The other option is \code{'trapezoidal'}, which means to give the 2 endpoint
wavelength values 1/2 the weight of the others.
Trapezoidal integration is provided mostly for compatibility with other software.}
}
}
}


    
\value{
\code{product()} returns either a \bold{colorSpec} object or a matrix, see \bold{Details}.

If \code{product()} returns a \bold{colorSpec} object, the \code{\link{organization}}
of the object is \code{'matrix'} or \code{'vector'};
any \code{\link{extradata}} is lost.
However, all terms in the product are saved in \code{attr(*,'sequence')}.
One can use \code{\link{str}()} to inspect this attribute. 

If \code{product()} returns a matrix,
this matrix can sometimes be ambiguous, see \bold{Note}.

All actinometric terms are converted to radiometric on-the-fly and the returned
\bold{colorSpec} object is also radiometric.

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
 \eqn{R_M} \tab \code{responsivity.material} \tab a material responder (e.g. a scanner)\cr 
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
 


4. \if{html}{\out{&nbsp;&nbsp;}} \eqn{L * M_1 * ... *} \bullet{} \eqn{* ... * M_m * R_L }  \mapsto{} \eqn{R_M'} \cr  
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
The numbers of spectra in the terms must be conformable as in sequence 1,
and the product is a material responder with \eqn{q} spectra.
\cr
In the function \code{product()} the location of the \bullet{} is marked
by any character string whatsoever - it's up to the user who might choose
something that describes the typical material (between the light source and camera).
For example one might choose:\cr
\code{scanner = product( A.1nm, 'photo', Flea2.RGB, wave='auto') }\cr
to model a scanner that is most commonly used to scan photographs.
Other possible strings could be \code{'artwork'}, \code{'crystal'}, \code{'varmaterial'},
or even \code{'slot'}.
See the \bold{gallery} vignette for a worked-out example.
\cr \cr
  
  
  
5. \if{html}{\out{&nbsp;&nbsp;}}  \eqn{ L * M_1 * M_2 * ... * M_m * R_L }  \mapsto{} \eqn{matrix} \cr  
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
The sequence of spectra counts is \code{c(1,24,3)}
which is splittable into \code{c(1,24)} and \code{c(3)}.
The product matrix is 24\times{}3.
See the \bold{gallery} vignette for a worked-out example.\cr
Note that is OK for there to be no materials in this product;
it is OK if \eqn{m=0}.
See the \bold{blueflame} vignette for a worked-out example.
\cr \cr  
  
  
6. \if{html}{\out{&nbsp;&nbsp;}}  \eqn{ M_1 * M_2 * ... * M_m * R_M } \mapsto{} \eqn{matrix} \cr  
The product of \eqn{m} materials followed by a material responder, is a matrix !
The sequence of numbers of spectra must be splittable into left and right
parts as in sequence 4, and the product matrix is formed the same way.
One reason for computing this matrix in 2 steps is that one can
\code{\link{calibrate}} the material responder separately in a customizable way.
See the \bold{gallery} vignette for a worked-out example with a flatbed scanner.
\cr \cr  

Note that sequences 5. and 6. are the only ones that 
use the usual matrix product \code{\%*\%}.
They may also use the Hadamard matrix product \code{*}, as in sequences 1 to 4.
  
  
The argument \code{wavelength} can also be \code{'auto'} or \code{NULL}.
In this case the intersection of all the wavelength ranges of the objects is computed.
If the intersection is empty, it is an ERROR and the function returns \code{NULL}.
The wavelength step \code{\link{step.wl}} 
is taken to be the smallest over all the object wavelength sequences.
If the minimum \code{\link{step.wl}} is less than 1 nanometer,
it is rounded off to the nearest power of 2 (e.g 1, 0.5, 0.25, ...).
}

\note{
The product for sequences 1, 2, and 3 is associative.
After all matrices are filled out to have  \eqn{q} columns,
the result is essentially a Hadamard product of matrices,
which is associative.
Also note that a subsequence of sequences 2 and 3 might be sequence 1.

The product for sequence 4 is never associative,
since subproducts that contain the variable \bullet{} are undefined.
However the  result is essentially a Hadamard product of matrices, and is unambiguous.

The product for sequence 5 is associative in special cases, but not in general.
The problem is that the left and right splitting point is not unique.
If all objects have only a single spectrum, then it *is* associative,
and therefore unambiguous.
If the left part has a different number of multiple spectra than the right part,
then it is not associative in general since some ways of grouping the 
product may be undefined.\cr
Moreover, in some cases the product can be ambiguous.
Suppose that the vector of spectrum counts is \code{c(1,3,1)};
this could come from a single light source, followed by 3 filters (e.g. RGB),
followed by a graylevel camera.
There are 2 ways to split this: \code{"1|3,1"} and \code{"1,3|1"}.
The first split is interpreted as the light source into a camera with 3 channels.
The second split is interpreted as 3 colored light sources into a graylevel camera.
In the first split the returned matrix is a \code{1x3} row vector.
In the second split the returned matrix is a \code{3x1} column vector.
For the vector \code{"1,3,1"}, one can show that
the computed components in the 2 vectors are equal, 
so the ambiguity is benign.
But consider the longer sequence  \code{"1,3,3,1"}.
There are 3 ways to split this, and the returned matrices are
\code{1x3}, \code{3x3}, and \code{3x1}.
So this ambiguity is obviously a problem.
Whenever there is an ambiguity, the function chooses a splitting
in which the left part is as long as possible, and issues a warning message.
The user should inspect the result carefully.
To avoid the ambiguity, the user should break the product into smaller pieces
and call \code{product()} multiple times.

The product for sequence 6 is essentially the same as sequence 5,
and the function issues a warning message when appropriate.
Note that a subproduct is defined only if it avoids the final multiplication with \eqn{R_M}.


}

\references{
Edward J. Giorgianni and Thomas E. Madden.
\bold{Digital Color Management: Encoding Solutions.} 2nd Edition
John Wiley. 2009.
Figure 10.11a. page 141.

Wikipedia.
\bold{Hadamard product (matrices)}.
\url{http://en.wikipedia.org/wiki/Hadamard_product_\%28matrices\%29}

ASTM E308-01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
(2001).
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
path = system.file( "extdata/objects/Midwest-SP700-2014.txt", package='colorSpec' )
blocker.IR = readSpectra( path )
product( blocker.IR, Hoya, wave='auto' )


#  sequence 2.
product( subset(solar.irradiance,1), atmosphere2003, blocker.IR, Hoya, wave='auto' )


#  sequence 3.
plumbicon = readSpectra( system.file( "extdata/cameras/plumbicon30mm.txt", package='colorSpec' ) )
product( blocker.IR, subset(Hoya,1:3), plumbicon, wave='auto' )


#   sequence 4.
#   make an RGB scanner
bluebalancer = subset(Hoya,'LB')
# combine tungsten light source A.1nm with blue light-balance filter
# use the string 'artwork' to mark the variable material slot
scanner = product( A.1nm, bluebalancer, 'artwork', Flea2.RGB, wave='auto' )


#  sequence 5.
product( D65.1nm, Flea2.RGB )  # a 1x3 matrix, no materials
product( D65.1nm, neutralMaterial(0.01), Flea2.RGB, wave='auto' ) # a 1x3 matrix, 1 material
path = system.file( "extdata/sources/Lumencor-SpectraX.txt", package='colorSpec' )
lumencor = readSpectra( path, wave=340:660 )
product( lumencor, Flea2.RGB, wave='auto' )   # a 7x3 matrix, no materials


#  sequence 6.
scanner = calibrate( scanner )
target = readSpectra( system.file( "extdata/targets/N130501.txt", package='colorSpec') )
product( target, scanner, wave='auto' )  #  a 288x3 matrix

}

\keyword{colorSpec}
