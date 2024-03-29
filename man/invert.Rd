% File estimate.Rd
\encoding{UTF-8}

\name{invert}
\title{estimate spectra from responses, effectively inverting the operator from spectrum to response}
\alias{invert}
\alias{invert.colorSpec}

\description{
Given a light responder (e.g. an eye or a camera),
two light spectra that produce the same response from the responder
are called \emph{metamers} for that responder.
Similarly, given a material responder (e.g. a scanner),
two reflectance spectra that produce the same response from the responder
are called \emph{metamers} for that responder.

For a given responder and response, 
there are typically infinitely many \emph{metamers}.
The set of all of them is often called the \emph{metameric suite}.
The goal of the function \code{invert()} is to calculate a "good" metamer
in the "suite".
\cite{Koenderink} calls this topic \emph{inverse colorimetry}.
In the case that the estimated spectrum is a reflectance spectrum,
the topic is often called \emph{reflectance estimation} or \emph{reflectance recovery},
see \cite{Bianco}.

The \emph{centroid method}, which is the default and the featured method
in this package, computes the centroid of the set of all the metamers (if any).
The centroid is computed in an infinite-dimensional context 
and is expounded further in
\cite{Davis}.

The \emph{Hawkyard method}, see \cite{Hawkyard} and \cite{Bianco},
has been around a long time.
The centroid and Hawkyard methods have similarities,
e.g. both are low-dimensional with the number of variables equal to
the number of responses (usually 3).
The Hawkyard method is very fast, but has a key problem, see below.

The \emph{Transformed Least Slope Squared} (TLSS) method was developed
by Scott Burns, see \bold{References}.
This is my name for it, not his.
What I call TLLS is actually is a combination of Burns' LHTSS and LLSS methods;
the one that \code{invert()} chooses depends on \code{type(x)}, see below.
Both of these are high-dimensional,
with the number of variables equal to #(responses) + #(wavelengths).

The first argument to \code{invert()} is the responder \code{x},
and the second is the matrix \code{response} of responses (e.g. XYZs or RGBs).

The goal is to return a "good" spectrum for each response so that:
\tabular{c}{
\code{product( invert(x,response), x )} \eqn{~\cong~}{~=~} \code{response}
}
The error is returned as column \code{estim.precis}, see below.


First consider the case where \code{x} has type \code{type='responsivity.material'}.
The goal is to compute a reflectance spectra.
All the methods will fail if the response is on the object-color boundary
(an \emph{optimal color}) or outside the boundary.
They \emph{may} also fail if the response is inside the object-color
solid (the \emph{Rösch Farbkörper}) and very close to the boundary.
\cr
The centroid method solves a non-linear system that contains a 
\emph{Langevin-function-based} squashing function, see \cite{Davis} for details.
When successful it always returns a feasible spectrum
with small \code{estim.precis}.
\cr
The Hawkyard method is linear and very fast,
but in raw form it may return a non-feasible reflectance spectrum.
In this case \code{invert()} simply clamps to the interval [0,1] and so
\code{estim.precis} can be large.
\cr
The TLSS method solves a non-linear system that contains
the squashing function \eqn{(\tanh(z) + 1)/2}, see \cite{Burns} for details.
When successful it always returns a feasible spectrum
with small \code{estim.precis}.

Now consider the case where \code{x} has \code{type='responsivity.light'}.
The goal is to compute the spectrum of a light source.
All the methods will fail if chromaticity of the response is on the boundary
of the inverted-U (assuming \code{x} models the human eye) or outside the boundary.
They \emph{may} also fail if the response is inside the inverted-U
and very close to the boundary.
\cr
The centroid method works on a relatively small range of chromaticities;
it will fail if the response is too far from the response to Illuminant E.
See \cite{Davis} for the details.
When successful it always returns an everywhere positive spectrum
with small \code{estim.precis}.
This method has the desirable property that if the response is multiplied by
a positive number, the computed spectrum is multiplied by that same number.
\cr
The Hawkyard method does not work in this case.
\cr
The TLSS method solves a non-linear system that contains
the squashing function \eqn{\exp(z)}, see \cite{Burns} for the details.
When successful it always returns an everywhere positive spectrum
with small \code{estim.precis}.
This method succeeds on a larger set of chromaticities than the centroid method.
It also has the desirable scaling multiplication property mentioned above.

The centroid and Hawkyard methods have an equalization option,
which is controlled by the argument \code{alpha}
and is enabled by default, see below.
When enabled, if the response comes from a constant spectrum
(a perfectly neutral gray material, or a multiple of Illuminant E),
then the computed spectrum is that same constant spectrum (up to numerical precision).
I call this the \emph{neutral-exact property}.
Equalization is a complicated mechanism, for details see \cite{Davis}.
For the TLSS method, the neutral-exact property is intrinsic,
and \code{alpha} is ignored.
\cr
NOTE:
If the responder has only one output channel (e.g. a monochrome camera)
and equalization is enabled,
then \emph{all} responses are inverted to a constant spectrum.
This may or may not be desirable.
}


\usage{
\S3method{invert}{colorSpec}( x, response, method='centroid', alpha=1 )

}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type(x)} = 
\code{'responsivity.material'} or \code{'responsivity.light'} and M responsivities.
The wavelengths must be regular (equidistant).}
  
\item{response}{ a numeric NxM matrix, or a numeric vector that can be converted
    to such matrix, by row. The N responses are contained in the rows.
    The \code{rownames(response)} are copied to the output \code{specnames}.}
    
\item{method}{either \code{'centroid'} or \code{'Hawkyard'} or \code{'TLSS'}.  
    \code{'Hawkyard'} is only valid when 
    \code{type(x)} is \code{'responsivity.material'}.
    Matching is partial and case-insensitive.}
    
\item{alpha}{a vector of M weighting coefficients, 
    or a single number that is replicated to length M. 
    When \code{method='centroid'}, \code{alpha} is used for \emph{equalizing} 
    the responsivities, which is recommended.
    For \code{alpha} to be valid, the linear combination of the M responsitivies,
    with coefficients \code{alpha}, must be positive.
    To disable equalization (not recommended) and use the original responsivities, 
    set \code{alpha=NULL}. 
    Similarly, when \code{method='Hawkyard'}, \code{alpha} is used for equalizing 
    the responsivities, which is also recommended.
    When \code{method='TLSS'}, \code{alpha} is ignored.}
}
    
\value{
If \code{type(x)='responsivity.material'} it returns a \bold{colorSpec} object
with \code{type} = \code{'material'} (\code{quantity} = \code{'reflectance'}).

If \code{type(x)='responsivity.light'} it returns a \bold{colorSpec} object
with \code{type} = \code{'light'}
(\code{quantity}=\code{'energy'} or \code{quantity}=\code{'photons'}
depending on \code{quantity(x)}).

In either case, the returned object has \code{organization} = \code{'df.row'} 
and the \code{extradata} is a \code{data.frame} with these columns:

  \item{response}{the input matrix of desired responses}
  
  \item{estim.precis}{the difference between the desired response and actual response.  It   is the mean of the absolute value of the differences.
  See \code{\link[rootSolve:multiroot]{rootSolve::multiroot}()} }  
  
  \item{time.msec}{the time to compute the spectrum, in msec. 
       When \code{method='Hawkyard'}, all N spectra are computed at once, 
       so all N spectra are assigned the same mean time.}
  
  \item{iters}{the number of iterations that were required to find the relevant root. 
              This is not present when \code{method='Hawkyard'}.}
    
  \item{clamped}{a logical indicating whether the reflectance was clamped to [0,1]. This is present only  when \code{method='Hawkyard'}.}
    

If a response could not be estimated, 
the row contains \code{NA} in appropriate columns,
and a warning is issued.\cr
  
In case of global error, the function returns \code{NULL}.
}


\details{
For \code{method='centroid'} the function calls the non-linear root-finder 
\code{\link[rootSolve:multiroot]{rootSolve::multiroot}()},
which is general purpose and "full Newton".
% \code{\link[rootSolve]{multiroot}()}.

For \code{method='Hawkyard'} the function solves a linear system by
inverting a small matrix (#[responses] x #[responses]).
The spectra are then clamped to [0,1].

For \code{method='TLSS'} the function solves a constrained least-squares problem
using Lagrange multipliers.
A critical point is found using a "full Newton" iteration.
The original MATLAB code is posted at \cite{Burns},
and was  ported from MATLAB to R with only trivial changes.
When computing a reflectance spectrum, the Hawkyard method is used for the
initial guess, after little extra clamping.
This improved guess cuts the number of iterations substantially, 
and the extra computation time is negligible.
}


\section{Known Issues}{
If \code{type(x)='responsivity.light'} the centroid method may fail
(not converge) if the response is too far from that of Illuminant E.
}


\seealso{
\code{\link{type}()},
\code{\link{quantity}()},
\code{\link{organization}()},
\code{\link{specnames}()},
\code{\link{product}()},
\code{\link{extradata}()},
\code{\link[rootSolve:multiroot]{rootSolve::multiroot}()},
vignette \href{../doc/inversion.pdf}{\bold{Estimating a Spectrum from its Response}}
}

\examples{
wave = 400:700
E.eye = product( illuminantE(1,wave), "material", xyz1931.1nm, wavelength=wave )
path = system.file( 'extdata/targets/CC_Avg30_spectrum_CGATS.txt', package='colorSpec' )
MacbethCC = readSpectra( path, wavelength=wave )
XYZ = product( MacbethCC, E.eye, wavelength=wave )
est.eq   = invert( E.eye, XYZ, method='centroid', alpha=1 )
extra   = extradata(est.eq)
range(extra$estim.precis)       # prints   0.000000e+00 3.191741e-08
}


\references{
Davis, Glenn.
A Centroid for Sections of a Cube in a Function Space, with Application to Colorimetry.
\url{https://arxiv.org/abs/1811.00990}.
[math.FA].
2018.

Bianco, Simone.
Reflectance spectra recovery from tristimulus values by adaptive estimation with metameric shape correction.
vol. 27, no 8.
\emph{Journal of the Optical Society of America A}.
pages 1868-1877.
2010
\url{https://opg.optica.org/josaa/abstract.cfm?uri=josaa-27-8-1868}.


Burns, Scott A. 
Generating Reflectance Curves from sRGB Triplets.
\url{http://scottburns.us/reflectance-curves-from-srgb/}.

Hawkyard, C. J.
Synthetic reflectance curves by additive mixing.
\emph{Journal of the Society of Dyers and Colourists}.
vol. 109. no. 10.
Blackwell Publishing Ltd.
pp. 323-329.
1993.

Koenderink, J.J.
Color for the Sciences.
MIT Press.
2010.
}

\keyword{colorSpec}
