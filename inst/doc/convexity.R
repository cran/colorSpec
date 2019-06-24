## ----setup, echo=FALSE, results="hide"--------------------------------------------------
require("knitr",quietly=TRUE)
opts_chunk$set(fig.path="figs/ag2-", fig.align="center",
  fig.width=7, fig.height=7, comment="")
knit_hooks$set(output = function(x, options) {
  paste('\\begin{Soutput}\n', x, '\\end{Soutput}\n', sep = '')
})
options(width=90)
if(!file.exists("figs")) dir.create("figs")

## ----packs, echo=TRUE, message=FALSE----------------------------------------------------
library( colorSpec )

## ----mets5, echo=TRUE, message=FALSE----------------------------------------------------
mets = responsivityMetrics( xyz1931.5nm )
mets$zeros

## ----mets6, echo=TRUE, message=FALSE----------------------------------------------------
mets$multiples

## ----mets7, echo=TRUE, message=FALSE----------------------------------------------------
mets$concavities

## ----mets8, echo=TRUE, message=FALSE----------------------------------------------------
wave  = wavelength(xyz1931.5nm)
E.eye = product( illuminantE(1,wave=wave), '*', xyz1931.5nm )
spec = canonicalOptimalColors( E.eye, c(580,585), spectral=TRUE ) 
bandRepresentation( spec )[[1]]

## ----fig1, echo=TRUE, fig.pos="H", fig.height=3.5, out.width='1.0\\linewidth', fig.cap='An example of a transmittance spectrum that is optimal, but has more than 2 transitions'----
par( omi=c(0,0,0,0), mai=c(0.5,0.6,0.2,0) )
plot( spec, main=FALSE, legend=FALSE )

## ----mets10, echo=TRUE, message=FALSE---------------------------------------------------
mets = responsivityMetrics( xyz1931.1nm )
mets$zeros
mets$multiples

## ----mets11, echo=TRUE, message=FALSE---------------------------------------------------
nrow( mets$concavities )

## ----mets12, echo=TRUE, message=FALSE---------------------------------------------------
fivenum( mets$concavities$extangle )
mets$concavities[ mets$concavities$extangle <= -0.01606, ]

## ----mets13, echo=TRUE, message=FALSE---------------------------------------------------
wave  = wavelength(xyz1931.1nm)
E.eye = product( illuminantE(1,wave=wave), '*', xyz1931.1nm )
spec = canonicalOptimalColors( E.eye, c(407,409), spectral=TRUE ) 
bandRepresentation( spec )[[1]]

## ----fig2, echo=TRUE, fig.pos="H", fig.height=3.5, out.width='1.0\\linewidth', fig.cap='An example of a transmittance spectrum that is optimal, but has more than 2 transitions'----
par( omi=c(0,0,0,0), mai=c(0.5,0.6,0.2,0) )
plot( spec, main=FALSE, legend=FALSE )

## ----finish, echo=FALSE, results="asis"-------------------------------------------------
knit_hooks$set(output = function(x, options) { x })
toLatex(sessionInfo(), locale=FALSE)

