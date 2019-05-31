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
source( "optimal-help.R" )

# make vctor of levels to be used for all the plots
Ylevel=c( seq( 0.10, 0.90, by=0.1 ), 0.95 ) 

## ----lee10, echo=TRUE, message=FALSE----------------------------------------------------
wave  = seq(380,800,by=5)
A.eye = product( A.1nm, "material", xyz1931.1nm, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), A.eye )
seclist = sectionOptimalColors( A.eye, normal=c(0,1,0), beta=white[2]*Ylevel )

## ----lee12, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth', fig.cap='MacAdam Limits for Illuminant A'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
plotOptimals( seclist, Ylevel, xyz1931.1nm, white )

## ----lee20, echo=TRUE, message=FALSE----------------------------------------------------
wave  = seq(380,800,by=5)
D65.eye = product( D65.1nm, "material", xyz1931.1nm, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), D65.eye )
seclist = sectionOptimalColors( D65.eye, normal=c(0,1,0), beta=white[2]*Ylevel )

## ----lee22, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth', fig.cap='MacAdam Limits for Illuminant D65'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
plotOptimals( seclist, Ylevel, xyz1931.1nm, white )

## ----lee40, echo=TRUE, message=FALSE----------------------------------------------------
wave  = seq(380,780,by=5)
C.eye = product( C.5nm, "material", xyz1931.1nm, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), C.eye )
seclist = sectionOptimalColors( C.eye, normal=c(0,1,0), beta=white[2]*Ylevel )

## ----lee42, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth', fig.cap='MacAdam Limits for Illuminant C'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
plotOptimals( seclist, Ylevel, xyz1931.1nm, white )

## ----lee30, echo=TRUE, message=FALSE----------------------------------------------------
wave = seq(420,680,by=5)
Flea2.scanner = product( A.1nm, "material", Flea2.RGB, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), Flea2.scanner )
seclist = sectionOptimalColors( Flea2.scanner, normal=c(0,1,0), beta=white[2]*Ylevel )

## ----lee32, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth', fig.cap='Approximate Output Limits for an RGB scanner'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
plotOptimals( seclist, Ylevel, Flea2.scanner, white )

## ----finish, echo=FALSE, results="asis"-------------------------------------------------
knit_hooks$set(output = function(x, options) { x })
toLatex(sessionInfo(), locale=FALSE)

