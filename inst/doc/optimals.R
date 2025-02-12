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

# make vector of levels to be used for the sections in all the plots
Ylevel=c( seq( 0.10, 0.90, by=0.1 ), 0.95 ) 

## ----lee10, echo=TRUE, message=FALSE----------------------------------------------------
wave  = seq(380,800,by=2)
A.eye = product( A.1nm, "material", xyz1931.1nm, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), A.eye )

## ----lee12, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth', fig.cap='MacAdam Limits for Illuminant A'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
seclist = sectionOptimalColors( A.eye, normal=c(0,1,0), beta=white[2]*Ylevel )
plotSections( seclist, Ylevel, xyz1931.1nm, white, col='red' )
seclist = sectionSchrodingerColors( A.eye, normal=c(0,1,0), beta=white[2]*Ylevel )
plotSections( seclist, Ylevel, xyz1931.1nm, white, add=TRUE )

## ----lee20, echo=TRUE, message=FALSE----------------------------------------------------
wave  = seq(380,800,by=2)
D65.eye = product( D65.1nm, "material", xyz1931.1nm, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), D65.eye )

## ----lee22, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth', fig.cap='MacAdam Limits for Illuminant D65'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
seclist = sectionOptimalColors( D65.eye, normal=c(0,1,0), beta=white[2]*Ylevel )
plotSections( seclist, Ylevel, xyz1931.1nm, white, col='red' )
seclist = sectionSchrodingerColors( D65.eye, normal=c(0,1,0), beta=white[2]*Ylevel )
plotSections( seclist, Ylevel, xyz1931.1nm, white, add=TRUE )

## ----lee25, echo=TRUE, message=FALSE----------------------------------------------------
wave  = seq(380,800,by=5)
D65.eye = product( D65.1nm, "material", lms2000.1nm, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), D65.eye )

## ----lee26, echo=TRUE, fig.pos="H", fig.height=4.6, out.width='1.0\\linewidth', fig.cap='MacAdam Limits for Illuminant D65, with updated cone fundamentals'----
par( omi=rep(0,4), mai=c(0.5,.6,0,0) )
normal = c(1,1,1)/3  ;  beta = sum(white*normal) * Ylevel
seclist = sectionOptimalColors( D65.eye, normal=normal, beta=beta )
plotSections( seclist, Ylevel, lms2000.1nm , white, col='red' )
seclist = sectionSchrodingerColors( D65.eye, normal=normal, beta=beta )
plotSections( seclist, Ylevel, lms2000.1nm , white, add=TRUE )

## ----lee40, echo=TRUE, message=FALSE----------------------------------------------------
wave  = seq(380,780,by=2)
C.eye = product( C.5nm, "material", xyz1931.1nm, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), C.eye )

## ----lee42, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth', fig.cap='MacAdam Limits for Illuminant C'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
seclist = sectionOptimalColors( C.eye, normal=c(0,1,0), beta=white[2]*Ylevel )
plotSections( seclist, Ylevel, xyz1931.1nm, white, col='red' )
seclist = sectionSchrodingerColors( C.eye, normal=c(0,1,0), beta=white[2]*Ylevel )
plotSections( seclist, Ylevel, xyz1931.1nm, white, add=TRUE )

## ----lee60, echo=TRUE, message=FALSE, fig.pos="H", fig.height=5, out.width='1.0\\linewidth'----
wave  = seq(320,600,by=2)
path  = system.file( 'extdata/eyes/BeeEye.txt', package='colorSpec' )
bee   = readSpectra( path, wavelength=wave )
plot( bee )

## ----lee62, echo=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth'----------
E.eye = product( illuminantE(1,wavelength=wave), "material", bee )
white = product( neutralMaterial(1,wave=wave), E.eye )
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
normal = c(1,1,1)/3  ;  beta = sum(white*normal) * Ylevel
seclist = sectionOptimalColors( E.eye, normal=normal, beta=beta )
plotSections( seclist, Ylevel, bee, white, col='red' )
seclist = sectionSchrodingerColors( E.eye, normal=normal, beta=beta )
plotSections( seclist, Ylevel, bee, white, add=TRUE )

## ----lee30, echo=TRUE, message=FALSE----------------------------------------------------
wave = seq(420,680,by=5)
Flea2.scanner = product( A.1nm, "material", Flea2.RGB, wavelength=wave )
white = product( neutralMaterial(1,wave=wave), Flea2.scanner )

## ----lee32, echo=TRUE, message=TRUE, fig.pos="H", fig.height=6.5, out.width='1.0\\linewidth'----
par( omi=rep(0,4), mai=c(0.5,0.6,0,0) )
normal = c(1,1,1)/3  ;  beta = sum(white*normal) * Ylevel
seclist = sectionOptimalColors( Flea2.scanner, normal=normal, beta=beta )
plotSections( seclist, Ylevel, Flea2.scanner, white, col='red' )
seclist = sectionSchrodingerColors( Flea2.scanner, normal=normal, beta=beta )
plotSections( seclist, Ylevel, Flea2.scanner, white, add=TRUE )

## ----finish, echo=FALSE, results="asis"-------------------------------------------------
knit_hooks$set(output = function(x, options) { x })
toLatex(sessionInfo(), locale=FALSE)

