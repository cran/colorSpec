
library( colorSpec )
options( width=144 )   
options( colorSpec.stoponerror=FALSE )

testReflectanceInversion <- function()
    {
    wave = 400:700
    
    E.eye = product( illuminantE(1,wave), "material", xyz1931.1nm, wavelength=wave )
    #   A.eye = product( A.1nm, "material", xyz1931.1nm, wavelength=wave )
    
    
    ##------------     MacbethCC     ----------##
    cat( "\n##------------     MacbethCC     ----------##\n" )
    path = system.file( 'extdata/targets/CC_Avg30_spectrum_CGATS.txt', package='colorSpec' )
    MacbethCC = readSpectra( path, wavelength=wave )
    #   MacbethCC = subset( MacbethCC, c( 21:24, 17:20, 13:16, 9:12, 5:8, 1:4 ) )
    
    XYZ = product( MacbethCC, E.eye, wavelength=wave )
    est.eq   = invert( E.eye, XYZ, method='centroid', alpha=1 )
    if( is.null(est.eq) ) return(FALSE)
    
    extra   = extradata(est.eq)
    #extra$response = NULL   # do not want to see it
    print( extra )
    if( any( is.na(extra$estim.precis) ) )    return(FALSE)

    cat( 'mean(iters) = ', mean(extra$iters), '\n' )
    cat( 'mean(estim.precis) = ', mean(extra$estim.precis), '\n' )
    
    
    ##------------     rectangular, i.e. Logvinenko     ----------##
    cat( "\n##------------     rectangular     ----------##\n" )    
    lambda = c(600,450, 650,500, 450,500, 500,550, 550,600, 600,650 )
    lambda = c( lambda, 500,525, 525,550, 550,575 )
    nearsemichrome = rectangularMaterial( lambda, 0.98, wave )
    XYZ = product( nearsemichrome, E.eye, wavelength=wave )
    est.c = invert( E.eye, XYZ, method='centroid' )
    if( is.null(est.c) ) return(FALSE)
        
    extra   = extradata(est.c)
    #extra$response = NULL   # do not want to see it
    print( extra )
    if( any( is.na(extra$estim.precis) ) )    return(FALSE)

    cat( 'mean(iters) = ', mean(extra$iters), '\n' )
    cat( 'mean(estim.precis) = ', mean(extra$estim.precis), '\n' )
    
    
    ##------------     camera + 3 illuminants     ----------##
    cat( "\n##------------     camera + 3 illuminants     ----------##\n" )        
    E.flea = product( illuminantE(1,wave), "material", Flea2.RGB, wavelength=wave )    
    A.flea = product( A.1nm, "material", Flea2.RGB, wavelength=wave )    
    P.flea = product( planckSpectra(9000), "material", Flea2.RGB, wavelength=wave )
    specnames( E.flea ) = c('rE','gE','bE')
    specnames( A.flea ) = c('rA','gA','bA')  
    specnames( P.flea ) = c('rP','gP','bP')
    PEA.flea = bind( P.flea, E.flea, A.flea )
    response = product( MacbethCC, PEA.flea, wavelength=wave )
    est.eq   = invert( PEA.flea, response, method='centroid', alpha=1 )
    if( is.null(est.eq) ) return(FALSE)
        
    extra   = extradata(est.eq)
    extra$response = NULL   # do not want to see it
    print( extra )
    if( any( is.na(extra$estim.precis) ) )    return(FALSE)

    cat( 'mean(iters) = ', mean(extra$iters), '\n' )
    cat( 'mean(estim.precis) = ', mean(extra$estim.precis), '\n' )
    
    return(TRUE)
    }
    
    
testSourceInversion <- function()
    {
    wave = 400:700
    
    eye = resample( xyz1931.1nm, wave )

    ##------------    light sources     ----------##
    cat( "\n##------------     light sources     ----------##\n" )    
    
    spec    = illuminantE( 1, wave )
    
    spec    = bind( spec, resample( D65.1nm, wave ) )
    
    pspec   = planckSpectra( c(2100,3000,4000,5000,6000,7000,8000,9000,10000,15000), wave )
    spec    = bind( spec, pspec )
    
    XYZ     = product( spec, eye )
    
    XYZ     = rbind( XYZ, c(0.3,0.4,0.3) ) ;  rownames(XYZ)[nrow(XYZ)] = 'inside1'
    #XYZ     = rbind( XYZ, c(0.3,0.2,0.5)*200 ) ;  rownames(XYZ)[nrow(XYZ)] = 'inside2'      #  singular jacobian !!
    XYZ     = rbind( XYZ, c(0.2,0.4,0.4) ) ;  rownames(XYZ)[nrow(XYZ)] = 'inside3'      
    #XYZ     = rbind( XYZ, c(0.3,0.5,0.2) ) ;  rownames(XYZ)[nrow(XYZ)] = 'outside1'
    #XYZ     = rbind( XYZ, c(0.5,0.3,0.2) ) ;  rownames(XYZ)[nrow(XYZ)] = 'outside2'
    #XYZ     = rbind( XYZ, c(0.2,0.7,0.1) ) ;  rownames(XYZ)[nrow(XYZ)] = 'outside3'
    
    est.c   = invert( eye, XYZ )    
    if( is.null(est.c) ) return(FALSE)
        
    extra   = extradata(est.c)
    #  extra$response = NULL   # do not want to see it
    print( extra )
    if( any( is.na(extra$estim.precis) ) )    return(FALSE)

    cat( 'mean(iters) = ', mean(extra$iters), '\n' )
    cat( 'mean(estim.precis) = ', mean(extra$estim.precis), '\n' )
    
    return(TRUE)
    }
    

    
if( ! testReflectanceInversion() )  stop( "testReflectanceInversion() failed !" )

if( ! testSourceInversion() )  stop( "testSourceInversion() failed !" )



cat( "Passed all inversion tests !\n", file=stderr() )
