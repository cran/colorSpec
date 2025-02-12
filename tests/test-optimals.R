
library( colorSpec )
options( width=180 )

testProbe <- function( tol=5.e-8 )
    {
    wave    = seq(400,700,by=5)
    
    D50.eye = product( D50.5nm, 'material', xyz1931.1nm, wavelength=wave )

    #   make a few random rectangular spectra
    set.seed( 0 )
    
    count   = 50
    alpha   = runif( count, min=0.01, max=0.99 )    
    #  alpha   = rep( 1, count )
    lambda  = matrix( runif(count*2,min=min(wave),max=max(wave)), count, 2 )
    colnames(lambda) = c('lambda.1','lambda.2')
    
    print( cbind(lambda=lambda, alpha=alpha) )
    
    rectspec = rectangularMaterial( lambda, alpha, wavelength=wave )
    
    #   compute XYZ
    XYZ = product( rectspec, D50.eye )  #; print(XYZ)
    
    #white.XYZ   = product( neutralMaterial(1,wavelength=wave), D50.eye )  #; print( white.XYZ/2 )
    
    #direction   = XYZ - matrix( white.XYZ/2, count, 3, byrow=TRUE )
    
    res = computeADL( D50.eye, XYZ )
    
    delta   = rowSums( abs(lambda - res$lambda) )
    
    print( cbind(res,delta=delta) )
    
    cat( sprintf( "max(abs(delta)) = %g\n", max(delta,na.rm=TRUE) ) )
    
    failures    = sum( is.na( delta ) )
    cat( sprintf( "inversion failures: %d of %d\n", failures, count ) )

    deltamax    = max( delta, na.rm=TRUE )  
    cat( sprintf( "deltamax = %g    tol = %g\n", deltamax, tol ) )
    
    return( failures==0  &&  deltamax<=tol )
    }
    

    
if( ! testProbe() )  stop( "testProbe() failed !" )

cat( "Passed all optimal color tests !\n", file=stderr() )
