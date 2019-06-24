
library( colorSpec )


#   n       dimension of n-cube
#   sd      standard-deviation of smoothing filter
#   scale   bigger means more 0's and 1's
#    
#   returns a random point in the n-cube, with emphasis on lots of 0's and 1's    
#
#   to get no 0's and 1's set sd=0 and scale=1
randomSpectrum  <-  function( n, sd=5, scale=100 )
    {
    out = stats::runif( n )
    
    if( 0 < sd )
        {
        kern    = stats::dnorm( seq(-3*sd,3*sd,len=n/2+1), sd=sd )
        kern    = kern / sum(kern)
        #   out = out[ is.finite(out) ] not needed when circular        
        out = stats::filter( out, kern, circular=TRUE ) #; print(out)
        out = as.numeric(out)
        }

    out = pmin( pmax( scale*(out-0.5) + 0.5, 0 ), 1 )
    
    return( out )
    }
    

testBands <- function( n, samples, sd=5, scale=100, tol=5.e-12 )
    {
    set.seed(0)
    
    wave    = 400:(400+n-1)

    mat = matrix( NA_real_, n, samples )
    
    for( j in 1:samples )
        {
        mat[ ,j] = randomSpectrum( n, sd=sd, scale=scale )  # ; print( spec )
        }
        
    obj = colorSpec( mat, wave, quantity='transmittance', specnames=sprintf( "random-%d", 1:samples ) )
    
    #  spectra to a list of matrices
    lambdalist  = bandRepresentation( obj )
    if( is.null(lambdalist) )   return(FALSE)
    #   print( lambdalist )
    
    #  list of matrices back to spectra
    obj.back    = bandMaterial( lambdalist, wave )
    if( is.null(obj.back) )   return(FALSE)
    
        
    #   now compare the two
    transmax    = 0
    transmin    = Inf
    delta       = numeric(samples)      
    
    mat.back    = as.matrix( obj.back )
    
    for( j in 1:samples )        
        {
        spec    = mat[ ,j]
        
        spec.back   = mat.back[ ,j]
       
        delta[j]    = max( abs(spec - spec.back) )
        
        lambdamat   = lambdalist[[j]]
                
        if( tol < delta[j] )
            {
            cat( "----------------\n" )
            print( spec )
            print( lambdamat )
            print( spec.back )
            }
            
        transmin    = min( length(lambdamat), transmin )
        transmax    = max( length(lambdamat), transmax )
        }
    
    count   = sum( tol < delta )
    mess    = sprintf( "%d violations of %d samples (delta > %g).    max(delta)=%g", 
                                count, samples, tol, max(delta) )
    cat( mess, '\n' )
    mess    = sprintf( "transition count range: %d to %d.", transmin, transmax )
    cat( mess, '\n' )
    
    return( count == 0 )
    }
    

    
if( ! testBands( 50, 1000 ) )  stop( "testBands() failed !" )

cat( "Passed all bands tests !\n", file=stderr() )
