

#   construct a zonotope from a material responder
#
#   x       the colorSpec object
#   ground  positive integers in increasing order with length=#rows in matrix
#           ground=NULL means to use the rownames of as.matrix(x), which are the wavelengths
#   ...     other args such as tolerances, passed to the constructor in package zonohedra

zonotope.colorSpec  <-  function( x, ground=NULL, ... )
    {
    if( ! requireNamespace( 'zonohedra', quietly=TRUE ) )
        {
        log_level( WARN, "Required package 'zonohedra' could not be loaded; returning NULL."  )
        return(NULL)
        }    

    if( type(x) != "responsivity.material" )
        {
        theName = deparse(substitute(x))
            
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(NULL)
        }

    n   = numSpectra(x)
    if( !( n %in% 1:3) )
        {
        log_level( ERROR, "n = %d  is invalid.", n )
        return(NULL)
        }
    
    wave    = wavelength(x)
    stepvec = breakandstep(wave)$stepvec
    
    W       = t( stepvec * as.matrix( x ) ) #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    #   ground  = 1:ncol(W)     ground moved to an argument
    
    if( n == 3 )
        out = zonohedra::zonohedron( W, ground=ground, ... )
    else if( n == 2 )
        out = zonohedra::zonogon( W, ground=ground, ... )
    else if( n == 1 )
        out = zonohedra::zonoseg( W, ground=ground, ... )       #   probably will never be used
        
    return( out )
    }
    

#--------       UseMethod() calls           --------------#

zonotope <- function( x, ground=NULL, ... )
    {
    UseMethod("zonotope")
    }    
    