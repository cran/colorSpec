
    
        
#   x      a colorSpec object with type 'light'   
#
#   returns CCT in Kelvin, or NA 
computeCCT.colorSpec   <- function( x )
    {
    n   = numSpectra( x )
    if( n == 0 )    return( numeric(0) )

    out = rep( NA_real_, n )
    names(out)  = specnames(x)
    
    if( type(x) != 'light' )
        {
        log.string( WARN, "The type of of '%s' is '%s', but it must be 'light'.",
                    deparse(substitute(x)), type(x) )
        return( out )
        }
    
    XYZ = product( radiometric(x), colorSpec::xyz1931.1nm, wavelength='auto' )

    for( k in 1:n )
        {
        out[k] = CCTfromXYZ( XYZ[k, ] )
        }
    
    return( out )
    }
        
            

    
#   XYZ     a 3-vector with XYZ
#   returns CCT, or NA if outside valid range
#
#   requires private data frame dataCCT,  which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()

CCTfromXYZ    <- function( XYZ )
    {
    if( any(is.na(XYZ)) )   return( NA_real_ )
    
    if( length(XYZ) != 3 )  return( NA_real_ )
    
    denom = XYZ[1]  +  15 * XYZ[2]  +  3 * XYZ[3]
    
    if( denom < 1.e-16 )    return( NA_real_ )
    us  = 4 * XYZ[1] / denom
    vs  = 6 * XYZ[2] / denom
    
    di = (vs - dataCCT$v) - dataCCT$t * (us - dataCCT$u)
    #   print( di )
    
    n = length(di)
    
    #   di should be decreasing, and with 1 zero crossing
    i   = which( di[1:(n-1)] * di[2:n] <= 0 ) + 1
   
    if( length(i) != 1  ||  i==1 )  return( NA_real_ )
    
    d0  = di[i]   / sqrt( 1 + dataCCT$t[i]^2 )
    dm  = di[i-1] / sqrt( 1 + dataCCT$t[i-1]^2 )
    p   = dm / (dm - d0)
    
    out = 1.e6 / ( (1-p)*dataCCT$mired[i-1]  +  p*dataCCT$mired[i] )
    
    return( out )
    }

#--------       UseMethod() calls           --------------#            
        
        
computeCCT <- function( x )
    {
    UseMethod("computeCCT")
    }
    