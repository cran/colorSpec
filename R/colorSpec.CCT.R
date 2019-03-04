
    
        
#   x      a colorSpec object with type 'light', and M spectra
#
#   returns an M-vector, containing CCT in Kelvin, or NA in case of error
computeCCT.colorSpec   <- function( x, isotherms='robertson', locus='robertson', strict=FALSE )
    {
    n   = numSpectra( x )
    if( n == 0 )    return( numeric(0) )

    out = rep(NA_real_,n) 
    names(out)  = specnames(x)
    
    if( type(x) != 'light' )
        {
        log.string( WARN, "The type of of '%s' is '%s', but it must be 'light'.",
                    deparse(substitute(x))[1], type(x) )
        return(out)
        }
        
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )
        {
        log.string( ERROR, "Required package 'spacesXYZ' could not be imported."  )
        return(out)
        }     
    
    XYZ = product( radiometric(x), colorSpec::xyz1931.1nm, wavelength='auto' )

    out = spacesXYZ::CCTfromXYZ( XYZ, isotherms=isotherms, locus=locus, strict=strict  )    

    return( out )
    }
        
            


   

#--------       UseMethod() calls           --------------#            
        
computeCCT <- function( x, isotherms='robertson', locus='robertson', strict=FALSE  )
    {
    UseMethod("computeCCT")
    }
    