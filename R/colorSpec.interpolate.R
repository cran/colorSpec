        
#   x       a colorSpec object, typically with multiple spectra
#   p       a numeric vector with length(p) == numSpectra(x)
#   pout    a numeric vector of parameter values at which interpolation of the spectra in x take place
#   pname   name of the parameter variable (p & pout).
#
#   Value
#           a colorSpec object with a spectrum for each value in pout
#           this has organization 'df.row' with 1 column of extradata() with name = pname

interpolate.colorSpec   <- function( x, p, pout, pname=deparse(substitute(p)) )
    {
    out = NULL

    ok  = is.numeric(p) && is.numeric(pout)
    if( ! ok )
        {
        log_string( ERROR, "One of argument p or pout is not numeric." )
        return( out )
        }    
        
    if( length(p) != numSpectra(x) )
        {
        log_string( ERROR, "length mismatch.  length(p)=%d  !=  %d=numSpectra(x)", length(p), numSpectra(x) )
        return( out )
        }    
        
    ok  = is.character(pname)  &&  length(pname)==1  #; print( ok )
    if( ! ok )
        {
        log_string( ERROR, "Argument pname must be a character vector with length 1." )
        return( out )
        }   
    
    myfun <- function( v )
        {       
        if( length(v) != length(p) )
            {
            log_string( FATAL, "%d != %d", length(v), length(p) )   # internal error
            return( out )
            }    
        return( stats::spline( p, v, xout=pout, method="natural" )$y )
        }
    
    mat = apply( as.matrix(x), 1, myfun )     #; print( str(mat) )
    
    if( is.null(dim(mat)) )
        #   this happens iff length(pout)==1
        dim(mat)    = c( length(mat), 1 )   
    else
        mat = t(mat)

    #   cname   = deparse( substitute(p) ) ; print(cname)
    #   print( sprintf( "%s%g", pname, pout ) )
    colnames(mat)   = sprintf( "%s=%g", pname, pout )       #; print( colnames(mat) )
    
    out = colorSpec( mat, wavelength=wavelength(x), quantity=quantity(x), organization='df.row' )
    
    #   add 1 column of extradata
    extra           = data.frame( pout=pout )
    colnames(extra) = pname
    extradata(out)  = extra
    
    return(out)
    }
    
    

        
#--------       UseMethod() calls           --------------#            
        
interpolate   <- function( x, p, pout, pname=deparse(substitute(p)) )        
    {
    UseMethod("interpolate")
    }

    