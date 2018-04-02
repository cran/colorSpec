
    
        
#   x      a colorSpec object with type 'light'   
#
#   returns CCT in Kelvin, or NA 
computeCCT.colorSpec   <- function( x, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    n   = numSpectra( x )
    if( n == 0 )    return( numeric(0) )

    out = rep( NA_real_, n )
    names(out)  = specnames(x)
    
    if( type(x) != 'light' )
        {
        log.string( WARN, "The type of of '%s' is '%s', but it must be 'light'.",
                    deparse(substitute(x))[1], type(x) )
        return( out )
        }
    
    XYZ = product( radiometric(x), colorSpec::xyz1931.1nm, wavelength='auto' )

    for( k in 1:n )
        {
        out[k] = CCTfromXYZ( XYZ[k, ], method=method, strict=strict, c2=c2 )
        }
    
    return( out )
    }
        
            

    
#   XYZ     a 3-vector with XYZ
#   returns CCT, or NA if outside valid range
#
#   requires private data frame dataCCT,  which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()

CCTfromXYZ    <- function( XYZ, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    if( any(is.na(XYZ)) )   return( NA_real_ )
    
    if( length(XYZ) != 3 )  return( NA_real_ )
    
    uv  = uv_from_XYZ( XYZ )
    
    if( is.na(uv[1]) )  return( uv[1] )
            
    return( CCTfromuv( uv, method=method, strict=strict, c2=c2 ) )
    }
    
    
CCTfromuv  <- function( uv, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    idx.method  = pmatch( method, c("robertson","lm"), nomatch=0 )
    if( idx.method == 0 )
        {
        log.string( ERROR, "method='%s' is invalid.", as.character(method) )
        return(  NA_real_ )
        }

    di = (uv[2] - dataCCT$v) - dataCCT$t * (uv[1] - dataCCT$u)
    #   print( di )
    
    n = length(di)
    
    #   di should be decreasing, and with 1 zero crossing
    i   = which( di[1:(n-1)] * di[2:n] <= 0 ) + 1
   
    if( length(i) != 1  ||  i==1 )
        {
        log.string( WARN, "For uv=%g,%g  cannot find unique zero crossing.", uv[1], uv[2] )
        return( NA_real_ )
        }
        
    d0  = di[i]   / sqrt( 1 + dataCCT$t[i]^2 )
    dm  = di[i-1] / sqrt( 1 + dataCCT$t[i-1]^2 )
    p   = dm / (dm - d0)
    
    mired   = (1-p)*dataCCT$mired[i-1]  +  p*dataCCT$mired[i]
    
    resid   = (1-p) * c(dataCCT$u[i-1],dataCCT$v[i-1])  +  p * c(dataCCT$u[i],dataCCT$v[i])  - uv
    
    if( idx.method == 2 )
        {
        #   use Robertson's mired as starting point for LM
        if( ! requireNamespace( 'minpack.lm', quietly=TRUE ) )
            {
            log.string( WARN, "method='lm' is unavailable, because package 'minpack.lm' cannot be loaded." )
            return(  NA_real_ )
            }
            
        wave    = wavelength( colorSpec::xyz1931.1nm )
        
        myresid <- function( mir )
            {
            XYZ = product( planckSpectra( 1.e6 / mir, wavelength=wave, normalize=FALSE, c2=c2 ), colorSpec::xyz1931.1nm )
            
            return( uv_from_XYZ( XYZ ) - uv )
            }
            
        res = minpack.lm::nls.lm( mired, lower=NULL, upper=NULL, fn=myresid ) #; print( str(res) )
        
        ok  = (1 <= res$info) && (res$info <= 4)
        if( ! ok )
            {
            log.string( WARN, "For uv=%g,%g Levenberg-Marquardt did not converge. info=%d", uv[1], uv[2], res$info )
            return( NA_real_ )
            }
            
        log.string( INFO, "LM polished mired from mired=%g to %g in %d iterations.", mired, res$par, res$niter )
        
        mired   = res$par    
        resid   = res$fvec
        }
        
    if( strict )
        {
        test    = sqrt( sum(resid*resid) )
        if( 0.05 < test )
            {
            log.string( WARN, "uv=%g,%g is invalid, because its distance to the Planckian locus = %g > 0.05.  (mired=%g)",
                                uv[1], uv[2], test, mired )
            return( NA_real_ )
            }
        }        
        
    out = 1.e6 / mired    
    
    return( out )
    }
    
uv_from_XYZ <- function( XYZ )
    {
    denom = XYZ[1]  +  15 * XYZ[2]  +  3 * XYZ[3]
    
    if( denom < 1.e-16 )
        {
        log.string( WARN, "XYZ=%g,%g,%g is invalid, because denom=%g.",
                            XYZ[1], XYZ[2], XYZ[3], denom )
        return( NA_real_ )
        }    
    
    return( c(4 * XYZ[1], 6 * XYZ[2]) / denom )
    }

#--------       UseMethod() calls           --------------#            
        
computeCCT <- function( x, method='robertson', strict=FALSE, c2=1.4388e7  )
    {
    UseMethod("computeCCT")
    }
    