
    
        
#   x      a colorSpec object with type 'light', and M spectra
#
#   returns an M-vector, containing CCT in Kelvin, or NA in case of error
computeCCT.colorSpec   <- function( x, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    n   = numSpectra( x )
    if( n == 0 )    return( numeric(0) )

    if( type(x) != 'light' )
        {
        log.string( WARN, "The type of of '%s' is '%s', but it must be 'light'.",
                    deparse(substitute(x))[1], type(x) )
        return( out )
        }
    
    XYZ = product( radiometric(x), colorSpec::xyz1931.1nm, wavelength='auto' )

    out = CCTfromXYZ( XYZ, method=method, strict=strict, c2=c2 )    

    return( out )
    }
        
            

    
#   XYZ     an Nx3 matrix, or a vector that can be converted to such a matrix
#   returns CCT, or NA if outside valid range

CCTfromXYZ    <- function( XYZ, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    XYZ = prepareNxM( XYZ, M=3 )
    if( is.null(XYZ) )   return(NULL)
    
    n   = nrow(XYZ)
    out = rep( NA_real_, n )
    names(out)  = rownames(XYZ)   
    
    if( pmatch( tolower(method), c('mccamy'), nomatch=0 ) != 0 )
        {
        # C. S. McCamy "Correlated color temperature as an explicit function of chromaticity coordinates" Color Research & Application Volume 17, Issue 2, pages 142-144, April 1992
        xyY = xyY_from_XYZ( XYZ )
        for( i in 1:n )
            out[i] = CCTfromxy_McCamy( xyY[i,1:2] )
        }
    else
        {
        for( i in 1:n )
            {
            uv      = uv_from_XYZ( XYZ[i, ] )
            out[i]  = CCTfromuv_single( uv, method=method, strict=strict, c2=c2 ) 
            }
        }
            
    return( out )
    }


#   xy  an Nx2 matrix, or a vector that can be converted to such a matrix
    
CCTfromxy  <- function( xy, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    xy  = prepareNxM( xy, M=2 )
    if( is.null(xy) )   return(NULL)
    
    n   = nrow(xy)
    out = rep( NA_real_, n )
    names(out)  = rownames(xy)    
    
    if( pmatch( tolower(method), c('mccamy'), nomatch=0 ) != 0 )
        {
        for( i in 1:n )
            out[i]  = CCTfromxy_McCamy( xy[i, ] )
        }
    else
        {
        uv  = uv_from_xy( xy, .year=1960 )        
        for( i in 1:n )
            out[i]  = CCTfromuv_single( uv[i, ], method=method, strict=strict, c2=c2 )
        }
        
        
    return(out)
    }

    
    
# C. S. McCamy "Correlated color temperature as an explicit function of chromaticity coordinates" Color Research & Application Volume 17, Issue 2, pages 142-144, April 1992

CCTfromxy_McCamy  <- function( xy )
    {
    n = (xy[1] - 0.3320)/(xy[2] - 0.1858)
    
    return( -449*n^3 + 3525*n^2 - 6823.3*n + 5520.33 )
    }
    
    
#   uv  an Nx2 matrix, or a vector that can be converted to such a matrix
    
CCTfromuv  <- function( uv, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    uv  = prepareNxM( uv, M=2 )
    if( is.null(uv) )   return(NULL)
    
    n   = nrow(uv)
    out = rep( NA_real_, n )
    names(out)  = rownames(uv)
        
    if( pmatch( tolower(method), c('mccamy'), nomatch=0 ) != 0 )
        {
        for( i in 1:n )
            {
            #   McCamy, so convert from uv to xy.  It is unlikely that a user will do this.
            denom   = uv[i,1] - 4*uv[i,2] + 2
            
            if( is.na(denom) ||  denom < 1.e-16 )
                {
                # log.string( WARN, "uv=%g,%g is invalid, because denom=%g.",   uv[1], uv[2], denom )
                }          
            else
                {
                xy      = c( 1.5*uv[i,1], uv[i,2] ) / denom            
                out[i]  = CCTfromxy_McCamy( xy )
                }
            }
        }
    else
        {        
        for( i in 1:n )
            out[i]  = CCTfromuv_single( uv[i, ], method=method, strict=strict, c2=c2 )
        }
        
        
    return(out)
    }


    
#   uv      a 2-vector with uv 1960
#   returns CCT, or NA if outside valid range
#
#   requires private data frame dataCCT,  which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
    
CCTfromuv_single  <- function( uv, method='robertson', strict=FALSE, c2=1.4388e7 )
    {
    if( any( is.na(uv) ) )  return(NA_real_)
    
    idx.method  = pmatch( tolower(method), c('robertson','lm'), nomatch=0 )
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
    
    if( is.na(denom)  ||  denom < 1.e-16 )
        {
        # log.string( WARN, "XYZ=%g,%g,%g is invalid, because denom=%g.",  XYZ[1], XYZ[2], XYZ[3], denom )
        return( c(NA_real_,NA_real_) )
        }    
    
    return( c(4 * XYZ[1], 6 * XYZ[2]) / denom )
    }
    
    

#--------       UseMethod() calls           --------------#            
        
computeCCT <- function( x, method='robertson', strict=FALSE, c2=1.4388e7  )
    {
    UseMethod("computeCCT")
    }
    