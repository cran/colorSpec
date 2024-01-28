
#   x           a (light or material) responder, to be modified to match y as closely as possible
#   y           the reference responder
#   filter      use a filter to modify x
#   matrix      use a matrix to modify x
#   preserve    if not NULL, then a light source whose coordinates must be preserved *exactly* (up to numerical truncation)
#               for this option, the number of spectra in x and y must be equal
#
#   x and y must have the same wavelength vector
#   and if preserve is a light source, then its wavelength vector must be the same as x and y
#
#   the value is a colorSpec object that is 'close' to y.  It is equal to product( filter, multiply( x, matrix ) )
#       the attribute "emulate" is attached, which is a list with items "filter" and "matrix"
#       filter          a colorSpec object with quantity='transmittance'.  Present when filter is TRUE
#       matrix          a numerical matrix.  Present when matrix is TRUE.


emulate.colorSpec  <-  function( x, y, filter=FALSE, matrix=TRUE )
    {
    if( ! requireNamespace( 'MASS', quietly=TRUE ) )
        {
        log_string( ERROR, "Required package 'MASS' could not be imported."  )
        return(NULL)
        }


    ok  = identical( wavelength(x), wavelength(y) )
    if( ! ok )
        {
        log_string( ERROR, "wavelengths of x and y are not identical." )
        return(NULL)
        }

    type.x  = type(x)
    type.y  = type(y)

    valid   = (type.x == 'responsivity.light'  ||  type.x == 'responsivity.material' )
    if( ! valid )
        {
        log_string( ERROR, "type(x) is '%s', which is invalid.\n", type.x )
        return(NULL)
        }

    if( type.y != type.x )
        {
        log_string( ERROR, "type(y) = '%s', which is not the same as type(x) = '%s'.\n",
                                    type.y, type.x )
        return(NULL)
        }

    M   = numSpectra(x)
    N   = numSpectra(y)

    if( M==0  ||  N==0 )
        {
        log_string( ERROR, "One of x or y has 0 spectra." )
        return(NULL)
        }

    if( ! filter  &&  ! matrix )
        {
        log_string( WARN, "Both 'filter' or 'matrix' are FALSE; returning x unmodified." )
        return(x)
        }


    #   ensure that *both* are radiometric
    x   = radiometric( x, warn=TRUE )
    y   = radiometric( y, warn=TRUE )

    wave = wavelength(x)

    X   = as.matrix( x )
    Y   = as.matrix( y )

    if( filter )
        {
        #   test X to see if this is doable
        test    = sqrt( rowSums(X*X) )
        mask    = (test < 1.e-4 * max(test))
        if( any(mask) )
            {
            log_object( ERROR, wave[mask] )
            log_string( ERROR, "x spectral values are too small when filter=TRUE, at the above %d wavelengths.  Try restricting the domain.",
                                sum(mask) )
            return(NULL)
            }
        }


    if( filter  &&  matrix )
        {
        if( N == 1  &&  2 <= M )
            {
            log_string( ERROR, "M=%d  and  N=%d; the equation to solve is under-determined. Try setting filter=FALSE.", M, N )
            return(NULL)
            }

        time_start  = as.double( Sys.time() )

        res = solveDAX( X, Y )
        if( is.null(res) )  return(NULL)


        #   rescale so that max(res$d) == 1
        dmax    = max( res$d )
        res$d   = res$d / dmax
        A       = dmax * res$X

        theFilter  = colorSpec( res$d, wavelength=wave, quantity='transmittance' )

        specnames(theFilter)    = 'filter'

        rownames( A )  = colnames(X)

        out = product( theFilter, multiply( x, A ) )

        attr( out, 'emulate' )  = list( filter=theFilter, A=A )

        time_elapsed    = as.double( Sys.time() ) - time_start

        log_string( TRACE, "Computed filter and matrix after %d iterations and %g seconds.",
                        res$iterations, time_elapsed )
        }
    else if( filter )
        {
        if( M != N )
            {
            log_string( ERROR, "matrix is FALSE and numSpectra of x and y are not equal; %d != %d.", M, N )
            return(NULL)
            }

        if( N == 1 )
            {
            #   a special case, X has already been "vetted"
            #   Y is a column vector, and so is X; because M==N
            #   the emulation is an exact match
            d   =  Y / X
            colnames(d) = "filter"
            }
        else
            {
            d   = solveDiagonal( X, Y )
            if( is.null(d) )    return(NULL)
            }

        ran = range(d)

        if( ran[1] < 0 )
            log_string( WARN, "The computed filter is not realizable, min(transmittance) = %g < 0.", ran[1] )

        if( 1 < ran[2] )
            log_string( WARN, "The computed filter is not realizable, max(transmittance) = %g > 1.", ran[2] )

        theFilter  = colorSpec( d, wavelength=wave, quantity='transmittance' )

        out = product( theFilter, x )

        attr( out, 'emulate' )  = list( filter=theFilter )
        }
    else if( matrix )
        {
        A  = MASS::ginv(X) %*% Y
        rownames( A )  = colnames(X)

        out = multiply( x, A )

        attr( out, 'emulate' )  = list( A=A )
        }


    #   append '.em' to spectrum names
    specnames( out )    = sprintf( "%s.em", specnames(y) )

    attr( out, "sequence" ) = NULL

    quantity( out )  = quantity( y )

    if( 0 )
    {
    #   out = multiply( out, t(res$X) )
    #   attr( out, "matchResponder" )   = list( d=res$d, X=res$X )

    attr( out$x.modified, "sequence")  = NULL


    #   compute residual for each channel in Y
    resid   = as.matrix( out$x.modified )  -  Y
    resid   = resid*resid
    resid   = colSums(resid)
    out$resid       = resid
    out$resid.sum   = sqrt( sum(resid) )
    out$resid       = sqrt(resid)
    }


    return( out )
    }


#   A   mxn matrix
#   Y   mxp matrix#
#
#   solve for D and X in
#           DAX = Y
#   where D is diagonal.  The best solution in the LS sense.
#
#   returns a list with vector d, and matrix X

solveDAX  <-  function( A, Y, iters=200, reltol=5.e-8 )
    {
    if( ! requireNamespace( 'MASS', quietly=TRUE ) )
        {
        log_string( ERROR, "Required package 'MASS' could not be imported."  )
        return(NULL)
        }

    m   = nrow(A)
    if( nrow(Y) != m )
        {
        log_string( ERROR, "Bad nrow(Y) = %d != %d.", nrow(Y), m )
        return(NULL)
        }

    n   = ncol(A)
    p   = ncol(Y)

    if( n == 1  &&  p == 1 )
        {
        #   special case, exact solution, no iterations required
        out = list( d=Y/A, X=as.matrix(1,1,1), iterations=0 )
        return( out )
        }

    #   initialize
    d   = rep( 1, m )

    resid_prev = Inf

    for( rep in 1:iters )
        {
        #   solve for X, for fixed D
        mat = matrix( d, nrow=m, ncol=n )
        B   = mat * A
        X   = MASS::ginv(B) %*% Y   #   ; print( X )

        #   compute resid
        resid   = B %*% X - Y
        resid   = sqrt( sum(resid*resid) )

        #if( rep == 1 || FALSE )  cat( "rep=", rep,  "    ", "resid=", resid, "\n" )

        bad =  reltol * (resid + reltol) < (resid_prev - resid)

        if( ! bad )    break

        #   solve for D, for fixed X
        d   = solveDiagonal( A %*% X, Y )
        if( is.null(d) )    return(NULL)

        resid_prev  = resid
        }

    if( rep == iters )
        {
        log_string( ERROR, "Failed to converge after %d iterations.", iters )
        return(NULL)
        }

    #cat( "rep=", rep,  "    ", "resid=", resid, "\n" )

    out = list( d=d, X=X, iterations=rep )

    return( out )
    }



#   B   mxp matrix
#   Y   mxp matrix#
#
#   solve for D in
#           DB = Y
#   where D is diagonal.  The best solution is in the LS sense.

solveDiagonal  <-  function( B, Y )
    {
    #   print( B )
    BB      = rowSums( B*B )

    test    = sqrt( BB )
    if( any( test < 1.e-6*max(test) ) )
        {
        log_string( ERROR, "Matrix B is deficient.  min(sqrt(BB))=%g, max(sqrt(BB))=%g.",
                            min(test), max(test) )
        return(NULL)
        }

    out = rowSums( Y*B ) / BB

    return(out)
    }


#--------       UseMethod() calls           --------------#

emulate <- function(  x, y, filter=FALSE, matrix=TRUE )
    {
    UseMethod("emulate")
    }

