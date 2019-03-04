
#   x           a colorSpec object, with type "responsivity.material".  The number of spectra must be 3.
#   gray        vector of numbers in (0,1) that define neutral gray on segment from black to white
#               this neutral gray point is the base of the probe ray
#   direction   a non-zero 3-vector defining the ray, or a matrix with 3 columns defining multiple rays
#   tol         convergence tolerance. Iterations continue until root-finding error < .tol.  
#   aux         return auxiliary performance data
#   spectral    return a colorSpec object with extradata()
#
#   return value
#   data.frame with a row for each direction and these columns:
#       gray        input
#       direction   input
#       s           position along the ray that intersects boundary
#       optimal     the optimal color on the boundary, where the ray intersects boundary
#       lambda      lambda.1 and lambda.2 of the ideal material producing the optimal color
#                       lambda.1 < lambda.2 => bandpass
#                       lambda.1 > lambda.2 => bandstop
#       dol         delta and omega, the Logvinenko parameters - analogous to latitude and longitude, plus lambda corresponding to omega
#
#
#   in case of ERROR returns NULL

probeOptimalColors.colorSpec <- function( x, gray, direction, aux=FALSE, spectral=FALSE, tol=1.e-6 )
    {
    theName = deparse(substitute(x))
    
    spectra = numSpectra(x)
    
    ok  = spectra %in% 2L:3L
    
    if( ! ok )
        {
        log.string( ERROR, "numSpectra(%s) = %d is invalid.", theName, spectra )     
        return(NULL)
        }        

    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )        
        return(NULL)
        }    
        
    ok  = is.numeric(gray)  &&  0<length(gray)
    if( ! ok )
        {
        log.string( ERROR, "Argument 'gray' is not a numeric vector with positive length." )
        return(NULL)
        }
                
    mask    =  0<gray & gray<1 
    if( ! all(mask) )
        {
        log.string( ERROR, "gray = %g is invalid.", gray[ which(!mask)[1] ] )
        return(NULL)
        }    
    
    direction   = prepareNxM( direction, M=spectra )
    if( is.null(direction) )    return(NULL)
 
    #   compute matrix W for zonotopes
    #   do not bother to optimize for regular wavelength sequence
    wave    = wavelength(x)
    #   n       = length(wave)    
    step    = breakandstep(wave)$stepvec
    #step    = diff( wave, lag=2 ) / 2
    #step    = c( wave[2]-wave[1], step, wave[n]-wave[n-1] )     #; print(step)
    W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    white   = colSums( W )
    
    #   funlist = makeReparamFunctionList( x, theName, 'equalize' )
    #    if( is.null(funlist) ) return(NULL)

    #   step.wl = step.wl(x)
   
    if( spectra == 3 )
        {
        zono    = zonohedron( W )
        }
    else if( spectra == 2 )
        {
        zono    = zonogon( W )          #;   print( zono )
        }

    if( is.null(zono) )  return(NULL)    

    out = NULL
    for( k in 1:length(gray) )
        {
        base    = gray[k] * white
        df      = raytrace( zono, base, direction )     # this does all the real work        
        if( is.null(df) )  return(NULL)

        #   add gray as initial column
        df  = cbind( gray=gray[k], df )
        
        #   append to bottom
        out = rbind( out, df )
        }
        
    # print( df )        

    #   compute the optimal spectra
    #   these spectra are all 0-1,
    #   except for a single alpha in 2D
    #   and a pair of alphas in 2D
    spectramat = invertboundary( zono, out )        #; print( str(spectramat) )
    if( is.null(spectramat) )   return(NULL)        
    
    lambda  = matrix( NA_real_, nrow(spectramat), 2 )    
    for( i in 1:nrow(spectramat) )
        lambda[i, ]  = compute_lambdas( spectramat[i, ], wave, step )
    
    out$lambda  = lambda
    
    #   remove columns no longer needed: 'base',  'sign' 
    out$base    = NULL          
    out$sign    = NULL
    out$idx     = NULL   
    
    if( ! aux )
        {
        #   remove even more columns
        out$alpha       = NULL         
        out$timetrace   = NULL      
        out$faces       = NULL   
        out$tested      = NULL        
        }
    
    
    cnames  = colnames(out)
    
    #   rename 'tmax' to 's'
    k   = which( cnames == 'tmax' )
    if( length(k) == 1 )    colnames(out)[k] = 's'
    
    #   rename 'boundary' to 'optimal'
    k   = which( cnames == 'boundary' )
    if( length(k) == 1 )    colnames(out)[k] = 'optimal'
    
    #   rename 'faces' to 'parallelograms'
    k   = which( cnames == 'faces' )
    if( length(k) == 1 )    colnames(out)[k] = 'parallelograms'
    
    #   add column dol. 
    #   funlist uses splinefun() with 'monoH.FC', which does not handle NAs properly
    funlist         = makeReparamFunctionList( x, theName, 'equalize' )
    valid           = is.finite( out$lambda[ ,1]  )
    omega1          = funlist$omega.from.lambda( out$lambda[valid,1] )
    omega2          = funlist$omega.from.lambda( out$lambda[valid,2] )
    bp              = omega1 < omega2       # BandPass
    dol             = matrix( NA_real_, nrow(out), 3 )
    colnames(dol)   = c( 'delta', 'omega', 'lambda' )
    dol[valid,1]    = ifelse( bp, omega2-omega1, 1 - (omega1-omega2) )
    dol[valid,2]    = ifelse( bp, 0.5*(omega1+omega2), ( 0.5*(omega1+omega2) + 0.5 ) %% 1 )
    dol[valid,3]    = funlist$lambda.from.omega( dol[valid,2] )
    out$dol         = dol
    
    

    if( FALSE )
        {
        time_elapsed    = as.double( Sys.time() ) - time_start
        log.string( INFO, "Processed %d rays in %g sec  (%g sec per rays)",
                            rays, time_elapsed, time_elapsed/rays ) 

        failures    = sum( is.na( out$s ) )
        if( 0 < failures )
            log.string( WARN, "There were %d failures out of %d rays.\n", failures, nrow(out) )
        
        if( aux )
            out = cbind( out, df.aux )
        }
    
    if( spectral )
        {
        extra   = out
        specnames   = as.character( 1:nrow(extra) )
        out = colorSpec( t(spectramat), wavelength=wavelength(x), quantity='reflectance', organization='df.row', specnames=specnames )
        extradata(out)  = extra
        }
    
    return(out)
    }
    
    
    
sectionOptimalColors.colorSpec <- function( x, normal, beta )
    {
    theName = deparse(substitute(x))
    
    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )        
        return(NULL)
        }        
    
    ok  = is.numeric(normal)   &&   (length(normal) %in% 2L:3L )  &&  all( is.finite(normal) )  &&   ! all( normal==0 )
    if( ! ok )
        {
        log.string( ERROR, "argument 'normal' is not numeric, with length 2 or 3. Or else it is zero."  )     
        return(NULL)
        }   
    dim(normal) = NULL

    if( numSpectra(x) != length(normal) )
        {
        log.string( ERROR, "numSpectra(x) = %d  !=  %d = length(normal).", numSpectra(x), length(normal) )
        return(NULL)
        }   
    
    ok  = is.numeric(beta)   &&   0<length(beta)
    if( ! ok )
        {
        log.string( ERROR, "argument 'beta' is not numeric, with positive length."  )
        return(NULL)
        }   
    dim(beta) = NULL
    
    
    #   compute matrix W for zonotopes
    #   do not bother to optimize for regular wavelength sequence
    wave    = wavelength(x)
    #   n       = length(wave)    
    step    = breakandstep(wave)$stepvec
    W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    white   = colSums( W )
    
    #   funlist = makeReparamFunctionList( x, theName, 'equalize' )
    #    if( is.null(funlist) ) return(NULL)

    #   step.wl = step.wl(x)
   
    if( length(normal) == 3 )
        {
        zono    = zonohedron( W )
        }
    else if( length(normal) == 2 )
        {
        zono    = zonogon( W )          #;   print( zono )
        }

    if( is.null(zono) )  return(NULL)   

    names(normal)   = specnames(x)  #; print( names(normal) )

    out = section( zono, normal, beta )
    if( is.null(out) )  return(NULL)

    return( invisible(out) )
    }
    
    
    
windingNumber <- function( iPoint, iX, iY )
    {
    n = length(iX)  #;    assert( length(iY)==n )
    
    #   look for intersections of edges with ray from iPoint
    
    x = iPoint[1]
    y = iPoint[2]
    
    winding = 0
    
    
    for( i in 1:n )
        {
        i_next = (i %% n) + 1

        y_diff = iY[i] - y
        
        y_diff_next = iY[i_next] - y
        
        if( 0 <= y_diff * y_diff_next ) next
        
        #  find intersection of edge and horizontal line through iPoint
        t = (y - iY[i]) / (iY[i_next] - iY[i])
        
        x_int = iX[i] + t * (iX[i_next] - iX[i])
        
        if( x < x_int )
            winding = winding + sign(y_diff_next)
        }
    
    winding
    }
    
    
 
               
               
       

plotOptimals3D.colorSpec <- function( x, size=50, type='w', both=TRUE )
    {
    theName = deparse( substitute(x) )
    
    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(FALSE)
        }        
    
    if( numSpectra(x) != 3 )
        {
        log.string( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )
        return(FALSE)
        }     
        
    typefull    = c( 'w', 'p' )
    
    idx = pmatch( type, typefull )
    if( is.na(idx) )
        {
        log.string( ERROR, "type='%s' is invalid", type )
        return(FALSE)
        }     
    type    = typefull[idx]
        
    if( is.finite(size) &&  0<size )
        {
        waverange    = range( wavelength(x) )
        x   = resample( x, wavelength = seq( waverange[1], waverange[2], len=size ) )
        }
        
    wave    = wavelength(x)
    step    = breakandstep(wave)$stepvec
    W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    
    zono    = zonohedron( W )
    if( is.null(zono) ) return(FALSE)
    
    plot( zono, type=type, both=both )
    
    return( invisible(TRUE) )
    }
    
    


plotOptimals2D.colorSpec <- function( x  )
    {
    theName = deparse( substitute(x) )
    
    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(FALSE)
        }        
    
    if( numSpectra(x) != 2 )
        {
        log.string( ERROR, "numSpectra(%s) = %d != 2", theName, numSpectra(x) )
        return(FALSE)
        }     
        
    wave    = wavelength(x)
    step    = breakandstep(wave)$stepvec
    W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    zono    = zonogon( W )
    if( is.null(zono) ) return(FALSE)
    
    plot( zono )
    
    title( main=theName )
    
    return( invisible(TRUE) )
    }
    
    



    
    
computeADL.colorSpec <- function( x, response )
    {
    theName = deparse( substitute(x) )
    
    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(NULL)
        }        
    
    if( numSpectra(x) != 3 )
        {
        log.string( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )
        return(NULL)
        }     
        
    response    = prepareNxM( response, M=3 )
    if( is.null(response) ) return(NULL)
           
    stepvec = breakandstep( wavelength(x) )$stepvec
    
    response.white  = colSums( stepvec * as.matrix(x) )     # vector stepvec is replicated to all columns
    
    gray    = 0.5 * response.white  #;    print(gray)
    
    direction   = response - matrix( gray, nrow(response), 3, byrow=T )
    
    data    = probeOptimalColors( x, 0.5, direction )
    
    if( is.null(data) ) return(NULL)
       
    ADL = cbind( 1/data$s, data$dol[ ,1], data$dol[ ,3] )
    rownames(ADL)   = NULL
    colnames(ADL)   = c('alpha','delta','lambda')
    #   class(ADL)      = "model.matrix"    
    
    colnames(response)  = toupper( specnames(x) )
    #   class(response)     = "model.matrix"
    
    rnames  = rownames(response)
    if( is.null(rnames) )   rnames = 1:nrow(response)
    
    out             = data.frame( row.names=rnames ) 
    out$response    = response
    out$ADL         = ADL
    out$omega       = data$dol[ ,2]
    out$lambda      = data$lambda
    
    #   out = data.frame( response=response, ADL=ADL, omega=data$dol[ ,2], lambda=data$lambda, row.names=rownames(response) )  # as.data.frame.model.matrix
    
    return( out )
    }
    
    
#--------       UseMethod() calls           --------------#            
              
probeOptimalColors <- function( x, gray, direction, aux=FALSE, spectral=FALSE, tol=1.e-6  )
    {
    UseMethod("probeOptimalColors")
    }
        
sectionOptimalColors <- function( x, normal, beta )
    {
    UseMethod("sectionOptimalColors")
    }
                  
computeADL <- function( x, response )
    {
    UseMethod("computeADL")
    }
 
 
plotOptimals3D <- function( x, size=50, type='s', both=TRUE )
    {
    UseMethod("plotOptimals3D")
    }
                 
plotOptimals2D <- function( x )
    {
    UseMethod("plotOptimals2D")
    }
                     
                 
                 