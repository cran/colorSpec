
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
        log_level( ERROR, "numSpectra(%s) = %d is invalid.", theName, spectra )     
        return(NULL)
        }        

    if( type(x) != "responsivity.material" )
        {
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )        
        return(NULL)
        }    
        
    ok  = is.numeric(gray)  &&  0<length(gray)
    if( ! ok )
        {
        log_level( ERROR, "Argument 'gray' is not a numeric vector with positive length." )
        return(NULL)
        }
                
    mask    =  0<gray & gray<1 
    if( ! all(mask) )
        {
        log_level( ERROR, "gray = %g is invalid.", gray[ which(!mask)[1] ] )
        return(NULL)
        }    
    
    direction   = prepareNxM( direction, M=spectra )
    if( is.null(direction) )    return(NULL)
 
    #   compute matrix W for zonotopes
    #   do not bother to optimize for regular wavelength sequence
    wave    = wavelength(x)
    #   n       = length(wave)    
    bslist  = breakandstep( wave )    
    step    = bslist$stepvec
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
    #   except for a single alpha in 1D, and a pair of alphas in 2D
    spectramat = invertboundary( zono, out$boundary, out )$source     #; print( str(spectramat) )
    if( is.null(spectramat) )   return(NULL)        
    
    lambda  = matrix( NA_real_, nrow(spectramat), 2 )    
    for( i in 1:nrow(spectramat) )
        {
        # lambda[i, ]  = compute_lambdas( spectramat[i, ], wave, step )        
        
        lambdamat   = lambdasFromSpectrum( spectramat[i, ], wave, bslist )
        if( ! is.null(lambdamat)  &&  nrow(lambdamat)==1 )
            lambda[i, ] = lambdamat
        }
        
        
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
        log_level( INFO, "Processed %d rays in %g sec  (%g sec per rays)",
                            rays, time_elapsed, time_elapsed/rays ) 

        failures    = sum( is.na( out$s ) )
        if( 0 < failures )
            log_level( WARN, "There were %d failures out of %d rays.\n", failures, nrow(out) )
        
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
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )        
        return(NULL)
        }        
    
    ok  = is.numeric(normal)   &&   (length(normal) %in% 2L:3L )  &&  all( is.finite(normal) )  &&   ! all( normal==0 )
    if( ! ok )
        {
        log_level( ERROR, "argument 'normal' is not numeric, with length 2 or 3. Or else it is zero."  )     
        return(NULL)
        }   
    dim(normal) = NULL

    if( numSpectra(x) != length(normal) )
        {
        log_level( ERROR, "numSpectra(x) = %d  !=  %d = length(normal).", numSpectra(x), length(normal) )
        return(NULL)
        }   
    
    ok  = is.numeric(beta)   &&   0<length(beta)
    if( ! ok )
        {
        log_level( ERROR, "argument 'beta' is not numeric, with positive length."  )
        return(NULL)
        }   
    dim(beta) = NULL
    
    
    x   = linearize(x)    
    
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
    
    
#   x           a colorSpec object with type "responsivity.material"
#   lambda      Mx2 matrix of wavelengths of x, each row defines a parallelogram on the boundary of the zonohedron
#   spectral    logical, whether to return just a data.frame, or a colorSpec object of type "material"
#               the spectra take value 1/2 at the 2 given wavelengths, and 0 or 1 elsewhere
    
canonicalOptimalColors.colorSpec <- function( x, lambda, spectral=FALSE )
    {
    theName = deparse(substitute(x))
    
    ok  = (numSpectra(x) == 3L)
    if( ! ok )
        {
        log_level( ERROR, "numSpectra(%s) = %d  != 3.", theName, numSpectra(x) )     
        return(NULL)
        }        

    if( type(x) != "responsivity.material" )
        {
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )        
        return(NULL)
        }    
        
    lambda  = prepareNxM( lambda, M=2 )
    if( is.null(lambda) )   return(NULL)
    
                
    wave    = wavelength(x)                
    mask    = lambda %in% wave
    ok  = all( mask )
    if( ! ok )
        {
        log_level( ERROR, "In argument lambda, %d of %d wavelengths are not in the wavelength sequence of '%s'.", 
                        sum(!mask), length(lambda), theName )     
        return(NULL)
        }        
        
    x   = linearize(x)    

    
    #   compute matrix W for zonotopes
    #   do not bother to optimize for regular wavelength sequence

    step    = breakandstep(wave)$stepvec
    W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    #   white   = colSums( W )
    
    #   the tolerance here is for collinearity, which can be larger than the face normal differences
    condlist    = condenseGenerators( W, tol=5.e-7 )
    if( is.null(condlist) )    return(NULL)
    
    ncond   = nrow(condlist$Wcond)    
    if( ncond < 3 )
        {
        log_level( ERROR, "Invalid number of condensed generators = %d < 3.", ncond )
        return(NULL)
        }
        
    n       = nrow(W)
    m       = nrow(lambda)
    
    #   optimal             = matrix( NA_real_, m, 3 )
    #   colnames(optimal)   = specnames(x)
    spectra     = matrix( NA_real_, m, n )
    transitions = rep( NA_integer_, m )
    
    tol2    = 5.e-10
    
    for( i in 1:m )
        {
        #   find the row indexes in W
        idxpair = match( lambda[i, ], wave )    
        if( any(is.na(idxpair)) )
            {
            log_level( FATAL, "bad logic for lambda=%g,%g.", lambda[i,1], lambda[i,2] )
            return(NULL)
            }
        
        #   print( idxpair ) ; print( condlist$groupidx[ idxpair ]  )
        
        #   find the row indexes in Wcond, in increasing order        
        kdxpair = sort( condlist$groupidx[ idxpair ],  na.last=TRUE )  
        
        if( any(is.na(kdxpair)) )
            #   the responsivity at one of these wavelengths is 0, so the parallelogram is degenerate
            next
            
        if( kdxpair[1] == kdxpair[2] )
            #   indexes are equal, so the parallelogram is degenerate
            next            

        #   print( kdxpair )
        
        #   collinear generators in W have already been identified, 
        #   so the parallelogram is non-degenerate and the normal is non-zero
        normal  = crossproduct( condlist$Wcond[kdxpair[1], ],  condlist$Wcond[kdxpair[2], ] )
        normal  = normal / sqrt( sum(normal*normal) )   # unitize
        #   print( normal )
        
        functional          = as.numeric( condlist$Wcond %*% normal )
        names(functional)   = rownames(condlist$Wcond)
        
        #   make useful logical masks
        between = kdxpair[1]<(1:ncond)  &  (1:ncond)<kdxpair[2]
        outside = ! between
        outside[ kdxpair ]  = FALSE
        
        #   coplanar generators certainly include kdxpair[1] and kdxpair[2], but maybe others
        coplanar    = abs(functional) < tol2
        
        #   all values in pcube are 0,1, or 1/2
        pcube   = (sign(functional) + 1)/2
        pcube[ kdxpair ] = 0.5  # override
        # print( pcube )
        
        if( sum(between) < sum(outside) )
            {
            #   count the number of 0s and 1s outside the band
            mask    = outside  &  !coplanar
            zeros   = sum( pcube[mask] == 0 )
            ones    = sum( pcube[mask] == 1 )
            bandstop    = (zeros < ones)    # this means that pcube is more like bandstop than a bandpass
            }
        else
            {
            #   count the number of 0s and 1s between the wavelengths (inside the band)
            mask    = between  &  !coplanar
            zeros   = sum( pcube[mask] == 0 )
            ones    = sum( pcube[mask] == 1 )
            bandstop    = (ones < zeros)    # this means that pcube is more like bandstop than a bandpass
            }
            
        if( bandstop )  
            # flip it so the spectrum is more like a bandpass
            pcube = 1 - pcube
            
        #   change the coplanars outside to 0
        pcube[ coplanar & outside ] = 0
        
        #   change the coplanars between to 1
        pcube[ coplanar & between ] = 1
                
        if( lambda[i,2] < lambda[i,1] )
            #   change bandpass to bandstop
            pcube   = 1 - pcube
                
        #   only 2 values in pcube are 1/2, and the rest are 0 or 1
        #  print( pcube )
        
        spectra[i, ]    = expandcanonical( condlist, idxpair, pcube  )
            
        transitions[i]  = counttransitions( spectra[i, ] )
        }
        
    rnames  = rownames(lambda)
    if( is.null(rnames) )   rnames = 1:m
    
    df  = data.frame( row.names=rnames )
    df$lambda       = lambda
    df$optimal      = spectra %*% W
    df$transitions  = transitions
    
    if( spectral )
        {
        out = colorSpec( t(spectra), wavelength=wave, quantity='transmittance', organization='df.row', specnames=rownames(df) )
        extradata(out)  = df
        }
    else
        out = df
        
        
    count   = sum( is.na(transitions) )
    if( 0 < count )
        {
        log_level( WARN, "%d of %d colors could not be computed, because the pair of wavelengths is invalid.",
                            count, length(transitions) )
        }
    
    return(out)
    }
    
    
    

#   x           a colorSpec object with type "responsivity.material"
#   lambda      2 wavelengths from x
plotfunctional <- function( x, lambda, gamma )
    {    
    wave        = wavelength(x)                
    
    idxpair   = match( lambda, wave )
    if( length(lambda)!=2  ||  any( is.na(idxpair) ) )
        {
        log_level( ERROR, "lambda='%s' is invalid.", as.character(lambda) )
        return(FALSE)
        }
    
    #   compute matrix W for zonotopes
    #   do not bother to optimize for regular wavelength sequence

    step    = breakandstep(wave)$stepvec
    W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    #   white   = colSums( W )
    
    #   the tolerance here is for collinearity, which can be larger than the face normal differences
    condlist    = condenseGenerators( W, tol=5.e-7 )
    if( is.null(condlist) )    return(FALSE)
    
    ncond   = nrow(condlist$Wcond)    
    if( ncond < 3 )
        {
        log_level( ERROR, "Invalid number of condensed generators = %d < 3.", ncond )
        return(FALSE)
        }

    kdxpair = sort( condlist$groupidx[ idxpair ],  na.last=TRUE )  
    
    if( any(is.na(kdxpair)) )
        #   the responsivity at one of these wavelengths is 0, so the parallelogram is degenerate
        return(FALSE)
        
    if( kdxpair[1] == kdxpair[2] )
        #   indexes are equal, so the parallelogram is degenerate
        return(FALSE)

    #   print( kdxpair )
    
    #   collinear generators in W have already been identified, 
    #   so the parallelogram is non-degenerate and the normal is non-zero
    normal  = crossproduct( condlist$Wcond[kdxpair[1], ],  condlist$Wcond[kdxpair[2], ] )
    normal  = normal / sqrt( sum(normal*normal) )   # unitize
    #   print( normal )
    
    functional  = as.numeric( condlist$Wcond %*% normal )
        
    functional[kdxpair] = 0     # make near 0 exactly 0
    
    if( nrow(condlist$Wcond) < nrow(W) )
        {
        #   expand functional
        funsaved    = functional
        functional  = rep( NA_real_, nrow(W) )
        for( k in 1:length(funsaved) )
            {
            group   = condlist$group[[k]]
            functional[ group ]   = funsaved[k]        
            }
        }
        
    names(functional)   = as.character(wave)
    # print( functional )
        
    #   ready to plot
    y   = powodd(functional,1/gamma)
    
    xlim    = range(wave)
    ylim    = range(y,na.rm=TRUE)
    ylab    = sprintf( "functional^(1/%g)", gamma )
    plot( xlim, ylim, type='n', xlab="wavelength (nm)", ylab=ylab )
    grid( lty=1 )
    abline( h=0 )
    points( wave, y )
    points( wave[idxpair], y[idxpair], pch=20 )
        
    return(TRUE)
    }
    
    
    
    
    
#   condlist    a list with W, Wcond, group, groupidx
#   idxpair     defining canonical optimal, idxpair[1] != idxpair[2] 
#   pcube       point in the condensed cube, thought of as a reflectance spectrum
#               all values are 0 or 1, except at condlist$groupidx[ idxpair ] where the value is 0.5

expandcanonical <- function( condlist, idxpair, pcube )
    {
    if( length(pcube) != nrow(condlist$Wcond) )
        {
        log_level( FATAL, "mismatch %d != %d", length(pcube), nrow(condlist$Wcond) )
        return(NULL)
        }
   
    n   = nrow(condlist$W)  #; print(n)
    if( nrow(condlist$Wcond) == n )    return( pcube )     # no condensation

    #   1<= kdxpair[1] < kdxpair[2] <= nrow(Wcond)        
    kdxpair = condlist$groupidx[ idxpair ]

    pdrop   = pcube[ -kdxpair ] 
    ok  = all( pdrop==0  |  pdrop==1 )
    if( ! ok )
        {
        log_level( FATAL, "pcube is invalid, because values not at %d and %d are not 0 or 1.", 
                                kdxpair[1], kdxpair[2] )
        return(NULL)
        }
        
    ok  = all( pcube[kdxpair] == 0.5 )
    if( ! ok )
        {
        log_level( FATAL, "pcube is invalid, because values at %d and %d are not 1/2.", 
                                kdxpair[1], kdxpair[2] )
        return(NULL)
        }

    knext   = c( 2:length(pcube), 1 )
    kprev   = c( length(pcube), 1:(length(pcube)-1) )
        
    out = numeric( n )
    
    
    #   examine the non-trivial groups
    for( k in 1:length(pcube) )
        {
        group   = condlist$group[[k]]
        alpha   = pcube[k]        
        
        if( length(group)==1  ||  alpha==0  ||  alpha==1 )
            {
            #   trivial cases
            out[ group ]   = alpha
            next
            }
       
        #   condlist$group[[k]] is a nontrivial splitting
        #   and since alpha is neither 0 nor 1, k should be one of kdxpair[1] or kdxpair[2]
        #   check this
        j   = match( k, kdxpair )
        if( is.na(j) )
            {
            log_level( FATAL, "k=%d is not in the pair (%d,%d).", k, kdxpair[1], kdxpair[2] )
            return(NULL)
            }
            
        #   another check
        ok  = idxpair[j] %in% group
        if( ! ok )
            {
            log_level( FATAL, "idxpair[%d]=%d is not in the proper group '%s'.", 
                                    j, idxpair[j], paste( group, collapse=' ' ) )
            return(NULL)
            }
            
        m           = length(group)        
        alphasplit  = numeric(m)            
        alphasplit[ group == idxpair[j] ]    = 0.5
            
        #   choose output form to minimize the number of transitions
        if( pcube[kprev[k]]  >  pcube[knext[k]] )
            #   output form is 1111...1/2...00000...
            alphasplit[ group < idxpair[j] ]    = 1
        else
            #   output form is 0000...1/2....11111...            
            alphasplit[ idxpair[j] < group ]  = 1            
            
        #   print( alphasplit )
            
        out[ group ]    = alphasplit
        }
        
    #    now examine the 0-generators.  
    #   They can be assigned any coefficient, but assign to minimize # of transitions.
    idx     = which( is.na(condlist$groupidx) )
    inext   = c( 2:n, 1 )
    iprev   = c( n, 1:(n-1) )
    for( i in idx )
        {
        if( out[iprev[i]]==1  ||  out[inext[i]]==1 )
            out[i]  = 1
        }
        
    #   print( out )

    return( out )
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
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(FALSE)
        }        
    
    if( numSpectra(x) != 3 )
        {
        log_level( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )
        return(FALSE)
        }     
        
    typefull    = c( 'w', 'p' )
    
    idx = pmatch( type, typefull )
    if( is.na(idx) )
        {
        log_level( ERROR, "type='%s' is invalid", type )
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
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(FALSE)
        }        
    
    if( numSpectra(x) != 2 )
        {
        log_level( ERROR, "numSpectra(%s) = %d != 2", theName, numSpectra(x) )
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
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(NULL)
        }        
    
    if( numSpectra(x) != 3 )
        {
        log_level( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )
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
    
canonicalOptimalColors <- function( x, lambda, spectral=FALSE )
    {
    UseMethod("canonicalOptimalColors")
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
                     
                 
                 