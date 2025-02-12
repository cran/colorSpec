
#   x           a colorSpec object, with type "responsivity.material".  The number of spectra must be 3.
#   gray        vector of numbers in (0,1) that define neutral gray on segment from black to white
#               this neutral gray point is the base of the probe ray
#   direction   a non-zero 3-vector defining the ray, or a matrix with 3 columns defining multiple rays
#   aux         return auxiliary performance data
#   spectral    return a colorSpec object with extradata()
#
#   return value
#   data.frame with a row for each direction and these columns:
#       gray        input
#       direction   input
#       s           position along the ray that intersects boundary
#       optimal     the optimal color on the boundary, where the ray intersects boundary
#
#   and if aux is TRUE
#       timetrace   time to trace the ray, in seconds
#       facetgens   # of generators of the zonogon facet that the ray intersects. Typically it is 2, which means the facet is a parallelogram.
#
#   and if spectral is TRUE
#       distance    the signed distance from optimal to the boundary; it should be very small.
#       transitions the number of transitions in the optimal spectrum; it is a non-negative even integer
#
#   in case of ERROR returns NULL

probeOptimalColors.colorSpec <- function( x, gray, direction, aux=FALSE, spectral=FALSE )
    {
    time_start  = as.double( Sys.time() )

    theName = deparse(substitute(x))

    spectra = numSpectra(x)

    #   ok  = spectra %in% 2L:3L

    if( spectra != 3 )
        {
        log_level( ERROR, "numSpectra(%s) = %d is invalid.  It must be 3.", theName, spectra )
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


    if( ! requireNamespace( 'zonohedra', quietly=TRUE ) )
        {
        log_level( WARN, "Required package 'zonohedra' could not be loaded; returning NULL."  )
        return(NULL)
        }


    direction   = prepareNxM( direction, M=spectra )
    if( is.null(direction) )    return(NULL)

    #   compute matrix W for zonotopes
    #   do not bother to optimize for regular wavelength sequence
    #wave    = wavelength(x)
    #   n       = length(wave)
    #bslist  = breakandstep( wave )
    #step    = bslist$stepvec
    #step    = diff( wave, lag=2 ) / 2
    #step    = c( wave[2]-wave[1], step, wave[n]-wave[n-1] )     #; print(step)
    #W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)

    #W   = t(W)
    
    #   zono    = zonohedra::zonohedron(W)
    
    zono    = zonotope( x )
    
    if( is.null(zono) )  return(NULL)
    
    white   = 2 * zonohedra::getcenter( zono )

    #   funlist = makeReparamFunctionList( x, theName, 'equalize' )
    #    if( is.null(funlist) ) return(NULL)

    #   step.wl = step.wl(x)


    out = NULL
    for( k in 1:length(gray) )
        {
        base    = gray[k] * white
        df      = zonohedra::raytrace( zono, base, direction, invert=spectral )     # zonohedra::raytrace() does all the real work
        if( is.null(df) )  return(NULL)

        #   add gray as initial column
        df  = cbind( gray=gray[k], df )

        #   append to bottom of out
        out = rbind( out, df )
        }


    if( aux )
        {
        #   add wanted columns

        #   facetgens is the number of generators of the zonogon facet, typically 2 meaning a trivial parallelogram facet
        out$facetgens   = lengths( zonohedra::gethyperplane( zonohedra::getsimplified(zono$matroid) )[ out$facetidx ] )
        }
    else
        {
        #   remove unwanted columns
        out$timetrace   = NULL
        #out$distance    = NULL      # present iff spectral is TRUE
        #out$transitions = NULL      # present iff spectral is TRUE
        }

    #   remove even more columns
    out$base        = NULL
    out$facetidx    = NULL
    out$sign        = NULL


    cnames  = colnames(out)

    #   rename 'tmax' to 's'
    k   = which( cnames == 'tmax' )
    if( length(k) == 1 )    colnames(out)[k] = 's'

    #   rename 'point' to 'optimal'
    k   = which( cnames == 'point' )
    if( length(k) == 1 )    colnames(out)[k] = 'optimal'


    rays        = nrow(out)
    failures    = sum( is.na( out$s ) )

    if( spectral )
        {
        #   the spectra we want are already in out$pcube, so just need to reorganize
        spectramat  = out$pcube
        out$pcube   = NULL          # erase from where it's not needed
        extra       = out
        specnames   = as.character( 1:nrow(extra) )
        out = colorSpec( t(spectramat), wavelength=wavelength(x), quantity='reflectance', organization='df.row', specnames=specnames )
        extradata(out)  = extra
        }


    if( TRUE )
        {
        time_elapsed    = as.double( Sys.time() ) - time_start

        log_level( INFO, "Processed %d rays in %g sec  (%g sec per ray)",
                            rays, time_elapsed, time_elapsed/rays )

        if( 0 < failures )
            log_level( WARN, "There were %d failures out of %d rays.\n", failures, rays )
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

    ok  = is.numeric(normal)   &&   (length(normal) %in% 2:3 )  &&  all( is.finite(normal) )  &&   ! all( normal==0 )
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

    m   = length(beta)

    ok  = is.numeric(beta)   &&   0 < m
    if( ! ok )
        {
        log_level( ERROR, "argument 'beta' is not numeric, with positive length."  )
        return(NULL)
        }
    dim(beta) = NULL



    x   = linearize(x)
    
    zono    = zonotope(x)
    
    if( is.null(zono) )  return(NULL)    
    
    #   white   = rowSums( zonohedra::getmatrix(zono) )

    names(normal)   = specnames(x)  #; print( names(normal) )

    dat = zonohedra::section( zono, normal, beta )

    if( is.null(dat) )  return(NULL)

    out = vector( m, mode='list' )

    for( k in 1:m )
        out[[k]]$beta       = beta[k]

    if( length(normal) == 3 )
        {
        for( k in 1:m )
            out[[k]]$section    = dat[[k]]$point
        }
    else if( length(normal) == 2 )
        {
        for( k in 1:m )
            {
            section     = matrix( 0, nrow=0, ncol=2 )

            if( all( is.finite(dat$boundary1[k, ]) ) )
                section = rbind( section, dat$boundary1[k, ] )

            if( all( is.finite(dat$boundary2[k, ]) ) )
                section = rbind( section, dat$boundary2[k, ] )

            out[[k]]$section    = section
            }
        }

    for( k in 1:m )
        colnames( out[[k]]$section )    = specnames(x)

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
    
    m   = nrow(lambda)
    
    if( m == 0 )
        {
        log_level( ERROR, "argument lambda is empty." )
        return(NULL)
        }
        
    x   = linearize(x)        

    wave    = wavelength(x)
    n       = length(wave)

    # make zonohedron with trivial ground set
    zono = zonotope( x, ground=1:n )
    if( is.null(zono) ) return(NULL)


    #   map from numeric lambda to indexes in the original matroid; this may generate NAs
    idxpair         = match( lambda, wave )
    dim(idxpair)    = dim(lambda)
    
    #   ground  = zonohedra::getground( zonohedra::getmatroid(zono) )
    
    #gndpair = gndsimp[ idxpair ]
    #dim(gndpair)    = dim(idxpair)


    #   in the next call, we can use idxpair, because the ground set is trivial
    data    = zonohedra::boundarypgramdata( zono, idxpair, cube=spectral )
    
    if( is.null(data) ) return(NULL)
    

    #   center  = zonohedra::getcenter( zono )
    
    rnames  = rownames(lambda)
    if( is.null(rnames) )   rnames = 1:m

    df  = data.frame( row.names=rnames )
    df$lambda       = lambda
    df$optimal      = data$center                   # +  matrix( center, nrow=m, ncol=3, byrow=TRUE )     # add the zonohedron center to the computed offsets
    df$transitions  = data$transitions

    if( spectral )
        {
        out = colorSpec( t(data$pcube), wavelength=wave, quantity='transmittance', organization='df.row', specnames=rownames(df) )
        extradata(out)  = df
        }
    else
        out = df

    count   = sum( is.na(data$hyperplaneidx) )
    if( 0 < count )
        {
        log_level( WARN, "%d of %d optimal colors could not be computed, because the pair of wavelengths is invalid.",
                            count, length(data$hyperplaneidx) )
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







plotOptimals3D.colorSpec <- function( x, ... )
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
        
    #if( is.finite(size) &&  0<size )
    #    {
    #    waverange    = range( wavelength(x) )
    #    x   = resample( x, wavelength = seq( waverange[1], waverange[2], len=size ) )
    #    }

    #wave    = wavelength(x)
    #step    = breakandstep(wave)$stepvec
    #W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)
    #zono    = zonohedra::zonohedron( t(W) )
    
    zono    = zonotope(x)
    
    if( is.null(zono) ) return(FALSE)

    out = zonohedra::plot.zonohedron( zono, ... )

    return( invisible(out) )
    }




plotOptimals2D.colorSpec <- function( x, ...  )
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


    #wave    = wavelength(x)
    #step    = breakandstep(wave)$stepvec
    #W       = step * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)
    #zono    = zonohedra::zonogon( t(W) )
    
    zono    = zonotope( x )
    
    if( is.null(zono) ) return(FALSE)

    out = zonohedra::plot.zonogon( zono, ... )

    title( main=theName )

    return( invisible(out) )
    }


insideOptimalColors.colorSpec  <-  function( x, p )
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
        
    p   = prepareNxM( p, M=3 )
    if( is.null(p) ) return(NULL)
    
    zono    = zonotope( x )
    if( is.null(zono) ) return(NULL)
    
    out = zonohedra::inside( zono, p )
    
    out$idxhyper    = NULL      # not meaninful here
    
    return( out )
    }




#   value   a dataframe with columns
#           response    the given response
#           lambda      the 2 transition wavelengths
#           ADL         computed ADL
#           omega       the reparameterized L in the interval [0,1]


computeADL.colorSpec <- function( x, response )
    {
    theName = deparse( substitute(x) )

    if( type(x) != "responsivity.material" )
        {
        log_level( ERROR, "The type is invalid.  type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(NULL)
        }

    if( numSpectra(x) != 3 )
        {
        log_level( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )
        return(NULL)
        }

    response    = prepareNxM( response, M=3 )
    if( is.null(response) ) return(NULL)
    
    zono    = zonotope( x, ground=1:length( wavelength(x) ) )

    if( is.null(zono) )  return(NULL)

    #wave    = wavelength(x)
    #bslist  = breakandstep( wave )
    #W       = bslist$stepvec * as.matrix( x )  #;  print(W)    # step   is replicated to all columns of as.matrix(x)
    #   rownames(W) = 1:nrow(W)     # erase the wavelengths, so they won't be used in the zonohedron
    #W   = t(W)
    #white   = rowSums( W ) #;  print( white )


    #   white   = rowSums( zonohedra::getmatrix(zono) )
    
    gray    = zonohedra::getcenter(zono)        # 0.5 * white  #;    print(gray)
    
    white   = 2 * gray

    direction   = response - matrix( gray, nrow(response), 3, byrow=T )


    data    = zonohedra::raytrace2trans( zono, gray, direction )

    if( is.null(data) ) return(NULL)

    #print( data$gndpair )
    #print( data$alpha )

    #   translate gdpair and alpha to wavelengths
    
    bslist  = breakandstep( wavelength(x) )

    valid           = is.finite( data$gndpair[ ,1]  )

    lambda  = matrix( NA_real_, nrow=nrow(data), ncol=2 )
    for( i in 1:nrow(data) )
        {
        if( ! valid[i] )    next

        k1  = data$gndpair[i,1]
        k2  = data$gndpair[i,2]

        a1  = data$alpha[i,1]
        a2  = data$alpha[i,2]

        #bp  = (k1 < k2)     # BandPass

        #if( ! bp )
        #    {
        #    #   bandstop, so reverse the sense
        #    a1 = 1 - a1 ;  a2 = 1 - a2
        #    }

        #   these equations are valid, for both bandpass and bandstop

        #   the 1st lambda
        lambda[i,1] = a1 * bslist$breakvec[k1]   +   (1 - a1) * bslist$breakvec[k1+1]

        #   the 2nd lambda
        lambda[i,2] = (1 - a2) * bslist$breakvec[k2]   +   a2 * bslist$breakvec[k2+1]
        }

    #print( lambda )

    #   compute dol[,]
    dol             = matrix( NA_real_, nrow(data), 3 )
    colnames(dol)   = c( 'delta', 'omega', 'lambda' )

    funlist         = makeReparamFunctionList( x, theName, 'equalize' )

    if( is.null(funlist) )
        return(NULL)


    omega1          = funlist$omega.from.lambda( lambda[valid,1] )
    omega2          = funlist$omega.from.lambda( lambda[valid,2] )
    bp              = omega1 < omega2       # BandPass

    dol[valid,1]    = ifelse( bp, omega2-omega1, 1 - (omega1-omega2) )
    dol[valid,2]    = ifelse( bp, 0.5*(omega1+omega2), ( 0.5*(omega1+omega2) + 0.5 ) %% 1 )
    dol[valid,3]    = funlist$lambda.from.omega( dol[valid,2] )


    ADL = cbind( 1/data$tmax, dol[ ,1], dol[ ,3] )

    #   override for the 3 "edge cases"
    for( i in 1:nrow(ADL) )
        {
        resp    = response[i, ]

        if( all( resp == 0 ) )
            {
            lambda[i, ] = NA_real_      # probably not necessary
            ADL[i, ]    = c( 1, 0, NA )
            dol[i,2]    = NA_real_
            }
        else if( all( resp == white ) )
            {
            lambda[i, ] = NA_real_      # probably not necessary
            ADL[i, ]    = c( 1, 1, NA )
            dol[i,2]    = NA_real_
            }
        else if( all( resp == gray ) )
            {
            lambda[i, ] = NA_real_      # probably necessary
            ADL[i, ]    = c( 0, NA, NA )
            dol[i,2]    = NA_real_
            }
        }


    rownames(ADL)   = NULL
    colnames(ADL)   = c('alpha','delta','lambda')

    colnames(response)  = toupper( specnames(x) )

    rnames  = rownames(response)
    if( is.null(rnames) )   rnames = 1:nrow(response)

    out             = data.frame( row.names=rnames )
    out$response    = response
    out$lambda      = lambda
    out$ADL         = ADL
    out$omega       = dol[ ,2]

    return( out )
    }





sectionSchrodingerColors.colorSpec <- function( x, normal, beta )
    {
    theName = deparse(substitute(x))

    if( type(x) != "responsivity.material" )
        {
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(NULL)
        }

    ok  = is.numeric(normal)   &&   (length(normal) == 3 )  &&  all( is.finite(normal) )  &&   ! all( normal==0 )
    if( ! ok )
        {
        log_level( ERROR, "argument 'normal' is not numeric, with length 3. Or else it is zero."  )
        return(NULL)
        }
    dim(normal) = NULL

    if( numSpectra(x) != length(normal) )
        {
        log_level( ERROR, "numSpectra(x) = %d  !=  %d = length(normal).", numSpectra(x), length(normal) )
        return(NULL)
        }

    m   = length(beta)

    ok  = is.numeric(beta)   &&   0 < m
    if( ! ok )
        {
        log_level( ERROR, "argument 'beta' is not numeric, with positive length."  )
        return(NULL)
        }
    dim(beta) = NULL


    x   = linearize(x)

    #   compute matrix W for zonotopes
    #   do not bother to optimize for regular wavelength sequence
    #wave    = wavelength(x)
    #step    = breakandstep(wave)$stepvec
    #W       = t( step * as.matrix( x ) ) #;  print(W)    # step   is replicated to all columns of as.matrix(x)
    #zono    = zonohedra::zonohedron( W )

    zono    = zonotope( x )
    
    if( is.null(zono) )  return(NULL)
    
    names(normal)   = specnames(x)  #; print( names(normal) )

    dat = zonohedra::section2trans( zono, normal, beta )

    if( is.null(dat) )  return(NULL)

    out = vector( m, mode='list' )

    for( k in 1:m )
        {
        out[[k]]$beta       = beta[k]
        out[[k]]$section    = dat[[k]]$point
        
        colnames( out[[k]]$section )    = specnames(x)        
        }

    return( invisible(out) )
    }


insideSchrodingerColors.colorSpec <- function( x, p )
    {
    theName = deparse(substitute(x))

    if( type(x) != "responsivity.material" )
        {
        log_level( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(NULL)
        }

    if( numSpectra(x) != 3 )
        {
        log_level( ERROR, "numSpectra(x) !=  3.", numSpectra(x) )
        return(NULL)
        }

    p   = prepareNxM( p, M=3 )
    if( is.null(p) ) return(NULL)
    
    zono    = zonotope( x )
    if( is.null(zono) ) return(NULL)
    
    out = zonohedra::inside2trans( zono, p )
    
    return( out )
    }
    
    
    
    
#--------       UseMethod() calls           --------------#

probeOptimalColors <- function( x, gray, direction, aux=FALSE, spectral=FALSE )
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
    
insideOptimalColors <- function( x, p )
    {
    UseMethod("insideOptimalColors")
    }    
    
plotOptimals3D <- function( x, ... )
    {
    UseMethod("plotOptimals3D")
    }

plotOptimals2D <- function( x, ... )
    {
    UseMethod("plotOptimals2D")
    }

    
    
    
    
computeADL <- function( x, response )
    {
    UseMethod("computeADL")
    }

sectionSchrodingerColors <- function( x, normal, beta )
    {
    UseMethod("sectionSchrodingerColors")
    }

insideSchrodingerColors <- function( x, p )
    {
    UseMethod("insideSchrodingerColors")
    }




