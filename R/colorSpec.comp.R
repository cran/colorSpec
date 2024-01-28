

#   temperature     a vector of N temperatures, in Kelvin
#   wavelength      at which to sample, in nm
#   normalize       logical to normalize to 1 at 560nm

#
#   returns colorSpec object with N spectra
#           in case of an individual temperature error, the spectrum is all NAs
#           in case of an argument error, the function returns NULL
#
#   if normalize is FALSE, then the output units are W*M^{-2}*nm^{-1}
#   see W&S pp. 12-13 and Figure 1(1.2.2)

planckSpectra <-  function( temperature, wavelength=300:830, normalize=TRUE, c2=1.4388e-2 )
    {
    #   fundamental constants - CODATA 2014
    h   = 6.62607004e-34    # Planck's Contant in, J*s or W*s^2
    c   = 2.99792458e8      # speed of light, in m/s
    k   = 1.38064852e-23    # Boltzmann constant, in J/K
    
    if( ! normalize )
        {
        c1  = 2*pi*h*c^2 * 1.e36    # the 1.e36 converts from W*nm^{-2}*nm^{-1} to W*M^{-2}*nm^{-1}    1.e18 nm^2 / M^2
        #   print( c1 )
        }
    
    if( pmatch( c2, 'calculate', nomatch=FALSE ) )
        c2  = h*c/k

    if( !(is.numeric(c2)  &&  length(c2)==1) )
        {
        log_string( ERROR, "c2 = '%s' is invalid.", as.character(c2) )
        return(NULL)
        }
        
    #   print( c2 )
    
    #   convert from (m*K) to (nm*K)
    c2nm    = c2 * 1.e9

    #   waveM   = wavelength * 1.e-9   #   convert from nm to m
        
    lambda0 = 560   # e-9
    
    pow0    = lambda0 ^ -5 
    pow1    = wavelength ^ -5 
        
    mat = matrix( NA_real_, length(wavelength), length(temperature) )
    
    for( j in 1:length(temperature) )
        {
        T   = temperature[j]
        
        if( is.na(T) || T<0 )   next    # leave full of NAs
        
        if( T == 0 )
            {
            if( ! normalize )   mat[ ,j] = 0
            next
            }
            
        if( T == Inf )
            {
            if( normalize )     mat[ ,j] = (wavelength / lambda0) ^ -4
            next
            }
        
        Power   = pow1 / ( exp( c2nm / (wavelength * T) ) - 1 )
    
        if( normalize )
            {
            #   normalize to 1 at 560 nm
            pow00   = pow0 / ( exp( c2nm / (lambda0 * T) )  -  1 )
            Power = Power / pow00
            }
        else
            # scale by c1 
            Power = c1 * Power        
            
        mat[ ,j]    = Power 
        }
        
    colnames(mat)   = sprintf( "P%g", round(temperature)  )    
    
    out = colorSpec( mat, wavelength, quantity="energy", organization="matrix" )

    metadata(out)   = list( description="Planck black-body power density" )
    
    return( out )
    }
    


erythemalSpectrum <- function( wavelength=250:400 )
    {
    mask.S  = wavelength <= 298
    mask.M  = wavelength <= 328
    mask.L  = wavelength <= 400
    
    mask.L  = mask.L  &  ! mask.M
    mask.M  = mask.M  &  ! mask.S
    
    out = numeric( length(wavelength) )
    
    out[ mask.S ]   = 1
    out[ mask.M ]   = 10 ^ (0.094*(298 - wavelength[mask.M] ) )
    out[ mask.L ]   = 10 ^ (0.015*(139 - wavelength[mask.L] ) )
    
    out = colorSpec( out, wavelength, quantity="energy->action" )
    
    specnames( out )    = "erythemal"
    
    metadata( out ) = list( description="A.F. McKinlay and B.L. Diffey (1987)" )
    
    return( out )
    }
    
    
illuminantE <- function( energy=1, wavelength=380:780 )
    {    
    n   = length(wavelength)
    
    theNames = sprintf( "E%g", energy )
    
    if( length(energy) == 1 )
        core    = rep(energy,n)
    else
        {
        core            = matrix( energy, n, length(energy), byrow=T )
        colnames(core)  = theNames
        }
        
    out = colorSpec( core, wavelength, quantity="energy" )
    
    if( length(energy) == 1 )    specnames( out ) = theNames
    
    metadata( out ) = list( description="Equal Energy White" )
        
    return( out )
    }
    
neutralMaterial <- function( gray=1, wavelength=380:780 )
    {    
    n   = length(wavelength)

    core    = matrix( gray, n, length(gray), byrow=T )
    
    colnames(core)  = sprintf( "Neutral%g", gray )

    out = colorSpec( core, wavelength, quantity="reflectance" )
    
    if( length(gray) == 1 )   organization( out ) = 'vector'
           
    metadata( out ) = list( description="Neutral Gray" )
                 
    return( out )
    }    
    
rectangularMaterial <- function( lambda, alpha=1, wavelength=380:780 )
    {     
    lambda  = prepareNxM( lambda, M=2 )
    if( is.null(lambda) )   return(NULL)
    
    if( ! is.numeric(alpha) )
        {
        log_string( ERROR, "alpha is invalid, because it is not numeric." )
        return(NULL)
        }
    
    n   = nrow(lambda)
    if( length(alpha) == 1 )    alpha = rep( alpha, n )
    
    if( length(alpha) != n )
        {
        log_string( ERROR, "length(alpha) = %d != %d = nrow(lambda).", length(alpha), nrow(lambda) )
        return(NULL)
        }
    
    if( ! all( -1<=alpha  &  alpha<=1 ) )
        {
        log_string( ERROR, "alpha is invalid; all values must be in [-1,1]." )
        return(NULL)
        }
        
        
    p   = length(wavelength) 
    
    #   compute the breaks and steps
    #   breaks are endpoints of the bins, and steps are widths of the bins
    bslist  = breakandstep( wavelength )
    lower   = bslist$breakvec[ 1:p ]        # lower breaks
    upper   = bslist$breakvec[ 2:(p+1) ]    # upper breaks
    stepvec = bslist$stepvec                # actually equal to upper - lower
    
    #print( lower )
    #print( upper )
    
    mat = matrix( NA_real_, p, n )
    
    gray    = 0.5 
    grayvec = rep( gray, p )
    
    for( i in 1:n )
        {
        if( alpha[i] == 0 )
            {
            #   special case; use this shortcut to save some time
            mat[ ,i]    = grayvec
            next
            }

        lam = lambda[i, ]
        
        if( lam[1] == lam[2] )
            {
            #   special case - undefined
            log_string( WARN, "material is undefined because lambda_min==lambda_max == %g.", lam[1] )
            next
            }

        #lambda_min  = min(lam)
        #lambda_max  = max(lam)
        

        
        #   compute spectrum for lam[1]
        spec1   = grayvec
        spec1[ upper <= lam[1] ] = 0
        spec1[ lam[1] <= lower ] = 1
        #print( spec1 )
        #   there might be 1 bin that needs interpolation
        k   = which( spec1 == gray )
        spec1[k]    = (upper[k] - lam[1]) / stepvec[k]
            
        #   compute spectrum for lam[2] in exactly the same way
        spec2   = grayvec
        spec2[ upper <= lam[2] ] = 0
        spec2[ lam[2] <= lower ] = 1
        #print( spec2 )
        #   there might be 1 bin that needs interpolation
        k   = which( spec2 == gray )
        spec2[k]    = (upper[k] - lam[2]) / stepvec[k]
            
        y   = spec1 - spec2
        
        if( lam[2] < lam[1] )   y   = y + 1     # raise it up
        
        #   y   = 0.25 * (sign(wavelength - lambda_min) + 1) * (sign(lambda_max - wavelength) + 1)
        
        a   = alpha[i]          # ifelse( lam[1] < lam[2], alpha[i], -alpha[i] )

        mat[ ,i]    = a*y + (1-a)*0.5
        }
        
    namevec = rownames(lambda)
    if( is.null(namevec) )
        {
        #   make suitable names
        namevec = character(n)
        for( i in 1:n )
            {
            lam = lambda[i, ]
            
            theDiff = diff( mat[ ,i] )
            pos     = any( 0 < theDiff )
            neg     = any( theDiff < 0 )
            
            if( any(is.na(theDiff)) )
                filter  = 'NA'
            else if( pos && neg )
                filter  = ifelse( lam[1] < lam[2], 'BP', 'BS' )
            else if( neg )
                filter  = 'SP'
            else if( pos )
                filter  = 'LP'
            else
                filter  = 'N'
             
            namevec[i] = sprintf( '%s_[%g,%g]', filter, lam[1], lam[2] )
            }
        }
        
    colnames(mat) = make.unique( namevec )
        
    out = colorSpec( mat, wavelength=wavelength, quantity='transmittance', organization='df.row' )
    
    extra   = data.frame( row.names=1:n )
    extra$lambda    = lambda
    extra$alpha     = alpha
    extradata(out)  = extra  
    
    return( out )
    }
    
    
#   lambda  an Mx2 matrix defining waveband intervals, or a list of such matrices
bandMaterial <- function( lambda, wavelength=380:780 )
    {     
    if( is.numeric(lambda) )
        lambdalist  = list(lambda=lambda)
    else if( is.list(lambda)  &&  0<length(lambda) )
        lambdalist  = lambda
    else
        {
        log_string( ERROR, "Argument 'lambda' is invalid." )
        return(NULL)
        }
        
    n   = length(lambdalist)        
    for( j in 1:n )
        {
        lambdalist[[j]] = prepareNxM( lambdalist[[j]], M=2, Nmin=0 )
        if( is.null(lambdalist[[j]]) )   return(NULL)
        }
        
        
    bslist  = breakandstep(wavelength)
    if( is.null(bslist) )   return(NULL)
    

    #   check that every item in the list is a valid matrix
    valid   = sapply( lambdalist, validLambdaMat, wavelength, bslist )
    idx     = which( ! valid )
    if( 0 < length(idx) )
        {
        log_string( ERROR, "List item %d is invalid.", idx[1] )
        return(NULL)
        }
        

    mat = matrix( NA_real_, length(wavelength), n )
    
    for( j in 1:n )
        {
        mat[ ,j]    = spectrumFromLambdas( lambdalist[[j]], wavelength, bslist )
        }
    
    specnames   = names(lambdalist)
    if( is.null(specnames) )  specnames = sprintf( "bandmat_%d", 1:n )
    
    colnames(mat) = specnames
        
    out = colorSpec( mat, wavelength=wavelength, quantity='transmittance', organization='matrix' )
    
    return( out )
    }
    


    
bandRepresentation.colorSpec  <-  function( x )
    {
    if( type(x) != "material" )
        {
        log_string( ERROR, "type(x)='%s' is invalid.", type(x) )
        return(NULL)
        }
    
    m       = numSpectra(x)

    out = vector( m, mode='list' )
    names(out)  = specnames(x)
    
    if( m == 0 )    return(out)
    
    mat     = as.matrix(x)
    
    wave    = wavelength(x)
    bslist  = breakandstep(wave)    
    
    for( k in 1:m )
        {
        out[[k]]    = lambdasFromSpectrum( mat[ ,k], wave, bslist )
        }
    
    return(out)
    }
    
    
#   lambdamat   Mx2 matrix defining M closed intervals of wavelengths, with at most being a "wrap-around"    
#               It checks that all intervals are inside the valid range, and are disjoint.
#
#   returns     TRUE or FALSE
 
validLambdaMat  <- function( lambdamat, wavelength, bslist=NULL )
    {
    ok  = is.matrix(lambdamat)  &&  ncol(lambdamat)==2
    if( ! ok )
        {
        log_string( ERROR, "lambdamat is not an Mx2 matrix."  )
        return(FALSE)
        }
    
    if( nrow(lambdamat) == 0 )  return(TRUE)
    
    if( is.null(bslist) )   bslist  = breakandstep(wavelength)
    
    lambdarange = range( bslist$breakvec )
    inside      = all( lambdarange[1] <= lambdamat  &  lambdamat <= lambdarange[2] )
    if( ! inside )
        {
        log_string( ERROR, "Not all lambda's are in the interval [%g,%g].", lambdarange[1], lambdarange[2] )
        return(FALSE)
        }
        
    bandpass <- function( lambda )
        {
        if( lambda[1] < lambda[2] )
            return(TRUE)
        else if( lambda[2] < lambda[1] )
            return(FALSE)
        else
            return(NA)
        }
        
    bpmask  = apply( lambdamat, 1, bandpass )
    if( any(is.na(bpmask)) )
        {
        k   = which( is.na(bpmask)[1] )
        log_string( ERROR, "Interval [%g,%g] is empty.", lambdamat[k,1], lambdamat[k,2] )
        return(FALSE)
        }
        
    #   count the number of bandstops
    idxbs   = which( ! bpmask )
    if( 1 < length(idxbs) )
        {
        log_string( ERROR, "There are %d > 1 bandstops.", length(idxbs) )
        return(FALSE)
        }
        
    lambdabs    = NULL
    if( length(idxbs) == 1 )
        {
        #   extract the single bandstop
        lambdabs    = lambdamat[idxbs, ]
        
        #   and now remove it from lambdamat
        lambdamat   = lambdamat[-idxbs, , drop=FALSE]
        }
        
    #   lambdamat now has only bandpass intervals
    #   check that they are disjoint
    m   = nrow(lambdamat) 
    if( 1 < m )
        {
        #   sort rows by the lower endpoints
        perm        = order( lambdamat[ ,1] )
        lambdamat   = lambdamat[perm, ]
        
        #   check that all lower endpoints are distinct
        dif     = diff( lambdamat[ ,1] )
        idx     = which( dif==0 )
        if( 1 <= length(idx) )
            {
            log_string( ERROR, "Two intervals have the same lower endpoint = %g.", lambdamat[idx[1],1] )
            return(FALSE)
            }
        
        #   check that each upper endpoint is less than the next lower endpoint
        dif = lambdamat[ 2:m, 1 ]  -  lambdamat[ 1:(m-1), 2 ]
        idx     = which( dif<0 )
        if( 1 <= length(idx) )
            {
            log_string( ERROR, "Interval [%g,%g] intersects another interval.", lambdamat[idx[1],1],  lambdamat[idx[1],2] )
            return(FALSE)
            }
        }
    
    if( ! is.null(lambdabs)  &&  0<m  )
        {
        #   lambdabs is a bandstop interval
        #   check that it does not intersect one of the others
        ok  = lambdamat[m,2] < lambdabs[1]  &&  lambdabs[2] < lambdamat[1,1]
        if( ! ok )
            {
            log_string( ERROR, "Bandstop interval [%g,%g] intersects another interval.", lambdabs[1], lambdabs[2] )
            return(FALSE)
            }
        }
        
    return(TRUE)
    }
    
    
        
#   lambdamat   Mx2 matrix defining bandpass filters, and at most 1 bandstop wrap-around filter
#               for bandpass lambda1 < lambda2 and for bandstop lambda1 > lambda2
#               for lambda1 the transition is from 0 to 1, and lambda2 from 1 to 0.
#               there must be no overlap, but this is NOT checked, it's up to the caller.
#               Also lambda1 != lambda2 is not checked.
#               the bandstop should come in row #1, by convention.
#
#   value   a transmittance spectrum parameterized by wavelength.
#           it is a point in the unit n-cube
    
spectrumFromLambdas <- function( lambdamat, wavelength, bslist=NULL )
    {
    n   = length(wavelength)
    
    if( is.null(bslist) )   bslist  = breakandstep(wavelength)
    
    breakvec    = bslist$breakvec       # strictly increasing
    stepvec     = bslist$stepvec        # these are all positive
        
    out = numeric(n)
    names(out)  = as.character(wavelength)
  
    if( is.matrix(lambdamat) &&  nrow(lambdamat)==0 )    return(out)
    
    lambdamat   = prepareNxM( lambdamat, 2 )
    if( is.null(lambdamat) )   return(NULL)

    m   = nrow(lambdamat)      
    
    for( i in 1:m )
        {
        lam     = lambdamat[i, ]
        lambda  = sort( lam )   # now lambda[1] < lambda[2]
        
        idx = findInterval( lambda, breakvec )  #; print(idx)
        
        spectrum = numeric(n)        
        
        if( idx[1] == idx[2] )
            {
            #   both transitions inside the same bin
            spectrum[ idx[1] ]  = (lambda[2] - lambda[1]) / stepvec[ idx[1] ]
            }
        else
            {
            #   transitions in 2 different bins
            spectrum[ idx[1] ]  = (breakvec[ idx[1]+1 ] - lambda[1]) / stepvec[ idx[1] ]
            
            if( idx[2] <= n )
                spectrum[ idx[2] ]  = (lambda[2] - breakvec[ idx[2] ]) / stepvec[ idx[2] ] 
            
            if( idx[1]+1 <= idx[2]-1 )
                spectrum[ (idx[1]+1) : (idx[2]-1) ] = 1
            }
            
        if( lam[2] < lam[1] )
            #   convert from bandpass to bandstop
            spectrum    = 1 - spectrum
        
        out = out + spectrum
        }
        
    return( out )
    }
    
    
# spectrum      a transmittance spectrum parameterized by wavelength.
#               it is a point in the unit n-cube.    
# wavelength    increasing wavelengths of length n
#
# returns   an Mx2 matrix defining a superposition of bandpass and bandstop filters
#           there is at most 1 bandstop, and if present it is put in the first row
#
# It is known that every point in the image of the cube is the image of a characteristic function on the wavelength interval.
# The returned matrix essentially defines such a characteristic function.
# In this case the responsivities are step functions and the set is a finite union of intervals,
# or considered circularly, as arcs on the circle.

lambdasFromSpectrum <- function( spectrum, wavelength, bslist=NULL )
    {
    n   = length(spectrum)
    if( n <= 1 )    return( matrix(0,0,2) )
    
    if( any( is.na(spectrum) ) )
        {
        out = matrix( NA_real_, 1, 2 )
        colnames(out)   = c('lambda1','lambda2')            
        return( out )
        }
    
    
    if( is.null(bslist) )
        bslist  = breakandstep(wavelength)

    # print( bslist )
    
    
    breakvec    = bslist$breakvec       # strictly increasing        
    #   stepvec     = bslist$stepvec        # these are all positive    
    

    interior    = 0<spectrum  &  spectrum<1
    if( all(interior) )
        {
        rows    = floor( n/2 )
        
        out     = matrix( NA_real_, rows, 2 )
        colnames(out)   = c('lambda1','lambda2')
        
        odd     = seq(1,n-1,by=2)
        even    = odd+1L
        alpha   = spectrum[ odd ]
        out[ ,1]    = alpha*breakvec[odd]  +  (1-alpha)*breakvec[even]
        alpha   = spectrum[ even ]
        out[ ,2]    = (1-alpha)*breakvec[even]  +  alpha*breakvec[even+1L]
            
        if( n %% 2 == 1 )
            {
            #   if n is odd, there is 1 more bandpass
            #   there is no unique way to do this, so just center it in the bin
            step    = breakvec[n+1]  -  breakvec[n]
            s       = step * (1 - spectrum[n])/2
            lambda1 = breakvec[n]   + s
            lambda2 = breakvec[n+1] - s
            
            out = rbind( out, c(lambda1,lambda2) )
            }
            
        return(out)
        }
    
    out = matrix( 0, 0, 2 )
    colnames(out)   = c('lambda1','lambda2')    
    
    if( all(spectrum==0) )
        #   special case
        return( out )
    else if( all(spectrum==1) )
        {
        #   2 transitions
        out = rbind( out, breakvec[ c(1,n+1) ] )
        return( out )
        }
        
    
    inverted    = all( 0<spectrum )
    if( inverted )  spectrum    = 1 - spectrum
    
    
    #   find runs of nonzero coords
    mat     = findRunsTRUE( 0<spectrum, periodic=TRUE )
    if( nrow(mat) == 0 )
        {
        log_string( FATAL, "Found 0 runs of non-zero coords, impossible!" )
        return(NULL)
        }
        

        
    for( k in 1:nrow(mat) )
        {
        start   = mat[k,1]
        stop    = mat[k,2]
        
        # make sequence from start to stop
        if( start <= stop )
            iseq    = start:stop
        else
            iseq    = c( start:n, 1:stop )
            
        if( length(iseq) == 1 )
            {
            #   special case, center this bandpass in the bin
            step    = breakvec[iseq+1]  -  breakvec[iseq]
            s       = step * (1 - spectrum[iseq])/2
            lambda1 = breakvec[iseq]   + s
            lambda2 = breakvec[iseq+1] - s
            
            out = rbind( out, c(lambda1,lambda2) )

            next
            }
            
        val = 0
        for( i in iseq )
            {
            #   cat( "------------- i=", i, "  val=", val, "--------\n" )
            
            alpha   = spectrum[i]            
            
            lambda2 = NA_real_
                        
            if( val == 0 )
                {                
                lambda1 = alpha*breakvec[i]  +  (1-alpha)*breakvec[i+1L] #; print( lambda1 )
                val = 1
                
                if( i == iseq[length(iseq)] )
                    {
                    #   the last i, so clean up
                    #print( i )
                    lambda2 = breakvec[i+1]
                    val = 0 # not really necessary
                    }                
                }
            else
                {
                #   val must be 1
                if( alpha == 1 )
                    {
                    #print( iseq )
                    if( i == iseq[length(iseq)] )
                        {
                        #   the last i
                        #   print( i )
                        lambda2 = breakvec[i+1]
                        val = 0 # not really necessary
                        }
                    }
                else
                    {
                    lambda2 = (1-alpha)*breakvec[i]  +  alpha*breakvec[i+1L]                    
                    val = 0
                    }
                }
                
            if( is.finite(lambda2) )            
                out = rbind( out, c(lambda1,lambda2) )                
            }
        }

    if( inverted )
        {
        #   fix output
        unshape = as.numeric( t(out) )
        unshape = unshape[ c( 2:length(unshape), 1 ) ]  # shift 1
        out     = matrix( unshape, length(unshape)/2, 2, byrow=TRUE )        
        colnames(out)   = c('lambda1','lambda2')    
        }
        
    #   if there is a bandstop, ensure that it is in the first row, and is unique
    bandstop    = which( out[ ,2] < out[ ,1] )
    count       = length(bandstop)
    if( 1 < count )
        {
        log_string( FATAL, "There are %d > 1 bandstops.", count )
        return(NULL)
        }    

    if( count==1  &&  1<bandstop  )
        {
        #   rotate into position
        perm    = c( bandstop:nrow(out), 1:(bandstop-1) )
        out     = out[ perm, , drop=FALSE]
        }
        
    rownames(out)   = 1:nrow(out)
        
    if( count==1 )
        {
        rownames(out)[1]    = "BS"
        }
        
    bandpass    = which( out[ ,1] < out[ ,2] )
    if( 0 < length(bandpass) )
        rownames(out)[bandpass] = sprintf( "BP%d", 1:length(bandpass) )
    
    return(out)
    }
    
    
#   test that spectrum -> lambdas -> spectrum is identity    
testSpecLamRoundTrip <- function( n, samples, sd=5, scale=100, tol=5.e-12 )
    {
    set.seed(0)
    
    wave    = 400:(400+n-1)
    
    delta   = numeric(samples)
    
    transmax    = 0
    transmin    = Inf
    
    for( i in 1:samples )
        {
        spec    = randomSpectrum( n, sd=sd, scale=scale )  # ; print( spec )
        
        lambdas = lambdasFromSpectrum( spec, wave )
        
        spec.back   = spectrumFromLambdas( lambdas, wave )
        
        delta[i]    = max( abs(spec - spec.back) )
        
        if( tol < delta[i] )
            {
            cat( "----------------\n" )
            print( spec )
            print( lambdas )
            print( spec.back )
            }
            
        transmin    = min( length(lambdas), transmin )
        transmax    = max( length(lambdas), transmax )
        }
    
    count   = sum( tol < delta )
    mess    = sprintf( "%d violations of %d samples (delta > %g).    max(delta)=%g", 
                                count, samples, tol, max(delta) )
    cat( mess, '\n' )
    mess    = sprintf( "transition range: %d to %d.", transmin, transmax )
    cat( mess, '\n' )
    
    return( count == 0 )
    }
    
#   n       dimension of n-cube
#   sd      standard-deviation of smoothing filter
#   scale   bigger means more 0's and 1's
#    
#   returns a random point in the n-cube, with emphasis on lots of 0's and 1's    
#
#   to get no 0's and 1's set sd=0 and scale=1
randomSpectrum  <-  function( n, sd=5, scale=100 )
    {
    out = runif( n )
    
    if( 0 < sd )
        {
        kern    = stats::dnorm( seq(-3*sd,3*sd,len=n/2+1), sd=sd )
        kern    = kern / sum(kern)
        #   out = out[ is.finite(out) ] not needed when circular        
        out = stats::filter( out, kern, circular=TRUE ) #; print(out)
        out = as.numeric(out)
        }

    out = pmin( pmax( scale*(out-0.5) + 0.5, 0 ), 1 )
    
    return( out )
    }
    

#   age     vector of ages, between 20 and Inf
#   return absorbance of human lens

lensAbsorbance  <-  function( age=32, wavelength=400:700 )
    {
    lens    =  LensAbsorbance1987       # colorSpec::LensAbsorbance1987

    #   we know lens[] is a matrix
    TL1     = lens[ , 'TL1' ]
    TL2     = lens[ , 'TL2' ]
    
    out = matrix( NA_real_, numWavelengths(lens), length(age) )
    
    for( k in seq_len( length(age) ) )
        {
        A = age[k]
        
        if( A < 20 )    next
        
        if( A < 60 )
            out[ ,k] = (1 + 0.02*(A - 32)) * TL1  +  TL2
        else
            out[ ,k] = (1.56 + 0.06667*(A - 60)) * TL1  +  TL2
        }
        
    colnames( out ) = sprintf( "age%g", age )
        
    out = colorSpec( out, wavelength(lens), quantity=quantity(lens) )
    
    out = resample( out, wavelength )
    
    return( out )
    }
    
    
#--------       UseMethod() calls           --------------#            
              
bandRepresentation <- function( x )
    {
    UseMethod("bandRepresentation")
    }        
    
    
    
#--------       obsolete below              --------------#            

#   input the spectrum
#   Use the the alphas to tweak the endpoint wavelengths from spectramat (usually integer values)
#   to fractional wavelengths that are returned
#   Sometimes the extension is in the positive direction, and sometimes negative.
    
#   spectrum        numeric N-vector with a spectrum
#                   all values 0-1 with at most 2 values in (0,1)
#                   and these non-trivial 'alphas' must be located between a 0 and a 1,
#                   or adjacent to each other
#   wave            N-vector of wavelengths
#   step            N-vector of wavelength steps (redundant)
#
#   return     c(lambda1,lambda2), lambda1 is step up, and lambda2 is step down   
#
#   this is a partial inverse for rectangularMaterial() in the case alpha=1  
   
compute_lambdas <- function( spectrum, wave, step=NULL )
    {
    lambda  = rep( NA_real_, 2 )
    
    dim( spectrum ) = NULL
    
    if( any( is.na(spectrum) ) ) return(lambda)
        
    n   = length(spectrum)

    if( is.null(step) )
        {
        #step    = diff( wave, lag=2 ) / 2
        #step    = c( wave[2]-wave[1], step, wave[n]-wave[n-1] )     #; print(step)
        step    = breakandstep(wave)$stepvec
        }
    
    if( length(wave) != n ) return(lambda)

    onemask     = spectrum == 1
    if( all(onemask) )   return(lambda)        # no transitions
    
    zeromask    = spectrum == 0
    if( all(zeromask) )   return(lambda)       # no transitions
    
    alphaidx    = which( ! (onemask | zeromask) )
    
    if( 2 < length(alphaidx) )    return(lambda)       # invalid spectrum, too many are not 0 or 1
    
    alpha   =   spectrum[ alphaidx ]
    
    if( any( alpha<0 | 1<alpha ) )   return(lambda)    # invalid spectrum, alpha out of range
    
    both01  = any(onemask)  &&  any(zeromask)
    
    if( both01 )
        {
        #   there are both 0s and 1s, which is the usual case
        if( 0 < length(alphaidx) )
            {
            # one or both of these must have 0 on one side and 1 on the other side
            inext   = c( 2:n, 1 )
            iprev   = c(n, 1:(n-1) )
            for( i in alphaidx )
                {
                ok  = (spectrum[ iprev[i] ] + spectrum[ inext[i] ] == 1)
                
                if( ! ok )  return(lambda)
                }
            }
        }
    else
        {
        #   spectrum is all 0s or all 1s, except at alphaidx
        if( length(alphaidx) == 2 )
            {
            #   they must be adjacent
            adjacent    = alphaidx[1]+1==alphaidx[2]  ||  (alphaidx[1]==n && alphaidx[2]==1)        
        
            if( ! adjacent )      return(lambda)
            }
        }
        
    #   prepare for transition-finding
    specopt = spectrum
    specopt[ alphaidx ] = 0
    
    
    #   specopt is now truly 0-1
    trans   = transitionMatrix( specopt )  #; print(trans) ; print( wave[t(trans)] )  # nrow(trans) must be even
    
    #print( specopt )    
    #print( trans )
    
    if( 2 < nrow(trans) )   return(lambda)  # too many transitions
    
    if( nrow(trans) == 0 ) 
        {
        #   there are 0 transitions
        #   specopt must be all 0s.  Color is almost black.
        #   there must be 1 or 2 alphas; if there are 0 then there could not be 0 transitions
        if( length(alphaidx) == 1 ) 
            {
            #   alphaidx has 0s on either side of it
            #   distribute alpha[i] evenly on both sides of this wavelength
            #   this is a bandpass
            halfwidth   = 0.5 * step[alphaidx] * alpha
            lambda[1]   = wave[alphaidx] - halfwidth       
            lambda[2]   = wave[alphaidx] + halfwidth     
            }
        else if( length(alphaidx) == 2 )
            {
            #   The 2 alphas must be adjacent, for otherwise there would have been more than 2 transitions.
            #   check this !
            adjacent    = alphaidx[1]+1==alphaidx[2]  ||  (alphaidx[1]==n && alphaidx[2]==1)
            if( ! adjacent )
                {
                log_string( ERROR, "For spectrum, alphaidx=%d,%d is invalid.", alphaidx[1], alphaidx[2] )
                return(lambda)
                }
            
            #   this is a bandpass
            j1  = alphaidx[1]
            lambda[1] = wave[j1] + (0.5 - alpha[1]) * step[j1]
            j2  = alphaidx[2]
            lambda[2] = wave[j2] + (alpha[2] - 0.5) * step[j2]
            }
        }

    if( nrow(trans) == 2 ) 
        {
        if( length(alphaidx)==1  &&  sum( trans == alphaidx )==2 )
            {
            #   alpha appears in *both* transitions, so if we use it as it is, its impact is doubled.
            #   this happens if an alpha is between two 1s
            #   cut the weight in half
            #cat( "halving alpha\n" )
            alpha   = alpha / 2 
            }
            
        #   examine both transitions, and compute lambda1 and lamba2 from them
        for( k in 1:2 )
            {
            i1  = trans[k,1]
            i2  = trans[k,2]
            
            #   transition is from i1 to i2
            if( specopt[i1]==0  && specopt[i2]==1 )
                {
                #   step-up from i1. assign lambda1
                jj  = match( i1, alphaidx )
                a   = ifelse( is.na(jj), 0, alpha[jj] ) 
                lambda[1] = wave[i1] + (0.5 - a) * step[i1]
                }
            else if( specopt[i1]==1  &&  specopt[i2]==0 )
                {
                #   step-down from i1. assign lambda2
                jj  = match( i2, alphaidx )     # jj is 1,2, or NA
                a   = ifelse( is.na(jj), 0, alpha[jj] ) 
                lambda[2] = wave[i2] + (a - 0.5) * step[i2]
                }
            else
                {
                log_string( ERROR, "transition matrix is invalid." )
                return( lambda )
                }
            }                      
        }
        
    #   by now, both lambda[1] and lambda[2] should be defined !  
    #   check it
    if( any( is.na(lambda) ) )
        {
        cat( "compute_lambdas()  internal error\n" )
        names(spectrum)  = wave
        print( spectrum )
        return( rep( NA_real_, 2 ) )
        }

    return( lambda )
    }
            