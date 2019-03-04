

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
        log.string( ERROR, "c2 = '%s' is invalid.", as.character(c2) )
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
        log.string( ERROR, "alpha is invalid, because it is not numeric." )
        return(NULL)
        }
    
    n   = nrow(lambda)
    if( length(alpha) == 1 )    alpha = rep( alpha, n )
    
    if( length(alpha) != n )
        {
        log.string( ERROR, "length(alpha) = %d != %d = nrow(lambda).", length(alpha), nrow(lambda) )
        return(NULL)
        }
    
    if( ! all( -1<=alpha  &  alpha<=1 ) )
        {
        log.string( ERROR, "alpha is invalid; all values must be in [-1,1]." )
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
            log.string( WARN, "material is undefined because lambda_min==lambda_max == %g.", lam[1] )
            next
            }

        #lambda_min  = min(lam)
        #lambda_max  = max(lam)
        
        #y   = as.numeric( lambda_min <= mid  &  mid < lambda_max )
        #y   = approx( mid, y, xout=wavelength, rule=2 )$y 
        
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
                log.string( ERROR, "For spectrum, alphaidx=%d,%d is invalid.", alphaidx[1], alphaidx[2] )
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
                log.string( ERROR, "transition matrix is invalid." )
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
        