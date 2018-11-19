

#   temperature     a vector of N temperatures, in Kelvin
#   wavelength      at which to sample, in nm
#   normalize       logical to normalize to 1 at 560nm

#
#   returns colorSpec object with N spectra
#
#   if normalize is FALSE, then the output units are W*M^{-2}*nm^{-1}
#   see W&S pp. 12-13 and Figure 1(1.2.2)

planckSpectra <-  function( temperature, wavelength=300:830, normalize=TRUE, c2=1.4388e7 )
    {
    #   fundamental constants - CODATA 2014
    h   = 6.62607004e-34    # Planck's Contant in, J*s or W*s^2
    c   = 2.99792458e17     # speed of light, in nm/s
    k   = 1.38064852e-23    # Boltzmann constant, in J/K
    
    if( ! normalize )
        {
        c1  = 2*pi*h*c^2 * 1.e18    # the 1.e18 converts from W*nm^{-2}*nm^{-1} to W*M^{-2}*nm^{-1}    1.e18 nm^2 / M^2
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

    #   waveM   = wavelength * 1.e-9   #   convert from nm to m
        
    lambda0 = 560   # e-9
    
    pow0    = lambda0 ^ -5 
    pow1    = wavelength ^ -5 
        
    mat = matrix( 0, length(wavelength), length(temperature) )
    
    for( j in 1:length(temperature) )
        {
        T       = temperature[j]
        
        Power   = pow1 / ( exp( c2 / (wavelength * T) ) - 1 )
    
        if( normalize )
            {
            #   normalize to 1 at 560 nm
            pow00   = pow0 / ( exp( c2 / (lambda0 * T) )  -  1 )
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
    
    #   compute the midpoints
    #mid = 0.5 * ( wavelength[1:(p-1)] + wavelength[2:p] )
    #mid = c( wavelength[1] - (mid[1]-wavelength[1]), mid, wavelength[p] + (wavelength[p]-mid[p-1]) )    # extend by 1 on either side
    
    mat = matrix( NA_real_, p, n )
    
    for( i in 1:n )
        {
        if( alpha[i] == 0 )
            {
            #   special case; use this shortcut to save some time
            mat[ ,i]    = 0.5
            next
            }

        lam = lambda[i, ]
        
        if( lam[1] == lam[2] )
            {
            #   special case - undefined
            log.string( WARN, "material is undefined because lambda_min==lambda_max == %g.", lam[1] )
            next
            }

        lambda_min  = min(lam)
        lambda_max  = max(lam)
        
        #y   = as.numeric( lambda_min <= mid  &  mid < lambda_max )
        #y   = approx( mid, y, xout=wavelength, rule=2 )$y 
        
        y   = 0.25 * (sign(wavelength - lambda_min) + 1) * (sign(lambda_max - wavelength) + 1)
        
        a   = ifelse( lam[1] < lam[2], alpha[i], -alpha[i] )
        
        #   if( lam[2] < lam[1] )   y = 1 - y   # flip it
        
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
        