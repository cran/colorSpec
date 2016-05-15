

#   temperature     a vector of N temperatures, in Kelvin
#   normalize       logical to normalize to 1 at 560nm
#   wavelength      at which to sample, in nm
#
#   returns colorSpec object with N spectra
#
#   see W&S pp. 12-13 and Figure 1(1.1.2)

planckSpectra <-  function( temperature, normalize=TRUE, wavelength=300:830 )
    {
    h   = 6.626176e-34      # Planck constant
    c   = 299792458         # speed of light
    k   = 1.380662e-23      # Boltzmann constant
    
    c1  = 2*pi*h*c^2        #; print( c1 )
    c2  = h*c/k             #; print( c2 )

    waveM   = wavelength * 1.e-9   #   convert from nm to m
        
    lambda0 = 560e-9
    
    pow0    = c1 * lambda0 ^ -5 

    pow1    = c1 * waveM ^ -5 
        
    mat = matrix( 0, length(waveM), length(temperature) )
    
    for( j in 1:length(temperature) )
        {
        T       = temperature[j]
        
        Power   = pow1 / ( exp( c2 / (waveM * T ) )  -  1 )
    
        if( normalize )
            {
            #   normalize to 1 at 560 nm
            pow00   = pow0 / ( exp( c2 / (lambda0 * T ) )  -  1 )
            Power = Power / pow00
            }
            
        mat[ ,j]    = Power 
        }
        
    if( ! normalize )
        # convert from 1/M to 1/nm
        mat = 1.e-9 * mat
        
    colnames(mat)   = sprintf( "P%g", round(temperature)  )    
    
    out = colorSpec( mat, wavelength, "power", "matrix" )

    #   if( normalize ) out = normalize( out, 560 )
    
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
    
    out = colorSpec( out, wavelength, "power->action" )
    
    specnames( out )    = "erythemal"
    
    metadata( out ) = list( description="A.F. McKinlay and B.L. Diffey (1987)" )
    
    return( out )
    }
    
    
illuminantE <- function( power=1, wavelength=380:780 )
    {    
    n   = length(wavelength)
    
    theNames = sprintf( "E%g", power )
    
    if( length(power) == 1 )
        core    = rep(power,n)
    else
        {
        core            = matrix( power, n, length(power), byrow=T )
        colnames(core)  = theNames
        }
        
    out = colorSpec( core, wavelength, "power" )
    
    if( length(power) == 1 )    specnames( out ) = theNames
    
    metadata( out ) = list( description="Equal Energy White" )
        
    return( out )
    }
    
neutralMaterial <- function( gray=1, wavelength=380:780 )
    {    
    n   = length(wavelength)

    core    = matrix( gray, n, length(gray), byrow=T )
    
    colnames(core)  = sprintf( "Neutral%g", gray )

    out = colorSpec( core, wavelength, "reflectance" )
    
    if( length(gray) == 1 )   organization( out ) = 'vector'
           
    metadata( out ) = list( description="Neutral Gray" )
                 
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
        
    out = colorSpec( out, wavelength(lens), quantity(lens) )
    
    out = resample( out, wavelength )
    
    return( out )
    }
        