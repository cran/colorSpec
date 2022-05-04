
        
#   temperature     a vector of temperatures, Kelvin
#   components      colorSpec object with components S0, S1, S2
#   wavelength      at which to resample.  If NULL then the wavelengths in components are used   
#   roundMs         round M1 and M2 to 3 decimal places
#
#   return value: colorSpec object with quantity  'energy'
#   Energy is normalized so it is 1 at 560 nm
#   depends on global variable - colorSpec::daylight1964
#
daylightSpectra <-  function( temperature, wavelength=NULL, components=colorSpec::daylight1964, roundMs=FALSE )
    {    
    locus = daylightLocus( temperature )
        
    m   = nrow(locus)   # same as length(temperature)
    
    mat3xm  = matrix( NA_real_, 3, m )
        
    for( k in 1:nrow(locus) )
        {
        x = locus$x[k]
        y = locus$y[k]
        
        if( is.na(x) ) next
    
        M   = 0.0241  +  0.2562 * x  -  0.7341 * y
        
        M1  = ( -1.3515 - 1.7703 * x  +  5.9114 * y ) / M
        
        M2  = ( 0.0300 - 31.4424 * x  +  30.0717 * y ) / M
        
        weight  =  c( 1, M1, M2 )
        
        if( roundMs )
            weight  = round( weight, 3 )
        
        mat3xm[ ,k]  =  weight
        
        #S  =  multiply( components, mat3x1 )                   #  data$S0  +  M1 * data$S1  +  M2 * data$S2
        #out[ ,k]  = coredata(S)
        #   change name of the last column
        #   colnames( out )[ ncol(out) ] = sprintf( "D%g", round(locus$temperature[k]) )
        }
        
    #   components always normalizes to 100 at 560nm        
    #   but we prefer energy=1 at 560nm
    mat3xm  = mat3xm / 100  #; print( mat3xm )
       
    colnames( mat3xm ) = sprintf( "D%g", round(temperature) )
              
    out = multiply( components, mat3xm )
        
    #  specnames( out ) = sprintf( "D%g", round(temperature) )   not necessary, since multiply() assigns from colnames(mat3xm)

    if( ! is.null(wavelength) )
        out = resample( out, wavelength )
    
    #   organization(out) = 'matrix'
    
    if( m == 1 )    organization(out) = 'vector'
    
    return( out )
    }
        
    
    
    
#   .temperature    a vector of temperatures in the interval [4000,25000]
#   .space          space to return 'xy', or 'XYZ'
#
#   returns a dataframe with 3 columns:  temperature, x, y
#   or a data.frame with 2 columns, temperature, XYZ
#   if any temperature is outside the valid range x,y are set to NA
#
#   the official lower limit is 4000 K, 
#   but to help with plotting I have extended it to 3500 K using the same polynomial.
#   the upper limit is 25000.

daylightLocus <-  function( .temperature, .space='xy'  )
    {
    t_inv = 10^3 / .temperature
    
    x = rep( as.numeric(NA),  length(.temperature) )
    
    idx = ( 3500 <= .temperature  &  .temperature <= 7000 )

    x[ idx ] =  0.244063  +  0.09911 * t_inv[idx]  +  2.9678 * t_inv[idx]^2  -  4.6070 * t_inv[idx]^3

    
    idx = ( 7000 <= .temperature  &  .temperature <= 25000 )
    
    x[ idx ] = 0.237040  +  0.24748 * t_inv[idx]  +  1.9018 * t_inv[idx]^2  -  2.0064 * t_inv[idx]^3

    y   = -3 * x^2  + 2.870 * x  - 0.275
        
    if( .space == 'xy' )
        out = data.frame( temperature=.temperature, x=x, y=y )
    else if( .space == 'XYZ' )
        {
        XYZ = XYZ_from_xyY( x, y, 1 )
        #   class(XYZ)  = 'model.matrix'
        
        out = data.frame( temperature=.temperature )    #, XYZ=XYZ )
        out$XYZ = XYZ
        }
    
    
    return( out )
    }
    