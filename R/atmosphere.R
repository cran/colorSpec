
#   transatmos()    calculate transmittance of atmosphere at sea level, along a horizontal path
#                   Since the path is horizontal, the air properties are constant along the path
#
#   distance        the distance of the optical path, in meters
#                   can also be a vector of distances
#
#   wavelength      at which to compute transmittance spectrum, in nanometers
#
#   molecules       list of parameters describing the molecules in the atmosphere
#       N   molecular density of atmosphere at sea level, molecules/m^3
#       n0  refractive index of air molecules
#   If molecules is NULL, then molecular transmittance is identically 1
#
#   aerosols        list of parameters describing the aerosols in the atmosphere
#       It uses the Angstrom model:  attenuation = beta * (lambda/lambda0)^(-alpha)
#       metrange    meteorological range, in meters. At this range transmittance=0.02 at lambda0=550nm. [Koschmieder]
#                       Used to calculate alpha and beta, but if metrange=NULL
#                       then alpha and beta can also given directly
#       alpha       Angstrom exponent, applies to the wavelength. dimensionless
#       beta        in m^-1
#   If aerosols is NULL, then aerosol transmittance is identically 1
#
#   Details on calculation of alpha and beta
#       If metrange is given, we use the Kruse function to calculate alpha.
#       alpha can also be a user-supplied custom function, taking meters to alpha,
#       and then this function is used in place of the Kruse model function.
#       beta is computed so that the product of
#       molecular and aerosol transmittance yields the desired metrange.
#
#
#   If both molecules and aerosols are NULL, then the returned transmittance is identically 1.
#   The atmosphere has become a vacuum.
#
#   value       a colorSpec object with quantity 'transmittance'
#               and number of spectra equal to length(distance) 

atmosTransmittance <- function( distance, wavelength=380:720, 
                                molecules=list(N=2.547305e25,n0=1.000293),
                                aerosols=list(metrange=25000,alpha=0.8,beta=0.0001) )
    {
    #   initialize attenuation vector to 0s
    mu  = numeric( length(wavelength) )
    
    if( ! is.null(molecules) )
        {
        #   Rayleigh scattering
        #   compute molecular attenuation.  units are m^-1
        if( ! is.numeric(molecules$N)  ||  ! is.numeric(molecules$n0) )
            {
            log_level( ERROR, "Either N or n0 is not numeric, or missing." )
            return(NULL)
            }
        
        mu_mol  = attenuation.molec( wavelength, molecules$N, molecules$n0 )
        
        if( is.null(mu_mol) )  return(NULL)
        
        mu = mu + mu_mol
        }
    
    if( ! is.null(aerosols) )
        {
        #   Mie scattering
        lambda0 = 550        
        
        if( is.null(aerosols$metrange) )
            {
            alpha   = aerosols$alpha
            beta    = aerosols$beta
            
            if( ! is.numeric(alpha)  ||  ! is.numeric(beta) )
                {
                log_level( ERROR, "Either alpha or beta is not numeric, or missing." )
                return(NULL)
                }            
            }
        else
            {
            if( ! is.numeric(aerosols$metrange)  ||  aerosols$metrange <= 0 )
                {
                log_level( ERROR, "metrange='%s' is not a positive number", as.character(aerosols$metrange) )
                return(NULL)
                }
            
            alpha   = alphaFromVr.Kruse( aerosols$metrange )
                
            if( ! is.null(molecules) )
                mu_mol0 = attenuation.molec( lambda0, molecules$N, molecules$n0 )
            else
                mu_mol0 = 0
                
            beta    = -mu_mol0 - log(0.02) / aerosols$metrange
            
            log_level( INFO, "From metrange=%g, computed alpha=%g, mu_mol0=%g, beta=%g.",
                                aerosols$metrange, alpha, mu_mol0, beta )
            #cat( 'alpha=', alpha, '\n' )
            #cat( 'mu_mol0=', mu_mol0, '\n' )            
            #cat( 'beta=', beta, '\n' )
            }

        #   compute aerosol attenuation.  units are m^-1            
        mu_aero = attenuation.aerosol( wavelength, alpha, beta, lambda0 )
        
        if( is.null(mu_aero) )  return(NULL)
                
        mu = mu + mu_aero 
        }
    
    #   compute optical depth = attenuation * distance.  It is dimensionless.        
    od  = mu %o% distance
        
    #   convert optical depth to transmittance
    out = exp( -od )
    colnames(out)   = sprintf( "dist=%gm", distance )
        
    out = colorSpec( out, wavelength=wavelength, quantity='transmittance' )
    
    return( out )    
    }
    
    
#   calculate molecular attenuation at each wavelength.
#   returned units are m^-1                
attenuation.molec   <- function( wavelength, N=2.547305e25, n0=1.000293 )
    {
    wavem   = wavelength * 1.e-9    # convert wavelength from nm to meters

    out = 8*pi^3 * (n0^2 - 1)^2 / (3 * N * wavem^4)
        
    return( out )
    }
    
    
#   calculate aerosol attenuation at each wavelength.  units are m^-1                
attenuation.aerosol   <- function( wavelength, alpha, beta, lambda0=550 )
    {
    ok  = is.numeric(alpha) && is.numeric(beta)  && length(alpha)==1  &&  length(beta)==1
    if( ! ok )
        {
        log_level( ERROR, "Either alpha or beta is not a numeric scalar." )
        return(NULL)
        }
        
    if( beta < 0 )
        {
        log_level( WARN, "The turbidity coefficient beta=%g < 0, and is not physically realizable.", beta )
        }        
            
    out = beta * (wavelength/lambda0)^(-alpha)
        
    return( out )
    }
    

#   Vr  visibility range in meters    
alphaFromVr.Kruse  <-  function( Vr )    
    {
    mask1   = Vr <= 6000
    mask2   = Vr <= 50000
    
    out     = Vr    #   just to get the size right
    
    out[ mask1 ]            = 0.585 * (Vr[mask1]/1000)^(1/3)
    out[ !mask1 & mask2 ]   = 1.3
    out[ !mask2 ]           = 1.6
        
    return(out)
    }
    
    
