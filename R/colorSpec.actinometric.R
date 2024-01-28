

#   actinometric.colorSpec()
#
#   force any spectrum to be actinometric-based, i.e. photon based
#   Colorimetry uses radiometric units and not actinometric units.
#   If a spectrum is already actinometric, it is returned unchanged.
#
#   for 'light':
#   the most commmon radiometric unit is W; convert to micromoles_of_photons/time
#
#   for 'responsivity.light':
#   the most commmon radiometric unit is A/W  or  J/C; convert to  Quantum Efficiency
#
#   multiplier  output multiplier; useful for unit conversion

actinometric.colorSpec <- function( x, multiplier=1, warn=FALSE )
    {
    if( ! is.radiometric(x) )   return(x)   # no change needed
        
    quantity    = quantity( x )

    type    = spectrumTypeFromQuantity( quantity )
    
    h   = 6.626070040e-34   # Planck's Contant in J*s or W*s^2
    c   = 2.99792458e17     # speed of light in nm/s
    N_A = 6.022140857e23    # Avogadro's Constant in mole^{-1}
    
    #   lambda  = 1.e-9 * wavelength(x) # convert lambda from nm to m
    lambda  = wavelength(x)

    if( type == "light" )
        {    
        #   the most common radiometric unit for "light" is joule
        #   the most commmon actinometric unit for "light" is micromoles_of_photons 
        #   a micromole of photons = 1.e-6*N_A photons
        phys    = 1.e6/(N_A*h*c)
        final   = multiplier * phys
        log_string( INFO, "final multiplier = multiplier * 1.e6/(N_A*h*c) * lambda = %g * %g * lambda = %g * lambda.",
                            multiplier, phys, final )        
        vec     = final * lambda
        quantity.actino  = "photons"
        }
    else if( grepl( "electrical$", quantity ) )
        {
        #   the most common input unit for responsivity is  J/C  or   A/W            
        #   the output actinometric unit for electrical "responsivity.light" is Quantum Efficiency
        e       = 1.6021766208e-19  # charge of electron, in C
        phys    = (h*c)/e
        final   = multiplier * phys
        log_string( INFO, "final multiplier = multiplier * ((h*c)/e) / lambda = %g * %g / lambda = %g / lambda.",
                            multiplier, phys, final )
        vec     = final / lambda
        quantity.actino  = sub( "^(energy|power)", "photons", quantity )    
        }
    else
        {
        #   the response type is either "neural" or "action"
        #   once again, assume micromoles of photons,
        #   this is the reciprocal of above "light"
        phys    = 1.e-6*N_A*h*c
        final   = multiplier * phys
        log_string( INFO, "final multiplier = multiplier * 10^{-6}*N_A*h*c / lambda = %g * %g / lambda = %g / lambda.",
                            multiplier, phys, final )
        vec     = final / lambda
        quantity.actino  = sub( "^(energy|power)", "photons", quantity )         
        }
        
    #   vec[] defined above has length(vec) == # of wavelengths
    #   use vec to define myfun()
    myfun <- function( y )  { vec * y }       
    
    #   apply the function to all spectra
    out = applyspec.colorSpec( x, myfun )
    
    #   and change the quantity to be actinometric
    quantity(out)   = quantity.actino
    
    if( warn )
        {
        log_string( WARN, "Object x has been converted from radiometric (quantity='%s') to actinometric (quantity='%s').",
                            quantity(x), quantity(out) )
        }
        
    return( out )
    }
    

is.actinometric.colorSpec <- function( x )
    {
    #   if( ! is.colorSpec(x) ) return( FALSE )     # very unlikely to happen
    
    return( grepl( "^photons", quantity(x) ) )
    }
    
    
    
#--------       UseMethod() calls           --------------#            

actinometric <- function( x, multiplier=1, warn=FALSE )
    {
    UseMethod("actinometric")
    }

is.actinometric <- function( x )
    {
    UseMethod("is.actinometric")
    }

    