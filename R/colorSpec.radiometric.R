
#   radiometric.colorSpec()
#
#   force any spectrum to be radiometric-based, i.e. energy based
#   Colorimetry uses radiometric units and not actinometric units.
#   If a spectrum is already radiometric, it is returned unchanged.
#
#   for conversion of 'light' a scaling factor of 10^18 is applied so that:
#       photons     -> attowatts
#       kilophotons -> femtowatts
#       megaphotons -> picowatts
#       gigaphotons -> nanowatts
#       teraphotons -> microwatts
#       petaphotons -> milliwatts
#       exaphotons  -> watts
#
#   for conversion of 'responsivity.light':
#   the most commmon actinometric unit for "responsivity.light" is Quantum Efficiency
#   convert to responsivity in A/W
#
#   x           a colorSpec object, typically actinometric
#   multiplier  scalar intended for unit conversion
#   warn        may be appropriate when called from another function

radiometric.colorSpec <- function( x, multiplier=1, warn=FALSE )
    {
    if( ! is.actinometric(x) )   return(x)   # no change needed
    
    quantity    = quantity( x )
        
    type    = spectrumTypeFromQuantity( quantity )
    
    h   = 6.626070040e-34   # Planck's Contant in J*s or W*s^2
    c   = 2.99792458e17     # speed of light in nm/s
    N_A = 6.022140857e23    # Avogadro's Constant
    
    #   lambda  = 1.e-9 * wavelength(x) # convert lambda from nm to m
    lambda  = wavelength(x)

    if( type == "light" )
        {    
        #   the most commmon actinometric unit for "light" is micromoles_of_photons 
        #   a micromole of photons = 1.e-6*N_A photons
        #   the output unit for the light is joule
        phys    = 1.e-6*N_A*h*c
        final   = multiplier * phys
        log.string( INFO, "final multiplier = multiplier * 10^{-6}*N_A*h*c / lambda = %g * %g / lambda = %g / lambda.",
                            multiplier, phys, final )
        vec     = final / lambda
        quantity.radio  = "energy"
        }
    else if( grepl( "electrical$", quantity ) )
        {
        #   the most commmon actinometric unit for electrical "responsivity.light" is Quantum Efficiency
        #   the output unit for responsivity is  J/C  or   A/W    
        e   = 1.6021766208e-19  # charge of electron, in C
        phys    = e/(h*c)
        final   = multiplier * phys
        log.string( INFO, "final multiplier = multiplier * e/(h*c) * lambda = %g * %g * lambda = %g * lambda.",
                            multiplier, phys, final )
        vec     = final * lambda
        quantity.radio  = sub( "photons", "energy", quantity )    
        }
    else
        {
        #   the response type is either "neural" or "action"
        #   once again, assume micromoles of photons,
        #   this is the reciprocal of above "light"
        phys    = 1.e6/(N_A*h*c)
        final   = multiplier * phys
        log.string( INFO, "final multiplier = multiplier * 1.e6/(N_A*h*c) * lambda = %g * %g * lambda = %g * lambda.",
                            multiplier, phys, final )
        vec     = final * lambda
        quantity.radio  = sub( "photons", "energy", quantity )         
        }
        
    #   vec[] defined above has length(vec) == # of wavelengths
    #   use vec to define myfun()
    myfun <- function( y )  { vec * y }   
    
    #   apply the function to all spectra
    out = applyspec.colorSpec( x, myfun )
    
    #   and change the quantity to be radiometric
    quantity(out)   = quantity.radio
    
    if( warn )
        {
        log.string( WARN, "Object x has been converted from actinometric (quantity='%s') to radiometric (quantity='%s').",
                            quantity(x), quantity(out) )
        }
    
    return( out )
    }
    
    
is.radiometric.colorSpec <- function( x )
    {
    #   if( ! is.colorSpec(x) ) return( FALSE )     # very unlikely to happen
    
    return( grepl( "^(energy|power)", quantity(x) ) )     # 'energy' is preferred, 'power is deprecated
    }
    
    
#--------       UseMethod() calls           --------------#            

radiometric <- function( x, multiplier=1, warn=FALSE )
    {
    UseMethod("radiometric")
    }

is.radiometric <- function( x )
    {
    UseMethod("is.radiometric")
    }    