

#   radiometric.colorSpec()
#
#   force any spectrum to be radiometric-based, i.e. power based
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


radiometric.colorSpec <- function( x )
    {
    quantity    = quantity( x )
    
    actinometric    = grepl( "^photons", quantity )
    
    if( ! actinometric )   return(x)   # no change needed
    
    type    = spectrumTypeFromQuantity( quantity )
    
    h   = 6.626070040e-34   # Planck's Contant in J*s or W*s^2
    c   = 2.99792458e8      # speed of light in m/s

    #   lambda  = 1.e-9 * wavelength(x) # convert lambda from nm to m
    lambda  = wavelength(x)
    c       = c * 1.e9                  # convert speed of light to nm/s

    if( type == "light" )
        {
        #   p10     = 1.e18     # 18 is a multiple of 3        
        #   the most commmon actinometric unit for "light" is micromoles_of_photons/sec 
        #   a micromole of photons = 1.e-6*N_A photons
        #   the output unit for the light is watt
        N_A = 6.022140857e23    # Avogadro's Constant
        myfun <- function( y )  { y * 1.e-6*N_A*h*c / lambda }
        quantity.radio  = "power"
        }
    else if( type == "responsivity.light" )
        {
        #   the most commmon actinometric unit for "responsivity.light" is Quantum Efficiency
        #   the output unit for responsivity is A/W  or  J/C      
        e   = 1.6021766208e-19  # charge of electron, in C
        #   print( e/(h*c) )
        myfun <- function( y )  { y * lambda * e/(h*c) }
        quantity.radio  = sub( "photons", "power", quantity )    
        }
    else
        {
        log.string( FATAL, "Internal error.  Bad type='%s'.", type )
        return(NULL)
        }
    
    #   apply the function to all spectra
    out = applyspec.colorSpec( x, myfun )
    
    #   and change the quantity to be radiometric
    quantity(out)   = quantity.radio
    
    return( out )
    }
    
#--------       UseMethod() calls           --------------#            

radiometric <- function( x )
    {
    UseMethod("radiometric")
    }

    