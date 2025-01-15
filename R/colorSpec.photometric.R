

#   photometric.colorSpec()
#
#   convert each spectrum of type 'light' to a single photometric number.
#   actinometric is converted to radiometric on-the-fly
#
#   x           colorSpec object with type 'light'
#   photopic    the Km conversion factor for photopic vision. 683 lumen/watt. 683.002 is also used.
#   scotopic    the Km' conversion factor for scotopic vision. 1700 lumen/watt.  1700.06 is also used.
#   multiplier  intended for unit conversion, and applies to both photopic and scotopic.

photometric.colorSpec <- function( x, photopic=683, scotopic=1700, multiplier=1 )
    {
    if( type(x) != 'light' )
        {
        log_level( ERROR, "type(x) = '%s, but it must be 'light'", type(x) )
        return(NULL)
        }
    
    #   x might be photon-based (actinometric)    
    x   = radiometric( x )
    
    out = product( x, colorSpec::luminsivity.1nm, wavelength='auto' )   #  partial argument match FIXED
    
    #   scale the columns
    phot    = grepl( "^photopic", colnames(out) )
    scot    = grepl( "^scotopic", colnames(out) )
    
    K   = photopic*phot  +  scotopic*scot
    
    if( multiplier != 1 )   K = multiplier * K
    
    K   = matrix( K, nrow(out), ncol(out), byrow=TRUE )
    #print(K)
    
    out = K * out 

    return( out )
    }
    
#--------       UseMethod() calls           --------------#            

photometric <- function( x, photopic=683, scotopic=1700, multiplier=1 )
    {
    UseMethod("photometric")
    }

    