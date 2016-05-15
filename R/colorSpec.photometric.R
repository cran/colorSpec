

#   photometric.colorSpec()
#
#   convert each spectrum of type 'light' to a single photometric number.
#   actinometric is converted to radiometric on-the-fly

photometric.colorSpec <- function( x )
    {
    if( type(x) != 'light' )
        {
        log.string( ERROR, "type(x) = '%s, but it must be 'light'", type(x) )
        return(NULL)
        }
    
    #   x might be photon-based (actinometric)    
    x   = radiometric( x )
    
    #   683 is the CIE-given photopic conversion factor (683.002 lumens/W is also used)
    out = 683 * product( x, subset(colorSpec::xyz1931.1nm,2), wave='auto' )
    
    return( out )
    }
    
#--------       UseMethod() calls           --------------#            

photometric <- function( x )
    {
    UseMethod("photometric")
    }

    