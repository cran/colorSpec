
#
#   x      colorSpec object
#

print.colorSpec  <-  function( x, ... )
    {
    return( summary.colorSpec( x, long=FALSE, ... ) )
    }


    
#   object       a colorSpec object
#   computes a pretty character string, prints it to console, and then returns it invisibly

summary.colorSpec  <-  function( object, long=TRUE, ... )
    {    
    if( ! is.colorSpec(object) )
        {
        log.string( ERROR, "object is not a valid colorSpec object." )
        return(FALSE)
        }
    
    #   theName = deparse( substitute(object) ) ; print( str(theName) )
    #   theName = as.character( substitute(object,.GlobalEnv) )    ; print( str(theName) )
    #   return(FALSE)
    #   if( 1 < length(theName) )   theName = 'unknown'
        
    theName = 'the object'
        
    out = ''
    
    mess    = sprintf( "colorSpec object.   The organization is '%s'.  Object size is %d bytes.", 
                        organization(object), object.size(object) )
    out = c( out, mess )                        
    
    wavelength  = wavelength(object)
    range.wl    = range( wavelength )    
    
    #   type
    type        = type(object)
    quantity    = quantity(object)
    spectra     = numSpectra(object)
    
    
    if( type == "light" )
        {
        if( quantity == "power" )
            aka = "power of photons, which is radiometric"
        else
            aka = "number of photons, which is actinometric"
            
        if( spectra == 1 )
            mess    = sprintf( "%s describes a single source of light, and the quantity is '%s' (%s).", 
                            theName, quantity, aka )
        else
            mess    = sprintf( "%s describes %d sources of light, and the quantity is '%s' (%s).", 
                            theName, spectra, quantity, aka  )
                            
        integral    = "Integral"
        }
    else if( type == "responsivity.light" )
        {        
        if( quantity == "power" )
            aka = "power of photons, or radiometric"
        else
            aka = "number of photons"
            
        mess    = sprintf( "%s describes a responder to light with %d output channels, and the quantity is '%s'.", 
                            theName, spectra, quantity )
                            
        integral    = "E.response"                            
        }
    else if( type == "material" )
        {
        if( quantity == "reflectance" )
            material = "opaque"
        else
            material = "transparent"
            
        mess    = sprintf( "%s describes %d %s materials, and the quantity is '%s'.", 
                            theName, spectra, material, quantity )
                            
        integral    = "Integral"                            
        }
    else if( type == "responsivity.material" )
        {           
        if( 0 < spectra )             
        
        mess    = sprintf( "%s describes a responder to a material with %d output channels; the quantity is '%s'.", 
                            theName, spectra, quantity )
                            
        integral    = "perfect.response"                            
        }
        
    out = c( out, mess )
        
    
    n   = length(wavelength)
    
    step.wl     = attr(object,"step.wl")   
    if( is.null( step.wl ) )
        {
        step.wl = step.wl.colorSpec(object)
        range.step  = range( diff(wavelength) )
        mess = sprintf( "Wavelength range: %g to %g nm.  Step size is irregular; range is %g to %g nm; mean is %g nm.",
                            range.wl[1], range.wl[2], range.step[1], range.step[2], step.wl )
        }
    else
        mess = sprintf( "Wavelength range: %g to %g nm.  Step size is %g nm.",
                            range.wl[1], range.wl[2], step.wl  )
                            
    out = c( out, mess )
    
    out = c( out, '' )    
    out = c( out, sprintf( "%d spectra", spectra ) )
    out = c( out, sprintf( "%d data points / spectrum", length(wavelength) ) )
    
    extra   = ncol( extradata(object) )
    if( 0 < extra )
        out = c( out, sprintf( "%d columns of extra data", extra ) )
    
    if( 0 < spectra )
        {
        data        = as.matrix( object )  #;    log.object( DEBUG, data  )
        
        #   compute matrix with some simple stats
        mat = matrix( as.numeric(NA), spectra, 4 )  #; print( dim(mat) )
        
        colnames(mat)           = c( "Min", "Max", "LambdaMax", integral )
        mat[ , "Min" ]          = base::apply( data, 2, min )
        mat[ , "Max" ]          = base::apply( data, 2, max )
        #   print( mat )
        mat[ , "LambdaMax" ]    = base::apply( data, 2, function(y) { ifelse( any(is.na(y)), NA_real_, mean(wavelength[y==max(y)]) ) } )
        mat[ , integral ]       = base::apply( data, 2, sum ) *  step.wl 
        #print( mat )            
        
        df  = cbind( Band=specnames(object), as.data.frame(mat) )

        
        extra   = extradata( object )   
        if( ! is.null(extra)  &&  0 < ncol(extra) )
            df  = cbind( df, extra )
        
        rownames(df)    = 1:nrow(df)
        
        if( type == "light" )
            colnames(df)[1] = "Source"        
        else if( type == "responsivity.light" )
            colnames(df)[1] = "Channel"  
        else if( type == "material" )            
            colnames(df)[1] = "Material" 
        else if( type == "responsivity.material" )            
            colnames(df)[1] = "Channel"   
            
        out = c( out, '' )
        out = c( out, capture.output( print( df ) ) )
        }

        
    if( long )
        {
        #   metadata
        out = c( out, charFromList("metadata:",metadata(object))  )  

        #   calibration
        out = c( out, charFromList("calibration:",attr(object,'calibration'))  )
        
        
        #   sequence
        mess    = character(0)
        if( ! is.null( attr(object,"sequence") ) )
            {
            theSequence = attr(object,"sequence") 
            names       = names( theSequence )
            
            for( k in 1:length(theSequence) )
                {
                if( is.character( theSequence[[k]] ) )
                    #   enclose in quotes
                    names[k] = sprintf( '"%s"', names[k] )
                }

            desc    = paste( names, collapse=', ' ) #;  print( desc )
            mess    = sprintf( "Product Terms:  %s is a product of %d terms: %s.", theName, length(names), desc )
            }
        else
            {
            #   mess    = sprintf( "Product Terms:  %s is not a product of other colorSpec objects.", theName )
            }
        out = c( out, '', mess, '' )
        }
        
    cat( out, sep='\n', file=stderr() )
    
    return( invisible(out) )
    }
    
    
#   .title  of .list    
#   .list   a named list
#   return  a pretty character vector
    
charFromList <- function( .title, .list )
    {
    if( length(.list) == 0 )    return( character(0) )
    
    out = c( sprintf( "--------------- %s start ------------------", .title ) )
    
    for( k in 1:length(.list) )
        {
        item    = .list[[k]]
        
        out = c( out, '' )

        if( is.character(item) )
            {
            if( length(item) == 1 )
                out = c( out, paste( names(.list)[k], ' = ', item, collapse=' ' ) )
            else
                {
                out = c( out, paste( names(.list)[k], ' = ', collapse='' ) )
                out = c( out, item )
                }
            #out = c( out,  item  )
            }
        else
            {
            #   out = c( out, '' )
            out = c( out, paste( names(.list)[k], ' = ' ) )                
            out = c( out, capture.output( print( item ) ) )
            }
        }    
        
    out = c( out, sprintf( "--------------- %s end ------------------", .title ) )
           
    return(out)
    }
    