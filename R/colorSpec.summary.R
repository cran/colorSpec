
#
#   x      colorSpec object
#

print.colorSpec  <-  function( x, ... )
    {
    return( summary.colorSpec( x, long=FALSE, ... ) )
    }


    
#   object  a colorSpec object
#   computes a pretty character string, prints it to stdout(), and then returns it invisibly

summary.colorSpec  <-  function( object, long=TRUE, ... )
    {    
    if( ! is.colorSpec(object) )
        {
        log_string( ERROR, "object is not a valid colorSpec object." )
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
        if( is.radiometric(object) )
            aka = "energy of photons, which is radiometric"
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
        data        = as.matrix( object )  #;    log_object( DEBUG, data  )
        
        #   compute matrix with some simple stats
        mat = matrix( as.numeric(NA), spectra, 4 )  #; print( dim(mat) )
        
        colnames(mat)           = c( "Min", "Max", "LambdaMax", integral )
        mat[ , "Min" ]          = apply( data, 2, min )
        mat[ , "Max" ]          = apply( data, 2, max )
        #   print( mat )
        mat[ , "LambdaMax" ]    = apply( data, 2, function(y) { ifelse( any(is.na(y)), NA_real_, mean(wavelength[y==max(y)]) ) } )
        
        mat[ , integral ]   = integral( object, method='rectangular' )
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

        #   calibrate
        out = c( out, charFromList("calibrate:",attr(object,'calibrate'))  )
           
        #   ptransform
        out = c( out, charFromList("ptransform:",attr(object,'ptransform'))  )
              
        #   sequence
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
            
            out = c( out, '', mess, '' )            
            }
        }
        
    cat( out, sep='\n' )
    
    return( invisible(out) )
    }
    
    
#   x       colorSpec object
#   method  'rectangular' or 'trapezoidal'
    
integral  <-  function( x, method='rectangular' )
    {
    wavelength  = wavelength(x)
    n   = length(wavelength)
    
    #   do not bother to optimize for regular wavelength sequence
#    weight  = diff( wavelength, lag=2 ) / 2
#    w.end   = c( wavelength[2] - wavelength[1], wavelength[n] - wavelength[n-1] )
#    if( tolower(method) == substr("trapezoidal",1,nchar(method)) )
#        w.end   = w.end / 2
#    weight  = c( w.end[1], weight, w.end[2] )     
 
    weight  = breakandstep( wavelength, method )$stepvec
 
    mat = as.matrix( x )
        
    #   weight  = matrix( weight, nrow=nrow(mat), ncol=ncol(mat) )
    
    out = colSums( weight * mat )    # weight is auto-replicated to all columns
    
    return( out )
    }
    
#   .title  of .list    
#   .list   a named list
#   return  a pretty character vector
    
charFromList <- function( .title, .list )
    {
    if( length(.list) == 0 )    return( character(0) )
    
    out = ''
    out = c( out, sprintf( "--------------- %s start ------------------", .title ) )
    
    for( k in 1:length(.list) )
        {
        if( 1 < k ) out = c( out, '' )     
                
        item    = .list[[k]]
        
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
    