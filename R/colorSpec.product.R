
#   ...     colorSpec objects that can be naturally multiplied
#           there are 5 different types of sequences
#
#   wavelength      'identical'     means check that all objects have the same wavelength (default)
#                   'auto'          resample all of them to a suitable computed new sequence
#                   numeric vector, new wavelengths to resample to
#
#
#   return:     a new colorSpec object,
#               or a matrix with a column for each output channel and a row for each input spectrum

product.colorSpec <- function( ... )            #,  wavelength='identical', white.response='auto' )
    {
    theList =  list(...) 
    
    n   = length(theList)
    if( n == 0 )
        {
        log.string( ERROR, "No arguments !" )
        return(NULL)
        }
        
    #   assign defaults
    wavelength      = 'identical'

    #   get the LHS names
    theNames    = names( theList )
    
    if( ! is.null(theNames) )
        {
        #   extract just the non-named arguments as the colorSpec objects
        mask.noname = nchar(theNames) == 0        
        
        #   and so the complement are named arguments, which we want to pass to resample()
        mask.named  = ! mask.noname
        
        sig = pmatch( theNames, c("wavelength") )     
        
        idx = which( sig == 1 )
        if( length(idx) == 1 )
            {
            wavelength      = theList[[idx]]
            mask.named[idx] = FALSE
            }
        
        #   print( mask.named )
        
        #   all named args except wavelength !
        resample.argv   = theList[ mask.named ]     #; print( argv_for_resample )
            
        log.object( TRACE, wavelength )
        }
    else
        {
        mask.noname     = rep( TRUE, n )
        resample.argv   = NULL        
        }

    if( sum(mask.noname) == 0 )
        {
        log.string( ERROR, "No colorSpec arguments." )
        return(NULL)
        }
        
    #   now get the actual RHS names
    theNames = as.character( substitute(list(...)) )    # ; print( theNames )
    if( length(theNames) == n+1 )
        theNames    = theNames[ 2:(n+1) ]
    else
        {
        log.string( WARN, "length(theNames) = %d != %d.  Using fake names.", length(theNames), n+1 )
        theNames    = sprintf( "Name%d", 1:n )
        }
        
    #   form subset for processing
    theList     = theList[ mask.noname ]
    theNames    = theNames[ mask.noname ]
        
    names(theList)  = theNames  #; print( str(theList) )
    
    n   = length( theList )
    log.string( TRACE, "Found %d objects in '...'", n )
    if( n == 0 )
        {
        log.string( ERROR, "No colorSpec arguments." )
        return(NULL)
        }
        
    type.return = returnTypeProduct( theList )    
    
    if( is.na(type.return) )
        {
        log.string( ERROR, "Cannot form product, because the argument sequence is invalid." )        
        return(NULL)
        }
        
    if( n == 1 )
        {
        #   a single colorSpec, a trivial case
        return( theList[[1]] )
        }

        
    #   gather type, spectra, and quantity
    type        = character(n)
    quantity    = character(n)    
    spectra     = integer(n)
    k.character  = 0
    for( k in 1:n )
        {
        if( is.character(theList[[k]]) )    
            {
            k.character = k
            type[k]     = "character"
            quantity[k] = "character"     
            spectra[k]  = 1L
            }
        else
            {
            type[k]     = type( theList[[k]] )      
            quantity[k] = quantity( theList[[k]] )  
            spectra[k]  = numSpectra( theList[[k]] )
            }
        }
        
    #   check # of spectra
    if( any(spectra == 0) )
        {
        log.object( ERROR, spectra )
        log.string( ERROR, "one or more of of the %d colorSpec objects has 0 spectra.", n )
        return(NULL)
        }    
    
    if( type.return != "matrix" )
        {
        #   easier case
        spectra.multi   = spectra[ 1 < spectra ]
        spectra.unique  = unique( spectra.multi )
        if( 2 <= length(spectra.unique) )
            {
            log.object( ERROR, spectra )
            log.string( ERROR, "# of spectra in list is invalid." )
            return(NULL)
            }       
        }
    else
        {
        #   more complex, try to split sequence into a begin and end (left and right) parts
        #   each part must meet conditions in previous if block
        res = splitSequence( spectra )
        if( is.null(res) )
            {
            log.object( ERROR, spectra )            
            log.string( ERROR, "# of spectra in list is invalid." )
            return(NULL)
            }       
        seq.beg = res$begin
        seq.end = res$end
        }
        
    

    #   examine wavelengths
    if( wavelength[1] == 'auto' )   wavelength = NULL
    
    ok  = is.null(wavelength)  ||  is.character(wavelength)  ||  is.numeric(wavelength)
    if( ! ok )
        {
        log.string( ERROR, "wavelength type '%s' is invalid", typeof(wavelength) )
        return(NULL)
        }        

    do.resample = TRUE
    
    if( is.character(wavelength) )
        {
        if( wavelength[1] == "identical" )
            {
            #   verify that all the objects have identical wavelengths
            wavelength  = wavelength.colorSpec( theList[[1]] )

            for( k in 2:n )
                {
                if( k == k.character )  next
                
                if( ! identical( wavelength.colorSpec(theList[[k]]), wavelength ) )
                    {
                    log.string( ERROR, "'%s' does not have the same wavelengths as '%s'. Consider option wavelength='auto'.",
                                        theNames[1], theNames[k] )
                    return(NULL)
                    }            
                }
                
            do.resample = FALSE                
            }
        else
            {
            log.string( ERROR, "wavelength option '%s' is invalid.", wavelength[1] )
            return(NULL)
            }
        }
    else if( is.null(wavelength) )
        {
        #   compute a new wavelength sequence
        mask    = sapply( theList, is.colorSpec )    #; print( mask )
        
        #   sublist = theList[mask] ; print( str(sublist) )
        
        wavelength  = commonWavelength( lapply( theList[mask], wavelength.colorSpec )  )

        log.string( INFO, "From %d objects, computed %d common wavelengths: %g to %g nm with step of %g nm",
                        sum(mask), length(wavelength), wavelength[1], wavelength[length(wavelength)], wavelength[2]-wavelength[1] )
        }        
        
        
    if( ! isRegularSequence(wavelength) )
        {
        log.string( ERROR, "Cannot form product, because resampling wavelength sequence is not regular." )
        return(NULL)
        }        
        
    if( do.resample )
        {
        #   resample all of them
        log.string( INFO, "Resampling all %d objects...", n )
        
        #print( resample.argv )
        
        #   make arguments for resample().  Always wavelength and maybe more in resample.argv.
        argv = c( list(wave=wavelength), resample.argv )
        
        #print( str(argv) )
        
        for( k in 1:n )
            {
            if( k == k.character )    next
            
            #print( "-------------" )
            #print( k )
            #print( str(theList[[k]]) )
            
            #   modifyList() will add argv to the colorSpec argument
            argv.resample   = modifyList( list(theList[[k]]), argv )    #; print( str(argv.resample) )
            
            theList[[k]] = do.call( resample, argv.resample )     #      resample( theList[[k]], wave=wavelength )        
            }
        }

    #if( ! is.regular( theList[[1]] ) )
    #    {
    #    log.string( ERROR, "Cannot form product, because wavelength sequence is not regular." )
    #    return(NULL)
    #    }
        
    step.wl = step.wl( theList[[1]] )
                    
    #   convert to proper quantities
    for( k in 1:n )
        {
        if( k == k.character )    next
                    
        pattern = "^absorb"
        if( grepl( pattern, quantity[k] ) )
            {
            # convert from absorbance to transmittance
            log.string( DEBUG, "converting spectrum %d from absorbance to transmittance.", k )
            theList[[k]]    = linearize(theList[[k]])
            }
            
        pattern = "^photons"
        if( grepl( pattern, quantity[k] ) )
            {
            # convert from actinometric to radiometric
            log.string( DEBUG, "converting spectrum %d from actinometric to radiometric.", k )
            theList[[k]]   = radiometric(theList[[k]]) 
            }
        }


    if( type.return == "matrix" )
        {
        #   output is *not* a colorSpec
        #   output is a matrix - a product of a left part and a right part
        
        mat.left    = simpleProduct( theList[ seq.beg ] )
        mat.right   = simpleProduct( theList[ seq.end ] )
        
        out = step.wl * t(mat.left) %*% mat.right  #; print( str(out) )
        
        #   assign row names from the left part
        idx = which( spectra == nrow(out) )
        idx = idx[ 1 ]
        rownames(out)   = specnames(theList[[idx]]) 
        
        #   assign col names from the right part
        idx = which( spectra == ncol(out) )
        idx = idx[ length(idx) ]
        colnames(out)   = specnames(theList[[idx]])

        if( all( nchar(colnames(out)) == 1 ) )
            #   change 'x' to 'X', and 'r' to 'R', etc.
            colnames(out)   = toupper(colnames(out))            
            
        #   attr( out, "sequence" ) = theList
        
        return( out )
        }        
        
        
    if( k.character == 0 )
        {
        #   the case without a variable material is easier
            
        #   compute simple product of matrices 1 to n
        mat.cum = simpleProduct( theList )

        #   channels    = max( spectra )
        #   for the column names, pick the first one with max number of spectra
        idx = which.max( spectra )
            
        colnames(mat.cum)   = specnames( theList[[ idx ]] )
            
        if( type[1] == "light" )
            #   starts with light and ends with a material
            quantity    = quantity( theList[[1]] )
        else if( type[n] == "responsivity.light" )
            #   starts with material and ends with a responder
            quantity    = quantity( theList[[n]] )
        else
            {
            #   just a sequence of materials !  could be transmittance or reflectance !
            #   use crude rule to select which one
            if( any( quantity == "transmittance" ) )
                quantity    = "transmittance" 
            else
                quantity    = "reflectance"
            }
            
        out = colorSpec( mat.cum, wavelength, quantity )
        
        attr( out, "sequence" ) = theList
                    
        return( out )
        }
        
    
    #   there is a variable material "slot".  
    #   since materials are not fluorescent, the position of the "slot" does not matter.
    #   the type of the returned colorSpec will be "responsivity.material"
    if( type.return != "responsivity.material" )
        {
        log.string( FATAL, "type.return='%s' != '%s'.", type.return, "responsivity.material" )
        return(NULL)
        }
    
    quantity    = quantity( theList[[n]] )
    pattern     = "^power"
    if( ! grepl( pattern, quantity ) )
        {
        log.string( FATAL, "Internal error. Last object quantity='%s' does not begin with 'power'.", quantity )
        return(NULL)
        }    
    quantity    = sub( pattern, "material", quantity )

    #   compute simple product of matrices 1 to n, but skip the variable slot
    mat.left    = simpleProduct( theList[ 1:(k.character-1) ] )
    mat.right   = simpleProduct( theList[ (k.character+1):n ] )
    
    channels    = max( ncol(mat.left), ncol(mat.right) )
    
    if( ncol(mat.left) < channels )
        mat.left    = matrix( mat.left, nrow(mat.left), channels )
        
    if( ncol(mat.right) < channels )
        mat.right   = matrix( mat.right, nrow(mat.right), channels )
                
    mat.cum = mat.left * mat.right    
        
    #   channels    = ncol( mat.cum ) 
            
    specnames   = specnames( theList[[n]] )
    if( length(specnames) < channels )
        specnames   = sprintf( "%s.%d", specnames[1], 1:channels )
    
        
    colnames(mat.cum)   = specnames    
           
    out = colorSpec( mat.cum, wavelength, quantity )        
                  
    attr( out, "sequence" ) = theList    
    
    return( out )
    }
    
    
#   .wavelist   a list of wavelength vectors
#   returns     a suitable intersection vector of wavelengths
#               always begins and ends on integral nm, and always regular
commonWavelength <- function( .wavelist )
    {
    if( length(.wavelist) == 0 )    return(NULL)
    
    if( length(.wavelist) == 1 )    return(.wavelist[[1]])
    
    wmin    = max( sapply( .wavelist, function(x) { x[1] } ) )
    wmax    = min( sapply( .wavelist, function(x) { x[length(x)] } ) )
    
    wmin    = round(wmin)
    wmax    = round(wmax)
    
    if( wmax <= wmin )
        {
        log.string( ERROR, "wavelength intersection [%g,%g] is empty.", wmin, wmax )
        return(NULL)
        }
    
    wstep   = min( sapply( .wavelist, function(x) { mean(diff(x)) } ) )
    
    #   round to nearest power of 2
    if( wstep < 1 )
        wstep   = 2 ^ round( log2(wstep) )
    
    out = seq( wmin, wmax, by=wstep )
    
    return( out )
    }
    

    
#   .list   a list of valid colorSpec objects    
#   returns product matrix, with possible replication of columns    
simpleProduct <- function( .list )
    { 
    n   = length( .list )
    
    out = coredata( .list[[1]], forcemat=T )
    
    if( n == 1 )    return(out)
    
    #   verify spectra counts
    spectra     = sapply( .list, numSpectra.colorSpec )
    
    channels    = max( spectra )

    mask    =  spectra==1  |  spectra==channels
    if( ! all(mask) )
        {
        log.object( FATAL, spectra )
        log.string( FATAL, "Internal error.  Invalid spectra counts." )
        return(NULL)
        }
        
    if( ncol(out) < channels )
        out = matrix( out, nrow(out), channels )

    for( k in 2:n )
        {
        mat = coredata( .list[[k]], forcemat=T )

        if( ncol(mat) < channels )
            mat = matrix( mat, nrow(mat), channels )
        
        out = out * mat
        }
    
    return( out )
    }
    
    
#   .list       as passed to product.colorSpec()
#
#   returns type of product colorSpec objects
#   can be "matrix", or NA in case of an invalid .list    
#   does not check the number of spectra, which might invalidate .list later
    
returnTypeProduct <- function( .list )
    {
    n           = length( .list )
    
    theNames    = names( .list )
    if( is.null(theNames) )
        {
        theNames    = as.character( substitute(.list) )     #; print( theNames )
        theNames    = theNames[ 2:(n+1) ]
        }

    out         = as.character( NA )
    
    type        = character(n)
    k.character = 0
    
    for( k in 1:n )
        {
        if( is.character(.list[[k]]) )    
            {
            if( 0 < k.character )
                {
                log.string( ERROR, "String '%s' is too many strings.  There can be at most 1 string.", .list[[k]] )
                return(out)
                }
                
            k.character = k
                        
            if( k==1 || k==n )
                {
                log.string( ERROR, "String '%s' must appear in the interior of the argument list.", .list[[k]] )
                return(out)
                }
                
            type[k] = "character"
            next
            }
            
        # print( sprintf( 'class %s', paste(clist[[k]],collapse=',') ) )

        if( ! "colorSpec" %in% class( .list[[k]] ) )
            {
            log.string( ERROR, "class(%s) = '%s', which is invalid.", theNames[k], paste( class(.list[[k]]) ,collapse=',') )
            return(out)
            }
            
        if( ! is.colorSpec( .list[[k]] ) )
            {
            log.string( ERROR, "%s is not a valid colorSpec object.", theNames[k] )
            return(out)
            }        
            
        type[k]     = type( .list[[k]] )      
        }
        
    if( n == 1 )
        #   trivial case
        return( type[1] )
        
        
    if( 3 <= n )
        {
        #   check interior
        for( k in 2:(n-1) )
            {
            if( k == k.character )  next

            if( type[k] != "material" )
                {
                log.string( ERROR, "type(%s) = '%s' which is invalid for list interior.  It must be 'material'.", 
                                        theNames[k], type[k] )
                return(out)
                }    
            }        
        }
        
    #   good so far, now determine the return type
    if( type[1]=="material"  &&  type[n]=="material" )
        out = "material"
    else if( type[1]=="light"  &&  type[n]=="material" )
        out = "light"
    else if( type[1]=="material"  &&  type[n]=="responsivity.light" )
        out = "responsivity.light"
    else if( type[1]=="light"  &&  type[n]=="responsivity.light" )
        {
        #   2 cases
        if( k.character == 0 )
            out = "matrix"
        else
            out = "responsivity.material"   # variable material slot
        }
    else if( type[1]=="material"  &&  type[n]=="responsivity.material" )
        out = "matrix"
    else
        {
        log.string( ERROR, "types of colorSpec objects are invalid" )
        }
        
    return( out )
    }
    
    
#   x   a sequence of 2 or more positive integers
#       split x into 2 parts, so each one has at most 1 unique value > 1    
#   returns a list with:
#       begin   a sequence from 1:k
#       end     a sequence from (k+1):n
#   or NULL in case that x cannot be split
splitSequence <- function( x )
    {
    n   = length(x)
    
    out = list()
    out$begin   = 1:(n-1)
    out$end     = n
    
    x.multi = x[ 1 < x ]
    
    if( length(x.multi) <= 1 )
        #   1 value > 1,  or none of them
        return(out)
        
    jump    = which( diff(x.multi) != 0 )
    
    if( length(jump) == 0 )
        #   no jumps, so x.multi is constant
        return( out )
    
    if( 1 < length(jump) )  
        #   cannot split
        return(NULL)
        
    #   only 1 jump
    
    #   find last occurence of x1
    k   =   which( x == x.multi[1] )
    k   =   k[ length(k) ]
    
    out$begin   = 1:k
    out$end     = (k+1):n
    
    return( out )
    }
    
#--------       UseMethod() calls           --------------#                    
        
product <- function( ... )      #,   wavelength="identical", white.response='auto' )
    {
    UseMethod("product")
    }    
    
    
    
    
#--------       testing           --------------#           
    
testVarArg <- function( ... )  #,   wavelength="identical", white.response='auto' )
    {
    theList = list(...)
    print( as.character( substitute(list(...)) ) )
    
    theNames = names(theList)
    print( theNames )
    
    #   split into unnamed and named parts
    if( is.null(theNames) )
        {
        list.noname = theList
        list.named  = NULL
        }
    else
        {
        mask    = nchar(theNames)==0
        
        list.noname = theList[ mask ]
        list.named  = theList[ ! mask ]
        }
    
    cat( "No names:\n" )
    print( str(list.noname) )
    
    cat( "\nNamed:\n" )
    print( str(list.named) )
    
        
    
    return(F)
    
    print( str(substitute(...)) )
    print( deparse(substitute(...)) )
    print( substitute(..2) )    
    print( str(...) )
    
    print( deparse(substitute(..1,parent.frame(n=2) ) ) )    
    print( deparse(substitute(..1) ) )    
        
    theList =  list(...) 
    print( str(theList) )
    #   print( str(pairlist(...)) )
    print( names(theList) )

    for( k in 1:length(theList) )
        print( deparse(substitute( theList[k], parent.frame(n = 2) )  ) )
    
    
    #   print( str( c(...) ) )
    
    #   print( white.response )
    print( wavelength )
    
    print( lapply( theList, class ) )
    
    return( invisible(TRUE) )
    }
        