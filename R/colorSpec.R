
#   data   vector or matrix  (univariate or multivariate spectra)
#           if data is a matrix that has column names, these are used as spectrum names
#           if data is a vector then deparse(substitute(data)) is used as the spectrum name

colorSpec  <-  function( data, wavelength, quantity='auto', organization='auto', specnames=NULL )
    {
    n = NROW(data)     # number of wavelengths
    
    #   log_object( INFO, .type )

    if( ! is.numeric(data)  )
        {
        log_string( ERROR, "data has incorrect type '%s'.", typeof(data) )
        return(NULL)
        }    
        
    if( ! is.numeric(wavelength)  )
        {
        log_string( ERROR, "wavelength has incorrect type '%s'.", typeof(wavelength) )
        return(NULL)
        }    
    
    if( length(wavelength) != n )
        {
        log_string( ERROR, "NROW(data) = %d != %d = length(wavelength) mismatch.", n, length(wavelength) )
        return(NULL)
        }     
        
    if( is.integer(wavelength) )
        #   avoid integers, make all wavelengths doubles
        wavelength = as.numeric(wavelength)

        
    if( ! isIncreasingSequence(wavelength) )
        {
        log_string( ERROR, "wavelength is not increasing." )
        return(NULL)
        } 

    #   specnames   = character(0)
    
    if( is.null( dim(data) ) )
        {
        #   data is a plain vector
        spectra     = 1

        if( is.null(specnames) )
            specnames = deparse(substitute(data))
        
        attr(data,'specname')  = NULL      # to prevent propagation below
        }
    else if( is.matrix(data) )
        {
        #   data is a matrix
        spectra     = ncol(data)      

        if( is.null(specnames) )        
            specnames   = colnames(data)   #; print( specnames )
        }
    else
        {
        log_string( ERROR, "data is numeric, but neither vector nor matrix." )
        #   log_object( ERROR, data, type='str' )
        return(NULL)
        }
                

    #   check validity of specnames
    dups    = sum(duplicated(specnames))
    ok      = is.character(specnames)  &&  length(specnames)==spectra  &&  dups==0
    
    if( ! ok )
        {
        if( 0 < spectra )
            {
            #   make up something simple, and issue a warning
            specnames   = sprintf( "S%d", 1:spectra )    
            combo       = paste(specnames,collapse=',')
            if( 80 < nchar(combo) ) combo = paste( substr(combo,1,76), '...', collapse='' )
            
            if( 0 < dups )
                log_string( WARN, "specnames are invalid (there are %d duplicates).  Using '%s' instead.", 
                            dups, combo )
            else
                log_string( WARN, "specnames are invalid.  Using '%s' instead.", combo )
            }
        else
            {
            #   assign the only possibility for specnames, silently
            specnames = character(0)
            }            
        }

    if( quantity == "auto" )
        {
        #   guess the quantity
        quantity = guessSpectrumQuantity( specnames, '' )
        if( is.na( quantity ) )
            {
            log_string( ERROR, "Cannot guess spectrum quantity from the specnames." )
            return(NULL)
            }    
        }

    if( is.na( spectrumTypeFromQuantity( quantity ) ) )
        {
        log_string( ERROR, "quantity='%s' is invalid.", quantity )
        return(NULL)
        }    
        
    
    if( organization == "auto" )
        organization   = ifelse( is.null( dim(data) ), "vector", "matrix" )
        
    ok  = organization %in% c("df.col","df.row","matrix","vector")       
    if( ! ok )
        {
        log_string( ERROR, "organization='%s' is invalid\n", organization )
        return(NULL)
        }
        
    if( organization == "vector" )
        {
        if( spectra != 1 )
            {
            log_string( ERROR, "organization='%s' is invalid, because there are %d spectra in data\n", organization, spectra )
            return(NULL)
            }        
        out         = as.numeric(data)
        #   names(out)  = as.character(wavelength)
        }
    else if( organization == "matrix" )
        {
        out             = data
        dim(out)        = c(n,spectra)
        #   rownames(out)   = as.character(wavelength)
        }
    else if( organization == "df.col" )
        {
        out = data.frame( Wavelength=wavelength )  #; print( out )   

        if( 0 < spectra )
            out = cbind( out, as.data.frame(data) )
            
        row.names(out)  = NULL        # or rownames ?
        }
    else if( organization == "df.row" )
        {
        tdata   = t(data)      #;      print( tdata )
        colnames(tdata) = wavelength
        #   class(tdata)    = "model.matrix"
        out = data.frame( row.names=specnames )    #; print( out )
        out$spectra = tdata
        #   row.names(out)   = sprintf( "S%d", 1:nrow(out) )        
        }
        
    class( out )    = c( "colorSpec", class(out) )
            
    #   assign wavelength, quantity, and specnames
    wavelength( out )   = wavelength
    quantity( out )     = quantity           
    specnames( out )    = specnames
        
    attr( out, "metadata" ) = list()    #   initially empty
    
    return( out )
    }
    
is.colorSpec <- function(x)
    {
    if( ! "colorSpec" %in% class(x) )   
        {
        log_string( DEBUG, "'%s' class '%s' is invalid.", deparse(substitute(x)), class(x) )
        return(FALSE)
        }

    if( is.na(organization(x)) )
        {
        log_string( DEBUG, "'%s' organization '%s' is invalid.", deparse(substitute(x)), organization(x) )
        return(FALSE)
        }
        
    if( ! type(x) %in% c("light","responsivity.light","material","responsivity.material")  )
        {
        log_string( DEBUG, "'%s' type '%s' is invalid.", deparse(substitute(x)), type(x) )
        return(FALSE)
        }
        
    wave    = wavelength(x) 
    if( ! is.double( wave ) )  
        {
        log_string( DEBUG, "'%s' wavelengths are not double-precision.", deparse(substitute(x)) )
        return(FALSE)
        }

    if( ! isIncreasingSequence(wave) )
        {
        log_string( DEBUG, "'%s' wavelengths are not increasing.", deparse(substitute(x)) )
        return(FALSE)
        }        
        
    #qvalid  = validQuantities( type )
    #if( ! quantity(x) %in% qvalid  )
    #    return(FALSE)
   
    return(TRUE)
    }
        
as.colorSpec.default <- function( ... )
    {
    log_string( WARN, "This function is designed to be called from other packages." )
    return(NULL)
    }
    
isValidQuantity <- function( .quantity )
    {
    return( ! is.na( spectrumTypeFromQuantity( .quantity ) ) )
    }
    
isValidType  <- function( .type )
    {
    return( .type %in% c("light","responsivity.light","material","responsivity.material") )
    }
    
    
guessSpectrumType  <-  function( .specnames, .header )    
    {
    pattern = "keyword[: \t]+scanner"
    if( any( grepl(pattern,.header,ignore.case=TRUE) ) )    return( "responsivity.material" )
       
    pattern = "power|radian|photon"
    if( all( grepl(pattern,.specnames,ignore.case=TRUE) ) )     return( "light" )
    
    #   print( .header )
    if( any( grepl( "^IT8", .header ) ) )   return( "material" )    # an IT8 target
    
    if( 1 < length(.specnames) )
        {
        pattern_vec = c("^R","^G","^B","^X","^Y","^Z","^L","^M","^S","^V")  #,"^S")
        
        if( allDistinctMatches( pattern_vec, .specnames, .ignore=T ) )    return( "responsivity.light" )
        }
    else
        {
        pattern = "quantum|QE|gray"
        if( all( grepl(pattern,.specnames,ignore.case=TRUE) ) )     return( "responsivity.light" )
        }
        
    pattern = "trans|reflect|absorb"
    if( all( grepl(pattern,.specnames,ignore.case=TRUE) ) )     return( "material" )
    if( any( grepl(pattern,.header,ignore.case=TRUE) ) )        return( "material" )
    
    #   give up !
    return( as.character(NA) )
    }
    

spectrumTypeFromQuantity <- function( .quantity )
    {
    #   first look for the responders
    pattern = "^(power|energy|photons|photons/sec)->(electrical|neural|action)$"
    if( grepl( pattern, .quantity ) )   return( "responsivity.light" )
    
    pattern = "^material->(electrical|neural|action)$"
    if( grepl( pattern, .quantity ) )   return( "responsivity.material" )
    
    pattern = "^(power|energy|photons|photons/sec)$"
    if( grepl( pattern, .quantity ) )   return( "light" )
    
    pattern = "^(reflectance|transmittance|absorbance)$"
    if( grepl( pattern, .quantity ) )   return( "material" )
        
    return( as.character(NA) )
    }
    
    
guessSpectrumQuantity <- function( .specnames, .header )
    {
    #   first look for a direct assignment of QUANTITY, not a guess
    pattern = '^#?[ \t]*QUANTITY:?[ \t]*"?([->A-Z]+)"?'
    out = sub( pattern, "\\1", .header, ignore.case=T )
    
    idx = which( nchar(out) < nchar(.header) )
    if( 0 < length(idx) )
        {
        #   found 1 or lines matching
        value   = tolower( out[idx[1]] )
        
        if( isValidQuantity(value) )    return(value)   # got it !
        
        #   try a little harder with photons
        #   it might be 'photon count', 'photon flux', 'photon-intensity', 'photon radiance' etc.
        pattern = "photon[- ]+|photon$"
        if( grepl( pattern, value ) )   return('photons')
        
        #   else just keep going
        #   return( tolower( out[idx[1]] ) )
        }
        
        
    #   first look for MEAS_TYPE, which appear in newer versions of ArgyllCMS .sp files
    pattern = '^MEAS_TYPE[ \t]*"?([A-Z]+)"?'
    out = sub( pattern, "\\1", .header, ignore.case=T )
    
    idx = which( nchar(out) < nchar(.header) )
    if( 0 < length(idx) )
        {
        #   found 1 or lines matching, only look at the first one
        value   = toupper( out[idx[1]] )
        
        if( value == "REFLECTIVE" )     return( 'reflectance' )
            
        if( value == "TRANSMISSIVE" )   return( 'transmittance' )
            
        if( value=="EMISSION"  ||  value=="AMBIENT" )   return( 'energy' )
            
        if( value == "SENSITIVITY" )   return( 'energy->neural' )
                        
        #   else just keep looking
        }
        
    #   scanners are electrical
    pattern = "keyword[: \t]+scanner"
    if( any( grepl(pattern,.header,ignore.case=TRUE) ) )    return( 'material->electrical' )    
    
    
    type    = guessSpectrumType( .specnames, .header )      ; log_string( DEBUG, "guessed type=%s", type )
    
    if( is.na(type)  ||  type == "light" )
        {
        #   photon count
        if( any( grepl( "photon", .specnames, ignore.case=T ) ) )   return( "photons" )
        
        #   energy
        if( any( grepl( "^(energy|power)", .specnames, ignore.case=T ) ) )   return( "energy" )
                
        #   photon count is rare, so assume "energy"
        return( "energy" )
        }
    else if( is.na(type)  || type == "responsivity.light" )
        {
        #   quantum
        if( all( grepl( "^(quant|QE)", .specnames, ignore.case=T ) ) )      return( 'photons->electrical' )
        
        if( any( grepl( "photon", .header, ignore.case=T ) ) )              return( "photons->action" )
        
        if( length(.specnames)==1  &&  grepl( "gray", .specnames, ignore.case=T ) )    return( 'energy->electrical' )
        
        #   RGB is electrical
        if( allDistinctMatches( c("^R","^G","^B"), .specnames, .ignore=T ) ) return( 'energy->electrical' )
        
        #   XYZ is neural
        if( allDistinctMatches( c("^x","^y","^z"), .specnames, .ignore=T ) ) return( 'energy->neural' )
                
        #   LMS is neural                
        if( allDistinctMatches( c("^L","^M","^S"), .specnames, .ignore=T ) ) return( 'energy->neural' )
                                
        #   V is neural - the definition of luminance               
        if( allDistinctMatches( "^V", .specnames, .ignore=T ) ) return( 'energy->neural' )
                                
        return( 'power->action' )        # grab-bag
        }
    else if( is.na(type)  ||  type == "material" )
        {
        #   look for "Transmittance" or "Reflectance" or "Absorbance"
        if( any( grepl( "reflect", .specnames, ignore.case=T  ) ) )     return( "reflectance" )
        
        if( any( grepl( "transmit", .specnames, ignore.case=T  ) ) )    return( "transmittance" )
        
        if( any( grepl( "absorb", .specnames, ignore.case=T  ) ) )      return( "absorbance" )     

        if( any( grepl( "IT8[.]7/1", .header ) ) )      return( "transmittance" )
        if( any( grepl( "IT8[.]7/2", .header ) ) )      return( "reflectance" )            
        
        if( any( grepl( 'QUANTITY[ \t]+"?TRANSMIT', .header, ignore.case=T  ) ) )   return( "transmittance" )
        if( any( grepl( 'QUANTITY[ \t]+"?REFLECT', .header, ignore.case=T  ) ) )    return( "reflectance" )
        }
    else if( is.na(type)  ||  type == "responsivity.material" )
        {
        #   RGB is electrical
        if( allDistinctMatches( c("^R","^G","^B"), .specnames, .ignore=T ) ) return( 'material->electrical' )
        
        #   XYZ is neural
        if( allDistinctMatches( c("^x","^y","^z"), .specnames, .ignore=T ) ) return( 'material->neural' )
                
        #   LMS is neural                
        if( allDistinctMatches( c("^L","^M","^S"), .specnames, .ignore=T ) ) return( 'material->neural' )

        return( 'material->action' )        # grab-bag
        }
    
    #   give up !
    #   log_object( DEBUG, type )
    log_object( DEBUG, .specnames )    
    #   log_object( DEBUG, .header )       
    log_string( DEBUG, "Cannot guess quantity from .specnames and .header." )
    
    return( as.character(NA) )
    }
    
    
validQuantities <- function( .type )
    {    
    qvalid  = as.character(NA) ; 
    
    if( .type == "light" )
        qvalid  = c("energy","power","photons")
    else if( .type == "responsivity.light" )
        qvalid  = c("energy->electrical","energy->neural","energy->action","power->electrical","power->neural","power->action","photons->electrical","photons->neural","photons->action")
    else if( .type == "material" )
        qvalid  = c("reflectance","transmittance","absorbance")
    else if( .type == "responsivity.material" )
        qvalid  = c("electrical","neural","action")
     
    return( qvalid )
    }

    
numSpectra.colorSpec <- function(x) 
    {
    org = organization.colorSpec(x)   #; print(org)
    
    if( org == "vector" )
        return(1L)
    else if( org == "matrix" )
        return( ncol(x) )
    else if( org == "df.col" )
        return( ncol(x)-1L )     # because the 1st column is Wavelength
    else if( org == "df.row" )
        return( nrow(x) )
    
    return( 0 )
    }
    
#   in this one, ignore the wavelengths attribute and use x itself    
numWavelengths.colorSpec <- function(x) 
    {
    org = organization.colorSpec(x)   #; print(org)
    
    if( org == "df.row" )
        return( ncol( x[[ ncol(x) ]] ) )
    else
        return( NROW(x) )
    }    
        
organization.colorSpec <- function(x) 
    {
    #   print( dim(x) )
    
    if( is.null(dim(x)) )  return( "vector" )
    
    if( is.matrix(x) )
        return( "matrix" )
        
    if( is.data.frame(x) )
        {
        #   cname   = colnames(x)
        #   cname   = colnames(x,do.NULL=FALSE)
        
        #   class.last  = class( x[[ncol(x)]] )
        
        if( is.matrix( x[[ncol(x)]] ) )     # "model.matrix"  %in%  class.last  ||  "matrix"  %in%  class.last    )
            #   the last column is a matrix, it has 2 dimensions
            return( "df.row" )
        else if( grepl( "^wave|^wl|nm$", colnames(x)[1], ignore.case=T ) )
            return( "df.col" )        
        }
        
    log_object( ERROR, x, 'str' )
    log_string( ERROR, "cannot determine organization of object x." )
    
    return( as.character(NA) )
    }
    
    
    
    
        
#   x       colorSpec object
#   value   the new organization, so x may be destroyed and recreated

"organization<-.colorSpec" <- function( x, value ) 
    {
    ok  = value %in% c("vector","matrix","df.row","df.col")        
    if( ! ok )
        {
        log_string( ERROR, "organization='%s' is invalid.", value )
        return(x)
        }    
    
    org = organization(x)   #; print(org)   #  the current organization

    if( org == value )
        return( x ) # nothing to do !
        
    m   = numSpectra(x)  #; print(m)
    
    if( value == "vector"  &&  m != 1 )
        {
        #   return x unchanged
        log_string( ERROR, "cannot reorganize to vector because #(spectra) = %d != 1.", m )
        return( x )
        }
            
    #   create a new object
    out = colorSpec( as.matrix(x), wavelength(x), quantity=quantity(x), organization=value )

    if( ! is.null(out) )
        {
        #   assign special attributes
        for( a in c('metadata','sequence','calibrate','emulate') )
            attr(out,a) = attr(x,a)
        }
        
    return(out)
    }        
    
    
#   return value:   a vector when a single spectrum, unless forcemat==TRUE
#                   a matrix when multiple spectra
coredata.colorSpec <- function( x, forcemat=FALSE )
    {
    wave        = wavelength(x)    
            
    if( numSpectra(x) == 0 )
        {
        #   a special case
        out = matrix( 0, length(wave), 0 )
        rownames(out)   = as.character(wave)
        return( out )
        }

    org         = organization.colorSpec(x)   #; print(org)

    if( org == "vector" )
        {
        out = as.double(x)         
        if( forcemat )
            dim(out)    = c(length(wave),1)     # out is now a matrix with 1 column, a column vector
        else
            names(out)  = as.character(wave)
        }
    else if( org == "matrix" )
        {
        out = x
        }
    else if( org == "df.col" )
        {
        out =   as.matrix.data.frame( x[ , 2:ncol(x), drop=F ] ) #; print( str(out) )     # because the 1st column is Wavelength
        }
    else if( org == "df.row" )
        {
        out = t( x[[ncol(x)]] )
        #   class( out ) = "matrix"
        }
    else
        {
        log_string( FATAL, "Internal error. organization='%s' is unknown.", org )
        return(NULL)
        }

    class( out ) = 'numeric'    # not a data.frame anymore
    
    #   erase all colorSpec attributes
    for( a in c( "wavelength","step.wl","quantity","metadata","sequence","calibrate","specname",'emulate') )
        attr(out,a) = NULL    
        
    if( is.matrix(out) )
        {
        #   assign the rownames and colnames
        colnames(out)   = specnames(x)          # this one is usually not necessary, since variable out already has these colnames
        rownames(out)   = as.character(wave)    # this one is new
            
        #   erase attribute meant for vectors only
        #   attr( out, "specname" ) = NULL
        }
        
    return(out)
    } 

as.matrix.colorSpec <- function( x, ... )    
    {
    #   log_string( TRACE, "as.matrix.colorSpec()" )
    return( coredata(x,forcemat=TRUE) )
    }
    
    
#   x               a colorSpec object
#   organization    'col', 'row', or 'auto'

as.data.frame.colorSpec  <-  function( x, row.names=NULL, optional=FALSE, organization='auto', ... )
    {
    if( organization == 'auto' )
        {
        organization    = ifelse( organization(x)=='df.row', 'row', 'col' )
        }

    out = x
    
    if( organization == 'col' )
        organization(out) = 'df.col'        
    else if( organization == 'row' )
        organization(out) = 'df.row'
    else
        {
        log_string( ERROR, "organization='%s' is invalid.", organization )
        return(NULL)
        }
        
    #   out is now a data.frame, but also a colorSpec.
    #   make it a data.frame only !
    class(out) = 'data.frame'

    return( out )
    }    
    
    
    
    
#   only applies if organization is "df.row"    
#   return value:  a data.frame = the left part before the column "spectra"
extradata.colorSpec <- function( x )
    {
    org = organization(x)   #; print(org)

    if( org != "df.row" )
        {    
        #   log_string( WARN, "organization = '%s' is invalid.", org )
        #   return data.frame with 0 columns
        return( data.frame( row.names=specnames(x) ) )
        }
        
    #   x is a data.frame
        
    m   = ncol(x)
    if( m == 1 )
        {    
        #   log_string( WARN, "Object '%s' has no extra data.", deparse(substitute(x)) )
        #   return data.frame with 0 columns        
        return( data.frame( row.names=row.names(x) ) )
        }
    
    out = x[ 1:(m-1) ]  # drop the spectra matrix itself and keep only the left part
    
    class(out)  = 'data.frame'      #   not a colorSpec anymore
    
    return( out )
    }
    
#   only applies if organization is "df.row"    
#   return value:  a data.frame = the left part before the column "spectra"
"extradata<-.colorSpec" <- function( x, add=FALSE, value )
    {
    #   print( str(value) )
    #   cat( sprintf( "type='%s'\n", typeof(value) ) )
    if( ! is.null(value)  &&  ! is.data.frame(value) )
        {    
        log_string( ERROR, "RHS object '%s' is not a data.frame, and not NULL.", deparse(substitute(value))[1] )
        return(x)
        }    
    
    if( is.null(value)  &&  add )
        {
        # no change to x, return silently
        return(x)   
        }
        
    if( is.data.frame(value)  &&  ncol(value)==0  )
        {
        # no data so no change to x, return silently
        return(x)   
        }
    
    org = organization(x)   #; print(org)

    if( org != "df.row" )
        {    
        log_string( ERROR, "organization(x) = '%s' is invalid.  Please change to 'df.row' first.", org )
        return(x)
        }

    if( is.data.frame(value) )
        {
        if( nrow(value) != nrow(x) )
            {    
            log_string( ERROR, "Row count mismatch. LHS %d != %d RHS.",  nrow(x), nrow(value) )
            return(x)
            }
            
        if( add )
            {
            #   check for duplicate column names
            namevec = unlist( sapply( list(value,x), colnames ) )
            idx     = which( duplicated(namevec) )
            if( 0 < length(idx) ) 
                {
                log_string( ERROR, "Cannot add new extradata, because one of its column names already appears in the current extradata, e.g. '%s'.",
                                namevec[idx[1]] )
                return(x)
                }       

            #   stick the extradata in front of the spectrum, and preserve any existing extradata  
            if( ncol(x) == 1 )
                out = cbind( value, x )
            else
                out = cbind( x[1:(ncol(x)-1)], value, x[ncol(x)] )  # stick new extradata between current and spectrum
            }
        else
            {
            #   stick the extradata in front of the spectrum, and erase any existing extradata
            out = cbind( value, x[ncol(x)] )
            }
        }
    else
        {
        #   value must be NULL, and add must be FALSE
        #   discard all extradata, and just keep the spectrum
        out = x[ncol(x)]
        }

    class(out)  = class(x)
    
    for( a in c('wavelength','step.wl','quantity','metadata','sequence','calibrate','emulate') )
        {
        attr( out, a ) = attr( x, a )
        }
        
    specnames(out)  = specnames(x)
        
    return( out )
    }
    
wavelength.colorSpec  <- function(x)    
    {
    wave    = attr(x,"wavelength")
    
    if( is.null(wave) && is.data.frame(x) )
        {
        #   use 1st column
        wave    = x[[1]]
        
        if( is.null(wave) )
            log_string( FATAL, "Internal Error. Cannot determine wavelength of object '%s'", deparse(substitute(x)) )
        }
        
    return( wave )
    }
    
"wavelength<-.colorSpec"  <- function( x, value )    
    {    
    n   = numWavelengths.colorSpec( x )
    
    if( length(value) != n )
        {
        log_string( ERROR, "length(wavelength)=%d is invalid.", length(value) )
        return(x)
        }
        
    if( ! isIncreasingSequence( value ) )
        {
        #   log_object( ERROR, value, 'str' )
        log_string( ERROR, "wavelength sequence is not increasing." )
        return(x)
        }
        
    #sym = deparse(substitute(x))        
    #if( bindingIsLocked(sym,environment(wavelength.colorSpec)) )
    #    {
    #    log_string( ERROR, "Cannot modify '%s', because it is locked.", sym )
    #    return(NULL)
    #    }
        
    org = organization.colorSpec( x ) 
       
    if( org == "df.col" )
        {
        #   wavelengths are stored in 1st column
        x[[1]]  = as.numeric(value)
        }
    else if( org == "df.row" )
        {
        #   character form of wavelengths are also stored in sub-colnames of the final column
        wavecurrent = as.numeric( colnames( x[[ncol(x)]] ) )
        if( any( 0.5 < abs(wavecurrent - value) ) )
            #   make the change
            colnames( x[[ncol(x)]] )  = as.character( value )              
        }       
       
    if( org != "df.col" )       
        attr( x, "wavelength" )   = as.numeric(value)

    #   check regularity and set "step.wl" as appropriate
    if( isRegularSequence(value) )
        step    = ifelse( 1 < n, (value[n] - value[1])/(n-1), 0 )
    else
        step = NULL
        
    attr( x,"step.wl")  = step
         
    return( x )
    }
        
type.colorSpec  <- function(x)    
    {
    return( spectrumTypeFromQuantity( quantity(x) ) )     # attr(x,"type") )
    }
    
quantity.colorSpec  <- function(x)    
    {
    return( attr(x,"quantity") )
    }
    
"quantity<-.colorSpec"  <- function( x, value )    
    {  
    if( ! isValidQuantity( value ) )
        {
        log_string( ERROR, "quantity=%s is invalid\n", value )
        return(x)
        }        
    
    attr(x,"quantity") = value
    
    return(x)
    }
    
specnames.colorSpec <- function(x) 
    {
    spectra = numSpectra.colorSpec(x)
    
    if( spectra == 0 )  return( character(0) )
    
    org = organization.colorSpec(x)   #; print(org)
           
    if( org == "vector"  )
        #   cannot use column names, so use special attr
        return( attr( x, "specname" ) )
        #   return( deparse(substitute(x)) )
    else if( org == "matrix" )
        return( colnames(x,do.NULL=F) )
    else if( org == "df.col" )
        {
        #   cname   = attr(x,"names")
        cname = colnames(x,do.NULL=F)
        return( cname[ 2:ncol(x) ] )     # because the 1st column is Wavelength
        }
    else if( org == "df.row" )
        {
        return( row.names(x) )
        }
        
    return( rep(as.character(NA),spectra) )
    }

"specnames<-.colorSpec" <- function(x, value) 
    {
    n   = numSpectra.colorSpec(x)
    
    if( length(value) != n )
        {
        log_string( ERROR, "length(value) = %d is invalid. It is not equal to numSpectra = %d", 
                            length(value), n )
        return(x)
        }
        
    if( n == 0 )    return(x)   # cannot assign any names
        
    #   check for dups
    count   = sum(duplicated(value))
    if( 0 < count )
        {
        log_string( ERROR, "The %d names for the spectra have %d duplicates.  Names are ignored.", 
                    n, count)        
        return(x)
        }
        
    org = organization.colorSpec(x)   #; print(org)        
        
    if( org == "vector" )
        {
        #   cannot use column names, so use special attr        
        attr( x, "specname" )   = value
        }
    else if( org == "matrix" )
        colnames(x) = value
    else if( org == "df.col" )
        colnames(x)[ 2:ncol(x) ] = value    # because the 1st column is Wavelength
    else if( org == "df.row" )
        row.names(x) = value
    
    return(x)        
    }    
    

is.regular.colorSpec <- function(x)
    {    
    return( ! is.null( attr(x,"step.wl") ) )
    }
    
step.wl.colorSpec  <-  function(x)
    {
    out = attr(x,"step.wl")    
    
    if( ! is.null(out) )    return(out)
    
    wave        = wavelength(x)
    n           = length(wave)
    range.wl    = range( wave )
    out =  ifelse( 2<=n, (range.wl[2]-range.wl[1])/(n-1), 1 )        

    return(out)
    }
    
metadata.colorSpec <- function( x, ... )
    {
    dots <- c(...) 
    
    metadata    = attr( x, "metadata" ) #   ; print(metadata)
    
    if (length (dots) == 0L)
        #   return the whole list
        metadata
    else if( length (dots) == 1L)
        #   return single item 
        metadata[dots][[1]]
    else
        #   return sublist
        metadata[dots]
    }    
    

"metadata<-.colorSpec" <- function( x, add=FALSE, value )
    {
    #   log_object( DEBUG, value)
    
    mask    = nzchar( names(value) )

    if( ! all (mask) )
        log_string( WARN, "options without name are discarded: %d", which(!mask) )

    if( add )
        {
        metadata    = attr( x, "metadata" ) #; log_object( DEBUG, metadata )
        attr(x,"metadata") <- modifyList( metadata, value[mask] )
        }
    else
        {
        attr(x,"metadata") <- value[mask]
        }
        
    return( x )
    }
    

    
if( 0 )
{    
#   in this one allow for monochromats and dichromats, but not tetrachromats    
isHuman.colorSpec <- function(x)
    {    
    m   = numSpectra(x)
    if( m<1  || 3<m )   return(FALSE)
    
    if( ! grepl( 'neural$', quantity(x) ) )  return(FALSE)
    
    #   finally check names
    theNames    = tolower(specnames(x))
    
    return( all( theNames %in% c('x','y','z') )  ||  all( theNames %in% c('l','m','s') ) )
    }
}        
    
    
    
#   force wavelengths to be equally spaced, keeping the endpoints fixed
regularizeWavelength  <-  function( .wavelength )
    {
    n   = length( .wavelength )
    
    if( n <= 1 )    return( .wavelength )
    
    return( seq( .wavelength[1], .wavelength[n], length.out=n ) )
    }

    
isRegularSequence <- function( x, epsilon=1.e-6 )
    {
    if( length(x) <= 1 )   return(TRUE)
    
    xr  = regularizeWavelength(x)
    
    return( max( abs(x - xr) ) <= epsilon )
    }
    
isIncreasingSequence <- function( x )
    {
    #   in this one, we enforce strictly increasing, 
    #   except at the endpoints where we allow 1 duplicate value
    n   = length(x)
    
    if( n <= 1 )    return(TRUE)
    
    sdiff   = diff(x)   #   so sdiff has length n-1
    
    if( n == 2 )    return( 0 < sdiff )
    
    if( any( sdiff < 0 ) )  return(FALSE)
    
    #   all sdiff are now 0 or +
    if( n == 3 )    
        #   at most one 0 allowed
        return( sum(sdiff==0) <= 1 )
    
    #   the middle part must be all +
    sdiff   = sdiff[ 2:(n-2) ]
    
    return( all( 0 < sdiff ) )
    }
    
isStrictlyIncreasingSequence <- function( x )
    {
    n   = length(x)
    
    if( n <= 1 )    return(TRUE)
    
    return( all( 0 < diff(x) ) )
    }
        
#--------       UseMethod() calls           --------------#    
    
is.regular <- function(x) 
    {
    UseMethod("is.regular")
    }    
    
step.wl <- function(x) 
    {
    UseMethod("step.wl")
    }       

numSpectra <- function(x) 
    {
    UseMethod("numSpectra")
    }   
    
numWavelengths <- function(x) 
    {
    UseMethod("numWavelengths")
    }   

    
coredata <- function( x, forcemat=FALSE )
    {
    UseMethod("coredata")
    }
    
extradata <- function( x )
    {
    UseMethod("extradata")
    }    
"extradata<-"  <- function( x, add=FALSE, value )        
    {
    UseMethod("extradata<-")    
    }  
    
if( FALSE )
{    
isValid <- function( x )
    {
    UseMethod("isValid")
    }    
}
   

wavelength <- function(x)
    {
    UseMethod("wavelength")
    }
    
"wavelength<-"  <- function( x, value )        
    {
    UseMethod("wavelength<-")    
    }    
    
type <- function(x)
    {
    UseMethod("type")
    }
    
quantity <- function(x)
    {
    UseMethod("quantity")
    }
    
"quantity<-" <- function( x, value )
    {
    UseMethod("quantity<-")
    }
    
organization <- function(x)
    {
    UseMethod("organization")
    }
    
"organization<-" <- function(x, value)
    {
    UseMethod("organization<-")
    }
  
specnames <- function(x)
    {
    UseMethod("specnames")
    }
    
"specnames<-" <- function(x, value)
    {
    UseMethod("specnames<-")
    }  
    
metadata <- function(x,...)
    {
    UseMethod("metadata")
    }
    
"metadata<-" <- function(x,add=FALSE,value)
    {
    UseMethod("metadata<-")
    }
    
as.colorSpec <- function(...)
    {
    UseMethod("as.colorSpec")
    }
    
