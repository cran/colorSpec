
#   pathvec         a vector of paths to files
#   ...             optional arguments to pass to resample.colorSpec
#
#   value:
#   a single colorSpec object, or NULL in case of error
#

readSpectra  <- function( pathvec, ... ) 
    {
    if( is.null(pathvec) || length(pathvec)==0  )    return(NULL)
    
    if( 1 < length(pathvec) )
        {
        #   call coroutine
        return( readAllSpectra( pathvec, ... ) )        
        }

    #   only a single path now
    path    = pathvec
        
    ftype = spectralFileType( path )
    
    if( is.na(ftype) )   
        {
        log.string( ERROR, "Cannot determine type of file '%s'.", path )
        return(NULL)
        }

    if( ftype == "Control"  ||   ftype == "CGATS" )
        {
        #   These 2 types are special because they return a *list* of objects.
        if( ftype == "Control" )
            out = readSpectraControl( path )
        else
            out = readSpectraCGATS( path )
        
        if( is.null(out) )  return(NULL)
        
        if( 0 < length( list(...) ) ) 
            {
            #   resample all of them
            for( k in 1:length(out) )
                out[[k]]    = resample( out[[k]], ... )
                
            log.string( INFO, "Resampled %d spectra from file '%s' to new wavelengths.", length(out), path )
            }    
        
        if( length(out) == 1 )
            {
            #   only 1 object, so just return i
            return( out[[1]] )
            }

        if( 0 )
        {
        if( ! areSpectraBindable(out) )
            {
            log.string( ERROR, "File '%s' has %d distinct spectra that cannot be combined.  Try assigning a new wavelength sequence for resampling.",
                                path, length(out) )
            return(NULL)
            }
        }

        out.bound   = bindSpectra( out )
                
        return( out.bound )
        }
        
        

    if( ftype == "scope" )
        {
        #   the simplest type
        out = readSpectrumScope( path )
        }
    else if( ftype == "XYY" )
        {
        #   the simplest type
        out = readSpectraXYY( path )
        }
    else if( ftype == "spreadsheet" )
        {
        #   this means tab-delimited text spreadsheet, e.g. the IT8 target files from Wolf-Faust
        out = readSpectraSpreadsheet( path )
        }        
    else if( ftype == "Excel" )
        {
        #   this means a real Excel spreadsheet
        out = NULL  # readSpectraExcel( path )
        }    
    else
        {
        return(NULL)
        }
        
    if( is.null(out) )  return(NULL)        
        
    # strip the extension off of .path, for nicer display
    #   attr( out, "path" ) = sub( "[.][a-zA-Z]+$", "", .path )
    #if( is.null(attr( out, "path" )) )  attr( out, "path" ) = .path
    #if( is.null(attr( out, "date" )) )  attr( out, "date" ) = file.info(.path)$mtime
    
    if( 0 < length( list(...) ) )
        {
        out = resample( out, ... )
        }
        
    return( out )
    }
    
    
    
    
#   pathvec         vector of  paths, always 2 or more
#   ...             optional arguments to pass to resample.colorSpec
#
#   value:
#   a colorSpec with organization 'df.row', or NULL in case of ERROR

readAllSpectra  <- function( pathvec, ... )
    {    
    #   read the 1st one
    out = readSpectra( pathvec[1], ... ) 
    
    if( is.null(out) )  return(NULL)
    
    organization(out)   = 'df.row'
    
    if( length(pathvec) == 1 ) return(out)

    #   specnames   = specnames(out)
    path        = rep( pathvec[1], numSpectra(out) )
    
    for( k in 2:length(pathvec) )
        {
        obj = readSpectra( pathvec[k], ... )
        
        if( is.null(obj) ) return(NULL)
        
        #   organization(obj)   = 'df.row'
            
        out = bind( out, obj )
        
        if( is.null(out) )  return(NULL)
        
        #   specnames   = c( specnames, specnames(obj) )
        path    = c( path, rep( pathvec[k], numSpectra(obj) ) )
        }
        
    extradata(out)   = cbind( data.frame( path=path, stringsAsFactors=F ), extradata(out) )
        
    return(out)
    }
        
    
#   readSpectraControl()
#
#   unlike most of these functions, this one returns a *list* of colorSpecs
#   Each colorSpec in the list may have a different wavelength sequence !    

readSpectraControl <- function( path)
    {
    line_vec = readLines( path)
    
    n = length(line_vec)
    
    if( n == 0 ) return( NULL )
    
    #   comments    = line_vec[ grepl( "^[ \t]*#", line_vec ) ]
    
    #   look for [Control]
    
    idx = which( grepl( "^\\[Control\\]", line_vec, ignore.case=TRUE ) )
    if( length(idx) != 1 )
        {
        log.string( ERROR, "Cannot find [Control] section in '%s'.", path)
        return( NULL )
        }
        
    header  = line_vec[ 1:(idx-1) ]     # every line above [Control]
        
    #   look for constant increment
    increment = 0
    pattern = "^#[ \t]+increment[ \t=]+([0-9.]+).*$"
    inc = which( grepl( pattern, line_vec ) )
    if( 0 < length(inc) )
        {
        increment   = sub( pattern, "\\1", line_vec[inc[1]] )
        increment   = as.numeric( increment )
        #   print( increment )
        }
        
        
    x_px_vec = numeric(0)
    y_px_vec = numeric(0)
    
    wavelength_vec = numeric(0)
    response_vec = numeric(0)
    
    for( k in (idx+1):n )
        {
        word = strsplit( line_vec[k], "[ \t]" )[[1]]
        #print(word)
        
        if( length(word) != 4 ) break
        
        vec = suppressWarnings( as.numeric(word) )
        
        #print( vec )

        if( all( ! is.na(vec) ) )
            {
            x_px_vec = c( x_px_vec, vec[1] )
            y_px_vec = c( y_px_vec, vec[2] )
            
            wavelength_vec  = c( wavelength_vec, vec[3] )
            response_vec    = c( response_vec, vec[4] )
            }
        }
        
    #print( x_px_vec )
    #print( y_px_vec )
    #print( wavelength_vec )
    #print( response_vec )
    
    lm.x = lm( wavelength_vec ~ x_px_vec + y_px_vec )
    lm.y = lm( response_vec   ~ y_px_vec + x_px_vec  )
    
    #   print( summary( lm.x ) )
    #   print( summary( lm.y ) )
    
    out = list()
    
    #   look for the data
    idx = grep( "^\\[[-._a-zA-Z0-9 ]+\\]", line_vec, value=F )
    
    #   print(idx)
    if( length(idx) <= 1 )
        {
        log.string( ERROR, "Cannot find any [data] sections, only a [Control] section in '%s'.", path)
        return( NULL )
        }
        
    #   rgb_vec = c( "Red", "Green", "Blue" )
    
    base = sub( "[.][a-zA-Z]+$", "", basename( path) )    #   take away the extension
    
    
    for( i in idx )
        {
        cname = gsub( "(\\[)|(\\])", "", line_vec[i] )
        
        #   print( cname )
        
        if( cname == "Control" ) next    # not a channel
            
        log.string( DEBUG, "Found [%s] channel on line %d.", cname, i )
        
        x_px_vec = numeric(0)
        y_px_vec = numeric(0)
    
        for( k in (i+1):n )
            {
            word = strsplit( line_vec[k], "[ \t]" )[[1]]
            #print(word)
            
            if( length(word) != 2 ) break
            
            vec = suppressWarnings( as.numeric(word) )
            
            #print( vec )

            if( all( ! is.na(vec) ) )
                {
                x_px_vec = c( x_px_vec, vec[1] )
                y_px_vec = c( y_px_vec, vec[2] )
                }
            }
        
        if( length(x_px_vec) < 2 )
            {
            log.string( WARN, "Found only %d points in [%s] channel", length(x_px_vec), cname )
            }
            
        data_new = data.frame( x_px_vec=x_px_vec, y_px_vec=y_px_vec )   # ; print( data_new )
        
        y = predict(lm.y,data_new)
        
        cextra  = paste( base, cname, sep='.', collapse=NULL )
        
        x_new   = predict(lm.x,data_new)
        
        if( 0 < increment )
            {
            #   round to nearest multiple of increment
            x_new   = increment * round( x_new / increment )
            }
            

        #   attr( y, "specname" ) = cname
        
        quant   = guessSpectrumQuantity( cname, header )            
        if( is.na(quant) )
            {
            quant = 'energy'
            log.string( WARN, "Cannot guess quantity from from contents of '%s', so assigning quantity='%s'.",
                                basename(path), quant )
            }      
            
        spec    = colorSpec( y, x_new, quantity=quant )   
        
        specnames( spec )   = cname
        metadata( spec )    = list( path=path, header=header )

        out[[ cname ]] = spec
        }
    
    if( length(out) == 0 )  return(NULL)

    return( out )
    }


readSpectraXYY <- function( path )
    {
    if( ! file.exists( path) )
        {
        log.string( ERROR, "File '%s' does not exist !", path)
        return(NULL)
        }
        
    #    peek to see whether tabs are used as separator    
    header  = readLines( path, 100 )
    
    if( length(header) == 0 ) return(NULL)
    
    pattern = "^(wave|wv?l)"
    idx = which( grepl( pattern, header , ignore.case=T ) )
    
    if( length(idx) == 0 ) 
        {
        log.string( ERROR, "Cannot find Wavelength column in '%s' !\n",  path)
        return(NULL)
        }

    #   comments = header[ grepl( "^[ \t]*#", header ) ]
    
    #   print( comments )
    
    line = header[ idx[1] ]
    
    if( grepl( "\t", line) )
        sep = '\t'
    else if( grepl( ",", line) )
        sep = ','        
    else
        sep = ''
        

    df = read.table( path, header=T, sep=sep, skip=idx[1]-1, stringsAsFactors=F )   #; print( str(df) )
    
    df[ is.na(df) ] = 0
    
    #   print( 2:ncol(df) )
    
    data    = as.matrix( df[  , 2:ncol(df) ] )  #; print( str(data) )
    
    colnames(data)  =  colnames(df)[ 2:ncol(df) ]
    
    m = ncol(data)

    if( sep == '' ) sep = ' '
    

    if( TRUE &&  m == 1 )
        {
        base = sub( "[.][a-zA-Z]+$", "", basename( path) )    #   take away the extension
        cname   = colnames(data)
        colnames(data) = paste( base, cname, sep='.', collapse=NULL )
        }
        
    header = header[ 1:(idx[1]-1) ]     #;   log.object( DEBUG, header )

    #   theType     = guessSpectrumType( colnames(data), header )   #; print( theType )
    theQuantity = guessSpectrumQuantity( colnames(data), header )
    if( is.na(theQuantity) )
        {
        theQuantity = 'energy'
        log.string( WARN, "Cannot guess quantity from from contents of '%s', so assigning quantity='%s'.",
                                basename(path), theQuantity )
        }
            
    out = colorSpec( data, df[[1]], quantity=theQuantity, organization='df.col' )
            
    metadata( out ) = list( path=path)            
            
    if( 0 < length(header) )
        metadata( out, add=TRUE ) = list(header=header)

    return(out)
    }
    
    
#   Spreadsheet here means one of:
#       the Wolf Faust files:  E131102.txt  etc.    

readSpectraSpreadsheet <- function( path )
    {
    if( ! file.exists( path) )
        {
        log.string( ERROR, "File '%s' does not exist !",  path)
        return(NULL)
        }
        
    #    peek to find first line with data
    line = readLines( path, 60 )
    
    if( length(line) == 0 ) return(NULL)
    
    idx_start = which( grepl( "^(ID|SAMPLE|Time)", line , ignore.case=F ) )
    
    if( length(idx_start) == 0 ) 
        {
        log.string( ERROR, "Cannot find header line in '%s' !\n",   path)
        return(NULL)
        }
        
        
    skip = idx_start[1] - 1
    header  = line[ 1:skip ]   #;   print( header )
    
    df  = read.table( path, skip=skip, sep='\t', header=T, quote='', stringsAsFactors=F )
    #   print( str(df) )
    
    pattern = "^[A-Z]+([0-9.]+)nm$"
    mask_nm = grepl( pattern, colnames(df) )
    count   = sum(mask_nm)
    
    log.string( INFO, "Found %d wavelengths in '%s'\n", count, path)
    if( count == 0 )
        {
        return(NULL)
        }
        
    #   print( colnames(df)[mask_nm] )
    
    #   convert wavelength from text to numeric
    wavelength  = as.numeric( sub( pattern, "\\1", colnames(df)[mask_nm] ) )
    #   print( wavelength )
    
    data    = t( as.matrix( df[  , mask_nm ] ) )
    #   print( str(val) )
    

    if( ! is.null( df$Name ) )
        {
        #   from Wolf Faust CD
        cnames  =   sub( "[ ]+$", '', df$Name )   # trim trailing space
        #   cnames  =   paste( cnames, ".Transmittance", sep='' )      
        }        
    else if( ! is.null( df$Sample ) )
        {
        if( is.null( df$Sequence ) )
            cnames = df$Sample
        else
            cnames  = paste( df$Sample, df$Sequence, sep='.' )
            
        #   cnames  =   paste( cnames, ".Power", sep='' )    # ;  print( cnames )        
        }            

    colnames(data)  = cnames

    #   type        = 'material'
    
    quantity    = guessSpectrumQuantity( cnames, header )
    if( is.na(quantity) )
        {
        quantity = 'reflectance'
        log.string( WARN, "Cannot guess quantity from from contents of '%s', so assigning quantity='%s'.",
                                basename(path), quantity )
        }
    
    out         = colorSpec( data, wavelength, quantity=quantity, organization="df.row" )
    
    metadata( out ) = list( path=path, header=header )    
    
    part_left   = df[  , ! mask_nm, drop=F  ]    
    
    if( 0 < ncol(part_left) )
        extradata(out)  = part_left
    
    for( w in c("date","originator","serial","white.point") )
        attr( out, "metadata" )[[w]]    = extractFieldFromHeader( header, w )    

    return( out )
    }
    
    

    
        
#   .path       to a CGATS file, e.g. .sp or .cal .cgt or .txt, or the Rosco files
#
#   returns a *list* of colorSpec objects, each with organization 'df.row'.   
#   WARNING:  does not return a single object.
    
readSpectraCGATS <-  function( path )
    {
    theData = readCGATS( path )
    
    if( is.null(theData) )  return(NULL)
    
    n = length(theData)     # n is the number of data.frames

    base = sub( "[.][a-zA-Z]+$", "", basename( path ) )
    
    preamble    = attr(theData,'preamble')
    
    #   attr_name = c( "Substrate", "RefractiveIndex", "Thickness" )
    
    #   iterate in reverse order, so we can drop the tables with non-spectral data
    for( i in n:1 )
        {
        obj = colorSpecFromDF( theData[[i]], path, preamble, i )                
        
        if( is.null(obj) )
            #   this is an error
            return(NULL)
            
        if( ! is.colorSpec(obj) )
            {
            #   table i does not have spectral data, so delete it
            theData = theData[-i]
            next
            }

        #   overwrite previous data.frame with new colorSpec object
        theData[[i]] = obj
        }
        
    if( length(theData) == 0 )
        {
        log.string( ERROR, "Found no tables with spectral data in '%s'.", path )
        return(NULL)
        }    
        
    attr( theData, 'path' ) = path
    
    return( theData )
    }
        
        
        
        
#   df          data.frame, from a table in a CGATS file
#   path        the file
#   preamble    from path.  might be NULL
#   idxtable    index of the table in path.  1-based
#
#   returns a colorSpec object
#   in case of ERROR, returns NULL
#   if there is no spectral data, returns as.logical(NA)

colorSpecFromDF <- function( df, path, preamble, idxtable )        
    {
    base = sub( "[.][a-zA-Z]+$", "", basename( path ) )
        
    colspec = spectralColumns( df, path, idxtable )  #; print( colspec )
        
    if( is.null(colspec) )
        #   this is an error
        return(NULL)
        
    if( ! is.list(colspec) )
        {
        #   table i does not have spectral data, so return any non-trivial thing
        return( as.logical(NA) )
        }
            
    idx_val     = colspec$idx_val
    idx_extra   = colspec$idx_extra
    wavelength  = colspec$wavelength
    
    #   idx     = which( mask )                            
    mat = as.matrix.data.frame( df[ ,idx_val, drop=F] ) / colspec$divisor
    
    mat = t( mat ) #;  print( dim(mat) ) 

    rownames( mat ) = as.character( wavelength )   #; print( str(mat) )      # colnames(df)[mask]
    
    #   print( mat )        
    
    #   specnames   = paste( base, as.character( 1:nrow(df)), sep='.' )     # just the defaults
    
    tableid = sprintf( "%s-%d", base, idxtable )
    
    if( nrow(df) == 1 )
        specnames   = tableid
    else
        {
        #specnames   = paste( tableid, as.character( 1:nrow(df)), sep='.' )   
        specnames   = sprintf( "%s.%d", tableid, 1:nrow(df) )
        }

    if( 0 < length(idx_extra) )
        {
        candidate   = c("SAMPLE_NAME","SAMPLE_ID","SampleID","Name")
        idx = which( candidate %in% colnames(df) )
        if( 1 <= length(idx) )
            specnames   = df[[ candidate[idx[1]] ]]
        }
        
    #   print(specnames)       
    
    colnames(mat)   = specnames

    header  = attr( df, "header" )
    
    #   look for SPECTRAL_NORM, and if present divide by it
    pattern = '^SPECTRAL_NORM[ \t]*"?([.0-9]+)"?'
    vec = sub( pattern, "\\1", header, ignore.case=T )
    
    idx = which( nchar(vec) < nchar(header) )
    if( 0 < length(idx) )
        {
        #   found 1 or lines matching, only look at the first one
        value   = as.numeric( vec[idx[1]] )
        
        if( is.numeric(value)  &&  is.finite(value)  &&  0 < value  &&  value != 1 )
            {
            log.string( INFO, "Dividing spectral values by SPECTRAL_NORM=%g.", value )
            mat = mat / value
            }
        }

    theQuantity = guessSpectrumQuantity( specnames, c(preamble,header) )
    
    if( is.na(theQuantity) )
        {
        theQuantity = 'energy'
        log.string( WARN, "Cannot guess quantity from from contents of '%s', so assigning quantity='%s'.",
                            basename(path), theQuantity )
        }
        
    out = colorSpec( mat, wavelength, quantity=theQuantity, organization="df.row" )
    
    #   ColorMunki wavelengths in hires mode are equally spaced at 3.333 nm, but rounded to the nearest integer
    #   check for this and fix it        
    wavereg = regularizeWavelength(wavelength)
    discrep = max( abs(wavereg - wavelength) )
    if( 0 < discrep  &&  discrep < 0.5 )
        {
        step.wl     = wavereg[2] - wavereg[1]
        log.string( WARN, "Perturbed wavelengths in '%s' to have equal increments of %g nm [%g to %g nm]", 
                            path, step.wl, wavereg[1], wavereg[ length(wavereg) ] )
        wavelength(out)  = wavereg                                
        }

    if( 0 < length(idx_extra) )
        {
        extra   = df[ , idx_extra, drop=F ]     # ;print( part_left )
    
        extradata(out)  = extra
        }
        
    metadata( out )    = list(  header=header, 
                                date=attr(df,"date"),
                                descriptor=attr(df,"descriptor") )    

    return( out )
    }
    
    
    
        
        
        
        
#   df          data.frame, from a table in a CGATS file
#   path        to the file
#   idxtable    index of the table in path
#
#   returns a list with items
#       idx_val     columns that correspond to a spectral values (but not to wavelength values)  
#       idx_extra   columns that correspond to extra data, maybe none
#       wavelength  vector of wavelengths, the same length as idx_val
#       divisor     for the values
#   in case of ERROR, returns NULL
#   if there is no spectral data, returns as.logical(NA)

spectralColumns <- function( df, path, idxtable )        
    {
    #cat( "spectralColumns()\n" )
    
    out = list()    
    
    cname   = colnames(df)      #; print(cname)
    
    #   search for standard convention for spectral data
    
    mask_dec    = ("SPECTRAL_DEC" == cname)   # grepl( "^SPECTRAL_DEC", cname )
    mask_pct    = ("SPECTRAL_PCT" == cname)   # grepl( "^SPECTRAL_PCT", cname )

    if( any( mask_dec | mask_pct ) )
        {
        #   the standard convention
        log.string( INFO, "In file '%s' and table %d, using standard convention for spectral data.", path, idxtable )            
        #cat( "Standard convention\n" )
                
        idx_dec = which( mask_dec )
        idx_pct = which( mask_pct )
                        
        #   pattern = "^SPECTRAL_NM"
        mask_wl = ("SPECTRAL_NM" == cname)    # grepl( pattern, cname )
        idx_wl  = which( mask_wl )
        n       = length(idx_wl)        # number of wavelengths

        if( n == 0 )
            {
            log.string( ERROR, "In file '%s' and table %d, there are %d spectral value fields, but there are no wavelengths.", 
                            path, idxtable, length(idx_dec)+length(idx_pct) )
            return(NULL)
            }

        divisor = 0
        if( length(idx_dec)==n  &&  all( idx_dec == idx_wl+1 ) )
            divisor = 1     # got a match
        else if( length(idx_pct)==n  &&  all( idx_pct == idx_wl+1 ) )
            divisor = 100   # got a match
        
        if( divisor == 0 )
            {
            log.string( ERROR, "In file '%s' and table %d, there are %d wavelength columns, but spectral value columns do not match up.", 
                            path, idxtable, n )      
            return(NULL)
            }
            
        wavelength  = lapply( df[ , idx_wl ], unique )  #; print(wavelength)
        waveunique  = sapply( wavelength, length )      #; print(wavelength)
        idx_max     = which.max( waveunique )
        if( 1 < waveunique[ idx_max ] )
            {
            #   not unique
            log.object( ERROR, wavelength[[idx_max]] )
            log.string( ERROR, "In file '%s', wavelength columns in table %d exist, but they are not unique.", 
                            path, idxtable )      
            return(NULL)
            }
            
        out$idx_val     = idx_wl + 1 
        out$idx_extra   = which( ! (mask_wl | mask_dec | mask_pct) )
        out$wavelength  = unlist( wavelength, use.names=FALSE )    #; print( out$wavelength )
        out$divisor     = divisor
        
        #   print( out )
        }
    else
        {
        #   try again with non-standard pattern
 
        pattern = "^(nm|SPEC_|SPECTRAL_)[A-Z_]*([0-9.]+)$"
        
        mask_val    = grepl( pattern, cname )  # ; print( idx )
        
        idx_val     = which( mask_val )         # ; print(idx_wl)
                
        n   = length(idx_val)   # number of wavelengths
        
        if( n == 0 )
            #   no spectral data found
            return( as.logical(NA) )
                            
        log.string( INFO, "In file '%s' and table %d, using non-standard convention for spectral data.", path, idxtable )
            
        contig = all( diff(idx_val) == 1 )
        
        if( ! contig )
            {
            log.string( ERROR, "In file '%s' and table %d, there are %d spectral columns, but they are not contiguous.", 
                            path, idxtable, n )
            #   log.object( WARN, idx_val )            
            return(NULL)
            }
            
        wavelength      = sub( pattern, "\\2", cname[idx_val] )     #; print( wavelength )
        
        out$wavelength  = type.convert( wavelength, as.is=TRUE )    #; print( wavelength )
            
        out$idx_val     = idx_val
        out$idx_extra   = which( ! mask_val )
        out$divisor = 1
        }
        
    #   verify that wavelength is valid
    ok  = is.numeric(out$wavelength)  &&  all( is.finite(out$wavelength) )
    if( ! ok )
        {
        log.string( ERROR, "In file '%s' and table %d, the %d wavelength values are not all numeric and finite.", 
                            path, idxtable, length(wavelength) )   
        return(NULL)
        }
        
    #   verify that spectral values are valid
    ok  = all( sapply( df[out$idx_val], function(y) { is.numeric(y)  &&  all( is.finite(y) ) } ) )
    if( ! ok )
        {
        log.string( ERROR, "In file '%s' and table %d, the %dx%d spectral values are not all numeric and finite.",
                            path, idxtable, nrow(df), length(idx_val) )
        return(NULL)
        }
     
    return(out)    
    }
    

#   read Ocean Optics format
#   we already know this is a light source of with quantity 'energy'
#
#   TODO:   extract integration time from header and organize as df.row
readSpectrumScope  <- function( path)
    {
    linevec = readLines( path)

    ok = (0 < length(linevec)) ;  #  assert( ok )

    if( ! ok )  { return(NULL) }

    header  = NULL
    idx = which( grepl( "^>+Begin", linevec[1:30] ) )
    if( length(idx) == 1 )
        header  = linevec[ 1:(idx-1) ]
    
    #   make a list of character vectors,  tab-delimited
    line_split = strsplit( linevec, '\t', fixed=T )

    #  form subset of those lines that have data - 2 tab-delimted words
    theData = base::subset( line_split, sapply(line_split,length) == 2 )
    
    #  unlist() is faster than sapply()
    theData = as.double( unlist(theData) ) ;   
    dim(theData) = c(2,length(theData)/2)

    y   = theData[ 2,  ]

    out = colorSpec( y, wavelength=theData[1, ], quantity="energy", organization="vector" )
    
    #   there is only 1 spectrum here, so use the path as the name
    specnames(out)  = stripExtension( basename(path) )

    metadata( out ) = list( path=path, header=header )
    
    return( out )
    }
    
   
    
    
#   plain text means:
#       the usual characters
#       TAB through CR    
#       © and µ
isPlainText <- function( .stringvec )
    {
    n   = length(.stringvec)
    
    out = logical(n)
    
    for( i in 1:n )
        {
        byte    = charToRaw( .stringvec[i] )       # utf8ToInt( .stringvec[i] )
        
        out[i]  = all(  (32 <= byte & byte <= 126) |  (9 <= byte & byte <= 13) | byte==169 | byte==181 )
        }

    return( out )
    }

is.UTF8 <- function( .stringvec )
    {
    n   = length(.stringvec)
    
    out = logical(n)
    
    for( i in 1:n )
        {
        ivec    = try( utf8ToInt( .stringvec[i] ), silent=TRUE )
        
        out[i]  = class(ivec) != "try-error"
        }

    return( out )
    }
    
    
#   does not look at file extension, except for .XLS(X)
spectralFileType <- function( .path )
    {
    out = as.character(NA)
    
    if( is.null(.path) )
        {
        log.string( ERROR, "File argument is NULL !" )
        return(out)
        }

    if( ! file.exists( .path) )
        {
        log.string( ERROR, "File '%s' does not exist !", .path)
        return(out)
        }
        

    #   read a few lines at the beginning of file
    line = readLines( .path, 32, warn=F )

    #   ignore comment lines   
    line = line[ ! grepl( "^#", line ) ]
    
    if( length(line) == 0 ) return(out)
    
    #   crude check for a text file
    if( ! all( is.UTF8(line) ) )
        {
        #   a binary file
        if( FALSE  &&  grepl( "[.]xlsx?$", basename(.path), ignore.case=T ) )
            return( "Excel" )
        else
            {
            log.string( WARN, "File type of '%s' unnknown. It appears to be binary. !", .path )
            return(out)
            }
        }
    
    pattern = c(  "^CGATS|^ISO28178|^NUMBER_OF_FIELDS|^KEYWORD" ,  "^\\[Control\\]", "^(wave|wv?l)", "^>+Begin", "^ID\tName", "^Time"  )
    
    type    = c( "CGATS", "Control", "XYY", "scope", "spreadsheet", "spreadsheet" )
    
    for( i in 1:length(pattern) )
        {
        #   print( type )
        #   print( pattern[type] )
        
        if( any( grepl(  pattern[i] , line , ignore.case=T ) ) )    return( type[i] )
        }
    
    return( out )
    }
    
stripExtension <- function( .path )
    {
    pattern = "[.][a-z0-9]+$"
    return( sub(pattern,'',.path,ignore.case=T) )
    }
    
    
    
##################      semi-deadwood below     ###############################

    
if( 0 )
{
    
#   Spreadsheet here means a real Excel spreadsheet, e.g.  ASTMG173.xls, sheet SMARTS2
#
#   path      full path to file
#   worksheet  name of the worksheet.  If NULL then search for the 1st one with data.
#
readSpectraExcel <- function( path, worksheet=NULL )
    {
    for( p in c('rJava','xlsxjars','xlsx') )
        {
        ok  = requireNamespace( p, quietly=TRUE )
        if( ! ok )
            {
            log.string( ERROR, "Package '%s' is required.  Please install it.", p )
            return(NULL)
            }
        }
        
        
    if( ! file.exists( path) )
        {
        log.string( ERROR, "File '%s' does not exist !\n", path)
        return(NULL)
        }    
        
    wb = xlsx::loadWorkbook( path)
    
    sheetlist   = xlsx::getSheets(wb)
    
    #   print( str(sheetlist) )
    
    idx_sheet   = 0
    
    if( is.null(worksheet) )
        {
        #   search for the first worksheet that has rows
        for( i in 1:length(sheetlist) )
            {
            #   print( str(getRows( sheetlist[[i]] ) ) )
            
            rowcount    = length( xlsx::getRows( sheetlist[[i]] ) )
            if( 0 < rowcount )
                {
                log.string( INFO, "Found sheet '%s'   rows=%d.", 
                                    names(sheetlist)[i], rowcount )
                idx_sheet = i
                break
                }
            }
            
        if( idx_sheet == 0 )
            {
            log.string( ERROR, "Cannot find a suitable worksheet, with rows, in '%s'.", path)
            return(NULL)
            }
        }
    else
        {
        #   search for 1st worksheet with matching name
        idx_sheet   = which( worksheet == names(sheetlist) )
        
        if( length(idx_sheet) != 1 )
            {
            log.string( ERROR, "Cannot find worksheet '%s' in '%s'.",
                                    worksheet, path )
            return(NULL)
            }
        }
        
    theSheet    = sheetlist[[ idx_sheet ]]
    
    
    #   look for Lambda in 1st column
    col1    = xlsx::readColumns(theSheet, 1, 1, 1, 16, header=F, as.data.frame=F )[[1]]
    #   print( str(col1) )
    
    iname   = sub( "^([A-Za-z0-9]+).*$", "\\1", col1[1] )
    
    pattern = "^lambda|nm$"
    
    idx = which( grepl( pattern, col1, ignore.case=T ) )
    #   print( idx )
    
    if( length(idx) != 1 )
        {
        log.string( ERROR, "Cannot find '%s' in worksheet %d in '%s'.",
                                pattern, idx_sheet, path)
        return(NULL)
        }
    
    row1    = xlsx::readRows( theSheet, 1, 1, 1  )
    #   print( row1 )
    
    spectra = length(row1)-1
    
    mat     = xlsx::readColumns( theSheet, 1, spectra+1, idx+1, header=F, colClasses='numeric' ) #; print( str(mat) )

    wave    = as.numeric( mat[ ,1] )
    mat     = as.matrix( mat[ ,2:ncol(mat)] )
    
    colnames(mat)   = sprintf( "%s%02d.Energy", iname, 1:spectra )
    
    out = colorSpec( mat, wave, quantity='energy' )

    metadata(out) = list( path=path)
    
    return( out )
    }
}        