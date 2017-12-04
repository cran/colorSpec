

#   path            path to a CGATS file
#   collapsesingle  if there is only one table, then return the single data.frame (instead of a list with 1 data.frame)
#    
#   returns a list of data.frames
#           the list is assigned attributes: "path", "preamble", and (if present) "date", "created", "originator", "file_descriptor", "url"
#           each data.frame is assigned attributes: "header", and (if present) "descriptor"
#
#   the column names are taken from the BEGIN_DATA_FORMAT line
#   no attempt is made to assign row.names to the data.frames

readCGATS <-  function( path, collapsesingle=FALSE )
    {  
    line_vec = readLines( path )
    
    out = parseCGATS( line_vec ) 
    
    if( is.null(out)  ) return(NULL)
    
    attr( out, 'path' ) = path
    
    log.string( INFO, "%d tables read from '%s'.", length(out), path )    
    
    if( length(out) == 1  &&  collapsesingle )
        {
        #   copy attributes from the list down to the 1st data.frame
        for( w in names( attributes(out) ) )
            {
            if( w %in% 'names' )    next
            
            attr( out[[1]], w ) = attr( out, w )
            }
            
        #   just return the 1st data.frame
        out = out[[1]]
        }
        
    return( out )
    }

    
#   .line_vec   a character vector of lines, as read from a CGATS file
#    
#   returns a list of data.frames
#       the list is assigned attributes: "preamble", and (if present) "date", "created", "originator", "file_descriptor"
#       each data.frame is assigned attributes: "header", and (if present) "descriptor"
#   or NULL in case of error
    
parseCGATS <-  function( .line_vec )
    {
    idx_end = which( grepl( "^END_DATA[ \t,]*$", .line_vec ) )
    
    tables  = length(idx_end)
    if( tables == 0 )
        {
        log.string( ERROR, "Cannot find any tables (^END_DATA lines)." )
        return(NULL)
        }
        
    idx_begin = c( 1, idx_end[1:(tables-1)] + 1 )
    
    out         = vector( tables, mode='list' )
    names(out)  = character(tables)
    
    for( k in 1:tables )
        {
        out[[k]] = extractTableCGATS( .line_vec[ idx_begin[k]:idx_end[k] ], k )
        
        if( is.null(out[[k]]) ) return(NULL)
        
        names(out)[k]            = attr( out[[k]], 'name' )
        attr( out[[k]], 'name' ) = NULL
        }
        
    #   move some some lines from header #1 to preamble
    header1     = attr( out[[1]], "header" )
    pattern     = "^[ \t]*ORIGINATOR|^[ \t]*CREATED|^[ \t]*FILE_DESCRIPTOR|^[ \t]*URL"
    preamble    = c( 1, which( grepl( pattern, header1 ) ) )
    preamble    = 1:max(preamble)

    #   preamble    = unique( c(1,preamble) )   #   add line #1
    
    attr( out[[1]], "header" )  = header1[ -preamble ]

    preamble    = .line_vec[ preamble ]
    attr( out, "preamble" ) = preamble
        
    for( w in c("date","created","originator","file_descriptor","url") )
        attr( out, w ) = extractFieldFromHeader( preamble, w )
        
    #   trim blank lines from the beginning of all headers
    for( k in 1:tables )
        {
        header  =     attr( out[[k]], "header" )
        idx = which( grepl( "[^ \t]+", header ) )
        if( length(idx) == 0 )
            header = character(0)   # ALL lines are blank
        else
            header = header[ idx[1]:length(header) ]
            
        attr( out[[k]], "header" )  = header
        }

    return( out )
    }
    
#   .line_vec   character vector of lines read from CGATS file, representing exactly 1 table
#               the last line is always "END_DATA"    
#
#   returns     a data.frame, with attributes: "header", and (if present) "descriptor"

extractTableCGATS <- function( .line_vec, table_idx )
    {
    #   find start of the table
    pattern =   "^NUMBER_OF_FIELDS[ \t]+([0-9]+).*$"
    
    idx_number = which( grepl( pattern, .line_vec ) )   #; print( idx_number )
    if( length(idx_number) != 1 )
        {
        log.string( ERROR, "In table %d, expected exactly 1 NUMBER_OF_FIELDS lines, but found %d.", 
                            table_idx, length(idx_number) )
        return(NULL)
        }
        
    number_of_fields    = sub( pattern, "\\1", .line_vec[idx_number] )
    number_of_fields    = as.integer(number_of_fields)
                  
       
    header  = as.character(0)       
       
    if( 1 < idx_number )
        header = .line_vec[ 1:(idx_number-1) ]

    table_name  = extractFieldFromHeader( header, 'TABLE_NAME' )    
    if( is.null(table_name) )
        table_name  = sprintf( "TABLE_%d", table_idx )
    
    #   find data format block
    idx_format_begin = which( grepl( "^BEGIN_DATA_FORMAT[ \t,]*$", .line_vec ) )
    if( length(idx_format_begin) != 1  )
        {
        log.string( ERROR, "In table %d, expected exactly 1 BEGIN_DATA_FORMAT lines, but found %d.", 
                                    table_idx, length(idx_format_begin) )
        return(NULL)
        }
        
    idx_format_end = which( grepl( "^END_DATA_FORMAT[ \t,]*$", .line_vec ) )
    if( length(idx_format_begin) != 1  )
        {
        log.string( ERROR, "In table %d, expected exactly 1 END_DATA_FORMAT lines, but found %d.", 
                                    table_idx, length(idx_format_end) )
        return(NULL)
        }
        
    if( idx_format_end <= idx_format_begin+1 )
        {
        log.string( ERROR, "In table %d, BEGIN_DATA_FORMAT and END_DATA_FORMAT lines do not enclose a valid block of lines.", 
                                    table_idx )
        return(NULL)
        }
   
    pattern = "^NUMBER_OF_SETS[ \t]+([0-9]+).*$"
    idx_number = which( grepl( pattern, .line_vec ) )   #; print( idx_number )
    if( length(idx_number) != 1 )
        {
        log.string( ERROR, "In table %d, expected exactly 1 NUMBER_OF_SETS lines, but found %d.", 
                            table_idx, length(idx_number) )
        return(NULL)
        }
    number_of_records   = sub( pattern, "\\1", .line_vec[idx_number] )
    number_of_records   = as.integer(number_of_records)

        
    #   find beginning of the table data
    idx_begin = which( grepl( "^BEGIN_DATA[ \t,]*$", .line_vec ) )
    if( length(idx_begin) != 1 )
        {
        log.string( ERROR, "In table %d, expected exactly 1 BEGIN_DATA lines, but found %d.", 
                                    table_idx, length(idx_begin) )        
        log.object( ERROR, idx_begin )
        return(NULL)
        }
        
    idx_end = which( grepl( "^END_DATA[ \t,]*$", .line_vec ) )
    if( length(idx_end) != 1 )
        {
        log.string( FATAL, "In table %d, expected exactly 1 END_DATA lines, but found %d.", 
                                    table_idx, length(idx_begin) )        
        log.object( FATAL, idx_begin )
        return(NULL)
        }
    if( idx_end != length(.line_vec) )
        {
        log.string( FATAL, "In table %d, expected END_DATA line at line %d, but found it at line %d.", 
                                    table_idx, length(idx_begin), idx_end )        
        return(NULL)
        }
    if( idx_end <= idx_begin+1 )
        {
        log.string( ERROR, "In table %d, BEGIN_DATA and END_DATA lines do not enclose a valid block of lines.", 
                                    table_idx )
        return(NULL)
        }        
        
    
    #  analyze the data format block
    block = .line_vec[ (idx_format_begin+1):(idx_format_end-1) ]

    #   first split using a _single_ tab
    cname = unlist( strsplit( block, '\t' ) )    # ; print( cname )
    
    if( length(cname) != number_of_fields )
        {
        log.string( INFO, "In table %d, using CGATS-standard white-space specification.", table_idx )

        #   try to split fields again, but this time using CGATS-standard white-space
        cname = strsplit( paste(block,collapse=''), '[ \t\n\r]+' )[[1]]     #; print( cname )
        if( length(cname) != number_of_fields )
            {
            log.string( ERROR, "In table %d, cannot split format lines so that NUMBER_OF_FIELDS == %d.",
                                table_idx, number_of_fields )
            return(NULL)
            }
            
        #   split data using standard CGATS convention
        
        #   extract quoted strings using scan()
        what = rep( list(character(0)), number_of_fields )
        #   names(what) = cname
        
        out = scan( text=.line_vec[ (idx_begin+1):(idx_end-1) ], what=what, flush=TRUE, quiet=TRUE )
        }
    else
        {
        log.string( INFO, "In table %d, using non-CGATS-standard single tabs for white-space.", table_idx )
        
        if( 0 )
        {
        #   check that field names do not contain a space and are non-empty
        mask1   = grepl( " ", cname )
        mask2   = (cname == '')
        idx = which( mask1 | mask2  )
        if( 0 < length(idx) )
            {
            log.object( ERROR, cname[idx] )
            log.string( ERROR, "In table %d, %d field names have embedded spaces, or are empty.",
                        table_idx, length(idx) )
            return(NULL)
            }
        }

        #   put the data into a matrix of strings
        data = strsplit( .line_vec[ (idx_begin+1):(idx_end-1) ], '\t' )
               
        #   quick check against NUMBER_OF_FIELDS
        for( i in 1:length(data) )
            {
            #	print( str(data[[i]]) )
            if( length( data[[i]] ) != number_of_fields )
                {
                log.object( ERROR, .line_vec[ idx_begin+i ] )
                log.string( ERROR, "In table %d and row %d, NUMBER_OF_FIELDS mismatch %d != %d.",
                                    table_idx, i, length( data[[i]] ), number_of_fields )
                return(NULL)
                }
            #   data_mat[ i, ] = data[[i]]
            }

        data_mat    = sapply( data, function(x) {x} )
        
        #   should be a matrix
        if( ! is.matrix(data_mat) )
            {
            log.string( FATAL, "Internal Error. Object returned from sapply should be a matrix." )
            return(NULL)
            }

        #   make data.frame with all strings.  Note that we transpose data_mat first
        out = as.data.frame( t(data_mat), stringsAsFactors=FALSE )
        }
        
    #   check that field names are valid
    pattern = "^[-_.:/A-Za-z0-9]+$"
    idx = which( ! grepl(pattern,cname) )
    if( 0 < length(idx) )
        {
        log.object( ERROR, cname[idx] )
        log.string( ERROR, "In table %d, %d field names are invalid.",
                    table_idx, length(idx) )
        return(NULL)
        }
        
        
    #   out is now a list of character vectors
    #   examine each character vector, and convert to vector of numbers if possible
    out = lapply( out, type.convert, as.is=TRUE )
    
    #   out is now a list, convert to a data.frame.  assigning class(out) = 'data.frame' does not work ?? !!
    out = as.data.frame( out, stringsAsFactors=FALSE )   #col.names=cname, 
    
    if( nrow(out) != number_of_records )
        {
        log.string( WARN, "In table %d, expected %d records, but found %d.", 
                                    table_idx, number_of_records, nrow(out) )        
        }
                                        
    colnames( out ) = cname    # duplicates in cname[] will stay that way !
    
    attr( out, "name" )     = table_name     
    attr( out, "header" )   = header        
    
    for( w in c("descriptor") )
        attr( out, w ) = extractFieldFromHeader( header, w )

    return( out )
    }
    
    

#   .data       a list of data.frames or colorSpec objects, or a single data.frame or colorSpec object.  
#               If .data is a list, the preamble can be stored as attr(df,"preamble").
#               For each object df in the list, the header data is stored as attr(df,"header")
#   .path       name of file to write to. Can also be 'clipboard'
#   .sep        either space ' ' or tab '\t'.  If space then character columns of df are quoted.
writeCGATS  <- function( .data, .path, .sep=' ' )
    {
    .sep    = .sep[1]
    
    if( .sep != ' '  &&  .sep != '\t' )
        {
        log.string( ERROR, "separator=''%s', but it must be a space or tab.\n", .sep )
        }
    
    if( is.data.frame(.data)  ||  is.colorSpec(.data) )
        {
        #   make a list with only 1 data.frame
        header  = attr(.data,'header')
        #   print( header )
        .data = list(.data)
        attr(.data[[1]],'header')   = header
        #   print( str(.data) )
        }    
        
    if( ! is.list(.data) )    
        {
        log.string( ERROR, ".data is invalid\n"  )
        return(FALSE)
        }    
        
    con = file( .path, open="wt" )
                   
    if( grepl( "[.]sp$", .path, ignore.case=T ) )
        {
        #   the user is saving to an Argyll .sp file
        writeLines( c( "SPECT", '' ), con )
        }
        
        
    preamble  = attr(.data,"preamble")
    if( ! is.null(preamble) )
        writeLines( preamble, con )
        
    for( df in .data )
        {
        header  = attr(df,"header")
    
        if( ! is.null(header) )
            writeLines( header, con )
            

        #   df  = makeCompatibleCGATS( df )
        
        if( is.null(df) )   next
        
        if( is.colorSpec(df) )
            {
            extra   = extradata( df )
            extradata( df ) = cbind( extra, data.frame( SAMPLE_NAME=specnames(df), stringsAsFactors=F ) )
            cols    = length( colnames(df) )               
            field   = colnames(df)[1:(cols-1)]                 
            field   = c( field, sprintf( "nm%g", wavelength(df) ) )  
            }
        else
            field   = colnames(df)  
            
        #   write description of the fields
        line_vec = character(0)
        #   line_vec = c( line_vec, '' )
        line_vec = c( line_vec, paste( "NUMBER_OF_FIELDS", length(field), sep=.sep ) )
        line_vec = c( line_vec, "BEGIN_DATA_FORMAT" )
        line_vec = c( line_vec, paste( field, collapse=.sep ) )
        line_vec = c( line_vec, "END_DATA_FORMAT" )
        writeLines( line_vec, con )
    
        #   write the data
        line_vec = character(0)   
        #   line_vec = c( line_vec, '' )        
        line_vec = c( line_vec, paste( "NUMBER_OF_SETS", nrow(df), sep=.sep ) )
        line_vec = c( line_vec, "BEGIN_DATA" )
        writeLines( line_vec, con )
        
        #   print( str(format(df,digits=9)) )
        if( .sep == ' ' )
            quote   = which( sapply( df, is.character ) )
        else
            quote   = FALSE
        
        #   convert all fields to character
        df  = format( df, digits=9, scientific=F )

        write.table( df, file=con, append=T, quote=quote, sep=.sep,
                        eol="\n", na="NA", dec=".", row.names=F, col.names=F )
            
        writeLines( "END_DATA", con )
        }
    
    close( con )

    return(TRUE)
    }


    
#   .data       a data.frame with color data, but possibly with columns that are model matrices
#
#   returns     a data.frame with 
#               'filter' or 'patch' or 'sample' converted to 'SAMPLE_ID'
#               RGB expanded to RGB_R, ...
#               XYZ expanded to XYZ_X, ...
#               Lab expanded to LAB_L, ...
    
makeCompatibleCGATS  <-  function( .data )
    {
    for( s in c( 'filter', 'patch', 'sample' ) )
        {
        idx = which( colnames(.data) == s )
        
        if( length(idx) == 1 )
            {
            colnames(.data)[idx]    = "SAMPLE_ID"
            break
            }
        }
    
    out = NULL

    colupper    = toupper( colnames(.data) )
    
    for( j in 1:ncol(.data) ) 
        {
        cname   = colupper[j]
        
        if( cname %in% c("RGB","XYZ","LAB","LCH") )
            {
            #   expand 1 column into 3
            mat = .data[ , j ]    #   a model.matrix
            
            if( ncol(mat) != nchar(cname) )    next    #   something is wrong
            
            if( cname  %in%  c("RGB","XYZ") )
                mat = 100 * mat     #   CGATS likes percentages for these
            
            for( k in 1:ncol(mat) )
                {
                if( is.null(out) )
                    out = data.frame( V=mat[ ,k]  )
                else
                    out = cbind( out, data.frame( V=mat[ ,k] ) )
                
                #   change the name of the last column               
                colnames( out )[ ncol(out) ] = sprintf( "%s_%s", cname, substr(cname,k,k) )
                }
            }
        else
            {
            #   copy with no change
            if( is.null(out) )
                out = .data[j]
            else
                out = cbind( out, .data[j], stringsAsFactors=F )
            }
        }
        
        
    header  = attr( .data, "header" )
        
    if( ! "CREATED" %in% header )
    header  = c( header, sprintf( 'CREATED "%s"', strftime( Sys.time(), format="%Y-%m-%dT%H:%M:%S" ) ) )
                
    attr( out, "header" ) = header
        
    return( out )
    }
    
    
    
#   .header     character vector
#   .name       name of field
#   returns NULL if not found    
extractFieldFromHeader  <-  function( .header, .name )
    {
    #   out     = NULL

    replacement = "\\1"
    
    if( .name == 'date' )   
        .name   = "CREATED"
    else if( .name == "descriptor" )
        {
        .name = "(DESCRIPTOR|TABLE_DESCRIPTOR)"
        replacement = "\\2"
        }
 
    keyword = toupper(.name)
    
    pattern = sprintf( '^[ \t]*%s[ \t]+"([^"]+)".*$', keyword )    
    
    idx = which( grepl( pattern, .header ) )
    
    if( length(idx) == 0 )  return(NULL)
    
    if( 1 < length(idx) )
        {
        log.string( WARN, "Found %d matches for '%s'; using the last one.", length(idx) )
        idx = idx[ length(idx) ]
        }

    out = sub( pattern, replacement, .header[idx] )  

    out = type.convert( out, as.is=TRUE )

    #   if( any(is.na(out))  ) out   = NULL
    
    return( out )
    }
                
                
    
#   .pathvec    to .cal files created by ArgyllCMS 
#               with RGB_I, RGB_R, etc.   

plotVideoLUTs  <-  function( .pathvec )
    {
    data_list   = list()
    
    for( i in 1:length(.pathvec) )
        {
        data = readCGATS( .pathvec[i] )
        if( is.null(data) ) return(FALSE)
        
         if( is.null(data[[1]]$RGB_I) )
            {
            log.string( ERROR, "%s(). RGB_I not found in '%s'\n", .pathvec[i] )
            return(FALSE)
            }
               
        
        #   print( str(data) )
        data_list[[i]]  = data[[1]]
        
        attr( data_list[[i]], 'path' )  = attr( data, 'path' )
        }
    

    par( omi=rep(0.25,4) )
    par( mai=c(0.8, 0.8, 0.5, 0.25 ) )

    plot( c(0,255), c(0,255), type='n', asp=1, xlab='RGB in', ylab='RGB out', las=1, lab=c(10,15,7) )
    grid( lty=1 )
    abline( h=0, v=0 )
    
    lines( c(0,255), c(0,255), col='gray', lty=2 )
    lines( c(255,255,0), c(0,255,255), col='gray', lty=2 )
    
    colvec = c('red','green','blue')
    
    legend  = character(0)
    
    for( i in 1:length(data_list) )
        {
        data    = data_list[[i]]
    
        x   = 255 * data$RGB_I
            
        for( k in 1:3 )
            {
            lines( x, 255 * data[[ k+1 ]], col=colvec[k], lty=i )
            }

        legend[i]   = basename(attr( data, 'path' ) )
        }
        
    if( length(data_list) == 1 )
        {
        path    = attr( data_list[[1]], 'path' ) 
        mess    = sprintf( "%s \n[%s]", basename(path), dirname(path) )        
        
        title( main=mess, cex.main=1 )
        }
    else
        {
        xleft   = par("usr")[1]  
        ytop    = par("usr")[4]
        
        legend( xleft, ytop, legend, lwd=2, bty='n', lty=1:length(data_list), cex=1.1, xpd=NA, xjust=0.07, yjust=0.12 )        
        }
        
        
    return( invisible(T) )
    }