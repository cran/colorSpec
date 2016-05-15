

    
#   returns a list of data.frames
#           all the data frames found in .path
importCGATS <-  function( .path )
    {  
    line_vec = readLines( .path )
    
    out = parseCGATS( line_vec ) 
    
    if( is.null(out)  ) return(NULL)
    
    attr( out, 'path' ) = .path
    
    return( out )
    }

#   returns a list of data.frames
    
parseCGATS <-  function( .line_vec )
    {
    pattern =   "^NUMBER_OF_FIELDS[ \t,]([0-9]+).*$"
    
    idx_number = which( grepl( pattern, .line_vec ) )
    if( length(idx_number) == 0 )
        {
        log.string( ERROR, "Cannot find any NUMBER_OF_FIELDS." )
        return(NULL)
        }

    n = length(idx_number)     

    number_of_fields    = sub( pattern, "\\1", .line_vec[idx_number] )
    number_of_fields    = as.integer(number_of_fields)
    
    
    #   find column names
    idx_format = which( grepl( "^BEGIN_DATA_FORMAT[ \t,]*$", .line_vec ) )
    if( length(idx_format) != n  )
        {
        log.string( ERROR, "Found %d BEGIN_DATA_FORMAT lines, but %d NUMBER_OF_FIELDS lines.", 
                                    length(idx_begin), n  )
        return(NULL)
        }
        
    idx_begin = which( grepl( "^BEGIN_DATA[ \t,]*$", .line_vec ) )
    if( length(idx_begin) != n )
        {
        log.string( ERROR, "Found %d BEGIN_DATA lines, but %d BEGIN_DATA_FORMAT lines.", 
                                    length(idx_begin), n )
        log.object( ERROR, idx_begin )
        return(NULL)
        }
        
        
    idx_end = which( grepl( "^END_DATA[ \t,]*$", .line_vec ) )
    if( length(idx_end) != n )
        {
        log.string( ERROR, "Found %d END_DATA lines, but %d BEGIN_DATA lines.", 
                                    length(idx_end), n )    
        log.object( idx_end )        
        return(NULL)
        }
        
    out = vector( n, mode='list' )
    
    for( k in 1:n )
        {
        i.min   = ifelse( k==1, 1,  idx_end[k-1] + 1 )
        i.max   = idx_number[k]-1    
        header  = .line_vec[ i.min:i.max ]         
        
        #  split the line after BEGIN_DATA_FORMAT
		line = .line_vec[ idx_format[k]+1 ]
		
		if( grepl( "\t", line) )
			sep = '\t'
		else if( grepl( ",", line) )
			sep = ','        
		else
			sep = ' '
		
        cname = strsplit( .line_vec[ idx_format[k]+1 ], sep )[[1]]      #; print( cname )
        
        if( length(cname) != number_of_fields[k] )
            {
            log.string( WARN, "NUMBER_OF_FIELDS mismatch %d != %d.",
                            length(cname), number_of_fields[k] )
            }
        
        data = strsplit( .line_vec[ (idx_begin[k]+1):(idx_end[k]-1) ], sep )
        
        #df = as.data.frame( data, stringsAsFactors=F )
        #print( str(df) )
        
        data_mat = matrix( "", length(data), length(cname) ) #	; print( str(data_mat) )
        
        for( i in 1:nrow(data_mat) )
            {
			#	print( str(data[[i]]) )
			if( length( data[[i]] ) != number_of_fields[k] )
                {
                log.object( ERROR, .line_vec[ idx_begin[k]+i ] )
                log.string( ERROR, "NUMBER_OF_FIELDS mismatch %d != %d.",
                                length( data[[i]] ), number_of_fields[k] )
                return(NULL)
                }
            data_mat[ i, ] = data[[i]]
            }
        #   print( data_mat )
        
        df = data.frame( type.convert(data_mat[ ,1],as.is=TRUE), stringsAsFactors=FALSE )
        
        if( 2 <= ncol(data_mat) )
            {
            for( j in 2:ncol(data_mat) )
                {
                df = cbind( df,   type.convert(data_mat[ ,j],as.is=TRUE) , stringsAsFactors=FALSE   )
                #   print( str(df) )
                }
            }
            
        #   simplify cname
        #   cname = gsub( "(XYZ_)|(LAB_)", "", cname )

        colnames( df ) = cname    
        
        out[[k]] = df
        
        attr( out[[k]], "header" ) = header        
        
        for( w in c("date","originator","serial","white.point") )
            attr( out[[k]], w ) = extractFieldFromHeader( header, w )
        }
        
    #   attr( out, "header" ) = header
        
    return( out )
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
    out     = NA

    if( .name == 'date' )
        {
        #   look for CREATED
        pattern = "^CREATED[ \t](.+).*$"
        idx = which( grepl( pattern, .header, ignore.case=T ) )
        
        if( length(idx) != 1 )  return(NULL)
        
        out = sub( pattern, "\\1", .header[idx] )   # ;       print(val)
        
        out = gsub( '\"', '', out )
        }
    else if( .name == "originator" )
        {
        pattern = "^ORIGINATOR[ \t,]([^,]+).*$"
        idx = which( grepl( pattern, .header ) )
        
        if( length(idx) != 1 )  return(NULL)
        
        out = sub( pattern, "\\1", .header[idx] )   # ;       print(val)
        
        out = gsub( '\"', '', out )
        }
        
    else if( .name == "serial" )
        {
        pattern = "^SERIAL[ \t,]([^,]+).*$"
        idx = which( grepl( pattern, .header ) )
        
        if( length(idx) != 1 )  return(NULL)
        
        out = sub( pattern, "\\1", .header[idx] )   # ;       print(val)
        }
        
    else if( .name == "white.point" )
        {
        pattern = "^WHITE_POINT[A-Za-z \t_,]+\"([^,]+)\".*$"
        idx = which( grepl( pattern, .header, ignore.case=T ) )     #    print( idx )
        
        if( length(idx) != 1 )  return(NULL)
        
        out = sub( pattern, "\\1", .header[idx] )   # ;       print(out)
        
        out = as.numeric( strsplit(out, ' ')[[1]] )
        }        

    if( any(is.na(out))  ) out   = NULL
    
    return( out )
    }
                
                
    
#   .pathvec    to .cal files created by ArgyllCMS 
#               with RGB_I, RGB_R, etc.   

plotVideoLUTs  <-  function( .pathvec )
    {
    data_list   = list()
    
    for( i in 1:length(.pathvec) )
        {
        data = importCGATS( .pathvec[i] )
        if( is.null(data) ) return(FALSE)
        
         if( is.null(data[[1]]$RGB_I) )
            {
            mess = sprintf( "%s(). RGB_I not found in '%s'\n", .pathvec[i] )
            cat(mess)
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