
#   .df1, .df2      data.frames with no duplicate rownames

#   always returns a data.frame with the right number of rows
#   and preserves all data from the 2 input data.frames
#   nrow(out) = nrow(.df1) + nrow(.df2)    
#
#   The columns of .df1 and .df2 are divided into 3 groups by their names
#       *) columns in .df1 but not in .df2
#       *) columns in .df2 but not in .df1
#       *) columns in common (the intersection of the names)


rbind.super <- function( .df1, .df2, .check=TRUE )   
    {
    if( .check )
        {
        out = rbind.super( .df1, .df2, .check=FALSE )
        if( is.null(out) )  return(NULL)
        
        if( nrow(out) != nrow(.df1) + nrow(.df2) )
            {
            mess    = sprintf( "rbind.super() FATAL  nrow(out)=%d != %d = nrow(.df1) + nrow(.df2).",
                                nrow(out), nrow(.df1),nrow(.df2)  )
            cat( mess )
            return(NULL)
            }
        return(out)
        }
    
    #   2 trivial cases
    if( nrow(.df1)==0 ) return(.df2)
    if( nrow(.df2)==0 ) return(.df1)    
        
    rnames1 = row.names( .df1 )
    rnames2 = row.names( .df2 )
    
    if( is.null(rnames1) )  rnames1 = rep( NA_character_, nrow(.df1) )
    if( is.null(rnames2) )  rnames2 = rep( NA_character_, nrow(.df2) )        
    
    rownames.out    =   c( rnames1, rnames2 )
    
    if( any( duplicated( rownames.out, incomparables=NA_character_ ) ) )
        {
        mess    = sprintf( "rbind.super() row names have duplicates.\n" )
        cat( mess )
        return(NULL)
        }
    
    
    
    common  = intersect( colnames(.df1), colnames(.df2) ) ;  print( common )
    
    if( 0 < length(common) )
        {
        df1.common  = .df1[ , common, drop=FALSE ]  ; print( df1.common )
        df2.common  = .df2[ , common, drop=FALSE ]  ; print( df2.common )
        
        ncol1   = lapply( df1.common, ncol )
        ncol2   = lapply( df2.common, ncol )
        if( ! identical( ncol1, ncol2 ) )
            {
            mess    = sprintf( "rbind.super() columns with same name in the data.frames do not have matching ncol." )
            cat( mess, '\n' )
            return(NULL)
            }

        out.common  = rbind( df1.common, df2.common )  #; print( out.common )
        }
    else
        out.common  = data.frame( row.names=rownames.out )
    
    
    mask1   = ! (colnames(.df1) %in% common )

    if( any(mask1) )
        {
        out.1   = .df1[ , mask1, drop=FALSE ]   #; print( out.1 )        
        mat     = matrix( NA, nrow(.df2), ncol(out.1), dimnames=list(rnames2,colnames(out.1)) )  
        #   colnames(mat)   = colnames(out.1)        
        out.1   = rbind(out.1,mat)  #; print( out.1 )
        }
    else
        out.1   = data.frame( row.names=rownames.out )
        
        
    mask2   = ! (colnames(.df2) %in% common )

    if( any(mask2) )
        {    
        out.2   = .df2[ , mask2, drop=FALSE ]        
        mat     = matrix( NA,  nrow(.df1), ncol(out.2), dimnames=list(rnames1,colnames(out.2)) )
        #   colnames(mat)   = colnames(out.2)        
        out.2   = rbind(mat,out.2)   #; print( out.2 )
        }
    else
        out.2   = data.frame( row.names=rownames.out )

    out     = cbind( out.1, out.2, out.common ) #; print(out)

    row.names(out)   = rownames.out
    
    return(out)
    }
    
    
    
    
    

#   ...    data.frames with no duplicate rownames

#   always returns a data.frame with the right number of rows
#   and preserves all data from the 2 input data.frames
#   nrow(out) = sum( nrow(...) )

#   if there is only 1 non-trivial schema, then the column order is preserved
#   if there are more than 1, then the columns are sorted alphabetically
#
#   in case of ERROR returns NULL
#
#   compare with plyr::rbind.fill

rbind.super.list <- function( ... )
    {    
    df.list = list(...)
    
    n   = length(df.list)       # of data.frames to bind
    
    #   trivial cases
    if( n == 0 )    return( data.frame() )  # 0 columns and 0 rows
    
    if( n == 1 )    return( df.list[[1]] )
    
    
    rownames.list   = lapply( df.list, rownames ) #; print( rownames.list )
    
    rownames.all    = unlist(rownames.list)      # ;  print( rownames.all )
    
    if( any( duplicated(rownames.all) ) )
        {
        log.string( ERROR, "Row names have duplicates." )
        return(NULL)
        }
        
        
    colnames.list   = lapply( df.list, colnames )       #; print( colnames.list )
    
    # drop the empty names
    colnames.len    = sapply( colnames.list, length )   #; print( colnames.len )
    if( all(colnames.len == 0) )
        {
        # trivial case.  all the data.frames have 0 columns
        out = data.frame( row.names=rownames.all )
        return( out )        
        }
    
    colnames.sub   = colnames.list[ 0 < colnames.len ]  #; print( colnames.sub )
    
    if( all.identical( colnames.sub ) )
        #   only 1 non-trivial schema
        field   = colnames.sub[[1]]
    else
        #   multiple non-trivial schemas.  Just sort all column names
        field   = sort( unique( unlist(colnames.sub) ) )   #; print(field)

    #   print(field)
    
    m    = length(field)    # of columns in output
    
    
    df.vert     = vector( n, mode='list' )
    mask.cname  = logical(n)    
    
    df.horz     = vector( m, mode='list' )

    for( j in 1:m )
        {
        cname   = field[j]  #;    print(cname)
        
        for( i in 1:n )
            {
            mask.cname[i]   = cname %in% colnames.list[[i]] 
            
            if( mask.cname[i] )
                {
                df.vert[[i]] = df.list[[i]] [ , cname, drop=FALSE ]     #; print(df.vert[[i]])
                }
            else
                {
                rnames = rownames.list[[i]]
                df.vert[[i]] = data.frame( rep(NA,length(rnames)), row.names=rnames )
                }
                
            colnames( df.vert[[i]] )    = cname     # these only have 1 column
            }
            
        if( 2 <= sum(mask.cname) )
            {
            #   check that all data.frames in df.vert, for which cname actually appears,
            #   are the same colnames, or NULL in case of a vector column
            colnames.mat   = lapply( df.vert[mask.cname], function(x) { colnames(x[[1]]) } )      #; print(colnames.mat)
            
            if( ! all.identical(colnames.mat) )
                {
                log.string( ERROR, "For column '%s', one data.frame has a matrix, and another has a vector (or a matrix with different # of columns or different column names).",
                                cname )
                return(NULL)
                }
            }
            
        df.horz[[j]] = do.call( rbind, df.vert )
        }
    
    out = do.call( cbind, df.horz )
    
    return(out)
    }
    
    
#   .list   list of objects
    
all.identical  <-  function( .list )
    {
    n   = length(.list)
    if( n <= 1 )    return(TRUE)
    
    for( k in 2:n )
        {
        if( ! identical(.list[[1]],.list[[k]]) )    return(FALSE)
        }

    return(TRUE)
    }
    