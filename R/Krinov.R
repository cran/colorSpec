

importKrinov <- function( .folder )
    {
    pathvec = list.files(  .folder, "[0-9]+", recursive=T, full.names=T )
    n       = length(pathvec)
    if( n == 0 )    return(NULL)
    
    wave    = seq( 400, 650, by=10 )
    
    mat     = matrix( 0, length(wave), n )
    idvec   = integer( n )
    namevec = character( n )
    
    pattern = "^([0-9]+)[  ]+([^ ].*$)"
        
    for( k in 1:n )
        {
        path    = pathvec[k]
        
        df  = read.table( path, header=F, skip=2 )
        #   print( df )
        #   break
        
        ok  = all( df$V1 == wave )
        if( ! ok )
            {
            mess    = sprintf( "bad wavelengths in '%s'\n", path )
            cat(mess)
            print(df)
            return(NULL)
            }
            
        if( any( 1 < df$V2 ) )
            #   probably lost the decimal point - Oak and Quay weird
            df$V2 = df$V2 / 1000
            
        mat[  ,k]   = df$V2
        
        line1       = readLines( path, n=1 )
        idvec[k]    = as.integer( sub( pattern, "\\1", line1 ) )
        namevec[k]  = sub( pattern, "\\2", line1 )
        }
        
        
    #   fix duplicated names
    name.unique = unique(namevec)
    
    for( k in 1:length(name.unique) )
        {
        thename = name.unique[k]
        
        idx = which( namevec == thename )
        if( length(idx) <= 1 )  next
        
        for( j in 1:length(idx) )
            namevec[ idx[j] ] = sprintf( "%s %d", thename, j )
        }
        
        
    
    
    colnames(mat)   = namevec
    out = colorSpec( mat, wavelength=wave, quantity='reflectance', organization='df.row' )
    
    extradata( out )    = data.frame( SAMPLE_ID=idvec, Group=findInterval( idvec, c( 1,50,169,176,230,319,333,355) ) )
    
    return( out )
    }
    
    