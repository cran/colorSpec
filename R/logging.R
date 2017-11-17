


g.word      = character(0)


setLogFormat <- function( .fmt )
    {
    #   print( "Entering setLogFormat()" )
    
    pattern = "%[a-zA-Z](\\{.+\\})?"
    
    res = gregexpr( pattern, .fmt )[[1]]    #; print( res )
    
    if( res[1] < 0 )
        {
        g.word  = .fmt[1]
        return(TRUE)
        }
    
    len = attr( res, "match.length" )
    
    code    = integer( nchar(.fmt) )
    
    for( k in 1:length(res) )
        {
        code[ res[k]:(res[k]+len[k]-1) ] = k
        }
    #print( code )
    
    rle = rle(code) #; print( rle )
    
    n   = length(rle[[1]])

    start   = c( 1, cumsum( rle$lengths ) + 1 )[1:n]        #;   print(start)
    stop    = start + rle$lengths-1                         #; print(stop)
    
    g.word  <<- character( n )  # get the size right, shrink it !
    for( k in 1:n )
        g.word[k]   <<- substr( .fmt, start[k], stop[k] )
        
    #   log.object( WARN, g.word )
    
    return(invisible(TRUE))    
    }
    
    

    
log.string <- function( level, msg, ... )
    {    
    conn    = stderr()
    
    if( ! is.integer(level) )
        {
        cat( "ERROR  log.string(). level is not an integer.\n", file=conn )
        return( invisible(FALSE) )
        }
    
    if( g.options$loglevel < level )
        #   do nothing
        return( invisible(FALSE) )
    
    msg = sprintf( msg, ... )
    
    #   print( g.word )
        
    idx = which( grepl("^%",g.word) )   #; print( idx )
    if( length(idx) == 0 )
        {
        cat( msg, '\n', file=conn )
        return( invisible(TRUE) )
        }
        
    word    = g.word
    
    #   where   = sys.parent(1) ; print(where)
            
    for( k in idx )
        {
        spec    = substr( g.word[k], 1, 2 )
        
        if( spec == "%l" )
            word[k] = sprintf( "%-5s", names(level) )
        else if( spec == "%t" )
            {
            fmt = ''
            m   = nchar(g.word[k]) 
            if( 5 <= m ) fmt = substr( g.word[k], 4, m-1 )
            word[k] = format( Sys.time(), fmt )
            }
        else if( spec == "%n" )
            word[k] = "colorSpec"
        else if( spec == "%f" )
            {
            where   = sys.parent(1)  # ; print(where)
          
            if( 0 < where )
                word[k] = tryCatch( deparse(sys.call(where)[[1L]]), error=function(e) "[console]" )
            else
                word[k] = "[console]"
            }
        else if( spec == "%m" )
            word[k] = msg
        }
    
    # print( word )
    #   cat( paste0(word,collapse=''), '\n' )   ; flush.console()
    
    mess    = paste0(word,collapse='')
    #   message( mess )   ; flush.console()
    
    #print( sys.parent(1) )
    #print( deparse(sys.call(-3L)) )        
    #print( deparse(sys.call(-2L)[[1L]]) )    
    #print( deparse(sys.call(-1L)[[1L]]) )    
    #print( deparse(sys.call(0L)[[1L]]) )    
    #print( deparse(sys.call(1L)[[1L]]) )    
    #print( deparse(sys.call(2L)[[1L]]) )    
    #print( deparse(sys.call(3L)[[1L]]) )    
 
   
    if( g.options$stoponerror  &&  level <= ERROR )
        stop( mess, '\n', "Stopping, because option stoponerror==TRUE", call.=FALSE )

    cat( mess, '\n', file=conn ) ;   flush(conn)

    return( invisible(TRUE) )
    }
    
log.object <- function( level, obj, type='whole', addname=TRUE )
    {    
    conn    = stderr()
    
    if( ! is.integer(level) )
        {
        cat( "ERROR  log.object(). level is not an integer.\n", file=conn )
        return( invisible(FALSE) )
        }
        
    if( g.options$loglevel < level )
        #   do nothing
        return( invisible(FALSE) )
        
    if( type == 'whole' )
        {
        line    = capture.output( print(obj) )
        if( addname )   
            {
            cat( deparse(substitute(obj)), file=conn )
            if( 1 < length(line) )
                cat( '\t=\n', file=conn )
            else
                cat( '\t=\t', file=conn )
            }
            
        for( k in 1:length(line) )
            cat( line[k], '\n', file=conn )
        }
    else if( type == 'str' )
        {
        if( addname )   
            {
            cat( deparse(substitute(obj)), '\n', file=conn )            
            }
        print( str(obj) )
        }
        
    flush.console()
        
    return( invisible(TRUE) )
    }
        
        
testLogging <- function()
    {
    log.string( ERROR, "Hello %s.  Matrix A is:", "World" )
    
    A   =  matrix( runif(3*3), 3, 3 )
    log.object( ERROR, A )
    
    return( invisible(TRUE) )
    }
    