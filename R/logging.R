

#   this function does not add level or timestamp or function call, or anything like that
#   it's meant be called right after a regular sprintf() logging event

log_object <- function( level, obj, type='whole', addname=TRUE )
    {
    if( ! is.integer(level) )
        {
        warning( "log_object(). level is not an integer."  )
        return( invisible(FALSE) )
        }

    if( logger::log_threshold(namespace='colorSpec') < level )
        #   do nothing
        return( invisible(FALSE) )
        
    conn    = stderr()

    if( type == 'whole' )
        {
        line    = capture.output( print(obj) )
                
        if( addname )
            {
            cat( deparse(substitute(obj)), file=conn )
            if( 1 < length(line) )
                cat( ' =\n', file=conn )
            else
                cat( ' =\t', file=conn )
            }
            
        logger::appender_console( line )
        }
    else if( type == 'str' )
        {
        line    = capture.output( print( str(obj) ) )
        
        if( addname )
            {
            cat( deparse(substitute(obj)), '\n', file=conn )
            }
            
        logger::appender_console( line )        
        }

    flush.console()

    return( invisible(TRUE) )
    }
    

#   this function is made to be called by me, from global env
testLogging <- function()
    {
    log_level( ERROR, "Hello %s.  Matrix A is:", "World" )

    A   =  matrix( runif(3*3), 3, 3 )
    log_object( ERROR, A )

    return( invisible(TRUE) )
    }

    

#####################       deadwood below      #########################################

setLogFormat <- function( .fmt )
    {
    #   print( "Entering setLogFormat()" )

    if( ! is.character(.fmt) )
        {
        if( is.null(.fmt) )
            mess = "setLogFormat() logformat=NULL  - ignored."
        else
            mess = sprintf( "setLogFormat() logformat='%s' is not a string - ignored.", as.character(.fmt) )

        warning(mess)
        return(NA_integer_)
        }

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

    #   log_object( WARN, g.word )

    #   now record it in g.options
    assignPrivateOption( 'logformat', .fmt )

    return(invisible(TRUE))
    }

