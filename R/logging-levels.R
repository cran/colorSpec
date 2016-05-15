
#   these numbers do not have to be exported
#   the user can use the unique initial letter instead

FATAL <- 1L
names(FATAL) <- "FATAL"
ERROR <- 2L
names(ERROR) <- "ERROR"
WARN <- 4L
names(WARN) <- "WARN"
INFO <- 6L
names(INFO) <- "INFO"
DEBUG <- 8L
names(DEBUG) <- "DEBUG"
TRACE <- 9L
names(TRACE) <- "TRACE"

logLevelFromString  <- function( .level )
    {    
    if( ! is.character(.level) )    return(FALSE)
    
    #   convert to the integer
    w   = toupper( substr(.level,1,1) ) #; print(w)
    
    if( w == "F" )
        return(FATAL)
    else if( w == "E" )
        return(ERROR)
    else if( w == "W" )
        return(WARN)
    else if( w == "I" )
        return(INFO)      
    else if( w == "D" )
        return(DEBUG)            
    else if( w == "T" )
        return(TRACE)
    else
        {
        warning( sprintf( "logLevelFromString() loglevel='%s' invalid, and ignored.", .level ) )
        return(FALSE)
        }
    }

    