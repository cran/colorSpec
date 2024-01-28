
library( colorSpec )

testSink <- function()
    {
    cs.options( stop=FALSE )

    path    = 'junk.log'

    #   start the sink
    con = file( path, open = "wt")
    sink( file=con, append=FALSE, type="message" )

    #   read 2 non-existent files
    readSpectra( "nonexistent.txt" )
    readSpectra( "nonexistent2.txt" )

    #   stop the sink
    sink( NULL, type='message' )
    close( con )

    cs.options( stop=TRUE )

    #   path should have 4 lines in it
    line    = readLines( path )   #;   cat( length(line), '\n' )

    return( length(line) == 4 )
    }


if( ! testSink() )  stop( "testSink() failed !" )

cat( "Passed all sink tests !\n", file=stderr() )
