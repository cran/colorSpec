
library( colorSpec )

cs.options( loglevel='TRACE' )

testCGATS <- function()
    {
    mess    = sprintf( "in testCGATS().  getwd() = '%s'", getwd() )
    cat( mess, '\n', file=stderr() )
    
    #   find the extdata folder
    if( grepl( "[..]Rcheck", getwd() ) )
        extdata = "../colorSpec/extdata"
    else
        extdata = "inst/extdata"
        
    #   extdata = system.file( "extdata", package="colorSpec" )    # override
    extdata  = '.'

    
    if( ! file.exists( extdata ) )
        {
        print( getwd() )    
        print( extdata )
        cat( "Cannot find the extdata folder !\n", file=stderr() )
        return(FALSE)
        }
        

        
    #   test-CGATS-2-2S.txt 2 tables, both are spectral
    path    =   file.path( extdata, "test-CGATS-2-2S.txt" )    
    junk    = readCGATS(path)   # should be 2 tables
    if( length(junk) != 2 )
        return(FALSE)
    junk    = readSpectraCGATS(path)  # should be 2 tables
    if( length(junk) != 2 )
        return(FALSE)
        
    #   test-CGATS-2-1S.txt 2 tables, only the 1st one is spectral
    path    =   file.path( extdata, "test-CGATS-2-1S.txt" )    
    junk    = readCGATS(path)           # should be 2 tables
    if( length(junk) != 2 )
        return(FALSE)
    junk    = readSpectraCGATS(path)    # should be 1 table 
    if( length(junk) != 1 )
        return(FALSE)
        
    #   test-CGATS-2-1S-swap.txt 2 tables, only the 2nd one is spectral
    path    =   file.path( extdata, "test-CGATS-2-1S-swap.txt" )    
    junk    = readCGATS(path)           # should be 2 tables
    if( length(junk) != 2 )
        return(FALSE)
    junk    = readSpectraCGATS(path)    # should be 1 table 
    if( length(junk) != 1 )
        return(FALSE)
        
        

    #   A70.ti3 has 2 tables, both non-spectral
    path    =   file.path( extdata, "A70.ti3" )
    junk    = readCGATS(path)   # should be 2 tables
    if( length(junk) != 2 )
        return(FALSE)
        
    #   the final call should generate an ERROR and return NULL, so disable stopping
    cs.options( stoponerror=FALSE )
    junk    = readSpectraCGATS(path) 
    if( ! is.null(junk) )
        return(FALSE)
        
    return(TRUE)
    }
    

    
if( ! testCGATS() )  stop( "testCGATS() failed !" )

cat( "Passed all CGATS read tests !\n", file=stderr() )
