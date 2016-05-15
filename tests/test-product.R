
require( colorSpec )

#   cs.options( loglevel='D' )

cs.options( stoponerror=FALSE )

testProducts <- function()
    {
    #   find the extdata folder
    if( grepl( "[..]Rcheck", getwd() ) )
        extdata = "../colorSpec/extdata"
    else
        extdata = "inst/extdata"
        
    extdata = system.file( "extdata", package="colorSpec" )    # override
                 
    if( ! file.exists( extdata ) )
        {
        print( getwd() )    
        print( extdata )
        cat( "Cannot find the extdata folder !\n" )
        return(FALSE)
        }
    
    #   read some materials
    IR.blocker  = readSpectra( file.path( extdata, "filters/Midwest-SP700-2014.txt" ) ) # 1 spectrum
    Hematoxylin = readSpectra( file.path( extdata, "stains/Hematoxylin.txt" ) )         # 1 spectrum
    Hoya        = readSpectra( file.path( extdata, "filters/Hoya.txt" ) )               # 4 spectra
    Rosco       = readSpectra( file.path( extdata, "filters/Rosco.txt" ) )              # 42 spectra
    
    #   read some sources
    Lumencor    = readSpectra( file.path( extdata, "sources/Lumencor-SpectraX.txt" ), 380:720 ) # 7 spectra
    
    print( getwd() )
    
    #   create an RGB camera from 3 filters and 1 sensor
    Zyla        = readSpectra( file.path( extdata, "cameras/Zyla_sCMOS.txt" ) )     # 1 spectrum
    cameraRGB   = product( subset(Hoya,1:3), Zyla, wave='auto' )                    # 3 spectra
    
    #   create an RGB scanner
    scannerRGB  = product( D65.1nm, 'MATERIAL', cameraRGB, wave='auto' )            # 3 spectra
    #   summary( scanner )
    
    #   4 product types return a colorSpec object
    
    cat( "--------------  M * ... * M  ------------------\n" )
    junk    = product( IR.blocker, Hoya, Hematoxylin, wave='auto' )
    if( is.null(junk) ) return(FALSE)
   
    junk    = product( Hoya, Rosco, wave='auto' )   # this should fail
    if( ! is.null(junk) ) return(FALSE)
   
   
    cat( "--------------  L * M * ... * M  ------------------\n" )
    junk    = product( D50.5nm, Hoya, IR.blocker, wave='auto' )
    if( is.null(junk) ) return(FALSE)

    junk    = product( Lumencor, Hoya, IR.blocker, wave='auto' )   # this should fail
    if( ! is.null(junk) ) return(FALSE)

    
    cat( "--------------  M * ... * M * R_L  ------------------\n" )
    junk    = product( Hoya, IR.blocker, Zyla, wave='auto' )
    if( is.null(junk) ) return(FALSE)

    junk    = product( subset(Hoya,1:3), IR.blocker, cameraRGB, wave='auto' )
    if( is.null(junk) ) return(FALSE)    
    
    junk    = product( Hoya, IR.blocker, cameraRGB, wave='auto' )   # this should fail
    if( ! is.null(junk) ) return(FALSE)


    cat( "--------------  L * M * 'VARMAT' * M * ... * M * R_L  ------------------\n" )
    junk    = product( Lumencor, IR.blocker, 'VARMAT', Zyla, wave='auto' )
    if( is.null(junk) ) return(FALSE)
    
    junk    = product( Lumencor, IR.blocker, 'VARMAT', cameraRGB, wave='auto' )   # this should fail
    if( ! is.null(junk) ) return(FALSE)
    
    
    
    #   2 product types return a matrix
    
    cat( "--------------  L * M * ... * M * R_L  ------------------\n" )
    junk    = product( D50.5nm, IR.blocker, Hoya, cameraRGB, wave='auto' )  # junk should be a 4 x 3 matrix
    if( is.null(junk) ) return(FALSE)
    if( ! all( dim(junk) == c(4,3) ) )    return(FALSE)
    
    junk    = product( Lumencor, IR.blocker, cameraRGB, wave='auto' )       # junk should be a 7 x 3 matrix
    if( is.null(junk) ) return(FALSE)
    if( ! all( dim(junk) == c(7,3) ) )    return(FALSE)
    
    junk    = product( Lumencor, IR.blocker, Hoya, cameraRGB, wave='auto' ) # this should fail
    if( ! is.null(junk) ) return(FALSE)
    
    
    cat( "--------------  M * ... * M * R_M  ------------------\n" )
    junk    = product( IR.blocker, Rosco, scannerRGB, wave='auto' )         # junk should be a 42 x 3 matrix
    if( is.null(junk) ) return(FALSE)
    if( ! all( dim(junk) == c(42,3) ) )    return(FALSE)
    
    junk    = product( IR.blocker, Hoya, IR.blocker, Rosco, scannerRGB, wave='auto' )         # this should fail
    if( ! is.null(junk) ) return(FALSE)
    
    
    
    #   some invalid sequences
    cat( "--------------  M * ... * M * L  ------------------\n" )
    junk    = product( IR.blocker, Hoya, Hematoxylin, D50.5nm, wave='auto' )         # this should fail
    if( ! is.null(junk) ) return(FALSE)
    
    cat( "--------------  M * ... * L * M   ------------------\n" )
    junk    = product( IR.blocker, Hoya, D50.5nm, Hematoxylin,  wave='auto' )         # this should fail
    if( ! is.null(junk) ) return(FALSE)
    
    cat( "--------------  L * L * R_L  ------------------\n" )
    junk    = product( D50.5nm, D65.1nm, cameraRGB,  wave='auto' )         # this should fail
    if( ! is.null(junk) ) return(FALSE)
    
    cat( "--------------  L * R_M  ------------------\n" )
    junk    = product( D50.5nm, scannerRGB,  wave='auto' )         # this should fail
    if( ! is.null(junk) ) return(FALSE)
    
    
    cat( "\nPassed all product tests !\n" )
    
    return( TRUE )
    }
    
    
if( ! testProducts() )  stop( "testProducts() failed !" )

        