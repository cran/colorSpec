
library( colorSpec )


testAllReads <- function()
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
        
    print( extdata )
    
    pathvec = list.files( extdata, recursive=TRUE, full=TRUE ) #; print( pathvec )
    
    #   ignore some file extensions
    pattern1    = "[.](png|jpg|xls|m|htm|html)$"
    mask1       = grepl( pattern1, pathvec, ignore.case=TRUE )
    
    pattern2    = "^dataCCT|^illuminant"
    mask2       = grepl( pattern2, basename(pathvec), ignore.case=TRUE )
    
    pathvec     = pathvec[ ! (mask1 | mask2) ]
    
    mess    = sprintf( "Found %d files, reading them all...\n", length(pathvec) )
    cat( mess )
    
    cs.options( stoponerror=TRUE )
        
    for( i in 1:length(pathvec) )
        {
        path    = pathvec[i]
        
        mess    = sprintf( "--------------  %s  ------------------\n", basename(path) )
        cat(mess)
        
        junk    = readSpectra( path, 400:700 )
        
        if( is.null(junk) ||  ! is.colorSpec(junk) )
            return(FALSE)
        }
    
    return(TRUE)
    }
    


checkQuantity <- function()
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
        
    testdata    = c(
        "illuminants/D65.1nm.txt",      "power",
        "illuminants/C.txt",            "power",     
        "illuminants/daylight1964.txt", "power",     
        "sources/BlueFlame.txt",        "power",     
        "sources/pos1-20x.scope",       "power",     
          
        
        "eyes/ciexyz31_1.csv",          "power->neural",
        "eyes/ciexyz64_1.csv",          "power->neural",
        "eyes/xyz2012.csv",             "power->neural",
                
#       "../inst/extdata/cameras/Flea2-spectral.txt",   "power->electrical",  
        "cameras/orthicon-5820-A.txt",  "power->electrical",  
        "cameras/Zyla_sCMOS.txt",       "photons->electrical",
        
        "action/BeanPhotosynthesis.txt",        "photons->action",
        
        "objects/Hoya.txt",                     "transmittance",
        "objects/Rosco.txt",                    "transmittance",        
        "objects/Midwest-SP700-2014.txt",       "transmittance",        

        "targets/CC_Avg20_spectrum_XYY.txt",    "reflectance",            
        "targets/N130501.txt",                  "transmittance",            
        "stains/Hematoxylin.txt",               "absorbance",            
        
        "scanners/SMPTE-ST-2065-2.txt", "material->electrical"        
        )
        
    ok  = TRUE
    
    testdata    = matrix( testdata, length(testdata)/2, 2, byrow=T )
    
    for( i in 1:nrow(testdata) )
        {
        path    =   file.path( extdata, testdata[i,1] )
        
        mess    = sprintf( "--------------  %s  ------------------\n", basename(path) )
        cat(mess)
        
        junk    = readSpectra(path)
        
        if( is.null(junk) ) 
            {
            ok  = FALSE
            break
            }
            
        quantity    = quantity(junk)
        mess    = sprintf( "quantity='%s'\n", quantity )
        cat(mess)
        
        if( quantity(junk) != testdata[i,2] )
            {
            mess    = sprintf( "'%s' != '%s'\n", testdata[i,2], quantity(junk) )
            cat(mess)
            ok  = FALSE
            break
            }
        }
        
    if( TRUE )
    {
    path    =  file.path( extdata, "cameras/Flea2-spectral.txt" )
    mess    = sprintf( "--------------  %s  ------------------\n", basename(path) )
    cat(mess)    
    junk    = readSpectra( path, 400:700, span=0.20 )
    
    if( quantity(junk) == "power->electrical" )
        {
        #   print( summary(junk) )
        }
    else
        {
        mess    = sprintf( "%s.   '%s' != '%s'\n", path, "power->electrical", quantity(junk) )
        cat(mess)
        ok  = FALSE
        }
    }

    
    cs.options( stoponerror=FALSE )
    
    for( path in c("test-combo1.txt","test-combo2.txt") )
        {
        mess    = sprintf( "--------------  %s  ------------------\n", basename(path) )
        cat(mess)    
        junk    = readSpectra( path )   # this should fail
        if( ! is.null(junk) )   return(FALSE)
        
        junk    = readSpectra( path, 400:700 )   # this should succeed
        if( is.null(junk) )     return(FALSE)
        }
    
    print( warnings() )
        

            
    return( ok )
    }
    
if( ! testAllReads() )  stop( "testAllReads() failed !" )

if( ! checkQuantity() ) stop( "checkQuantity() failed !" )

cat( "Passed all read tests !\n" )
