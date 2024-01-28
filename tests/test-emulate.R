
library( colorSpec )

testEmulate <- function()
    {
    wave = 400:700
    
    # read the 2 cameras
    path = system.file( 'extdata/cameras/Plumbicon30mm.txt', package='colorSpec') 
    plumbicon =  readSpectra( path, wavelength=wave )
    path = system.file( 'extdata/cameras/Red-Epic-Dragon.txt', package='colorSpec') 
    dragon = readSpectra( path, wavelength=wave )
    dragon.mod = emulate( dragon, plumbicon, filter=FALSE, matrix=TRUE ) 
    
    if( ! all( 0 < dragon.mod ) )   return(FALSE)
    
    return(TRUE)
    }
    

    
if( ! testEmulate() )  stop( "testEmulate() failed !" )

cat( "Passed all emulation tests !\n", file=stderr() )
