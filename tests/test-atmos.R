
library( colorSpec )

testAtmos <- function()
    {
    wave    = 400:700
    
    junk = atmosTransmittance( c(500,1000,2000), wavelength=wave, aerosols=list(metrange=1000) )
        
    #   for special sequence wave[], reference wavelength 550nm is index 151
    err = junk[ 550 - wave[1] + 1,  2 ] - 0.02
    
    if( 1.e-10 < abs(err) ) return(FALSE)
    
    return(TRUE)
    }
    

    
if( ! testAtmos() )  stop( "testAtmos() failed !" )

cat( "Passed all atmosphere tests !\n", file=stderr() )
