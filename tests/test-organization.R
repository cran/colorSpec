


library( colorSpec )

rotateOrganization <- function()
    {
    #   start with matrix and return again
    xyz1931.5nm = colorSpec::xyz1931.5nm
    
    organization(xyz1931.5nm)    = 'matrix'
    xyz1931.5nm.mat    = xyz1931.5nm  # save for later
    organization( xyz1931.5nm ) = 'df.col'
    organization( xyz1931.5nm ) = 'df.row'
    organization( xyz1931.5nm ) = 'matrix'
    stopifnot( identical(xyz1931.5nm,xyz1931.5nm.mat) )
    
    #   start with vector and return again 
    C.5nm       = colorSpec::C.5nm
    C.5nm.vec   = colorSpec::C.5nm
    organization( C.5nm ) = 'matrix'
    organization( C.5nm ) = 'df.col'
    organization( C.5nm ) = 'df.row'
    organization( C.5nm ) = 'vector'
    stopifnot( identical(C.5nm,C.5nm.vec) )
    
    cat( "\n", "Passed all organization tests.", "\n", sep='' )
    
    return( TRUE )
    }
    
rotateOrganization()
