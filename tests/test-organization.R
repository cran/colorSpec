


library( colorSpec )

rotateOrganization <- function()
    {
    #   start with matrix and return again
    testobj = colorSpec::xyz1931.5nm

    organization(testobj)    = 'matrix'
    testobj.mat = testobj  # save for later
    organization( testobj ) = 'df.col'
    organization( testobj ) = 'df.row'
    organization( testobj ) = 'matrix'

    # override class value, so identical() effectively ignores it
    # this is necessary for R v 4.0, which changes class(a matrix object) to c('matrix','array')
    #print( class(testobj) )
    class(testobj)      = class(testobj.mat)
    dimnames(testobj)   = dimnames(testobj.mat)     # added Oct 28, 2023
    
    if( ! identical(testobj,testobj.mat) )
        {
        print( str(testobj) )
        print( str(testobj.mat) )
        stopifnot( identical(testobj,testobj.mat) )
        }

    #   start with vector and return again
    testobj     = colorSpec::C.5nm
    testobj.vec = colorSpec::C.5nm
    organization( testobj ) = 'matrix'
    organization( testobj ) = 'df.col'
    organization( testobj ) = 'df.row'
    organization( testobj ) = 'vector'
    
    class(testobj)      = class(testobj.vec)
    dimnames(testobj)   = dimnames(testobj.vec)     # added Oct 28, 2023
        
    if( ! identical(testobj,testobj.vec) )
        {
        print( str(testobj) )
        print( str(testobj.vec) )
        stopifnot( identical(testobj,testobj.vec) )
        }

    cat( "\n", "Passed all organization tests.", "\n", sep='' )

    return( TRUE )
    }

rotateOrganization()
