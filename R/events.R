
g.microbenchmark    = FALSE     # logical value, whether the package microbenchmark is loaded.  It must be unlocked.


    
.onLoad <- function( libname, pkgname )
    {            
    #   requireNamespace( "utils" )
    #   desc    = packageDescription( pkgname )
    #   mess = sprintf( "This is %s %s.  %s.  Author: %s  Built: %s\n", 
    #                    pkgname, desc$Version, desc$Title, desc$Author, desc$Built )
    #   print( pkgname )
    
    #   unlockBinding() fails here in .onLoad(), so use .onAttach() instead
    

    #mess    = environmentName( environment(.onLoad) )
    #mess = sprintf( "Environment: '%s'\n", mess )
    #cat( mess )
    }
    

.onAttach <- function( libname, pkgname )
    {
    #print( libname )
    #print( pkgname )
    
    info    = library( help='colorSpec' )        #eval(pkgname) 
    info    = format( info )
    mask    = grepl( "^(Version|Author|Built)", info )     #Title
    info    = gsub( "[ ]+", ' ', info[mask] )
    mess    = sprintf( "This is %s", pkgname )
    mess    = paste( c( mess, info ), collapse='.  ' )   #; cat(mess)
    packageStartupMessage( mess )

    unlockBinding( "g.options", asNamespace('colorSpec') )          # asNamespace(pkgname) here generates a NOTE !
    unlockBinding( "g.word", asNamespace('colorSpec') )             # asNamespace(pkgname) here generates a NOTE !
    unlockBinding( "g.loglevel", asNamespace('colorSpec') )         # asNamespace(pkgname) here generates a NOTE !
    unlockBinding( "g.microbenchmark", asNamespace('colorSpec') )   # asNamespace(pkgname) here generates a NOTE ! 
    # unlockBinding( "g.Trap30x301", asNamespace('colorSpec') )       # asNamespace(pkgname) here generates a NOTE ! 
    
    g.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "g.microbenchmark=", g.microbenchmark, '\n' )
        
    
    #   print( g.options )  
    
    #   initialize the colorSpec options managed by base package    
    for( opt in names(g.options) )
        {
        fopt    = sprintf( "%s.%s", pkgname, opt )

        value   = getOption( fopt )
        
        if( is.null( value ) )
            {
            #   fopt is not set, so set it now from g.options
            argv        = list( g.options[[ opt ]] )
            names(argv) = fopt
            options( argv ) 
            }
        else
            {
            #   option has been set in the base package, probably from Profile.site
            assignPrivateOption( opt, value )
            }
        }
        
    #   force update of derived variables
    setLogLevel( .Options$colorSpec.loglevel )      # updates g.loglevel, an integer
    setLogFormat( .Options$colorSpec.logformat )    # updates g.word[], a character array
        
    checkBaseOptions()    
    
    #   the options are now in synch !
    }
    

    