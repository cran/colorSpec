
g.microbenchmark    = FALSE     # logical value, whether the package microbenchmark is loaded, and therefore available.

#
#   the 3 hooks below are called in this order
#
#   1)  .onLoad()   synchronizes the global R colorSpec.* options with the private options (3 of them)
#   2)  myonLoad()  unlocks the private options
#   3)  .onAttach() prints a startup message.  It may not be called at all.
#

.onLoad <- function( libname, pkgname )
    {
    #   requireNamespace( "utils" )
    #   desc    = packageDescription( pkgname )
    #   mess = sprintf( "This is %s %s.  %s.  Author: %s  Built: %s\n",
    #                    pkgname, desc$Version, desc$Title, desc$Author, desc$Built )
    #   print( pkgname )

    #   unlockBinding( "g.microbenchmark", asNamespace('colorSpec') )   # asNamespace(pkgname) here generates a NOTE !

    #   unlockBinding() is not needed for g.microbenchmark, because it is currently unlocked in .onLoad()
    g.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "g.microbenchmark=", g.microbenchmark, '\n' )

    #mess    = environmentName( environment(.onLoad) )
    #mess = sprintf( "Environment: '%s'\n", mess )
    #cat( mess )

    initOptions( pkgname )

    #   base::options( warnPartialMatchArgs=TRUE )      #   suggested by Pedro Aphalo

    #   after exiting this function, private options are then locked
    #   set hook to unlock private options
    base::setHook( base::packageEvent(pkgname,"onLoad"), myonLoad )
    }



myonLoad <- function( pkgname, pkgpath )
    {
    #  packageStartupMessage( "myonLoad()" )

    if( ! unlockPrivateOptions() )
        {
        #   should not happen, this is FATAL
        mess = sprintf( "FATAL  Cannot unlock private options." )
        stop( mess, '\n', call.=FALSE )
        }
    }






.onAttach <- function( libname, pkgname )
    {
    #packageStartupMessage( libname )
    #packageStartupMessage( pkgname )

    info    = library( help='colorSpec' )        #eval(pkgname)
    info    = format( info )
    mask    = grepl( "^(Version|Author|Built)", info )     #Title
    info    = gsub( "[ ]+", ' ', info[mask] )
    mess    = sprintf( "Attaching %s", pkgname )
    mess    = paste( c( mess, info ), collapse='.  ' )   #; cat(mess)
    packageStartupMessage( mess )

    #initOptions()
    }


