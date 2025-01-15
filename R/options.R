
#   g.* are global variables for colorSpec that are "unlocked" in .onAttach()

#   g.options is a list used to initialize the colorSpec.* options managed by the base package
#   after initOptions() is called, g.options is not used again (starting in 2025 that is),
#   except that it is used to determine the correct variable type in cs.options()
g.options   <- list(
                    #loglevel = 'WARN',                                 #   must be character
                    #logformat = "%t{%H:%M:%OS3} %l %n::%f(). %m",      #   must be character
                    stoponerror = TRUE                                  #   must be logical
                    )

#   g.loglevel  is an integer derived from g.options$loglevel which is a string
#               it is managed in logging.R
#g.loglevel  <- 0L               # 0L is invalid

#   g.word[]    is a vector of strings derived from g.options$logformat
#               it is managed in logging.R
#g.word      <- character(0)     # character(0) is invalid




#   set the R global .Options in base
#   this is only called once, from .onLoad()

initOptions <- function( pkgname )
    {
    #unlockBinding( "g.options", asNamespace('colorSpec') )          # asNamespace(pkgname) here generates a NOTE !
    #unlockBinding( "g.word", asNamespace('colorSpec') )             # asNamespace(pkgname) here generates a NOTE !
    #unlockBinding( "g.loglevel", asNamespace('colorSpec') )         # asNamespace(pkgname) here generates a NOTE !


    #   print( g.options )

    #   initialize the colorSpec options managed by base package
    for( opt in names(g.options) )
        {
        fopt    = sprintf( "%s.%s", pkgname, opt )

        value   = base::getOption( fopt )

        if( is.null( value ) )
            {
            #   fopt is not set, so set it now from g.options
            argv        = list( g.options[[ opt ]] )
            names(argv) = fopt
            base::options( argv )         # this changes the global object .Options
            }
        else
            {
            #   the option has already been set in the base package, probably from Profile.site. so do nothing
            # assignPrivateOption( opt, value )
            }
        }

    #   force update of derived variables
    #setLogLevel( .Options$colorSpec.loglevel )      # updates g.loglevel, an integer
    #setLogFormat( .Options$colorSpec.logformat )    # updates g.word[], a character array

    # checkBaseOptions()

    
    return(TRUE)
    }




#   cs.options()
#   retrieve and assign items in g.options AND in the global options()
#
#   ...     items with a valid name are assigned
#           items with no name must be character and if valid are returned returned
cs.options  <-  function(...)
    {
    myfun   = function() { .Options[ grepl( "^colorSpec[.]", names(.Options) ) ] }

    theList = list(...)         #; print(theList)

    n   = length(theList)

    if( n==0 )  return( do.call(myfun,list()) )

    namevec = names(theList)    #; print( namevec )

    if( is.null(namevec) )
        {
        cat( sprintf( "WARN  cs.options() All of the %d arguments are unnamed, and ignored.\n", n ), file=stderr() )
        return( do.call(myfun,list()) )
        }

    #   extract just the args with names and ignore the rest
    mask    = nzchar( namevec )

    idx.unnamed    = which( ! mask )
    if( 0 < length(idx.unnamed) )
        {
        mess    = paste( idx.unnamed, sep=', ', collapse=' ' )
        mess    = sprintf( "WARN  cs.options() Arguments indexed '%s' are unnamed, and ignored.\n", mess )
        cat( mess, file=stderr() )
        }

    theList = theList[ mask ]       #; print(theList)
    namevec = names( theList )      #; print( namevec )


    #   extract just the args with names that partially match the full names, and ignore the rest
    fullname    = names(g.options)

    idx     = pmatch( namevec, fullname, duplicates.ok=TRUE )

    idx.na  = which( is.na(idx) )

    if( 0 < length(idx.na) )
        {
        mess    = paste( namevec[idx.na], collapse=' ' )
        mess    = sprintf( "WARN  cs.options() Arguments named '%s' cannot be matched, and are ignored.\n", mess )
        cat( mess, file=stderr() )

        if( length(idx.na) == length(idx) ) return( do.call(myfun,list())  )
        }

    mask    = ! is.na(idx)

    #print( idx )    ;   print( mask )

    theList     = theList[ mask ]  #    ; print(theList)
    idx         = idx[ mask ]           # to be valid this must have positive length
    fullname    = fullname[idx]

    if( length(theList) == 0 )
        {
        mess    = sprintf( "ERROR  cs.options(). no options are recognized.\n" )
        cat( mess, file=stderr() )
        return( do.call(myfun,list())  )
        }

    names(theList)  = sprintf( "colorSpec.%s", fullname )    # prepend 'colorSpec.'

    #   check the type
    for( k in length(theList):1 )
        {
        value   = theList[[k]]

        thetype = typeof( g.options[[ fullname[k] ]] )

        if( thetype == "logical" )
            ok  = is.logical(value)
        else if( thetype == "integer" )
            ok  = is.integer(value)
        else if( thetype == "double" )
            ok  = is.double(value)
        else if( thetype == "character" )
            ok  = is.character(value)

        if( ! ok )
            {
            mess    = sprintf( "ERROR cs.options()  Value of option %s is '%s', which has the wrong type.\n",
                                             names(theList)[k], as.character(value) )
            cat( mess, file=stderr() )

            #   drop this one
            theList = theList[-k]
            }
        }


    if( 0 < length(theList) )
        {
        #   finally ready to assign global options
        base::options( theList )  # this changes the global object .Options

        #print( theList )
        #cat( "cs.options().   stop=", .Options$colorSpec.stoponerror, '\n'  )

        #   checkBaseOptions()

        # updatePrivateOptions()
        }
    else
        {
        mess    = sprintf( "WARN  cs.options(). no options were assigned.\n" )
        cat( mess, file=stderr() )
        }

    return( do.call(myfun,list())  )
    }





#####################       deadwood below      #########################################


checkBaseOptions    <- function()
    {
    for( name in names(g.options) )
        {
        baseopt    = sprintf( "colorSpec.%s", name )

        value   = base::getOption( baseopt )

        if( is.null(value) )
            {
            mess    = sprintf( "ERROR checkBaseOptions() Option '%s' is unassigned.", baseopt  )
            cat( mess, '\n', file=stderr() )
            next
            }

        if( name == 'stoponerror' )
            ok  = is.logical(value)
        else
            ok  = is.character(value)

        if( ! ok )
            {
            mess    = sprintf( "ERROR checkBaseOptions()  Value of option %s is '%s', which has the wrong type.", baseopt, as.character(value) )
            cat( mess, '\n', file=stderr() )
            }
        }
    }




#   this is only called once, from myonLoad()

unlockPrivateOptions <- function()
    {
    env =  asNamespace('colorSpec')

    #if( ! base::bindingIsLocked( "g.options", env ) )
    #    #   all 3 variables have already been unlocked, nothing to do
     #   return(TRUE)

    #   this is the first time the function has been called
    for( sym in c("g.options","g.word","g.loglevel") )
        {
        base::unlockBinding( sym, env )

        #   now test it !
        if( base::bindingIsLocked( sym, env ) )
            return(FALSE)
        }

    return(TRUE)
    }


#   name        one of 'loglevel', 'logformat', 'stoponerror'
#   value       an appropriate value
assignPrivateOption <- function( name, value )
    {
    g.options[[ name ]]   <<- value

    #   colorSpec::g.options$stoponerror =  .Options$colorSpec.stoponerror      does not work, does not understand the namespace 'colorSpec'

    #   this one works, but takes 3 lines
    #opts    = g.options
    #opts[[ name ]]  = value
    #assign( "g.options", opts, envir=asNamespace('colorSpec') )

    #   test assignment - should be an assert()
    if( ! identical( g.options[[name]], value ) )   cat( "ERROR assignPrivateOption() failed.\n", file=stderr() )

    #   cat( name, sprintf( '%s\n', as.character(g.options[[name]]) ) )
    }

