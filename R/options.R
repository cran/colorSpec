
#   g.* are global variables for colorSpec that are "unlocked" in .onAttach()

#   g.options is a list with private copies of the colorSpec.* options managed by the base package
#   these copies are used because g.loglevel and g.word[] are derived variables,
#   and we want to see whether the user has changed a global colorSpec.* option
g.options   <- list(loglevel = 'WARN',                                  #   must be character
                    logformat = "%t{%H:%M:%OS3} %l %n::%f(). %m",       #   must be character
                    stoponerror = TRUE                                  #   must be logical
                    )

#   g.loglevel  is an integer derived from g.options$loglevel which is a string
#               it is managed in logging.R
g.loglevel  <- 0L               # 0L is invalid

#   g.word[]    is a vector of strings derived from g.options$logformat
#               it is managed in logging.R
g.word      <- character(0)     # character(0) is invalid



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
        cat( sprintf( "WARN  cs.options() All of the %d arguments are unnamed, and ignored\n", n ), file=stderr() )
        return( do.call(myfun,list()) )
        }

    #   extract just the args with names and ignore the rest
    mask    = nzchar( namevec )

    idx.unnamed    = which( ! mask )
    if( 0 < length(idx.unnamed) )
        {
        mess    = paste( idx.unnamed, sep=', ', collapse=' ' )
        mess    = sprintf( "WARN  cs.options() Arguments indexed '%s' are unnamed, and ignored.", mess )
        cat( mess, '\n', file=stderr() )
        }

    theList = theList[ mask ]       #; print(theList)
    namevec = names( theList )      #; print( namevec )


    #   extract just the args with names that partially match the full names, and ignore the rest
    fullname    = names(g.options)
    #fullname    = c( "loglevel", "logformat", "stoponerror" )

    idx     = pmatch( namevec, fullname )

    idx.na  = which( is.na(idx) )

    if( 0 < length(idx.na) )
        {
        mess    = paste( namevec[idx.na], collapse=' ' )
        mess    = sprintf( "WARN  cs.options() Arguments named '%s' cannot be matched, and are ignored.", mess )
        cat( mess, '\n', file=stderr() )

        if( length(idx.na) == length(idx) ) return( do.call(myfun,list())  )
        }

    mask    = ! is.na(idx)

    #print( idx )    ;   print( mask )

    theList = theList[ mask ]  #    ; print(theList)
    idx     = idx[ mask ]   # this must have positive length

    names(theList)  = sprintf( "colorSpec.%s", fullname[idx] )    # prepend

    #print( theList )

    #   finally ready to assign global options
    options( theList )  # this changes the global object .Options

    checkBaseOptions()

    updatePrivateOptions()

    return( do.call(myfun,list())  )
    }


#   updatePrivateOptions    copy values from public base package .Options to private copies in g.options
#                           and update derived globals (g.word and g.loglevel) if necessary
updatePrivateOptions <- function()
    {
    #if( ! unlockPrivateOptions() )
    #    {
    #    #   should not happen, this is FATAL
    #    mess = sprintf( "FATAL  Cannot unlock private options." )
    #    stop( mess, '\n', call.=FALSE )
    #    }
        

    if( ! identical( g.options$loglevel, .Options$colorSpec.loglevel ) )
        {
        setLogLevel( .Options$colorSpec.loglevel )
        }

    if( ! identical( g.options$logformat, .Options$colorSpec.logformat ) )
        {
        setLogFormat( .Options$colorSpec.logformat )
        }

    if( is.logical(.Options$colorSpec.stoponerror) )
        {
        if( ! identical( g.options$stoponerror, .Options$colorSpec.stoponerror ) )
            {
            assignPrivateOption( 'stoponerror', .Options$colorSpec.stoponerror )
            }
        }
    else
        {
        mess = sprintf( "WARN  updatePrivateOptions() colorSpec.stoponerror='%s' is not logical - ignored.",
                            as.character(.Options$colorSpec.stoponerror) )
        cat( mess, '\n', file=stderr() )
        }

    return(TRUE)
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


checkBaseOptions    <- function()
    {
    for( name in names(g.options) )
        {
        baseopt    = sprintf( "colorSpec.%s", name )

        value   = getOption( baseopt )

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
    
    
    

#   synchronize R global .Options in base,  and g.options in colorSpec
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

        value   = getOption( fopt )

        if( is.null( value ) )
            {
            #   fopt is not set, so set it now from g.options
            argv        = list( g.options[[ opt ]] )
            names(argv) = fopt
            base::options( argv )         # this changes the global object .Options
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
    return(TRUE)
    }
    
    
    