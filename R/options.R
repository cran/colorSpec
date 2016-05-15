
#   g.options is a global variable for colorSpec that is "unlocked" in .onAttach()
g.options <- list (loglevel = WARN,
                  logformat = "%t{%H:%M:%OS3} %l %n::%f(). %m",                   
                  stoponerror = TRUE
                  )

#   cs.options()
#   retrieve and assign items in g.options
#
#   ...     items with a valid name are assigned
#           items with no name must be character and if valid are returned returned

cs.options  <-  function(...)
    {
    theList = list(...)         #; print(theList)
    
    n       = length(theList)
    namevec = names(theList)    #; print( namevec )
    
    if( n == 0 )    return( g.options )     # return them all
    
    if( is.null(namevec) )
        {
        #   return a sublist of g.options
        #print( "no names !  return sublist !" )
        namevec = theList[ sapply( theList, is.character) ]
        namevec = unlist( namevec )    #; print(out)
        
        if( length(namevec) == 0 )  return(NULL)
        
        out = g.options[ namevec ]
        
        return( out )
        }
    

    #   split theList into a list with names and one without
    mask    = nzchar( namevec )
    
    list.named  = theList[ mask ]   
    list.noname = theList[ ! mask ]   
    
    #   iterate over the list of options to assign
    m   = length(list.named)
    if( 0 < m )
        {
        opts    = g.options     # make copy for modification
        
        for( k in seq_len(m) )
            {
            theName     = names(list.named)[k]
            theValue    = list.named[[k]]
            
            if( theName == 'loglevel' && is.character(theValue) )
                {
                #   convert to the integer
                theValue = logLevelFromString( theValue )
                
                if( ! is.integer(theValue) )    next    # ignore unknown string                     
                }            
            
            idx = which( theName == names(opts) )
            
            if( length(idx) == 1 )
                {
                # got a match !
                opts[[idx]] = theValue
                }
            else
                log.string( ERROR, "Unknown name '%s'.", theName )
            }
            
            
        env = asNamespace('colorSpec')      # ; print( env )

        assign( "g.options", opts, envir=env )  #; print( g.options )
            
        if( "logformat" %in% names(list.named) )
            {
            setLogFormat( g.options$logformat )        
            #   print( g.word )
            }                 
        }
             
    #   return all values mentioned in theList
    out = names(list.named)
    
    if( 0 < length(list.noname) )
        {
        namevec = list.noname[ sapply(list.noname, is.character) ]
        out = c( out, unlist( namevec ) )
        }
        
    out = g.options[ out ]
    
    return( out )
    }
                  
                  
#############################  deadwood below     ############################
                  
##' Options for package colorSpec
##' Functions to access and set colorSpec's options.
##' 

cs.getOptions <- function (...)
    {
    dots <- c(...)
    if (length (dots) == 0L)
        #   return the whole list
        g.options
    else if( length (dots) == 1L)
        #   return single item 
        g.options[dots][[1]]
    else
        #   return sublist
        g.options[dots]
    }



cs.setOptions <- function(...)
    {
    new <- list (...) #; print(new)
    names <- nzchar (names (new)) 

    if (! all (names))
        warning ("options without name are discarded: ", which (! names))

    #   print( new[names] )
    
    loglevel.prev   =   g.options$loglevel    
    
    opts <- modifyList (g.options, new [names])

    #if (sys.parent() == 0  &&  isNamespace("colorSpec")) 
    #    env <- asNamespace("colorSpec")
    #else
    #    env <- parent.frame()
     

    #print( str(g.options[["loglevel"]]) )
    
    if( is.character(opts$loglevel) )
        {
        #   convert to the integer
        level = logLevelFromString( opts$loglevel )
        
        if( is.integer(level) ) 
            opts$loglevel   = level
            #assign( "g.options", , envir=env )     # success
            #g.options[["loglevel"]] <<- level
        else
            opts$loglevel   = loglevel.prev        
            #assign( "g.options$loglevel", loglevel.prev, envir=env )     # revert to previous value        
            #g.options[["loglevel"]] <<- loglevel.prev
        }
        
    env = environment(cs.setOptions)    # ;  print( env )

    env = asNamespace('colorSpec')      # ; print( env )
    
    # unlockBinding("g.options", env )
        
    assign( "g.options", opts, envir=env )  #; print( g.options )
        
    if( "logformat" %in% names(new) )
        {
        setLogFormat( g.options$logformat )        
        #   print( g.word )
        }

    invisible(g.options)
    }

