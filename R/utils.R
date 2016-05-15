

.onAttach <- function (libname, pkgname)
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


    unlockBinding( "g.options", asNamespace('colorSpec') )      # asNamespace(pkgname) here generates a NOTE !
    unlockBinding( "g.word", asNamespace('colorSpec') )         # asNamespace(pkgname) here generates a NOTE !
    
    #   force initialization of variable g.word, which is in logging.R
    cs.options( logformat=g.options$logformat )              

    #   unlockBinding("g.options[['loglevel']]", asNamespace("colorSpec")) ; cat( "unlocked g.options$loglevel\n" )
    }
    
    
.onLoad <- function( libname, pkgname )
    {            
    #   requireNamespace( "utils" )
    #   desc    = packageDescription( pkgname )
    #   mess = sprintf( "This is %s %s.  %s.  Author: %s  Built: %s\n", 
    #                    pkgname, desc$Version, desc$Title, desc$Author, desc$Built )
    #   print( pkgname )
    

    #mess    = environmentName( environment(.onLoad) )
    #mess = sprintf( "Environment: '%s'\n", mess )
    #cat( mess )
    }
    
    
    
#   package.file() is a simple wrapper around system.file()
#   but gets the current package automatically
package.file <-  function( .filename )
    {
    pname   = environmentName( environment(package.file) )      #   pname is 'colortools'
    
    if( pname == "R_GlobalEnv" )
        {
        #   not in any package, so must be in development mode
        if( grepl( "^exec", .filename ) )
            folder = ".."
        else
            folder = "../inst"
        
        return( file.path( folder, .filename ) )
        }
        
        
    return( system.file(.filename,package=pname) )
    }
    
#   always returns a datframe with the right number of rows
#   and preserves all data from the 2 parts
#   nrow(out) = nrow(.df1) + nrow(.df2)    
    
rbind.super <- function( .df1, .df2 )   
    {
    #   2 trivial cases
    if( 0==nrow(.df1) ) return(.df2)
    if( 0==nrow(.df2) ) return(.df1)    
    
    if( 0<ncol(.df1)  &&  0<ncol(.df2) )
        #   the typical case
        return( merge( .df1, .df2, all=T ) )
        
    rnames1 = row.names( .df1 )
    rnames2 = row.names( .df2 )
    
    if( is.null(rnames1) )  rnames1 = rep( NA_character_, nrow(.df1) )
    if( is.null(rnames2) )  rnames2 = rep( NA_character_, nrow(.df2) )    
    
    if( 0==ncol(.df1)  &&  0==ncol(.df2) )
        #   special case
        return( data.frame( row.names=c(rnames1,rnames2) ) )
        
    #   now the 2 hard cases
    if( 0==ncol(.df1)  &&  0<ncol(.df2) )
        {
        mat = matrix( NA, nrow(.df1), ncol(.df2) )
        rownames(mat)   = rnames1
        colnames(mat)   = colnames(.df2)
        
        return( rbind(mat,.df2) )
        }
    
    if( 0<ncol(.df1)  &&  0==ncol(.df2) )
        {
        mat = matrix( NA, nrow(.df2), ncol(.df1) )
        colnames(mat)   = colnames(.df1)        
        rownames(mat)   = rnames2
        
        return( rbind(.df1,mat) )
        }    
    }
    
    
listDepth <- function (x) 
    {
    if (is.list(x)) {
        maxdepth <- 1
        for (lindex in 1:length(x)) {
            newdepth <- listDepth(x[[lindex]]) + 1
            if (newdepth > maxdepth) 
                maxdepth <- newdepth
        }
    }
    else maxdepth <- 0
    return(maxdepth)
    }    
    

hogs <- function( iPos=1 )
    {
    theNames    = ls(iPos)
    
    theList = sapply( theNames, function(x) get(x, pos=iPos )  )

    n = length(theList)    
    
    if( n == 0 ) { return(NULL) }    
    
    class1  <- function( x )    { return( class(x)[1] ) }
    
    out     = data.frame( name=theNames, size=sapply( theList, object.size ), 
                                mode=sapply(theList,mode),  class=sapply(theList,class1),  
                                stringsAsFactors=F )

    perm    = order( out$size )
    
    out     = out[ perm, ]
    
    out     = rbind( out, data.frame(name="Total:", size=sum(out$size), mode=NA, class=NA ) )
    
    #   row.names( out )    = theNames
    
    #z[ n+1 ] = sum(z)
    
    #names(z)[n+1] = "Total:"
    
    return( out )
    }
    
removeFunctions  <-  function( iPos=1, .exceptions=c("removeFunctions","hogs") )
    {
    theHogs = hogs()
    
    if( is.null(theHogs) )  return(FALSE)
    
    #   print( theHogs )
     
    df_sub  = subset( theHogs, mode=='function' )
    
    df_sub  = df_sub[ order(df_sub$name), ]
    
    #   print( df_sub )
    
    idx = match( .exceptions, df_sub$name )
    if( 0 < length(idx) )
        #   do not remove these
        df_sub = df_sub[ -idx, ]
        
    #	print( df_sub )

    n   = nrow( df_sub )
    
    if( n == 0 )
        {
        cat( "No functions to remove !\n" )
        return(FALSE)
        }
        
    mess    = sprintf( "About to remove %d functions...\n", n )
    cat( mess )
    
    print( df_sub$name )
    
    keydown <- function(key) { return( key )}

    
    key = readline( prompt="Proceed with removal ?  [Y for Yes]" )
    
    #   print(key)
    
    ok  = toupper(key) == "Y"
    
    if( ! ok )  return(FALSE)
    
    rm( list=df_sub$name, pos=iPos )
    
    return(TRUE)
    }
    
            
#   .pattern    a character vector of patterns
#   .string     a vector of strings
#
#   return value: a matrix of logicals
#   a row for each pattern and a column for each string 
#   value of each entry is whether the corresponding string matches the corresponding pattern            
multiPatternMatch <- function( .pattern, .string, .ignore=FALSE )
    {
    out = matrix( FALSE, length(.pattern), length(.string) )
    
    for( i in 1:length(.pattern) )
        out[i, ]    = grepl( .pattern[i], .string, ignore.case=.ignore )
    
    rownames(out)   = .pattern
    colnames(out)   = .string
    
    return(out)
    }
    
#   returns a logical - does every string match exactly one pattern ?    
allDistinctMatches  <- function( .pattern, .string, .ignore=FALSE )
    {
    mask    = multiPatternMatch( .pattern, .string, .ignore )
    
    return( all( colSums(mask) == 1 )  &&  max(rowSums(mask)) == 1 )
    }
            
    
compileText  <-  function( .text )
    {
    return( eval( parse( text=.text ) ) )    
    }
    
#   .y  a numerical vector - thought of as a spectrum    
#   .interval   integer vector with 2 values - giving the blending interval for lo and hi
#
#   returns:    a matrix with 2 columns:  y.lo, y.hi
#               where .y = y.lo + y.hi
splitSpectrum <- function( .y, .interval, adj=0.5 )    
    {
    .interval   = as.integer( .interval )
    
    n   = length(.y)
    ok  = all( 1 <= .interval  &  .interval <= n )
    if( ! ok )  return(NULL)    
    
    i1  = .interval[1]
    i2  = .interval[2]
    
    m   =  i2 - i1
    if( m < 2 )   return(NULL)

    s   = 1:(m-1) / m
    
    ramp.lo = (1-s)*.y[i1]
    ramp.hi = s*.y[i2]
    bridge  = ramp.lo + ramp.hi
        
    residual    = .y[ (i1+1):(i2-1) ] - bridge
    
    y.lo    = numeric(n)
    y.lo[1:i1]  = .y[1:i1]
    y.lo[ (i1+1):(i2-1) ]   = ramp.lo + adj*(residual)
    
    y.hi    = numeric(n)
    y.hi[i2:n]  = .y[i2:n]
    y.hi[ (i1+1):(i2-1) ]   = ramp.hi + (1-adj)*(residual)
    
    out = cbind( y.lo, y.hi )
    
    colnames(out)   = c( "y.lo", "y.hi" )
    
    return( out )
    }
    
    
#   .dim        dimensions of a matrix
#   .start      starting coordinates of a matrix entry

#   returns
#       a matrix with 2 columns, and prod(.dim) rows
#       the first row is .start, following by ring around .start, etc.
#       the rows form a suitable search pattern, starting at .start    
    
ringOrder <- function( .dim, .start )
    {
    stopifnot( length(.dim)==2  &&  length(.start)==2 )
    stopifnot( all( 1 <= .start  &  .start <= .dim ) )
    
    mat         = expand.grid( 1:.dim[1], 1:.dim[2] )
    
    #mat.diff    = abs( mat - matrix( .start, nrow(mat), 2, byrow=T ) )
    
    #   mat.diff    = sweep( mat, 2, .start, "-" )

    col.max = pmax( abs(mat[ ,1] - .start[1]),  abs(mat[ ,2] - .start[2]) )     # fastest by far
    
    mat     = cbind( mat, col.max )
    
    perm    = order( col.max )
    
    return( mat[perm, ] )
    }
    
time.ringOrder <- function( .dim=c(32,32), .reps=10 )
    {
    out = numeric( .reps )
    
    for( i in 1:.reps )
        {
        start   = round( runif(2,1,min(.dim)) )
        time_start  = as.double( Sys.time() )
        ringOrder( .dim, start )
        out[i]  = as.double( Sys.time() ) - time_start
        }
    
    return( out )
    }
    
    
#   .x          vector of numbers
#   .range      min and max    
CLAMP   <- function( .x, .range )
    {
    out = .x
    
    out[ .x < .range[1] ] = .range[1]
    
    out[ .range[2] < .x ] = .range[2]
    
    return( out )
    }
    
    
 