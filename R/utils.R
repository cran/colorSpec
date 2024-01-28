
    
    
#   returns time in seconds, from an arbitrary origin
gettime <- function()
    {
    if( g.microbenchmark )
        return( microbenchmark::get_nanotime() * 1.e-9 )
    else
        return( as.double( base::Sys.time() ) )
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
    

#   vec     numeric vector, regarded as periodic.  no NAs
#
#   returns an Nx2 matrix with indices of the transitions, including possible first and last    
transitionMatrix <- function( vec )
    {
    n   = length(vec)
    if( n == 0 )    return( matrix( 0L, 0, 2 ) )
    
    idx = which( diff( vec ) != 0 )
    if( length(idx) == 0 )  return( matrix( 0L, 0, 2 ) )
    
    out = cbind( idx, idx+1 )
    
    if( vec[n] != vec[1] )  out = rbind( out, c(n,1) )
    
    return( out )
    }
    
#   wavelength      numeric vector of length n
#    
#   computes parameters for n bins, roughly centered at the wavelengths
#
#   returns a list with numeric vectors
#       breakvec    breaks between n bins, the length is n+1
#       stepvec     width of the bins, the length is n
#
breakandstep <- function( wavelength, method='rectangular' )
    {    
    n   = length(wavelength)
    
    out = list()  
    
    breakvec        = 0.5 * (wavelength[1:(n-1)] + wavelength[2:n])
    
    if( tolower(method) == substr("trapezoidal",1,nchar(method)) )     
        out$breakvec    = c( wavelength[1], breakvec, wavelength[n]  )   # cutoff first and last bins
    else
        out$breakvec    = c( 2*wavelength[1] - breakvec[1], breakvec, 2*wavelength[n] - breakvec[n-1] )   # extend first and last bins symmetric
    
    out$stepvec     = diff( out$breakvec )
    
    return(out)
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
    
    
#   .x      a numeric vector whose sum is 1 - accurate to .digits digits
#   .digits the number of fractional decimal digits to round.
#
#   returns .x but with components rounded to .digits digits
#           and with sum still 1
#           in case of error, returns NULL

roundAffine  <- function( .x, .digits )
    {
    ok  = (1 <= .digits)  &&  (.digits <= 10 )  &&  (round(.digits) == .digits)
    if( ! ok )
        {
        log_string( ERROR, ".digits=%g is invalid\n", .digits )
        return(NULL)
        }
    
    n   = 10^.digits
    
    isum    = round( n * sum(.x) )
    
    if( isum != n )
        {
        log_string( ERROR, "sum(.x) = %g is not accurate to %g fractional digits\n", 
                                sum(.x), .digits )
        return(NULL)
        }
    
    out     = round( n * .x )
    
    delta   = sum(out) - n #;     print( delta ) ;
    
    if( delta == 0 )
        #   easy case
        return( out / n )
        
    if( length(.x) < abs(delta) )
        {
        log_string( ERROR, "abs(delta) = %g is too large.  This should not happen.", abs(delta) )
        return(NULL)
        }
    
    #   find the delta largest values of abs(.x)
    perm    = order( abs(.x) , decreasing=TRUE )    #; print(perm)
    
    idx     = perm[ 1:abs(delta) ] #;   print( idx)
    
    out[ idx ] = out[ idx ] - sign(delta)
    
    return( out / n )
    }

#   .vec1 and .vec2     non-zero vectors of the same dimension    
#
angleBetween  <-  function( .vec1, .vec2, eps=5.e-14 )
    {
    len1    = sqrt( sum(.vec1^2) )
    len2    = sqrt( sum(.vec2^2) )     #;    print( denom )
    
    denom   = len1 * len2
    
    if( abs(denom) < eps )    return( NA_real_ )
    
    q   = sum( .vec1*.vec2 ) / denom  #; print(q)
        
    if( abs(q) < 0.99 )
        {
        #   the usual case uses acos
        out = acos(q)
        }    
    else
        {
        #   use asin instead
        .vec1   = .vec1 / len1
        .vec2   = .vec2 / len2
        
        if( q < 0 ) .vec2 = -.vec2
        
        d   = .vec1 - .vec2
        d   = sqrt( sum(d*d) )
        
        out = 2 * asin( d/2 )
        
        if( q < 0 ) out = pi - out
        }

    return(out)
    }
    
is.identity <- function( x )
    {
    if( ! is.matrix(x) )    return(FALSE)
    
    m   = nrow(x)
    if( m != ncol(x) )  return(FALSE)
    
    return( all( x==diag(m) ) )       # was identical()
    }

    
#   x   numeric vector
#   y   positive number
#
#   returns x^y, with extension for negative x to make an odd function
powodd <- function( x, y )
    {
    ok  = is.numeric(x)  &&  is.numeric(y)  &&  length(y)==1  &&  0<y  
    if( ! ok )  return(NULL)
    
    s   = sign(x)
    
    return( s * (s*x)^y )
    }
    
###########     argument processing     ##############
#
#   A   a non-empty numeric NxM matrix, or something that can be converted to be one
#
#   Nmin    the minimum allowed number of rows
#
#   returns such a matrix, or NULL in case of error
#
prepareNxM  <-  function( A, M=3, Nmin=1 )
    {
    ok  = is.numeric(A) &&  M*Nmin<=length(A)  &&  (length(dim(A))<=2)  # &&  (0<M) 
    
    ok  = ok  &&  ifelse( is.matrix(A), ncol(A)==M, ((length(A) %% M)==0)  )
    
    if( ! ok )
        {
        #print( "prepareNx3" )
        #print( sys.frames() )
        mess    = substr( paste0(as.character(A),collapse=','), 1, 10 )
        #arglist = list( ERROR, "A must be a non-empty numeric Nx3 matrix (with N>=%d). A='%s...'", mess )
        #do.call( log_string, arglist, envir=parent.frame(n=3) )
        #myfun   = log_string
        #environment(myfun) = parent.frame(3)
        
        Aname = deparse(substitute(A))        
        
        #   notice hack with 2L to make log_string() print name of parent function
        log_string( c(ERROR,2L), "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>=%d). %s='%s...'", 
                                    Aname, M, Nmin, Aname, mess )
        return(NULL)
        }
    
    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )
        
    return( A )
    }
         