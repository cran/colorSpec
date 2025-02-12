


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
        log_level( ERROR, ".digits=%g is invalid\n", .digits )
        return(NULL)
        }

    n   = 10^.digits

    isum    = round( n * sum(.x) )

    if( isum != n )
        {
        log_level( ERROR, "sum(.x) = %g is not accurate to %g fractional digits\n",
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
        log_level( ERROR, "abs(delta) = %g is too large.  This should not happen.", abs(delta) )
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
    
    

    
findRunsTRUE <- function( mask, periodic=FALSE )
    {
    #   put sentinels on either end, to make things far simpler
    dif = diff( c(FALSE,mask,FALSE) )
    
    start   = which( dif ==  1 )
    stop    = which( dif == -1 )
    
    if( length(start) != length(stop) )
        {
        log_level( FATAL, "Internal error.  %d != %d", length(start), length(stop) )
        return(NULL)
        }
        
    stop    = stop - 1L
    
    if( periodic  &&  2<=length(start) )
        {
        m   = length(start)
        
        if( start[1]==1  &&  stop[m]==length(mask) )
            {
            #   merge first and last
            start[1]    = start[m]
            start   = start[ 1:(m-1) ]
            stop    = stop[ 1:(m-1) ]
            }
        }
        
    return( cbind( start=start, stop=stop ) )
    }
    
    
#   mask    logical vector, mostly with TRUEs.  No NAs allowed.
#   return value:   Nx2 matrix with columns
#       start   start position of a run of FALSEs
#       stop    stop position of that run of FALSEs
#   the number rows == number of runs in mask
findRunsFALSE <- function( mask )
    {
    #   put sentinels on either end, to make things far simpler
    dif = diff( c(TRUE,mask,TRUE) )
    
    start   = which( dif == -1 )
    stop    = which( dif ==  1 )
    
    if( length(start) != length(stop) )
        {
        log_level( FATAL, "Internal error.  %d != %d", length(start), length(stop) )
        return(NULL)
        }
        
    return( cbind( start=start, stop=stop-1L ) )
    }
    
        
    
    
crossproduct  <-  function( v, w )
    {
    out = c( v[2]*w[3] - v[3]*w[2], -v[1]*w[3] + v[3]*w[1], v[1]*w[2] - v[2]*w[1] )
    
    return( out )
    }
        
    
    
#   W       Nx3 matrix with generators in the rows.  Some rows may be all 0s, and these are omitted from Wcond.
#   tol     tolerance for collinear rows.  
#           This is distance between unitized vectors, treated one dimension at a time, thus L^inf-norm.
#           Also, any row with Linf-norm < tol will be set to all 0s, and thus omitted.
#
#   The purpose of this function is to make a simpler 'condensed' matrix in which
#       *) rows that are all 0s are omitted.
#       *) identify rows that are positive multiples of each other and replace such a group by their sum
#
#   If two rows are negative multiples of each other (opposite directions) then W is invalid and the function returns NULL.

#   returns a list with
#       W           the original generators
#       Wcond       Gx3 matrix with G rows, a row for each group of generators that are multiples of each other
#       group       a list of length G.  group[[k]] is an integer vector of the indexes of the rows of W.
#                   usually this is a vector of length 1.  Rows of W that are all 0s do not appear here.
#       spread      numeric vector of length G.  spread[k] is the spread of the rows of Wunit in group k.
#       groupidx    integer vector of length N.  groupidx[i] is the index of the group to which this row belongs,
#                   or NA_integer_ if the row is all 0s.
#
#   In case of ERROR the function returns NULL.
#
#   These groups are based on the given generators unitized, i.e. the point on the unit sphere.  
#

    
condenseGenerators <- function( W, tol )
    {
    #   unitize all the rows of W.
    #   any 0 rows become NaNs.
    W2  = .rowSums( W*W, nrow(W), ncol(W) )      
    
    zeromask    = (W2 == 0)
        
    if( all(zeromask) )
        {
        log_level( ERROR, "W is invalid; all rows are 0." )
        return(NULL)
        }
    
    #   in computing Wunit, the vector sqrt(W2) is replicated to all columns of W
    #   rows of all 0s in W become rows of all NaN in Wunit
    Wunit   = W / sqrt(W2)
    
    #   print( Wunit ) 
    
    res = findRowClusters( Wunit, tol=tol, projective=TRUE )
    if( is.null(res) )  return(NULL)
    
    n   = nrow(W)

    out = list()        
        
    #   res$cluster has the non-trivial groups        
    if( length(res$cluster)==0  &&  all( !zeromask ) )      # all(is.finite(res$clusteridx))  && 
        {
        #   no condensation, so easy
        out$Wcond       = W
        out$groupidx    = 1:n
        out$group       = as.list( 1:n )
        return( out )
        }
        
    groupidx    = res$clusteridx
    
    #  this line forces 0 rows to be omitted from Wcond
    groupidx[ zeromask ]   = NA_integer_
        
    if( 0 < length(res$cluster) )
        {
        #   put the non-trivial clusters in row order
        perm        = order( base::sapply( res$cluster, function(x) {x[1]} ) )
        res$cluster = res$cluster[perm]     
        
        #print( res$cluster )
        
        #   reassign groupidx[]
        for( k in 1:length(res$cluster) )
            groupidx[ res$cluster[[k]] ]    = k
            
        #print( groupidx )
        }
        
    #   count the number of trivial and non-trivial groups
    groups  = sum( groupidx==0, na.rm=TRUE )  +  length(res$cluster)
    
    Wcond   = matrix( NA_real_, groups, 3 )
    

    newidx    = groupidx
    
    visited = logical( length(res$cluster) )
    k   = 0     #   k is the new group index
    for( i in 1:n )
        {
        if( W2[i] == 0 )    next    # ignore 0 rows in W, so they belong to no group
        
        if( groupidx[i] == 0 )
            {
            #   trivial group
            k   = k + 1
            Wcond[k, ]  = W[i, ]
            newidx[i] = k            
            }
        else if( ! visited[ groupidx[i] ] )
            {
            group   =  res$cluster[[ groupidx[i] ]]
            Wgroup  = W[group, ] 
            # print(Wgroup)    #   all rows are multiples of each other
            
            #   check that all generators have consistent direction
            test    = Wgroup %*% Wgroup[1, ]
            if( any( test < 0 ) )
                {
                log_level( ERROR, "Generators are invalid. One generator is a negative multiple of another one." )
                return(NULL)
                }
                
            k   = k + 1
            Wcond[k, ]  = colSums( Wgroup )
            visited[ groupidx[i] ]  = TRUE
            newidx[ group ] = k             
            }
        }
        
    groupidx    = newidx
    
    #   now compute group from groupidx
    group   = vector( groups, mode='list' )
    spread  = numeric( groups )
    
    for( k in 1:length(group) )
        {
        group[[k]]  = which( groupidx == k )
    
        if( 1 < length(group[[k]]) )
            {
            edge        = bbedge( Wunit[ group[[k]], ] )
            spread[k]   = sqrt( sum( edge^2 ) )
            
            #   compare with un-unitized
            #   sp  = sqrt( sum( bbedge( W[ group[[k]], ] )^2 ) ) ; print( sp )
            #print( spread[k] )            
            #print( edge )
            #print( Wunit[ group[[k]], ] )
            }
        }
        
    #   assign names to Wcond
    colnames(Wcond) = colnames(W)
    
    singleton   = (sapply( group, length ) == 1)
    rownames(Wcond)[singleton]  = rownames(W)[ as.integer(group[singleton]) ]
    
    multiple    = group[ !singleton ]
    #first   = sapply( multiple, function(y) { y[1] } )
    #last    = sapply( multiple, function(y) { y[length(y)] } )
    #count   = sapply( multiple, function(y) { length(y) } )
    #rownames(Wcond)[!singleton] = sprintf( "%s...%s (%d)", rownames(W)[ first ],  rownames(W)[ last ], count )
    rownames(Wcond)[!singleton] = sapply( multiple, function(y) { sprintf( "%s...%s (%d)", rownames(W)[ y[1] ],  rownames(W)[ y[length(y)] ], length(y) ) } )   #; print(namevec)
    
    out$W           = W 
    out$Wcond       = Wcond
    out$group       = group
    out$spread      = spread
    out$groupidx    = groupidx
    
    return( invisible(out) )
    }
    
    
    
    
#   p   a point in the n-cube, which we can think of as a transmittance spectrum
#
#   returns the min number of transitions between 0 and 1 necessary to achieve such a p, including interpolation
#   it always returns an even integer    
counttransitions <- function( p )    
    {
    n   = length(p)
    if( n == 0 )    return(0L)
    
    #   find all runs of points in the interior of [0,1], in a periodic way
    interior    = 0<p  &  p<1
    if( all(interior) )
        #   special case
        return( 2 * floor( (n+1)/2 ) )

    mat     = findRunsTRUE( interior, periodic=TRUE )
    
    transitions = 0
    if( 0 < nrow(mat) )
        {
        inext   = c(2:n,1)
        iprev   = c(n,1:(n-1))
        
        for( i in 1:nrow(mat) )
            {
            start   = mat[i,1]
            stop    = mat[i,2]
            m       = stop - start + 1
            
            if( m < 0 ) m = m + n
                
            same    = p[ iprev[start] ] == p[ inext[stop] ]
            inc     = ifelse( same, 2*floor( (m+1)/2 ), 2*floor( m/2 ) )
            transitions = transitions + inc
            }
        }
        
    #   now remove all the interior coordinates
    ppure   = p[ ! interior ]
    
    #   add usual 0-1 transitions
    transitions = transitions  +  sum( diff(ppure) != 0 )
    
    #   add wrap-around transition, if present
    if( ppure[1] != ppure[ length(ppure) ] )    transitions = transitions + 1

    return( as.integer(transitions) )
    }
        
    
    
frame3x2fun <- function( normal )
    {
    ok  = is.numeric(normal)  &&  length(normal)==3  &&  0<sum(abs(normal))
    if( ! ok )
        {
        log_level( ERROR, "argument normal is invalid." )
        return(NULL)
        }
        
    out = base::svd( normal, nu=3 )$u[ , 2:3 ]     #; print( out )
    
    test    = crossproduct( out[ ,1], out[ ,2] )
    
    if( sum(test*normal) < 0 )
        {
        #   swap columns
        out = out[ , 2L:1L ] #; cat( 'frame3x2fun().  columns swapped !\n' )
        }
        
    if( isaxis( -out[ ,1] ) )
        {
        out[ ,1]    = -out[ ,1]
        out         = out[ , 2L:1L ]    # swap
        }
    else if( isaxis( -out[ ,2] ) )
        {
        out[ ,2]    = -out[ ,2]    
        out         = out[ , 2L:1L ]    # swap
        }
        
    if( FALSE )
        {
        #   multiplication check
        test    = normal %*% out
        if( 1.e-14 < max(abs(test)) )
            {
            log_level( ERROR, "frame3x2fun() failed orthogonal test = %g!", max(abs(test))  )
            }
        
        test    = t(out) %*% out  -  diag(2)
        if( 1.e-14 < max(abs(test)) )
            {
            log_level( ERROR, "frame3x2fun() failed product test = %g!", max(abs(test))  )
            }
        }
        
    return(out)
    }
    
    
#   idx     a vector of integers in 1:m, usually unique    
#   m       the upper limit for values in idx
#
#   returns TRUE iff the values in idx are contiguous, including period wraparound
contiguous  <-  function( idx, m )
    {
    if( length(idx) <= 1 )  return(TRUE)
    
    mask        = logical(m)
    mask[idx]   = TRUE
    
    mat = findRunsTRUE( mask, periodic=TRUE )
    
    out = nrow(mat)==1
    
    #   if( ! out ) print(idx)
    
    return( out )
    }

#   A           a matrix, possibly with NAs or NaNs
#   tol         difference tolerance, used to 'collapse' one column at a time
#   projective  if TRUE, then 2 rows that differ only in sign are considered the same.
#               Also, any row with L-inf norm less than tol is set to 0; probably a separate tolerance should be used for this.
#
#   returns a list with
#       clusteridx  an integer vector with length(cluster) = nrow(A)
#                   an index into cluster.  0 means this row is a trivial singleton cluster (most common).
#       cluster     list of integer vectors, the indexes of the rows of A forming the non-trivial clusters
#
#   a row with an NA is always in its own singleton cluster
#

findRowClusters <- function( A, tol, projective )
    {
    if( ! (ncol(A) %in% c(2,3) ) )
        {
        log_level( FATAL, "ncol(A)=%d is invalid.", ncol(A) )
        return(NULL)
        }    
    
    if( projective )
        {
        #   set near 0 entries to exactly 0.  There should probably be a separate tolerance for this.
        A[ abs(A) < tol ]   = 0
        
        #   extract the first non-zero entry in each row
        first   = apply( A, 1, function(r) { r[ which(r!=0)[1] ] } )
        first[ ! is.finite(first) ] = 0     # change NAs to 0
        A   = sign(first) * A               # vector sign(first) is auto-replicated to all columns
        }

    #   cluster each column
    Acollapsed  = matrix( NA_real_, nrow(A), ncol(A) )
    for( j in 1:ncol(A) )
        {
        Acollapsed[ ,j]  = collapseClusters1D( A[ ,j], tol=tol )
        }

    #   order the columns lexicographically.  A row with an NA will be put last.
    if( ncol(A) == 3 )
        #   most common
        perm    = base::order( Acollapsed[ ,1], Acollapsed[ ,2], Acollapsed[ ,3] )
    else if( ncol(A) == 2 )
        perm    = base::order( Acollapsed[ ,1], Acollapsed[ ,2] )  # ; print(perm)


    Asorted = Acollapsed[ perm, ]    #; print( Asorted )
    
    Adiff   = diff( Asorted )                                   #; print(Adiff)  # one fewer rows than A
    jump    = .rowSums( abs(Adiff), nrow(Adiff), ncol(Adiff) )  #; print( sort(jump) )
    
    #   now look for positive jumps
    jump[ ! is.finite(jump) ] = 1
    mask    = 0 < jump      # tol < jump
    runidx  = findRunsFALSE( mask )
    
    out = list()      
    out$clusteridx  = integer( nrow(A) )    
    out$cluster     = vector( nrow(runidx), mode='list' )
    out$spread      = numeric( nrow(runidx) )
    
    if( 0 < nrow(runidx) )
        {
        # perminv = order( perm )
        #   sort by run size, in decreasing order    
        #cat( "------------\n" )
        #print( runidx )
        perm2    = order( as.integer( diff( t(runidx) ) ), decreasing=TRUE )    #  ; print( as.integer( diff( t(runidx) ) ) )
        runidx  = runidx[ perm2,  ,drop=FALSE]
        #print( runidx )

        for( k in 1:nrow(runidx) )
            {
            run     = runidx[k,1]:(runidx[k,2]+1)   # add back the +1  !
            clust   = perm[ run ]  
            out$clusteridx[ clust ] = k
            out$cluster[[k]]        = clust
            
            #   set spread to the L2-length of the diagonal of the bounding box of the points
            edgevec         = bbedge( A[clust, ] )
            out$spread[k]   = sqrt( sum(edgevec^2) )
            }
        }
        
    return(out)
    }
    
    
#   mat     an NxM matrix, with points in R^M in the rows.
#
#   returns an M-vector, with the M lengths of the edges of the bounding box around these N point.
#   i'th entry = total range of i'th column of mat    
#
bbedge <- function( mat )
    {
    as.numeric( diff( apply( mat, 2, range, na.rm=T ) ) )   
    }
    

    
#   vec     numeric vector, possibly with NAs
#   tol     numerical tolerance for cluster separation
#
#   find gaps larger than tol, and change the values of each remaining 'cluster' to
#       *) the mean of that 'cluster'
#       *) but if the cluster contains a single integer, then change to that integer
#   return resulting vector of the same length    
collapseClusters1D <- function( vec, tol  )
    {
    perm        = base::order( vec )    # any NAs are put *last*
    
    vecsorted   = vec[ perm ]    #; print( Asorted )
    
    jump   = diff( vecsorted )
    
    jump[ ! is.finite(jump) ] = tol + 1
    mask    = tol < jump
    runidx  = findRunsFALSE( mask )
    
    if( nrow(runidx) == 0 )
        #   no clusters and no change
        return( vec )

    # perminv = order( perm )
    #   sort by run size, in decreasing order    
    #cat( "------------\n" )
    #print( runidx )
    #perm2    = order( as.integer( diff( t(runidx) ) ), decreasing=TRUE )    #  ; print( as.integer( diff( t(runidx) ) ) )
    #runidx  = runidx[ perm2,  ,drop=FALSE]
    #print( runidx )
    
    for( k in 1:nrow(runidx) )
        {
        run                 = runidx[k,1]:(runidx[k,2]+1)   # add back the +1  !
        cluster             = vecsorted[run]
        idx                 = which( floor(cluster) == cluster )
        if( length(idx) == 1 )
            value   = cluster[idx]
        else
            value   = mean( cluster )
            
        vec[ perm[run] ]    = value
        }
    
    return( vec )
    }
    
    
isaxis <- function( vec )
    {
    n   = length(vec)
    
    return( sum(vec==0)==n-1  &&  sum(vec==1)==1 )
    }
        
    
    
#  Computes the signed area of a polygon.
#  The area is positive if the vertices are clockwise
#  in a right-handed (y increases going up) 2D coordinate system.
#
#   Formula taken from p. 213 of

#   M. I. Shamos
#   Computational Geometry (his PhD thesis)    

areaOfPolygon <- function( iXY )
    {
    n = nrow(iXY)
    
    area = 0

    for( i in 1:n )
        {
        i_prev = ((i-2) %% n) + 1
        i_next = (i %% n) + 1

        area = area + iXY[i,2] * ( iXY[i_next,1] - iXY[i_prev,1] )
        }

    return( 0.5 * area )
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
        #do.call( log_level, arglist, envir=parent.frame(n=3) )
        #myfun   = log_level
        #environment(myfun) = parent.frame(3)

        Aname = deparse(substitute(A))

        #   notice hack with .topcall to make log_level() print name of parent function, and *NOT* prepareNxM()
        log_level( ERROR, "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>=%d). %s='%s...'",
                                    Aname, M, Nmin, Aname, mess, .topcall=sys.call(-1L) )       # or -2L ?  NO !
        return(NULL)
        }

    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )

    return( A )
    }
