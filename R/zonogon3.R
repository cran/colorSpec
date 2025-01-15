
#   'zonogon3' is an S3 class
#
#   It only looks at the boundary, and is designed for 'raising' to a zonohedron.
#   Each vertex is raised to a parabaloid.
#
#   It is stored as a list with
#       W               the original given nx2 matrix, defining n line segments in R^2  (with the other endpoint at 0).  
#                       These rows are the generators.  Must be rank 2, and no NAs.
#       center          middle gray = white/2
#       nonnegative     logical.  All entries of W are non-negative
#       zonohedron      the zonohedron obtained by 'raising' the generators to a cone
#
#   zonogon3() constructor for a zonogon object
#       W       nx2 matrix with rank 2
#       tol     relative tolerance for degenerate faces
#
#   returns: a list as above
#

zonogon3  <- function( W, tol=1.e-9 )
    {
    ok  = is.numeric(W)  &&  is.matrix(W)  &&  2<=nrow(W)  &&  ncol(W)==2
    if( ! ok )
        {
        log_level( ERROR, "argument W is not an nx2 numeric matrix, with n>=2." )
        return(NULL)
        }
        
    if( ! all( is.finite(W) ) )
        {
        log_level( ERROR, "matrix W is invalid because it has %d entries that are not finite.", sum(! is.finite(W) ) )
        return(NULL)
        }
    
        
    normal  = cbind( W[ ,2], -W[ ,1] )
    
    #   set rows close to 0 to all NAs
    # WW2 = .rowSums( W*W, nrow(W), ncol(W) )
    
    normal2 = .rowSums( normal*normal, nrow(normal), ncol(normal) )
    bad     = normal2 <= tol*tol * mean(normal2)
    if( all(bad) )
        {
        log_level( ERROR, "argument W does not have rank 2, with relative tol=%g.", tol )
        return(NULL)
        }
    
    normal[ bad, ]  = NA_real_
    
    #   now unitize
    normal  = normal / sqrt(normal2)    # sqrt(normal2) is replicated to all 2 columns

    n   = nrow(W)
    

    
    out = list()
    
    out$W   = W

    out$center  = 0.5 * .colSums( W, nrow(W), ncol(W) )
    out$nonnegative = all( 0 <= W )


    
    class(out)  = c( 'zonogon3', class(out) )    
    
    #   W3  = cbind( W, rowSums(W*W) )
    
    #   'raise' generators to a cone
    W3  = cbind( W, sqrt(rowSums(W*W)) )
    
    #   compute zonohedron, but without any clustering by setting tol2=NA
    out$zonohedron  = zonohedron( W3, tol2=NA )
        
    return( out )
    }
    
    

plot.zonogon3  <- function( x, normal=c(0,0,1), side=-1,  ... )
    {
    plotprojection( x$zonohedron, normal=normal, side=side )
    
    dat = getcenters( x )
    
    center  = dat$center + matrix( x$center, nrow(dat), 2, byrow=TRUE )
    
    points( center[ ,1], center[ ,2], col='red', pch=20 )
    
    return( invisible(TRUE) )
    }
    
#   a zonogon3 has an implicit partition into parallelograms, the bottom surface of the zonohedron    
#   return a data.frame with a row for each parallelogram and these columns
#       idx     pairs of indexes for the parallelogram
#       center  center of p-gram in 2D, in centered zonogon
getcenters.zonogon3  <-  function( x )
    {
    face    = x$zonohedron$face  
    
    out = data.frame( row.names=rownames(face) )
    
    out$idx = face$idx
    
    signvec = as.numeric( sign( -face$normal[ ,3] ) )

    out$center  = signvec * face$center[ , 1:2 ]  #+  matrix( x$center, nrow(out), 2, byrow=TRUE )

    return( out )
    }
    
    
#   x       a zonogon3, with generators counterclockwise
#   count   number of points to test    
#
#   for each test, 2-transition coefficients are generated at random,
#   which are then mapped to a point in the zonogon3.
#   This point is passed to raytrace2() to try to recover the original coefficients.
testsorted.zonogon3 <- function( x, count=1, comp=FALSE )    
    {
    n   = nrow( x$W )
    
    set.seed(0)
    
    tol = 5.e-12
    
    out = TRUE
    
    for( i in 1:count )
        {
        #   cat( "------------------\n" )
        
        p   = random2trans( n, comp=comp )     #; print( p )
        
        xy  = p %*% x$W 
        
        res = raytrace2( x$zonohedron, c(xy,0), c(0,0,1) )
        
        #print( res )
        
        delta   = max( abs(res$source - p) )
        if( tol < delta )
            {
            cat( i, "source delta: ", delta, '\n' )
            out = FALSE
            }
        }
        
    return( out )
    }
    
    
#   returns n 2-vectors in quadrant #1 sorted counter-clockwise by angle;  lengths are also random.
#   the returned  nx2 matrix is suitable for input to zonogon3().
randomsorted <- function( n )
    {
    theta = sort( runif( n, min=0, max=pi/2 ) )
    
    out = cbind( cos(theta), sin(theta) )
    
    out = runif( n, min=0.1,max=1) * out
    
    return( out )
    }
    
    
    
#   returns a random n-vector
#   2 coordinates are uniform in [0,1]
#   all coords between are 1, and rest are 0    
random2trans <- function( n, comp=FALSE )
    {
    cmat    = utils::combn( n, 2 )
    m       = ncol( cmat )
    k       = m * runif(1) + 1
    pair    = cmat[ ,k]
    
    out     = numeric( n )  # all 0s
    
    out[pair]   = runif(2)
    
    if( pair[1]+1 <= pair[2]-1 )
        out[ (pair[1]+1) : (pair[2]-1) ]    = 1
    
    #if( runif(1) < 0.5 )    
    if( comp )
        out    = 1 - out  # complement
    
    return( out )
    }
    
    
#   the most simple-minded I can think of - just choose n points on the circle !    
randomcxpoly <- function( n )    
    {
    theta   = sort( runif( n, min=-pi, max=pi ) )   ; print(theta)
    
    return( cbind( cos(theta), sin(theta) ) )
    }
    
    
#   return random generators, on the boundary of a convex cone in cyclic order
#   n           number of starting points, n >=3, all with z=1
#   between     number to put between 2 starting points, all on the plane spanned by the 2 points
#   range       random range of z, used separately for each output

#   value       an nx3 matrix, with generators in cyclic order

randomgenerators <- function( n=3, between=0, range=c(1,1) )
    {
    if( n < 3 ) return(NULL)
    
    if( length(between) == 0 )  return(NULL)
    
    if( length(between) < n )
        {
        between = rep( between, ceiling( n/length(between) ) )
        between = between[1:n]
        #   print( between )
        }
    
    out = cbind( randomcxpoly(n), 1 )
    
    total   = sum(between)
    if( 0 < total )
        {
        outsave = out
        out = matrix( NA_real_, n+total, 3 )
        
        k   = 0
        for( i in 1:n )
            {
            k   = k+1
            out[k, ]    = outsave[i, ]
            
            if( 0 < between[i] )
                {
                inext   = ifelse( i<n, i+1, 1 )
                
                w   = sort( runif( between[i] ) )
                W   = cbind( 1-w, w )   #; print(W)
                out[ (k+1):(k+between[i]), ] = W  %*%  outsave[c(i,inext), ]
                k   = k + between[i]
                }
            }
        }
    
    if( range[1] < range[2] )
        {
        s   = runif( nrow(out), min=range[1], max=range[2] )
        out = s * out       # multiplies the rows
        }
        
    return( out )
    }
    
    
#--------       UseMethod() calls           --------------#    
    
testsorted <- function( x, count=1, comp=FALSE )
    {
    UseMethod("testsorted")
    }    
    
getcenters  <-  function( x )
    {
    UseMethod("getcenters")
    }    