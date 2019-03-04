
#   'zonogon' is an S3 class
#
#   It presents the zonogon as an intersection of 'slabs'.
#   Each slab is the intersection of 2 halfplanes with the same normal vector.
#   The equation of a slab is:
#       -beta <=   <x,normal>   <= beta        (beta is always positive)
#
#   It is stored as a list with
#       W               the original given nx2 matrix, defining n line segments in R^2  (with the other endpoint at 0).  
#                       Must be rank 2, and no NAs.
#       center          middle gray = white/2
#       nonnegative     logical.  All entries of W are non-negative
#       face            a data frame with n rows.  Because of central symmetry, only 1/2 of the faces need to be stored.
#                           normal  to a face of the zonohedron, then unitized
#                           center  of the face, a segment, *after* centering the zonogon
#                           beta    the plane constant, *after* centering the zonogon.  All are positive.
#       clusteridx      an integer vector with length(cluster) = nrow(W) = the faces of the zonogon.
#                       an index into cluster list.  0 means this row is a trivial singleton cluster (most common).
#       cluster         list of integer vectors, the non-trivial clusters.
#       parallelogram   data frame with n(n-1)/2 rows
#                           idx     a pair 1<= i < j <=n
#                           center  center of parallelogram
#       vertex          a matrix with 2 columns with all vertices in the rows.
#                       These are the vertices of the original and *non-centered* zonogon.
#                       There may be "inline" vertices.

#   zonogon() constructor for a zonogon object
#       W       nx2 matrix with rank 2
#       tol     relative tolerance for degenerate faces
#
#   returns: a list as above
#
#   Issues: if normal vectors are antipodal, then the centers of faces are incorrect.  Fix this.

zonogon  <- function( W, tol=1.e-9, partition=FALSE )
    {
    ok  = is.numeric(W)  &&  is.matrix(W)  &&  2<=nrow(W)  &&  ncol(W)==2
    if( ! ok )
        {
        log.string( ERROR, "argument W is not an nx2 numeric matrix, with n>=2." )
        return(NULL)
        }
        
    if( ! all( is.finite(W) ) )
        {
        log.string( ERROR, "matrix W is invalid because it has %d entries that are not finite.", sum(! is.finite(W) ) )
        return(NULL)
        }
    
    n   = nrow(W)
    

    normal  = cbind( W[ ,2], -W[ ,1] )
    
    #   set rows close to 0 to all NAs
    # WW2 = .rowSums( W*W, nrow(W), ncol(W) )
    
    normal2 = .rowSums( normal*normal, nrow(normal), ncol(normal) )
    bad     = normal2 <= tol*tol * mean(normal2)
    if( all(bad) )
        {
        log.string( ERROR, "argument W does not have rank 2, with relative tol=%g.", tol )
        return(NULL)
        }
    
    normal[ bad, ]  = NA_real_
    
    #   now unitize
    normal  = normal / sqrt(normal2)    # sqrt(normal2) is replicated to all 2 columns

    out = list()
    
    out$W   = W
    
    out$center  = 0.5 * .colSums( W, nrow(W), ncol(W) )
    
    out$nonnegative = all( 0 <= W )
    
    out$face        = data.frame( row.names=1:n )
    #   out$face$idx    = 1:n    
    out$face$normal = normal
    
    if( TRUE )
        {
        #   cluster the unit normal vectors  [reminds me of something similar at Link]
        res = findRowClusters( normal )
        if( is.null(res) )  return(NULL)
        
        out$clusteridx  = res$clusteridx
        out$cluster     = res$cluster
        }

    #   the columns of functional are linear combinations of columns of W, and the points in the n-cube
    #                         nx2  nx2    
    functional  = tcrossprod( W, normal )   #; print( functional )      This is large:  nxn

    


    #   functional is n x n 
    #   the diagonal entries should be exactly 0, 
    #   but floating-point multiplication (after unitizing the normal) is *not* commutative, so some diagonal entries are e-17 etc. !        
    
    #dia     = cbind( 1:n, 1:n )
    delta   = functional[ cbind(1:n,1:n) ]      #; print(check)      # all these are very small, and often exactly 0
    
    epsilon = max( 10*abs(delta), 1.e-11, na.rm=TRUE )       #; print( epsilon )
    
    out$epsilon = epsilon   # used later in invertboundary()
    
    if( partition )
        {
        #   The columns of matrix vertex are vertices in the unit n-cube
        vertex  = 0.5 * (sign(functional) + 1)
        # force all diagonals to exactly 0.5.  This is not necessary 
        #   vertex[ cbind(1:n,1:n) ]    = 0.5    #; print( vertex )
        
        if( requireNamespace( 'arrangements', quietly=TRUE ) )
            idx = arrangements::combinations(n,2)       # faster
        else
            idx = t( utils::combn(n,2) )   # matrix of pairs.  slower
        
        center  = matrix( NA_real_, nrow(idx), 2 )
        #base    = matrix( NA_real_, nrow(idx), n )
        for( k in 1:nrow(idx) )
            {
            i   = idx[k,1]
            j   = idx[k,2]   # so i <  j
            
            #cat( "---------",  i,j, "-----------------\n" )
            #   calculate the 'basepoint' of parallelogram i,j
            #   this is the starting vertex of edge #i that goes counterclockwise around the zonogon
            basek   = vertex[ ,i]        
            
            if( functional[j,i] < 0 )
                {
                #   expansion is from the 'clockwise' side of face j,
                #   so invert the cube vertex to get the clockwise edge instead
                #cat( 'inverting...\n' )
                basek   = 1 - basek
                }
            basek[ j:n ] = 0         # trim away the end, since this is an expansion of the subzonogon generated by 1:(j-1)            

            basek[ c(i,j) ] = 0.5    # for the center of parallelogram
            
            #cat( "basek=", basek, '\n' )

            #   linearly map from n-cube to the plane
            center[k, ] =  basek %*% out$W
            #base[k, ]   =  basek
            }
            
        out$parallelogram           = data.frame( row.names=1:nrow(idx) )
        out$parallelogram$idx       = idx
        out$parallelogram$center    = center
        #out$parallelogram$base      = base
        }
    
            
    #   now find all entries within epsilon of 0, and change to 0, +1, or -1
    idxsmall    = which( abs(functional) <= epsilon, arr.ind=TRUE )   #; print( idxsmall )
    
    #   change small entries above diagonal to +1, below to -1, and on the diagonal to 0 (near 0 already)
    functional[idxsmall]    = sign( idxsmall[ ,2] - idxsmall[ ,1] )

   
    #   calculate the n face centers of the unit cube,
    #   but then all translated by -0.5 and therefore centered
    #   sign(0) = 0 exactly, so 2 (or more) of the 'spectrum' coords are 0
    face_center = sign(functional)  #;    print(face_center)  

    face_center = 0.5 * face_center
    
    out$face$center = crossprod( face_center, W )
    
    #   calculate the n plane constants beta
    #   these plane constants are for the centered zonogon
    out$face$beta   = .rowSums( normal * out$face$center, nrow(normal), ncol(normal) )
    
    #   since W has rank 2, 0 is in the interior of the _centered_ zonogon,
    #   and this implies that all n beta's are positive, except the NAs.
    #   verify this
    betamin = min( out$face$beta, na.rm=TRUE )
    if( betamin <= 0 )
        {
        log.string( FATAL, "Internal Error.  min(beta)=%g <= 0.", tol, betamin )
        return(NULL)
        }
        
    #   find the vertices in counter-clockwise order
    #   re-order the faces by angle of the normal vector     
               
    vertex  = face_center
    vertex[ vertex == 0 ] = -0.5  # back from center to beginning of edge.    


    # vertex  = vertex + 0.5         # translate back to [0,1] cube
    #vertex[ vertex == 0.5 ] = 0             # back from center to beginning of edge.  The columns of vertex are now vertices in the n-cube.
    #   print( face_center )
    
    #   map from vertices of centered n-cube to the vertices of centered zonogon in the plane
    vertex = crossprod( vertex, W )
    

    #   remove any rows that are NA
    mask    = .rowSums( vertex, nrow(vertex), ncol(vertex) )
    vertex  = vertex[ is.finite(mask), ]
    
    vertex  = rbind( vertex, -vertex )   #  reflect about the origin
        
    theta   = atan2( vertex[ ,2], vertex[ ,1] )     #; print( duplicated(theta) )
    #out$face$theta  = theta 

    perm    = base::order( theta )     
    
    vertex  = vertex[perm,  ]    # put vertices in proper order 

    out$vertex = vertex + matrix( out$center, nrow(vertex), 2, byrow=TRUE )     # back to original

    #   these are only half the vertices, now add the rest by symmetry around out$center
    #   out$vertex  = rbind( out$vertex, matrix( 2*out$center, nrow(out$vertex), 2, byrow=TRUE ) - out$vertex )
        
    #   out$cyclicposition  = cyclicposition( W, idx, normal )   TODO:  needs more work on tolerances


    class(out)  = c( 'zonogon', class(out) )
        
    return(out)
    }
    
    
##----------        zonogon methods    -------------##
    
#   x   a zonogon object
#   p   an Nx3 matrix, etc.
#
#   value   a dataframe with columns
#           inside      logical
#           distance    numeric signed distance to zonogon boundary, non-positive means inside.  
#                       If positive then this is only approximate.
inside.zonogon <- function( x, p )
    {
    p   = prepareNxM( p, 2 )
    if( is.null(p) )    return(NULL)
    
    #   translate p to the centered zonogon
    gcentered   = p - matrix( x$center, nrow(p), 2, byrow=TRUE ) #; print(gcentered)
    
    hg  = tcrossprod( x$face$normal, gcentered )    #; print( str(hg) )
    
    distance    = abs(hg) - matrix( x$face$beta, nrow(hg), ncol(hg) )
    
    distance    = base::apply( distance, 2, function(z) {suppressWarnings( max(z,na.rm=TRUE) ) } )  #;   print(distance)
    
    distance[ ! is.finite(distance) ]   = NA_real_
    
    if( x$nonnegative )
        {
        #   special override for black
        black   = apply( p, 1, function(v) { isTRUE(all(v==0)) } )  #; print(black)
        if( any(black) )
            distance[black] = 0
            
        #   special override for white.  Fortunately multiplication by 0.5 and 2 preserves all precision.
        white   = apply( p, 1, function(v) { isTRUE(all(v==2*x$center)) } )  #; print(white)
        if( any(white) )
            distance[white] = 0
        }
        
    rnames  = rownames(p)
    if( is.null(rnames) )   rnames = 1:nrow(p)
        
    out = data.frame( row.names=rnames )
    
    out$p           = p
    out$inside      = distance <= 0
    out$distance    = distance

    return(out)
    }
    
    
#   x           a zonogon object
#   base        a numeric vector of length 2, the basepoint of all the rays
#               base must be in the interior of x,
#               or if x is non-negative, base can also be the black or white point on the boundary(x)
#   direction   an Nx2 matrix with directions in the rows
#
#   value   a dataframe with columns
#           base        given basepoint of the ray (all the same)
#           direction   given direction of the ray
#           idx         of the face (edge) where ray exits the zonohedron
#           tmax        ray parameter of intersection with face
#           boundary    point of intersection with the intersection with face
#           alpha       coordinate of boundary in the face coords.  in [0,1]

raytrace.zonogon <- function( x, base, direction )
    {
    ok  = is.numeric(base)  &&  length(base)==2  &&  all( is.finite(base) )
    if( ! ok )
        {
        log.string( ERROR, "base is invalid. It must be a numeric vector of length 2, and all entries finite." )
        return(NULL)
        }
        
    direction   = prepareNxM( direction, 2 )
    if( is.null(direction) )    return(NULL)
    

    base    = as.numeric(base)
    
    #   translate base to the centered zonohedron
    gcentered   = base - x$center  #; print(gcentered)
    
    #   test whether base is black or white point, no tolerance here
    blackwhite  = ifelse( x$nonnegative, all(gcentered == x$center) || all(gcentered == -x$center), FALSE )  #; print( blackwhite )

    hg  = as.numeric( x$face$normal  %*%  gcentered )      #; print( str(hg) )
    
    distance    = abs(hg) - x$face$beta 
    
    rtol    = 1.e-14    
    if( blackwhite )
        {
        #   change distance that is very close to 0
        #   print( abs(distance) )
        bound    = abs(distance) < rtol * mean(abs(x$center))     #; print(bound)
        distance[bound]  = 0                                     #; print(distance)
        distance    = max( distance, na.rm=TRUE )   #;   print(distance)
        ok  = distance <= 0
        }
    else
        {
        distance    = max( distance, na.rm=TRUE )   #;   print(distance)
        ok  = distance < 0
        }
        
    if( ! ok )
        {
        log.string( ERROR, "point base=(%g,%g) is not in the interior of the zonogon.  distance=%g >= 0.", 
                                base[1], base[2], distance )
        return(NULL)
        }
        
    n   = nrow(direction)
    

    tmax        = rep(NA_real_,n)
    idx         = rep(NA_integer_,n)
    sign        = rep(NA_integer_,n)    
    boundary    = matrix(NA_real_,n,2)
    alpha       = rep(NA_real_,n)
    timetrace   = rep(NA_real_,n)
    
    for( k in 1:n )
        {
        time_start  = gettime()
        
        v   = direction[k, ]
        
        if( any( is.na(v) ) )   next
        
        if( sum(v*v) == 0 ) next    # 0-vector
        
        hv      = x$face$normal  %*%  v
        
        numerator   = x$face$beta - sign(hv)*hg
        tvec    = numerator / abs(hv) 
            
        if( blackwhite )
            {
            #   find those planes that contain black or white
            bound   = abs(numerator) < rtol * sum(abs(x$center))    #; print( sum(is.na(boundary)) )
            
            bound[ is.na(bound) ] = FALSE
            #numerator[boundary]  = 0

            tvec[bound]  = ifelse( 0 < hv[bound], 0, Inf )    # 0 here will later make the corresponding boundary = NA 
            }

        tvec[ ! is.finite(tvec) ]   = Inf
        
        j   = which.min( tvec )     # this ignores ties and Infs
        
        tmax[k] = tvec[j]       # tmax[k] is not negative, and might be 0 if base is black or white

        idx[k]  = j  # x$face$idx[j] 
        
        if( tmax[k] <= 0 )  next    # failed to intersect properly
        
            
        optcentered     = gcentered  +  tmax[k] * v     #; print( optcentered )  
        boundary[k, ]   = optcentered + x$center        
        
        cidx  =  x$clusteridx[j]    
        if( cidx != 0 )
            {
            cluster = x$cluster[[ cidx ]]   # always more than 1
            cat( sprintf( "Searching in cluster %d, with %d edges.\n", cidx, length(cluster) ) )
            }
        else
            cluster = j     # a singleton cluster

        #   print( cluster )
        
        for( j in cluster )
            {
            #   test each edge in the cluster, to find alpha in [0,1]
            theSign         = sign(hv[j])
            facecenter      = theSign * x$face$center[j, ]              #; print( facecenter )            
            edge            = t( x$W[ j, ,drop=F] )    #; print( edge )     # 1 edge of the polygon side, as 2x1 matrix
            M               = cbind( edge, x$face$normal[j, ] )             #; print(M) # M is 2x2
            y               = base::solve( M, optcentered - facecenter )    #; print(y)
            #   alpha[k]        = pmin( pmax( y[1] + 0.5, 0 ), 1 )          # translate from [-0.5,0.5] to [0,1] and clamp it
            if( abs(y[1]) <= 0.5 )
                {
                idx[k]      = j             # override above
                sign[k]     = theSign
                alpha[k]    = y[1] + 0.5    # alpha is in [0,1]
                break
                }
            }
            
            
        timetrace[k]   = gettime() - time_start
        }
    
    out = data.frame( row.names=1:n )
    
    out$base        = matrix( base, n, 2, byrow=TRUE )  # replicate base to all rows
    out$direction   = direction
    out$idx         = idx
    out$sign        = sign    
    out$tmax        = tmax
    out$boundary    = boundary
    out$alpha       = alpha
    out$timetrace   = timetrace
    
    #if( blackwhite )    tmax[ tmax==0 ] = NA_real_      # so boundary becomes NA too !
    #out$boundary     = out$base  +  tmax * direction   # tmax is replicated to 3 columns

    cnames  = colnames(base)
    if( is.null(cnames) )   cnames = colnames(direction)    
        
    colnames(out$boundary)   = cnames
    
    return( out )
    }
    

#   given a point on the boundary, return a point in the unit n-cube that maps to it, under W  
#       x       the zonogon
#       data    a data.frame as returned from raytrace.zonogon().  Columns used are:
#               boundary    coordinates of the point on the boundary
#               idx         index of the face
#               sign        sign of the face = +-1
#               alpha       coordinate along the face, in [0,1]
#       tol     tolerance for verification.  Set to NULL for RELEASE.
#
#   returns:
#       an mxn matrix, where m=nrow(data)  n=number of rows in x$W 
#       each row of the matrix is a point in the n-cube
#       that maps to the point on the zonogon

invertboundary.zonogon <- function( x, data, tol=NULL )     #1.e-9 )
    {
    ok  = is.data.frame(data)  &&  0<nrow(data)
    ok  = ok  &&  !is.null(data$boundary)  &&   !is.null(data$idx)  &&  !is.null(data$sign)  &&  !is.null(data$alpha)
    if( ! ok )
        {
        log.string( ERROR, "argument data is invalid." )
        return(NULL)
        }
    
    n   = nrow(x$W)
    m   = nrow(data)
    
    #   Recompute a submatrix of functional in the constructor, which was n x n
    #   This functional is n x m    
    #   This matrix could have been stored in the zonogon list,
    #   but recomputing it saves memory.  A trade-off.
    #                        n x 2   m x 2
    functional  = tcrossprod( x$W, x$face$normal[data$idx, ,drop=F] )   #; print( functional ) #  n x m

    #check   = cbind( data$idx, 1:m )    
    #check   = functional[ check ]   ; print(check)      # all these are very small, and often 0
    
    #   now find all entries within epsilon of 0, and change to 0, +1, or -1
    idxsmall    = which( abs(functional) <= x$epsilon, arr.ind=TRUE )   #; print( idxsmall )
    
    functional[idxsmall]    = sign( data$idx[ idxsmall[ ,2] ]  -  idxsmall[ ,1] )

    #   calculate the m face centers of the unit cube,
    #   but then all translated by -0.5 and therefore in the centered cube
    #   sign(0) = 0 exactly 
    face_center = matrix( 0.5*data$sign, n, m, byrow=T ) * sign(functional)  #;    print(face_center)  

    vertex      = face_center + 0.5         # translate back to [0,1] cube
    #vertex[ vertex == 0.5 ] = 0             # back from center to beginning of edge.  The columns of vertex are now vertices in the n-cube.
    
    #   now move distance alpha along the proper edge of the n-cube
    vertex[ cbind( data$idx, 1:m ) ]    = data$alpha
    
    if( is.numeric(tol) && is.finite(tol) && 0<tol )
        {
        #   map from vertices of n-cube to the vertices of zonogon in the plane
        bpoint  = crossprod( vertex, x$W )
        delta   = abs(bpoint - data$boundary)
        
        if( tol <= max(delta,na.rm=TRUE) )
            log.string( WARN, "Verification failed. max(delta)=%g > %g.", max(delta), tol )        
        }
        
    vertex  = t(vertex)
        
    colnames(vertex)    = rownames( x$W )
    rownames(vertex)    = as.character( data$idx )
        
    return( vertex )
    }
    
    
    
    

#   section() compute intersection of zonogon and line(s)
#    
#   x           a zonogon object
#   normal      a non-zero numeric vector of length 2, the normal of all the lines
#   beta        a vector of line-constants.  The equation of line k is: <x,normal> = beta[k]
#
#   value   a list of length = length(beta).  Each list has the items:
#           beta        given line constant of the line
#           section     Mx2 matrix of points on the section
#                       usually M=2 (lines that intersect interior - the usual case) or M=1 for support lines
#                       If there is no intersection, then M=0 rows.

section.zonogon <- function( x, normal, beta )
    {    
    ok  = is.numeric(normal)  &&  length(normal)==2  &&  all( is.finite(normal) )
    if( ! ok )
        {
        log.string( ERROR, "normal is invalid. It must be a numeric vector of length 2, and all entries finite." )
        return(NULL)
        }
        
    ok  = is.numeric(beta)  &&  0<length(beta)   &&  all( is.finite(beta) )
    if( ! ok )
        {
        log.string( ERROR, "beta is invalid. It must be a numeric vector of positive length, and all entries finite." )
        return(NULL)
        }
    
    normal  = as.numeric(normal)

    #   compute functional in R^n
    functional  = as.numeric( x$W %*% normal )   #; print( functional ) 
    
    #   find vertex where functional is maximized
    vertex_max  = 0.5 * sign( functional )      #; print( vertex_max )  # this a vertex of the n-cube, translated by -1/2
    #   vertex  = crossprod( vertex_max, x$W )
    betamax = sum( functional * vertex_max )    # the maximum of <x,normal> for x in the zonotope
    betamin = -betamax  # by symmetry
    
    #   the face centers are much easier to work with when the antipodal ones are added to the originals
    center      = rbind( x$face$center, -x$face$center )
    W2          = rbind( x$W, x$W )    
    functional  = rep( functional, 2 )     
    delta       = 0.5 * abs(functional)  

    
    cn  = center %*% normal

    #   compute range of normal over each edge/face
    cnpos   = cn + delta
    cnneg   = cn - delta

    #   translate beta to the centered zonogon
    betacentered   = as.numeric(beta) - sum( x$center * normal )  #; print(betacentered)

    out         = vector( length(beta), mode='list' )
    names(out)  = sprintf( "normal=%g,%g. beta=%g", normal[1], normal[2],  beta )
    
    
    tol = 1.e-10
    
    for( k in 1:length(beta) )
        {
        beta_k  = betacentered[k]
        if( beta_k < betamin-tol  ||  betamax+tol < beta_k )
            {
            # line does not intersect the zonogon. there is no section
            out[[k]]    = list( beta=beta[k],  section=matrix( 0, 0, 2 ) )              
            next    
            }
            
        if( abs(beta_k - betamax) < tol )
            {
            #   special case - only one point of intersection
            out[[k]]    = list( beta=beta[k],  section=matrix( vertex_max %*% x$W  + x$center, 1, 2 ) )           
            next
            }
        if( abs(beta_k - betamin) < tol )
            {
            #   special case - only one point of intersection
            out[[k]]    = list( beta=beta[k],  section=matrix( -vertex_max %*% x$W  + x$center, 1, 2 ) )                  
            next
            }            
        
        idx = which( cnneg <= beta_k  &  beta_k < cnpos )  #; print( idx )
        
        count   = length(idx)
        ok  = count %in% 1:2
        if( ! ok )  
            {    
            # should not happen
            # log.string( WARN, .... )
            next
            }

        section = matrix( NA_real_, count, 2 )

        for( i in 1:count )
            {
            j   = idx[i]
            
            alpha   = (beta_k - cn[j])/ functional[j]   #; print(alpha)
            
            p   = center[j, ] + alpha * W2[j, ] + x$center  # translate from centered back to original
            
            section[i, ] = p
            }
            
        out[[k]]            = list()      
        out[[k]]$beta       = beta[k]
        out[[k]]$section    = section   

        colnames(out[[k]]$section)  = names(normal)         
        }
    
    return( invisible(out) )
    }
    
    
print.zonogon  <-  function( x, ... )
    {
    cat( "number of faces:     ", nrow(x$face), ' [only half due to central symmetry]\n' )
    
    count   = sum( is.na( rowSums(x$face$normal) ) )
    cat( "degenerate faces:    ", count, ' [length too small]\n' )
    
    count   = ifelse( is.null(x$parallelogram), 0, nrow(x$parallelogram) )
    cat( "parallelograms:      ", count, '\n' )   
    
    cat( "vertices:            ", nrow(x$vertex), ' [complete]\n' )    
    cat( "center:              ", x$center, '\n' )
    cat( "non-trivial clusters:", length(x$cluster), '\n' )
    
    if( 0 < length(x$cluster) )
        {
        for( k in 1:length(x$cluster) )
            {
            cluster = sort( x$cluster[[k]] )
            cat( "    #", k, "  normal=", x$face$normal[ cluster[1], ], "   ", cluster, '\n' )
            }
        }
        
    return( invisible(TRUE) )
    }

plot.zonogon  <- function( x, interior=TRUE, ... )
    {
    xlim    = range( x$vertex[ ,1], na.rm=T )
    ylim    = range( x$vertex[ ,2], na.rm=T )
    
    plot( xlim, ylim, type='n', las=1, asp=1, lab=c(10,10,7) )
    grid( lty=1 )
    abline( h=0, v=0 )
    polygon( x$vertex[ ,1], x$vertex[ ,2], col='white', border='red' )
    
    n   = nrow(x$W)
    
    
    if( interior  &&  ! is.null(x$parallelogram) )
        {
        #   plot parallelogram partition
        cvec    = grDevices::rainbow( n-1, s=0.3 )
        
        par4    = matrix( 0, 4, 2, byrow=TRUE )
        par4[1, ] = c(-0.5,-0.5)                
        par4[2, ] = c(0.5,-0.5)
        par4[3, ] = c(0.5,0.5)
        par4[4, ] = c(-0.5,0.5)
        #print( par4 )        
        
        for( k in 1:nrow(x$parallelogram) )
            {
            i   = x$parallelogram$idx[k,1]
            j   = x$parallelogram$idx[k,2]

            #cat( "---------",  i,j, "-----------------\n" )
            z   = par4 %*% x$W[ c(i,j), ]
            
            z   = z + matrix( x$parallelogram$center[k, ], 4, 2, byrow=TRUE )
            
            polygon( z[ ,1], z[ ,2], border='black', col=cvec[j-1] )
            }
        }
        
    if( nrow(x$vertex) <= 32 )    
        #   mark vertices nicely
        points( x$vertex[ ,1], x$vertex[ ,2], pch=19, cex=0.8 )
            
    #   plot the center
    points( x$center[1], x$center[2], cex=1, pch=21, bg='white' )
    
    # title( main=sprintf( "zonogon with %d vertices, partitioned into %d parallelograms", n, n*(n-1)/2 ) )
    
    return( invisible(TRUE) )
    }
    
##-------------------   helpers below       -------------------##

#   A           a matrix, possibly with NAs
#   tol         difference tolerance, consecutive in nature
#   projective  2 rows that differ only in sign are considered the same
#
#   returns a list with
#       clusteridx  an integer vector with length(cluster) = nrow(A)
#                   an index into cluster.  0 means this row is a trivial singleton cluster (most common).
#       cluster     list of integer vectors, the non-trivial clusters
#
#   a row with NAs is ignored, and always in a trivial cluster
#
findRowClusters <- function( A, tol=1.e-9, projective=TRUE )
    {
    if( projective )
        {
        #   extract the first non-zero entry in each row
        first   = apply( A, 1, function(r) { r[ which(r!=0)[1] ] } )
        first[ is.na(first) ] = 0   # change NAs to 0
        A   = sign(first) * A     # vector sign(first) is auto-replicated to all columns
        }
    
    if( ncol(A) == 3 )
        #   most common
        perm    = base::order( A[ ,1], A[ ,2], A[ ,3] )         
        # perm    = statnet.common::order( A )  # ; print(perm)    
    else if( ncol(A) == 2 )
        perm    = base::order( A[ ,1], A[ ,2] )  # ; print(perm)
    else
        {
        log.string( FATAL, "ncol(A)=%d is invalid.", ncol(A) )
        return(NULL)
        }
        
    Asorted = A[ perm, ] 
    
    Adiff   = diff( Asorted )   #; print(Adiff)  # one fewer rows than A
    jump    = .rowSums( abs(Adiff), nrow(Adiff), ncol(Adiff) )  #; print( sort(jump) )
    
    
    jump[ is.na(jump) ] = tol + 1
    mask    = tol < jump
    runidx  = findRunsFALSE( mask )
    
    out = list()      
    out$clusteridx  = integer( nrow(A) )    
    out$cluster     = vector( nrow(runidx), mode='list' )
    out$diameter    = numeric( nrow(runidx) )
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
            out$diameter[k]         = max( spread( Asorted[run, ] ) )
            }
            

        }
        
    return(out)
    }
    
    
#   mat     an NxM matrix
#   returns an M-vector.  i'th entry = spread of i'th column of mat    
spread <- function( mat )
    {
    as.numeric( diff( apply( mat, 2, range, na.rm=T ) ) )   
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
        log.string( FATAL, "Internal error.  %d != %d", length(start), length(stop) )
        return(NULL)
        }
        
    return( cbind( start=start, stop=stop-1L ) )
    }