
#   'zonohedron' is an S3 class
#   It presents the zonohedron as an intersection of 'slabs'.
#   Each slab is the intersection of 2 halfspaces with the same normal vector.
#   The equation of a slab is:
#       -beta <=   <x,normal>   <= beta      (beta is always positive)
#   All the slabs are centered, so their intersection is centered too.
#
#   It is stored as a list with
#       W           the original given nx3 matrix, defining n line segments in R^3  (with one endpoint at 0).  Must be rank 3.
#       center      middle gray = white/2
#       nonnegative logical.  All entries of W are non-negative
#       face        a data.frame with n*(n-1)/2 rows.  Because of central symmetry, only 1/2 of the faces need to be stored.
#                   idx     2 columns with 2 indices to distinct rows of W
#                   normal  to a face of the zonohedron, the cross-product of the 2 rows, then unitized
#                   center  of the face, a parallelogram, *after* centering the zonohedron
#                   beta    the plane constant, *after* centering the zonohedron.  All are positive.
#
#   If there are C clustered (compound) faces, then there are also:
#       cluster     a list of length C.  cluster[[k]] is an integer vector of the face indexes of the cluster
#       clusteridx  integer vector of length n(n-1)/2.  It points back to cluster.
#                   for face i, clusteridx[i] is the index of the cluster to which it belongs, or 0 if the face is a singleton.
#       diameter    a numeric vector of length C.  diameter[k] is the diameter of the unit face normals in cluster k.
#                   This is controlled by a numerical threshold,  and is usually < 1.e-9.  It can be 0 too.
#                   The diameter is useful for diagnostics.
#   These clusters are based on the face normals, i.e. the point on the unit sphere.  Antipodal points are considered the same.


#   zonohedron() constructor for a zonohedron object
#       W       nx3 matrix with rank 3
#       tol     sin(angle) tolerance for degenerate parallelogram faces
#       perf    if TRUE, the print performance (timing) data
#   returns: a list as above

zonohedron  <- function( W, tol=1.e-6, perf=FALSE )
    {
    if( perf )  time0   = gettime()
    
    ok  = is.numeric(W)  &&  is.matrix(W)  &&  3<=nrow(W)  &&  ncol(W)==3
    if( ! ok )
        {
        log.string( ERROR, "argument W is not an nx3 numeric matrix, with n>=3." )
        return(NULL)
        }
        
    if( any( is.na(W) ) )
        {
        log.string( ERROR, "matrix W is invalid because it has %d entries that are NA.", sum(is.na(W)) )
        return(NULL)
        }
        

    n   = nrow(W)    
    
    if( perf )
        {
        cat( sprintf( "computing %d face normals...", n*(n-1)/2 ) )
        flush.console()
        time_start  = gettime()
        }
        
        
    p12 = tcrossprod( W[ ,1,drop=F], W[ ,2,drop=F] )
    p13 = tcrossprod( W[ ,1,drop=F], W[ ,3,drop=F] )
    p23 = tcrossprod( W[ ,2,drop=F], W[ ,3,drop=F] )
    
    d12 = p12 - t(p12)
    d13 = p13 - t(p13)
    d23 = p23 - t(p23)
    

    
    if( requireNamespace( 'arrangements', quietly=TRUE ) )
        idx = arrangements::combinations(n,2)       # faster
    else
        idx = t( utils::combn(n,2) )   # matrix of pairs.  slower
    
    normal  = cbind( d23[idx], -d13[idx], d12[idx] )
    
    #   set rows close to 0 to all NAs
    W1      = W[ idx[,1], ]
    W2      = W[ idx[,2], ]
    W1W2    = .rowSums( W1*W1, nrow(W1), ncol(W1) ) * .rowSums( W2*W2, nrow(W2), ncol(W2) )
    
    normal2 = .rowSums( normal*normal, nrow(normal), ncol(normal) )
    
    #   sin(angle) = sqrt(normal2 / W1W2)
    bad     = normal2 <= tol*tol * W1W2     #; print( bad )
    if( all(bad) )
        {
        log.string( ERROR, "matrix W does not have rank 3, with relative tol=%g.", tol )
        return(NULL)
        }
        
    area    = sqrt(normal2)
    angle   = asin( sqrt( normal2 / W1W2 ) )
    
    if( any(bad) )
        {    
        #mess    = sprintf( "%d bad face normals, out of %d.\n", sum(bad), length(bad) )
        #cat(mess)
        #   log.string( INFO, "%d normals flagged as too small, out of %d.", sum(bad), length(bad) )
        normal[ bad, ]  = NA_real_        
        #   angle[ bad ]    = NA_real_            
        }
    
    #   now unitize
    normal  = normal / area    # sqrt(normal2) is replicated to all 3 columns

    if( perf ) cat( sprintf( "[%g]\n", gettime()-time_start ) )

    

    out = list()
    
    out$W   = W
    
    out$center  = 0.5 * .colSums( W, nrow(W), ncol(W) )
    
    out$nonnegative = all( 0 <= W )
    
    out$face        = data.frame( row.names=1:nrow(idx) )
    out$face$idx    = idx
    out$face$area   = area
    out$face$angle  = angle    
    out$face$normal = normal

    
    if( TRUE )
        {
        if( perf )
            {
            cat( sprintf( "clustering %d face normals...", n*(n-1)/2 ) )
            flush.console()            
            time_start = gettime()
            }

        
        #   cluster the unit normal vectors  [reminds me of something similar at Link]
        res = findRowClusters( out$face$normal, tol=1.e-10 )
        if( is.null(res) )  return(NULL)
        
        out$cluster     = res$cluster       #; print( out$cluster )        
        out$clusteridx  = res$clusteridx    
        out$diameter    = res$diameter
        
        out$clusteridx[bad] = NA_integer_        #; print( out$clusteridx )
        
        if( all( out$clusteridx==1  |  is.na(out$clusteridx) ) )
            {
            log.string( ERROR, "matrix W does not have rank 3." )
            return(NULL)
            }
        
        if( 0 < length(out$cluster) )
            {
            #   normal vectors within each cluster may differ in sign,
            #   so make them consistent !
            
            for( k in 1:length(out$cluster) )
                {
                cluster = out$cluster[[k]]
                normalk = out$face$normal[cluster[1], ]            # take the first one's normal to represent all of them
                
                s   = as.numeric( out$face$normal[cluster, ] %*% normalk )  #; print(s)
                
                #   s should already be +1 or -1, but may be a little off because of truncation, so use sign(s)
                out$face$normal[cluster, ]  = sign(s) * out$face$normal[cluster, ]
                }
            }
            
        if( perf ) cat( sprintf( "[%g]\n", gettime()-time_start ) )    
        }

        

    
    #   the next line is the biggest bottleneck
    #   the size of the output is on the order n^3
    #   functional is n x n(n-1)/2
    #   each column corresponds to an i<j pair, from array idx[]
    #   the corresponding entries in that column, in rows i and j, should be very close to 0    
    
    if( perf )
        {
        cat( sprintf( "computing %d face centers...", n*(n-1)/2 ) )
        flush.console()            
        time_start = gettime()
        }
        
    functional  = tcrossprod( W, out$face$normal )  # output size ~ n^3 / 2

    if( FALSE  &&   0 < length(out$cluster) )
        {
        W2  = W*W
        zeromask    = ( .rowSums( W2, nrow(W2), ncol(W2) ) == 0 )
            
        #   use functional to check clusters
        for( k in 1:length(out$cluster) )
            {
            cat( "-----------------\n" )
            cluster = out$cluster[[k]]  # indexes into the pairs
            i   = cluster[1]
            
            print( normal[i, ] )
            inprod  = functional[ ,i]
            
            j   = which( abs(inprod) < 1.e-9  &  ! zeromask )
            print( length(j) )
            print(j)
            
            #   repeat
            j2  = unique( c(idx[cluster,1],idx[cluster,2]) )
            print( length(j2) )
            print(j2)
            }
        }        
    

    if( TRUE )
        {
        #time_start = gettime()
        
        #colnames( functional )  = apply( idx, 1, function(r){ paste(r,collapse=',') } )   ;  print( functional )
        
        across  = 1:nrow(idx)
        idx2    = rbind( cbind(idx[,1],across),  cbind(idx[,2],across) )
        
        #check   = functional[idx2]  ; print(check)      # all these are very small, or actually 0
        
        #   without this next line the computed face$center could be one of 9 different points in the parallelogram.
        #   it could be center, or a vertex, or a midpoint of an edge
        #   By forcing this to 0, it forces it to the center
        functional[idx2] = 0
        
        #time_elapsed    = gettime() - time_start
        #cat( "time_elapsed =", time_elapsed, '\n' )    # takes  <1% to 5% of total time 
        }
    
    
    #   calculate the n(n-1)/2 face centers of the unit cube,
    #   but then all translated by -0.5 and therefore centered
    #   sign(0) = 0 exactly, so 2 (or more) of the 'spectrum' coords are 0
    face_center = 0.5 * sign(functional)      #   exact 0 maps to exact 0

    out$face$center = crossprod( face_center, W )
    
    if( perf ) cat( sprintf( "[%g]\n", gettime()-time_start ) )        
    
    if( 0 < length(out$cluster) )
        {
        #   if a simple p-face is in a cluster of p-faces with the same normal vector, 
        #   then face$center is probably not correct.
        #   Fix these by transforming to 2D, computing zonogon, and then back to 3D again
        
        if( perf )
            {
            cat( sprintf( "partitioning %d compound faces...", length(out$cluster) ) )
            flush.console()            
            time_start = gettime()
            }
        
        for( k in 1:length(out$cluster) )
            {
            cluster     = out$cluster[[k]]
            
            #   take the first one's normal to represent all of them
            #   we have already made them consistent in sign
            normalk     = out$face$normal[ cluster[1], ]            
            frame3x2    = base::svd( normalk, nu=3 )$u[ , 2:3 ]   #; print(frame3x2)
            
            idx = out$face$idx[ cluster, ]
            segmentvec  = sort( unique( as.integer(idx) ) ) #; print(segmentvec)
            
            #  project the segments from 3D to 2D
            W2      = out$W[ segmentvec, ] %*% frame3x2     #; print(W2)
            zono    = zonogon( W2, partition=TRUE )         #;   plot( zono )
            if( is.null(zono) ) return(NULL)
            
            #   find center of the compound face in 3D
            functional  = as.numeric(out$W %*% normalk)     #; print(functional)
            functional[ segmentvec ] = 0    # exactly
            face_center = 0.5 * sign(functional)            #; print(face_center) # in the unit cube
            centerface  = crossprod( face_center, out$W )   #; print(centerface)  # in the centered zonohedron
            
            # transform parallelogram centers from zono to 3D
            pgrams  = nrow(zono$parallelogram$center)
            
            #   center2 contains pgram centers, in the *centered* zonogon
            center2 = zono$parallelogram$center - matrix( zono$center, pgrams, 2, byrow=TRUE )
            
            #   center3 contains the pgram centers, in the centered zonohedron
            center3 = tcrossprod( center2, frame3x2 )  +  matrix( centerface, pgrams, 3, byrow=TRUE )
            #  print( center3 )
            
            #   copy the centers from center3 to face$center, while tracking the pair indexes
            for( kk in 1:pgrams )
                {
                idx = zono$parallelogram$idx[kk, ]
                
                #   map from zonogon to zonohedron indices
                idx = segmentvec[idx]   #; print(idx)
                
                #   find the p-face number !
                i   = which( out$face$idx[ ,1] == idx[1]   &   out$face$idx[ ,2] == idx[2] )
                
                if( i %in% cluster )    out$face$center[i, ] = center3[kk, ]
                }
            }
            
        if( perf ) cat( sprintf( "[%g]\n", gettime()-time_start ) )                
        }

    
    if( perf )
        {
        cat( sprintf( "computing %d plane constants...", n*(n-1)/2 ) )
        flush.console()            
        time_start = gettime()
        }
            
    #   calculate the n(n-1)/2 plane constants beta
    #   these plane constants are for the centered zonohedron
    out$face$beta   = .rowSums( out$face$normal * out$face$center, nrow(out$face$normal), 3 )
    
    #   since W has rank 3, 0 is in the interior of the _centered_ zonohedron,
    #   and this implies that all n(n-1)/2 beta's are positive, except the NAs.
    #   verify this
    betamin = min( out$face$beta, na.rm=TRUE )
    if( betamin <= 0 )
        {
        log.string( FATAL, "Internal Error.  min(beta)=%g <= 0.", tol, betamin )
        return(NULL)
        }
        
    if( FALSE  &&  0 < length(out$cluster) )  # Set to FALSE for RELEASE.   
        {
        #   within each cluster, all the faces should have the same plane constant (up to a tiny epsilon)
        betarange   = sapply( out$cluster, function(ivec) { diff( range( out$face$beta[ivec] ) ) }  )   # print( betarange )
        
        tiny    = 5.e-8
        if( tiny < max(betarange) )
            {
            log.string( WARN, "For cluster %d, plane constant range = %g > %g.", 
                                which.max(betarange), max(betarange), tiny )
            }
        }

    if( perf ) cat( sprintf( "[%g]\n", gettime()-time_start ) ) 
        
    class(out)  = c( 'zonohedron', class(out) )
        
    if( perf )  cat( sprintf( "Total constructor time: %g.\n",  gettime() - time0 ) )
    
    return(out)
    }
    
    
##----------        zonohedron methods    -------------##
    
#   x   a zonohedron object
#   g   an Nx3 matrix, etc.
#
#   value   a dataframe with columns
#           inside      logical
#           distance    numeric, non-positive means inside
inside.zonohedron <- function( x, g )
    {
    g   = prepareNxM( g, 3 )
    if( is.null(g) )    return(NULL)
    
    #   translate g to the centered zonohedron
    gcentered   = g - matrix( x$center, nrow(g), 3, byrow=TRUE ) #; print(gcentered)
    
    hg  = tcrossprod( x$face$normal, gcentered )    #; print( str(hg) )
    
    distance    = abs(hg) - matrix( x$face$beta, nrow(hg), ncol(hg) )
    
    distance    = base::apply( distance, 2, function(z) {suppressWarnings( max(z,na.rm=TRUE) ) } )  #;   print(distance)
    
    distance[ ! is.finite(distance) ]   = NA_real_
    
    if( x$nonnegative )
        {
        #   special override for black
        black   = apply( g, 1, function(v) { isTRUE(all(v==0)) } )  #; print(black)
        if( any(black) )
            distance[black] = 0
            
        #   special override for white.  Fortunately multiplication by 0.5 and 2 preserves all precision.
        white   = apply( g, 1, function(v) { isTRUE(all(v==2*x$center)) } )  #; print(white)
        if( any(white) )
            distance[white] = 0
        }
        
    out = data.frame( inside=(distance<=0) )
    rownames(out)   = rownames(g)
    
    out$g           = g
    out$distance    = distance

    return(out)
    }
    
    
#   x           a zonohedron object
#   base        a numeric vector of length 3, the basepoint of all the rays
#               base must be in the interior of x,
#               or if x is non-negative, base can also be the black or white point on the boundary(x)
#   direction   an Nx3 matrix with directions in the rows
#
#   value   a dataframe with columns
#           base        given basepoint of the ray (all the same)
#           direction   given direction of the ray
#           idx         of the parallelogram face where ray exits the zonohedron
#           tmax        ray parameter of intersection with face
#           boundary    the intersection with face
#           alpha       2 coordinates of boundary in the parallelogram coords.  both in [0,1]

raytrace.zonohedron <- function( x, base, direction )
    {
    ok  = is.numeric(base)  &&  length(base)==3  &&  all( is.finite(base) )
    if( ! ok )
        {
        log.string( ERROR, "base is invalid. It must be a numeric vector of length 3, and all entries finite." )
        return(NULL)
        }
        
    direction   = prepareNxM( direction, 3 )
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
        boundary    = abs(distance) < rtol * mean(abs(x$center))     #; print(boundary)
        distance[boundary]  = 0                                     #; print(distance)
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
        log.string( ERROR, "point base=(%g,%g,%g) is not in the interior of the zonohedron.  distance=%g >= 0.", 
                                base[1], base[2], base[3], distance )
        return(NULL)
        }
        
    n   = nrow(direction)
    

    tmax        = rep(NA_real_,n)
    idx         = matrix(NA_integer_,n,2)
    sign        = rep(NA_integer_,n)        
    boundary    = matrix(NA_real_,n,3)
    alpha       = matrix(NA_real_,n,2)
    timetrace   = rep(NA_real_,n)
    faces       = rep(NA_integer_,n)
    tested      = rep(NA_integer_,n)
    
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
            bound   = abs(numerator) < rtol * sum(abs(x$center))    #; print( sum(is.na(bound)) )
            
            bound[ is.na(bound) ] = FALSE
            #numerator[bound]  = 0

            tvec[bound]  = ifelse( 0 < hv[bound], 0, Inf )    # 0 here will later make the corresponding boundary = NA 
            }

        tvec[ ! is.finite(tvec) ]   = Inf
        
        j   = which.min( tvec )     # this ignores Infs
        
        tmax[k]     = tvec[j]       # tmax[k] is not negative, and might be 0 if base is black or white

        idx[k, ]    = x$face$idx[j, ] 
        
        if( tmax[k] <= 0 )  next    # failed to intersect properly
        

        #   cat( "--------------", v, "---------------------------\n" )
        
        optcentered     = gcentered  +  tmax[k] * v             #;  print( optcentered )         
        boundary[k, ]   = optcentered + x$center
        
        cidx  =  x$clusteridx[j]    
        
        if( is.na(cidx) )   next    # should not happen
        
        if( cidx != 0 )
            {
            cluster = x$cluster[[ cidx ]]   # always more than 1
            #   cat( sprintf( "Searching in cluster %d, with %d subfaces.\n", cidx, length(cluster) ) )
            }
        else
            cluster = j     # a singleton cluster
        
        faces[k]    = length(cluster)
        
        for( i in 1:length(cluster) )
            {     
            j               = cluster[i]
            
            theSign         = sign(hv[j])            

            facecenter      = theSign * x$face$center[j, ]          #;  print( facecenter ) 
            
            edges           = t( x$W[ x$face$idx[j, ], ] )         # 2 edges of the parallelogram, as 3x2 matrix
            M               = cbind( edges, x$face$normal[j, ] )    #; print( M )    # M is 3x3
            y               = base::solve( M, optcentered - facecenter )    #; print(y)
            
            #   test for inside parallelogram, with a tolerance
            if( all(abs(y[1:2]) <= 0.5 + 5.e-7 ) )
                {
                #   found it
                idx[k, ]    = x$face$idx[j, ]   # override above
                sign[k]     = theSign                
                alpha[k, ]  = pmin( pmax( y[1:2] + 0.5, 0), 1 )     # translate from [-0.5,0.5] to [0,1] and clamp
                tested[k]   = i
                break
                }
            }

        timetrace[k]    = gettime() - time_start
        }
    
    out = data.frame( row.names=1:n )
    
    out$base        = matrix( base, n, 3, byrow=TRUE )  # replicate base to all rows
    out$direction   = direction
    out$idx         = idx    
    out$sign        = sign
    out$tmax        = tmax
    out$boundary    = boundary
    out$alpha       = alpha
    out$timetrace   = timetrace
    out$faces       = faces
    out$tested      = tested
    
    
    
    #if( blackwhite )    tmax[ tmax==0 ] = NA_real_      # so boundary becomes NA too !
    #out$boundary     = out$base  +  tmax * direction   # tmax is replicated to 3 columns

    cnames  = colnames(base)
    if( is.null(cnames) )   cnames = colnames(direction)    
        
    colnames(out$boundary)   = cnames
    
    return( out )
    }
    
    

#   given a point on the boundary, return a point in the unit n-cube that maps to it, under W  
#       x       the zonogon
#       data    a data.frame as returned from raytrace.zonohedron().  Columns used are:
#               boundary    coordinates of the point on the boundary
#               idx         index of the face
#               sign        sign of the face = +-1
#               alpha       coordinate along the face, in [0,1]
#               faces       in case the faces is a member of a compound face
#               tested      the number of faces actually tested during the search
#       tol     tolerance for verification.  Set to NULL for RELEASE.
#
#   returns:
#       an mxn matrix, where m=nrow(data)  n=number of rows in x$W 
#       each row of the matrix is a point in the n-cube
#       that maps to the point on the zonogon

invertboundary.zonohedron <- function( x, data, tol=NULL )      #1.e-9 )
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
    
    faceindex   = apply( data$idx, 1, function(pair) {  which( x$face$idx[ ,1] == pair[1]  &  x$face$idx[ ,2] == pair[2] ) } )
    
    vertex  = matrix( NA_real_, m, n )
    
    for( i in 1:m )
        {
        if( 1 < data$faces[i] ) next    # compound faces not working yet
        
        k   = faceindex[i]
                
        normalk             = x$face$normal[k, ]
        
        functional          = as.numeric( x$W %*% normalk )       #; print(functional)
        pair                = as.integer( x$face$idx[k, ] )
        functional[pair]    = 0    # exactly
        vertex[i, ]         = 0.5 * ( data$sign[i] * sign(functional) + 1 )      #; print(face_center) # in the n-cube         
        vertex[i,pair]      = data$alpha[i, ]                   # overwrite 0.5 and 0.5
        }
        
    if( is.numeric(tol) && is.finite(tol) && 0<tol )
        {
        #   map from boundary of n-cube to the boundary of zonohedron
        bpoint  = vertex %*% x$W
        delta   = abs(bpoint - data$boundary)   #; print( delta )
        
        attr( vertex, 'delta' ) = delta
        
        if( tol <= max(delta,na.rm=TRUE) )
            log.string( WARN, "Verification failed. max(delta)=%g > %g.", max(delta), tol )        
        }
           
    colnames(vertex)    = rownames( x$W )
    #   rownames(vertex)    = as.character( data$idx )
        
    return( vertex )
    }
    
    
    

#   section() compute intersection of zonohedron and plane(s)
#    
#   x           a zonohedron object
#   normal      a non-zero numeric vector of length 3, the normal of all the planes
#   beta        a vector of plane-constants.  The equation of plane k is: <x,normal> = beta[k]
#
#   value   a list of length = length(beta).  Each list has the items:
#           beta        given plane constant of the plane
#           section     Mx3 matrix of points on the section, in order around the boundary
#                       M=1 for support lines, and if there is no intersection, then M=0 rows.

section.zonohedron <- function( x, normal, beta )
    {    
    ok  = is.numeric(normal)  &&  length(normal)==3  &&  all( is.finite(normal) )
    if( ! ok )
        {
        log.string( ERROR, "normal is invalid. It must be a numeric vector of length 3, and all entries finite." )
        return(NULL)
        }
        
    ok  = is.numeric(beta)  &&  0<length(beta)   &&  all( is.finite(beta) )
    if( ! ok )
        {
        log.string( ERROR, "beta is invalid. It must be a numeric vector of positive length, and all entries finite." )
        return(NULL)
        }
        
    cnames  = names(normal) # save these
    dim(normal) = NULL
    
    n   = nrow( x$W )

    #   compute functional in R^n
    functional  = as.numeric( x$W %*% normal )   #; print( functional ) 
    
    #   find vertex of centered zonohedron where functional is maximized
    vertex_max  = 0.5 * sign( functional )      #; print( vertex_max )  # this a vertex of the n-cube, translated by -1/2
    #   vertex  = crossprod( vertex_max, x$W )
    betamax = sum( functional * vertex_max )    # the maximum of <x,normal> for x in the zonotope
    betamin = -betamax  # by symmetry
    
    # print( betamax )
    
    
    #   find points on boundary of centered zonohedron where <x,normal> is maximized and minimized
    bp_max  = as.numeric( crossprod( x$W, vertex_max ) )
    bp_min  = -bp_max        # by symmetry
    
    #   the face centers are much easier to work with when the antipodal ones are added to the originals
    
    
    #   make matrix of all parallelogram centers with n*(n-1) rows
    center      = rbind( x$face$center, -x$face$center )
    
    #   make matching products of these centers and normal vector
    cn  = as.numeric( x$face$center %*% normal )
    cn  = c( cn, -cn )
    
    #   make delta1 and delta2, these are signed and with length n*(n-1)/2
    delta1      = functional[ x$face$idx[ ,1] ]
    delta2      = functional[ x$face$idx[ ,2] ]
    
    #   make delta with length n*(n-1), this is the max delta of the inner product with normal, and we apply it on both sides
    delta       = pmax( abs(delta1), abs(delta2) )
    delta       = rep( delta, 2 )


    #   compute range of normal over each parallelogram face
    cnneg   = cn - delta/2    
    cnpos   = cn + delta/2

    #   print( data.frame( idx=x$face$idx, cn=cn, delta=delta ), max=300 )

    #   translate beta to the centered zonogon
    betacentered   = as.numeric(beta) - sum( x$center * normal )  #; print(betacentered)

    out         = vector( length(beta), mode='list' )
    names(out)  = sprintf( "normal=%g,%g,%g. beta=%g", normal[1], normal[2], normal[3], beta )
    
    
    tol = 1.e-10
    
    #   make tiny lookup table for wraparound
    inext   = c( 2L, 3L, 4L, 1L )
    
    #   make 4x2 weights for the vertices of the centered parallelogram, in order around circumference
    vertmat =  matrix( c( -0.5,-0.5, -0.5,0.5, 0.5,0.5, 0.5,-0.5), 4, 2, byrow=TRUE  )
    
    #   find 3x2 matrix for projection to 2D
    frame3x2    = base::svd( normal, nu=3 )$u[ , 2:3 ]   #; print(frame3x2)
    
    
    for( k in 1:length(beta) )
        {
        beta_k  = betacentered[k]
        if( beta_k < betamin-tol  ||  betamax+tol < beta_k )
            {
            #   plane does not intersect the zonohedron, there is no section
            out[[k]]    = list( beta=beta[k],  section=matrix( 0, 0, 3 ) )      
            next    
            }
            
            
        if( abs(beta_k - betamax) < tol )
            {
            #   special case - only one point of intersection
            out[[k]]    = list( beta=beta[k],  section= matrix( vertex_max %*% x$W  + x$center, 1, 3 ) )
            next
            }
        if( abs(beta_k - betamin) < tol )
            {
            #   special case - only one point of intersection
            out[[k]]    = list( beta=beta[k],  section= matrix( -vertex_max %*% x$W  + x$center, 1, 3 ) )       
            next
            }            
        
        #   find indexes of all parallelograms that intersect this plane
        indexvec = which( cnneg <= beta_k  &  beta_k < cnpos )  #; print( indexvec )

        if( length(indexvec) == 0 )  
            {    
            # should not happen
            # log.string( WARN, .... )
            next
            }

        # visit all the parallelograms
        # cat( sprintf( "In section %d, found %d pgrams, beta=%g.\n", k, length(indexvec), beta[k] ) )

        section = matrix( NA_real_, length(indexvec), 3 )

        for( j in 1:length(indexvec) )
            {
            #   cat( "------------", j, "-------------\n" )            
            
            #   the boundary of each parallelogram intersects the plane in 2 points
            #   compute one of them and add to the matrix
            #   i is the index of the pgram in the doubled arrays: center, cn, delta, cnneg, cnpos
            i   = indexvec[j]   
            
            #   imod is the index into x$face
            imod    = ( (i-1) %% nrow(x$face) ) + 1
            
            if( is.na(x$face$beta[imod]) )    next    # pgram is degenerate sliver, so ignore it
            
            pcenter = center[i, ]
            
            if( FALSE )
            {
            print( i )
            print( x$face$idx[imod, ] )
            print( pcenter + x$center )
            print( cn[i] )
            pdelta  = delta[i]            
            print( pdelta )
            }
            
            #   compute signed weights of the 4 vertices, matching the order in vertmat
            #   there should be exactly 1 transition from negative to positive
            weight  =  cn[i] - beta_k   +  as.numeric( vertmat %*% c( delta1[imod], delta2[imod] ) )   #; print( weight )
            
            sdiff   = diff(  sign( c( weight, weight[1] ) ) )
            itrans  = which( sdiff == 2 )   #;    print( itrans )
            
            if( length(itrans) != 1 )   
                {
                # failure to intersect properly, so ignore it
                cat( "Did not intersect pgram !\n" )
                print( x$face$idx[imod, ] )
                next    
                }
                
            itransnext  = inext[ itrans ]
            #itransnext  = itrans + 1
            #if( itransnext == 5 )   itransnext = 1
            
            #   compute interpolation parameter s, in the inteval [0,1]
            s   = -weight[itrans] / ( weight[itransnext] - weight[itrans] )     #; cat( "s=", s, '\n' )
            
            edgemat = (1-s) * vertmat[itrans, ]  +  s * vertmat[itransnext, ]
            
            
            #   extract 2x3 matrix of sides of this parallelogram
            W2x3    = x$W[ x$face$idx[imod, ], ]    #;print( W2x3 )                        
            
            #   finally compute intersection of boundary of pgram and the plane, in the centered zonohedron
            section[j, ] = pcenter  +  edgemat %*% W2x3      #; print( section[j, ] )
            
            #if( 2 < j ) break            
            }
            
        #   remove any NAs, which may have come from degenerate pgrams
        mask = is.finite( .rowSums( section, nrow(section), ncol(section) ) )
        section = section[ mask, ]            
            
        #   project all 3D points in section to the plane
        #   find suitable point in the interior of this section, using bp_min and bp_max
        s   = (beta_k - betamin) / (betamax - betamin)
        
        center_section  = (1-s)*bp_min + s*bp_max
        
        #print( center_section + x$center )

        p2D =  ( section - matrix(center_section,nrow(section),3,byrow=TRUE) )  %*%  frame3x2
        
        #   not a polygon yet, the points must be ordered by angle using atan2()
        perm    = order( atan2(p2D[ ,2],p2D[ ,1]) )
        
        #   reorder and translate from centered zonohedron to the original
        section = section[ perm, ] + matrix( x$center, nrow(section), 3, byrow=TRUE )
            
        out[[k]]            = list()      
        out[[k]]$beta       = beta[k]
        out[[k]]$section    = section
        
        colnames(out[[k]]$section)  = cnames  
        }

    return( invisible(out) )
    }
        
    
    
    
    
    
    
    
    
        
    
dumpface.zonohedron <- function( x, i1, i2=0 )
    {
    n   = nrow(x$face)
    if( i2 == 0 )
        {
        k = i1
        ok  = 1<=k  &&  k<=n
        if( ! ok )  return(FALSE)
        i   = x$face$idx[k,1]
        j   = x$face$idx[k,2]
        }
    else
        {
        k   = which( x$face$idx[ ,1]==i1  &  x$face$idx[ ,2]==i2 ) 
        if( length(k) != 1 )    return(FALSE)     
        i   = i1
        j   = i2
        }        
        
    #ok  = 1<=i1  &&  i1<i2  &&  i2<=nrow(x$W)
    #if( ! ok )  return(FALSE)
    
    
    print( x$face[k, ] )
    
    w   = x$W[ x$face$idx[k, ], ]
    rownames(w) = x$face$idx[k, ]
    
    #vec = c( angle=angleBetween(w[1, ],w[2, ]) )
    #cat( '\n' )
    #print( vec )    
    
    lens    = sqrt( rowSums(w*w) )
    w   = cbind( w, length=lens)
    
    cat( '\n' )
    print( w )
    
    cat( "\n" )    
    if( is.na(x$face$beta[k]) )
        {
        mess    = sprintf( "Face %d(%d,%d) is degenerate.\n", k, i, j )
        cat( mess )
        return(FALSE)
        }
    

    normalk     = x$face$normal[k, ]
    functional  = as.numeric( x$W %*% normalk )
    names(functional)    = 1:length(functional)
    func    = functional
    func[ c(i,j) ] = 0
    vertex  = 0.5*( sign(func) + 1 )
    functional  = rbind( functional, vertex )
    print(functional)    
    
    #   examine cluster situation
    cat( "\n" )    
    if( x$clusteridx[k] == 0 )
        {
        mess    = sprintf( "Face %d(%d,%d) is in no cluster.\n", k, i, j )
        cat( mess )
        return(  invisible(TRUE) )
        }
        
    cidx    = x$clusteridx[k]
    cluster = sort( x$cluster[[ cidx ]] )
    
    idx = x$face$idx[ cluster, ]
    segmentvec  = sort( unique( as.integer(idx) ) )    
    mess    = sprintf( "Face %d(%d,%d) is in cluster %d, with %d faces and %d segments.\n",  
                        k, i, j, cidx, length(cluster), length(segmentvec) )
    cat( '\n' )                        
    cat( mess )
    
    cat( '\n' )    
    df  = x$face[cluster, ]
    print( df )
    
    
    cat( "\nSpread of normal vectors:\n" )
    print( as.numeric( diff( apply( df$normal, 2, range, na.rm=T ) ) ) )

    cat( "\nSegments:", segmentvec, '\n' )
    
    return( invisible(TRUE) )
    }
    
dumpcluster.zonohedron <- function( x, k, maxrows=50 )
    {
    ok  = 1<=k  &&  k<=length(x$cluster)
    if( ! ok )
        {
        log.string( ERROR, "argument k=%d is invalid; it is not in [1,%d].", k, length(x$cluster) )
        return(FALSE)
        }
        
    cluster = sort( x$cluster[[ k ]] )
    
    idx = x$face$idx[ cluster, ]
    segmentvec  = sort( unique( as.integer(idx) ) )    
    mess    = sprintf( "Cluster %d has %d faces and %d segments.\n",  
                        k, length(cluster), length(segmentvec) )
    cat( '\n' )                        
    cat( mess )
    
    cat( '\n' )    
    df  = x$face[cluster, ]
    print( df, max=maxrows*ncol(df) )
    
    
    cat( "\nSpread of normal vectors:\n" )
    print( as.numeric( diff( apply( df$normal, 2, range, na.rm=T ) ) ) )

    cat( "\nSegments:", segmentvec, '\n' )        
        
    return( invisible(TRUE) )    
    }
    
    
plotcluster.zonohedron <- function( x, k )
    {
    if( length(x$cluster) == 0 )
        {
        log.string( ERROR, "The zonohedron has no clusters." )
        return(FALSE)
        }
    
    ok  = 1<=k  &&  k<=length(x$cluster)
    if( ! ok )
        {
        log.string( ERROR, "k=%d is outside [1,%d].", length(x$cluster) )
        return(FALSE)
        }    
    
    
    cluster     = x$cluster[[k]]
    

    #   take the first one's normal to represent all of them
    #   we have already made them consistent in sign
    normalk     = x$face$normal[ cluster[1], ]            
    frame3x2    = base::svd( normalk, nu=3 )$u[ , 2:3 ]   #; print(frame3x2)
    
    idx = x$face$idx[ cluster, ]
    segmentvec  = sort( unique( as.integer(idx) ) ) #; print(segmentvec)
    
    #   find center of the compound face in 3D
    functional  = as.numeric(x$W %*% normalk)       #; print(functional)
    functional[ segmentvec ] = 0    # exactly
    face_center = 0.5 * sign(functional)            #; print(face_center) # in the unit cube
    centerface  = crossprod( face_center, x$W )     #; print(centerface)  # in the centered zonohedron
        
        
    faces   = length(cluster)   # parallelograms
            
    pgramlist   = vector( faces, mode='list' )
    
    xlim    = NULL
    ylim    = NULL
    
    pgram   = matrix( 0, 4, 2 )    
    for( i in 1:faces )
        {
        #   get center of parallelogram in 2D
        center  = (x$face$center[ cluster[i], ] - centerface) %*% frame3x2      #; print(center)
        
        idx     = x$face$idx[ cluster[i], ]
        
        # cat( i, "area=", x$face$area[ cluster[i] ], '\n' )
        
        #  project the segments from 3D to 2D
        edge    = x$W[ idx, ] %*% frame3x2     #; print(edge)

        pgram[1, ] = center - 0.5 * edge[1, ] - 0.5*edge[2, ]
        pgram[2, ] = center - 0.5 * edge[1, ] + 0.5*edge[2, ]
        pgram[3, ] = center + 0.5 * edge[1, ] + 0.5*edge[2, ]
        pgram[4, ] = center + 0.5 * edge[1, ] - 0.5*edge[2, ]
        
        pgramlist[[i]]  = pgram
        
        xlim    = range( pgram[ ,1], xlim )
        ylim    = range( pgram[ ,2], ylim )
        }
     
    plot( xlim, ylim, type='n', las=1, asp=1, lab=c(10,10,7) )
    grid( lty=1 )
    abline( h=0, v=0 )
    
    for( i in 1:faces )
        {    
        pgram   = pgramlist[[i]]
        polygon( pgram[ ,1], pgram[ ,2], col='white', border='red' )
        
        center  = (x$face$center[ cluster[i], ] - centerface) %*% frame3x2    
        idx     = x$face$idx[ cluster[i], ]
        text( center[1], center[2], sprintf( "%d,%d",idx[1],idx[2]  ) )
        
        points( pgram[ ,1], pgram[ ,2], pch=19 )
        }
        
    title( main=sprintf( "cluster #%d, faces=%d, center=%g,%g,%g", 
                    k, faces, centerface[1], centerface[2], centerface[3] ) )

    return( invisible(TRUE) )
    }
    
    
#   x       a zonohedron object    
#   type    'w' for wireframe, 'p' for points
#   axes    c(1,2), or 1, or 2
plot.zonohedron <- function( x, type='w', both=TRUE, axes=c(1,2), cidx=0, faceidx=0, labels=FALSE,  ... )
    {
    if( ! requireNamespace( 'rgl', quietly=TRUE ) ) 
        {
        log.string( ERROR, "Package 'rgl' is required.  Please install it." )        
        return(FALSE)
        }
            
    center  = x$center
    white   = 2 * center
    
    white4x3    = matrix( white, 4, 3, byrow=TRUE )
    
    #   start 3D drawing
    rgl::bg3d("gray50")
    # rgl::light3d()
    
    cube    = rgl::scale3d( rgl::cube3d(col="white"), center[1], center[2], center[3] )
    cube    = rgl::translate3d( cube, center[1], center[2], center[3]  )
    rgl::wire3d( cube, lit=FALSE )    
    
    #   exact diagonal of box        
    rgl::lines3d( c(0,white[1]), c(0,white[2]), c(0,white[3]), col=c('black','white'), lwd=3, lit=FALSE )       
    
    if( 1 <= cidx  &&  cidx <= length(x$cluster) )
        cluster = x$cluster[[cidx]]
    else
        cluster = NULL

        
        
    if( type == 'p' )
        {
        #   draw first half in 'black'
        offset  = matrix( center, nrow(x$face), 3, byrow=TRUE )
            
        xyz = x$face$center + offset
        rgl::points3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3] )
        
        if( both )
            {
            #   draw 2nd half in 'red'
            xyz = -x$face$center + offset
            rgl::points3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col='red' )
            }
        }
        
        
    if( type == 'w' )
        {
        faces   = nrow(x$face)
        
        step    = 2 * length(axes)  # 2 * the number of segments per face
        
           
        mat     = matrix( 0, step*faces, 3 )
        
        # quad    = matrix( 0, 4, 3 )
        
        for( i in 1:faces )
            {
            center  = x$face$center[i, ] 

            edge    = x$W[ x$face$idx[i, ], ]
            
            k       = step*(i-1)
            
            mat[k+1, ] = center - 0.5 * edge[1, ] - 0.5*edge[2, ]            
            
            if( all( axes == c(1,2) ) )
                {
                mat[k+2, ] = center - 0.5 * edge[1, ] + 0.5*edge[2, ]
                
                mat[k+3, ] = mat[k+2, ]
                mat[k+4, ] = center + 0.5 * edge[1, ] + 0.5*edge[2, ]
                }
            else if( axes == 1 )
                mat[k+2, ] = center + 0.5 * edge[1, ] - 0.5*edge[2, ]
            else if( axes == 2 )
                mat[k+2, ] = center - 0.5 * edge[1, ] + 0.5*edge[2, ]
                
                
            #fmatch  = i %in% cluster ||  i %in% faceidx
            
            #col = ifelse( fmatch, 'red', 'black' )
            }
            
        offset  = matrix( x$center, step*faces, 3, byrow=TRUE )
        xyz = offset + mat
        
        rgl::segments3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col='black' )    #, front=polymode, back=polymode, col=col, lit=FALSE )
              

        if( both )
            {
            #   draw 2nd half in 'red'
            xyz = offset - mat
            
            rgl::segments3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col='red' )
            }
              
        # polymode    = "lines"  #ifelse( fmatch, "filled", "lines" )
        
        # rgl::quads3d( quad[ ,1], quad[ ,2], quad[ ,3], front=polymode, back=polymode, col=col, lit=FALSE )
        

        
            #if( both )
            #    {
            #    quad    = white4x3 - quad
            #    #rgl::quads3d( quad[ ,1], quad[ ,2], quad[ ,3], front=polymode, back=polymode, col=col, lit=FALSE )
            #    }
            
            
            #if( i %in% cluster )
            #    rgl::points3d( center[1], center[2], center[3], col=col )
            
            #if( labels || fmatch )
            #    {
            #    label   = sprintf( "%d,%d", x$face$idx[i,1], x$face$idx[i,2] )
            #    rgl::text3d( center[1], center[2], center[3], label, col='black', cex=0.75 )
            #    }
        }
        
    
    return( invisible(TRUE) )
    }
    
print.zonohedron  <-  function( x, maxrows=40, ... )
    {
    cat( "number of segments:  ", nrow(x$W), '\n' )
        
    cat( "number of faces:     ", nrow(x$face), ' [parallelogram faces. only half the total, due to central symmetry]\n' )
    
    mask_degen  = is.na( x$face$beta )
    count   = sum( mask_degen )
    cat( "degenerate faces:    ", count, ' [angle too small, or undefined]\n' )
   
    area_total  = 2*sum( x$face$area )
    area_degen  = 2*sum( x$face$area[ mask_degen ] ) 
    mess    = sprintf( "total area:          %g  (degenerate face area: %g.    %g of total area).\n",
                        area_total, area_degen, area_degen/area_total )
    cat( mess )    
    
    cat( "center:              ", x$center, "    [ white:", 2*x$center, ']\n' )
              
    if( length(x$cluster) == 0 )
        {
        cat( "No clusters (compound faces) were computed.\n" )
        return( invisible(TRUE) )        
        }

    cat( '\n' )
    cat( "non-trivial clusters:", length(x$cluster), '   [clustered by face normals]\n' )
    
    n       = length(x$cluster) 

    faces       = integer(n)
    segments    = integer(n)
    area        = numeric(n)
    normal      = matrix( NA_real_, n, 3 )
    
    for( k in 1:length(x$cluster) )
        {
        cluster     = sort( x$cluster[[k]] )
        faces[k]    = length(cluster)
        
        idx = x$face$idx[ cluster, ]
        segments[k] = length( unique( as.integer(idx) ) )
        
        area[k]     = sum( x$face$area[cluster] )   
        
        normal[k, ] = x$face$normal[ cluster[1], ]
        
        #mess1 = sprintf( "#%d (faces=%d, segments=%d, area=%g, diameter=%g)", 
        #            k, length(cluster), length(segments), area_cluster, x$diameter[k] )
        #mess2 = sprintf( "%d(%d,%d)", cluster, idx[ ,1], idx[ ,2] )
        #cat( "    ", mess1, "  normal=", x$face$normal[ cluster[1], ], "   ", mess2, '\n' )
        }
        
    data    = data.frame( row.names=1:n )
    data$faces          = faces   
    data$segments       = segments  
    data$area           = area  
    data$normal         = normal  
    data$diameter       = x$diameter
    
    cat( '\n' )        
    print( data, max=maxrows*ncol(data) )
    
    area_compound =  2 * sum( data$area )
    cat( '\n' )    
    cat( "area of clustered (compound faces): ", area_compound, '     ', area_compound/area_total, 'of total area\n' )        

    cat( '\n' )
    mess = sprintf( "maximum cluster normal diameter: %g  (for cluster %d).\n", 
                        max(data$diameter), which.max(data$diameter) )
    cat( mess )
    
    return( invisible(TRUE) )
    }    
    
    
#   W           nx3 matrix with zonohedron generators; may have NAs
#
#   returns TRUE or FALSE, as defined by Paul Centore
cyclicposition <- function( W  )
    {
    n   = nrow(W)
    
    inext   = c( 2:n, 1 )
    iprev   = c( n, 1:(n-1) )
    
    #   unitize all the rows
    W2      = .rowSums( W*W, nrow(W), ncol(W) )  
    Wnorm   = W / sqrt( W2 )
    
    myfun   <- function( i )
        {
        Wsub    = Wnorm[ c(iprev[i],i,inext[i]), ]
        out     = determinant( Wsub, logarithm=FALSE )
        out     = out$sign * out$modulus
        return( out )
        }
    
    dets    = sapply( 1:n, myfun )   # ; print( out )
    
    signdets    = sign( dets )
    
    maskpos     = ( signdets == 1 )
    maskneg     = ( signdets == -1 )
    mask0       = ! maskpos  &  ! maskneg
    
    countpos    = sum( maskpos )
    countneg    = sum( maskneg )
    count0      = sum( mask0 )
    
    cat( "rows:     ", n, '\n' )
    
    rangepos    = rep( NA_real_, 2 )
    rangeneg    = rep( NA_real_, 2 )
    
    if( 0 < countpos )  rangepos    = range( dets[ maskpos ] )
        
    mess    = sprintf( "positive:  %d   range = [ %g, %g ]\n", countpos, rangepos[1], rangepos[2] )
    cat( mess )
    
    if( 0 < countneg )  rangeneg    = range( dets[ maskneg ] )
           
    mess    = sprintf( "negative:  %d   range = [ %g, %g ]\n", countneg,  rangeneg[1], rangeneg[2] )
    cat( mess )
    cat( "zero:     ", count0, '\n' )
    
    if( countneg < countpos )
        {
        #   print the worst case exception
        cat( "matrix with det=", rangeneg[2], '\n' )
        i   = which( dets == rangeneg[2] )[1]
        Wsub    = W[ c(iprev[i],i,inext[i]), ]
        print( Wsub )
        }
    else
        {
        #   positive are the exceptions, print some of them
        cat( '\n' )
        dets[ maskneg | mask0 ] = Inf
        perm    = order( dets )
        
        for( k in 1:min(countpos,50) )
            {
            i   = perm[k]
            cat( "----------------------\n" )
            cat( "matrix with det=", dets[i], '\n' )
            Wsub    = W[ c(iprev[i],i,inext[i]), ]
            print( Wsub )
            }
        }
    
    
    
    out = diff( range(signdets) ) <= 1
    
    return( out )
    }
    
    
#--------       UseMethod() calls           --------------#    
    
inside <- function( x, g )
    {
    UseMethod("inside")
    }    
    
raytrace <- function( x, base, direction )
    {
    UseMethod("raytrace")
    }    
    
section <- function( x, normal, beta )
    {    
    UseMethod("section")
    }    
    
invertboundary <- function( x, data, tol=NULL )     #1.e-9 )    
    {
    UseMethod("invertboundary")
    }    
    
dumpface <- function( x, i1, i2=0 )    
    {
    UseMethod("dumpface")
    }    
    
dumpcluster <- function( x, k, maxrows=50 )    
    {
    UseMethod("dumpcluster")
    }    
    
plotcluster  <- function( x, k )    
    {
    UseMethod("plotcluster")
    }    

    
gettime <- function()
    {
    return( microbenchmark::get_nanotime() * 1.e-9 )
    }
    