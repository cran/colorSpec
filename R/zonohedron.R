
#   'zonohedron' is an S3 class
#   It presents the zonohedron as an intersection of 'slabs'.
#   Each slab is the intersection of 2 halfspaces with the same normal vector.
#   The equation of a slab is:
#       -beta <=   <x,normal>   <= beta      (beta is always positive)
#   All the slabs are centered, so their intersection is centered too.
#
#   It is stored as a list with
#       W           the original given Nx3 matrix, defining N line segments in R^3  (with one endpoint at 0).  Must be rank 3.
#                   If a row is a negative multiple of another (both non-zero), then W is invalid.
#
#       Wcond       Gx3 matrix with G rows, a row for each group of generators that are multiples of each other; the 'condensed' generators
#       group       a list of length G.  group[[k]] is an integer vector of the indexes of the rows of W.
#                   usually this is a vector of length 1.
#       groupidx    integer vector of length N.  groupidx[i] is the index of the group to which this row belongs,
#                   or NA_integer_ if the row is 0.
#
#       center      middle gray = white/2
#       nonnegative logical.  All entries of W are non-negative
#       face        a data.frame with G(G-1)/2 rows.  Because of central symmetry, only 1/2 of the faces need to be stored.
#                   idx     a pair of integers.  indices to distinct rows of Wcond
#                   area    of the parallelogram.  All are positive
#                   angle   at origin of the parallelogram.  All are positive
#                   normal  a numeric 3-vector.  The normal to the face of the zonohedron, the cross-product of the 2 rows, then unitized
#                   center  a numeric 3-vector.  center of the face, a parallelogram, relative to the center of the zonohedron.  *after* centering the zonohedron
#                   beta    the plane constant, *after* centering the zonohedron.  All are positive.
#
#   If there are C>0 clustered (compound) faces, then there are also:
#       cluster     a list of length C.  cluster[[k]] is an integer vector of the face indexes of cluster k.
#       center3D    a Cx3 matrix, with the centers of the compound faces in the rows.  This is relative to the center of the zonohedron.
#       clusteridx  integer vector of length G(G-1)/2.  It points back to cluster[[]].
#                   for face i, clusteridx[i] is the index of the cluster to which it belongs, or 0 if the face is a singleton, not part of a compound face.
#       spread      a numeric vector of length C.  spread[k] is the spread of the unit face normals in cluster k.
#                   This is controlled by a numerical threshold.  It can be 0 too.
#                   The spread is useful for diagnostics.
#       frame3x2    a list of length C.  frame3x2[[k]] is a 3x2 matrix transforming from the plane of the cluster to R^2
#       zonogon3    a list of length C.  zonogon3[[k]] is the zonogon3 for cluster k
#       lookupgen   a list of length C.   lookup zonogon3 generator -> zonohedron generator
#       lookupface  a list of length C.   lookup zonogon3 faceidx -> zonohedron faceidx
        
#   These clusters are based on the face normals, i.e. the point on the unit sphere.  Antipodal points are considered the same.


#   zonohedron() constructor for a zonohedron object
#       W       nx3 matrix with rank 3
#       tol1    sin(angle) tolerance for clustering the generators
#       tol2    sin(angle) tolerance for clustering face normals.  If NULL or NA or negative, then do not cluster.
#       perf    if TRUE, then print performance (timing) and diagnostic data
#
#   returns: a list as above, or NULL in case of global error.

zonohedron  <- function( W, tol1=5.e-7, tol2=5.e-10, perf=FALSE )
    {
    if( perf )  time0   = gettime()
    
    ok  = is.numeric(W)  &&  is.matrix(W)  &&  3<=nrow(W)  &&  ncol(W)==3
    if( ! ok )
        {
        log.string( ERROR, "argument W is not an nx3 numeric matrix, with n>=3." )
        return(NULL)
        }
        
    mask    = is.finite(W)
    if( ! all(mask) )
        {
        log.string( ERROR, "matrix W is invalid because it has %d entries that are not finite.", sum(! mask) )
        return(NULL)
        }
        

    if( perf )
        {
        cat( sprintf( "Condensing %d generators, with tol1=%g...\n", nrow(W), tol1 ) )
        flush.console()
        time_start  = gettime()
        }

    #   the tolerance here is for collinearity, which can be larger than the face normal differences
    res = condenseGenerators( W, tol=tol1 )
    if( is.null(res) )    return(NULL)
        
    if( perf )
        {
        cat( sprintf( "non-trivial groups: %d.   maxspread=%g  [%g sec]\n", 
                    sum(0 < res$spread), max(res$spread), gettime()-time_start ) )
        time_start  = gettime()
        }

        
    Wcond   = res$Wcond

    n   = nrow(Wcond)    
    
    if( perf )
        {
        cat( sprintf( "For %d condensed generators, computing %d face normals...", n, n*(n-1)/2 ) )
        flush.console()
        time_start  = gettime()
        }
        
        
    p12 = tcrossprod( Wcond[ ,1,drop=F], Wcond[ ,2,drop=F] )
    p13 = tcrossprod( Wcond[ ,1,drop=F], Wcond[ ,3,drop=F] )
    p23 = tcrossprod( Wcond[ ,2,drop=F], Wcond[ ,3,drop=F] )
    
    d12 = p12 - t(p12)
    d13 = p13 - t(p13)
    d23 = p23 - t(p23)
    

    
    if( requireNamespace( 'arrangements', quietly=TRUE ) )
        idx = arrangements::combinations(n,2)       # faster
    else
        idx = t( utils::combn(n,2) )   # matrix of pairs.  slower
    
    normal  = cbind( d23[idx], -d13[idx], d12[idx] )
    

    W1      = Wcond[ idx[,1], ]
    W2      = Wcond[ idx[,2], ]
    W1W2    = .rowSums( W1*W1, nrow(W1), ncol(W1) ) * .rowSums( W2*W2, nrow(W2), ncol(W2) )     #;  print(range(W1W2))
    
    normal2 = .rowSums( normal*normal, nrow(normal), ncol(normal) )
    
    area    = sqrt(normal2)
    angle   = asin( pmin( sqrt( normal2 / W1W2 ), 1 ) )    #   sin(angle) = sqrt(normal2 / W1W2)
        
    if( FALSE )
        {
        #   set rows close to 0 to all NAs    
        bad     = normal2 <= tol2*tol2 * W1W2     #; print( bad )
        if( all(bad) )
            {
            log.string( ERROR, "matrix W does not have rank 2, with relative tol2=%g.", tol2 )
            return(NULL)
            }
            
        if( any(bad) )
            {    
            #mess    = sprintf( "%d bad face normals, out of %d.\n", sum(bad), length(bad) )
            #cat(mess)
            #   log.string( INFO, "%d normals flagged as too small, out of %d.", sum(bad), length(bad) )
            normal[ bad, ]  = NA_real_        
            #   angle[ bad ]    = NA_real_            
            }
        }
        
        
    #   now unitize
    normal  = normal / area    # sqrt(normal2) is replicated to all 3 columns

    if( perf )
        {
        cat( sprintf( "[%g]\n", gettime()-time_start ) )
        cat( "collinearity tolerance: ", tol1, '\n' )
        i   = which.min(angle)
        cat( "minimum angle between condensed generators: ", angle[i], '\n' )
        Wmin    = Wcond[ idx[i, ],  ]
        rownames(Wmin)  = idx[i, ]
        print( Wmin )
        }
    

    out = list()
    
    out$W       = W
    
    out$Wcond       = Wcond
    out$groupidx    = res$groupidx
    out$group       = res$group
    
    out$center  = 0.5 * .colSums( Wcond, nrow(Wcond), ncol(Wcond) )
    
    out$nonnegative = all( 0 <= Wcond )
    
    out$face        = data.frame( row.names=1:nrow(idx) )
    out$face$idx    = idx
    out$face$area   = area
    out$face$angle  = angle    
    out$face$normal = normal

    
    if( ! is.null(tol2)  &&  is.finite(tol2)  &&  0<=tol2 )
        {
        if( perf )
            {
            cat( sprintf( "clustering %d face normals...", n*(n-1)/2 ) )
            flush.console()            
            time_start = gettime()
            }

        
        #   cluster the unit normal vectors  [reminds me of something similar at Link]
        res = findRowClusters( out$face$normal, tol=tol2, projective=TRUE )
        if( is.null(res) )  return(NULL)
        
        out$cluster     = res$cluster       #; print( out$cluster )        
        out$clusteridx  = res$clusteridx    
        out$spread      = res$spread
        
        #   out$clusteridx[bad] = NA_integer_        #; print( out$clusteridx )
        
        if( all( out$clusteridx==1  |  is.na(out$clusteridx) ) )
            {
            log.string( ERROR, "matrix W does not have rank 3." )
            return(NULL)
            }
        
        if( 0 < length(out$cluster) )
            {
            #   normal vectors within each cluster may differ in sign,
            #   so make them consistent !
            out$frame3x2    = vector( length(out$cluster), mode='list' )
            
            for( k in 1:length(out$cluster) )
                {
                cluster = out$cluster[[k]]
                normalk = out$face$normal[cluster[1], ]            # take the first face's normal to represent all of them
                
                s   = as.numeric( out$face$normal[cluster, ] %*% normalk )  #; print(s)
                
                #   s should already be +1 or -1, but may be a little off because of truncation, so use sign(s)
                out$face$normal[cluster, ]  = sign(s) * out$face$normal[cluster, ]
                
                out$frame3x2[[k]]   = frame3x2fun( normalk )
                }
            }
            
        if( perf )
            {
            cat( sprintf( "[%g]\n", gettime()-time_start ) )    
            
            clusteridx  = out$clusteridx
            
            dumpclusters( cbind( out$face, clusteridx ) )
            }
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
        
    functional  = tcrossprod( Wcond, out$face$normal )  # output size ~ n^3 / 2

    if( FALSE  &&   0 < length(out$cluster) )
        {
        #   use functional[,] to check clusters        
        W2  = Wcond*Wcond
        zeromask    = ( .rowSums( W2, nrow(W2), ncol(W2) ) == 0 )
            
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

    out$face$center = crossprod( face_center, Wcond )
    
    if( perf ) cat( sprintf( "[%g]\n", gettime()-time_start ) )        
    
    if( 0 < length(out$cluster) )
        {
        #   if a p-face is in a cluster of p-faces with the same normal vector, 
        #   then face$center is probably not correct.
        #   Such a p-face is just a part of a compound face.
        #   Fix these by transforming to 2D, computing zonogon3, and then back to 3D again.
        #   Store the zonogon3 for later use when raytracing.
        
        clusters    = length(out$cluster) 
        
        if( perf )
            {
            cat( sprintf( "partitioning %d compound faces...", clusters ) )
            flush.console()            
            time_start = gettime()
            }
            
            
        #   make 2D lookup table from zonohedron p-face indexes, to zonohedron p-face number
        lookup2D    = matrix( 0L, nrow(Wcond), nrow(Wcond) )
        lookup2D[ out$face$idx ]    =  1:nrow(out$face$idx)
            
        out$zonogon3    = vector( clusters, mode='list' )
        out$center3D    = matrix( NA_real_, clusters, 3 )
        out$lookupgen   = vector( clusters, mode='list' )   # lookup zonogon3 generator -> zonohedron generator
        out$lookupface  = vector( clusters, mode='list' )   # lookup zonogon3 faceidx -> zonohedron faceidx
        
        
        for( k in 1:clusters )
            {
            cluster     = out$cluster[[k]]
            
            #   take the first one's normal to represent all of them
            #   we have already made them consistent in sign
            normalk     = out$face$normal[ cluster[1], ]            
            #   frame3x2    = base::svd( normalk, nu=3 )$u[ , 2:3 ]   #; print(frame3x2)
            frame3x2    = out$frame3x2[[k]]
            
            idx = out$face$idx[ cluster, ]
            generatorvec    = makeconsecutive( unique( as.integer(idx) ), n ) #; print(generatorvec)
            
            #   find center of the compound face in 3D.  center3D
            functional  = as.numeric(out$Wcond %*% normalk)     #; print(functional)
            functional[ generatorvec ] = 0    # exactly
            face_center = 0.5 * sign(functional)            #; print(face_center) # in the unit cube
            out$center3D[k, ]   = crossprod( face_center, out$Wcond )   #; print(out$center3D[k, ])  # in the centered zonohedron

            
            #  project the generators from 3D to 2D
            W2      = out$Wcond[ generatorvec, ] %*% frame3x2       #; print(W2)
            zonogon = zonogon3( W2 )              
            if( is.null(zonogon) ) return(NULL)

            out$zonogon3[[k]]   = zonogon
            
            out$lookupgen[[k]]  = generatorvec      # lookup from zonogon generator index to zonohedron generator index (row number of Wcond)

            df  = getcenters( zonogon )
            
            #   facecenter2D is in the plane and not centered
            facecenter2D    = df$center

            # transform parallelogram centers from 2D to 3D
            pgrams  = nrow(df)
            
            facecenter3D    = facecenter2D %*% t(frame3x2)  +  matrix( out$center3D[k, ], pgrams, 3, byrow=TRUE ) 
            #print(facecenter3D)  # in the centered zonohedron for the compound face

            
            #   center2 contains pgram centers, in the *centered* zonogon
            #center2 = zono$parallelogram$center - matrix( zono$center, pgrams, 2, byrow=TRUE )
            
            #   center3 contains the pgram centers, in the centered zonohedron
            #center3 = tcrossprod( center2, frame3x2 )  +  matrix( centerface, pgrams, 3, byrow=TRUE )
            #  print( center3 )
            
            
            #   copy the centers from center3 to face$center, while tracking the index pairs
            lookupface  = integer(pgrams)
            
            for( kk in 1:pgrams )
                {
                #   get the zonogon indexes for this p-gram
                idx = df$idx[kk, ]
                
                #   map from zonogon to zonohedron indexes
                idx = sort( generatorvec[idx] ) #; print(idx)
                
                #   find the p-face index in the zonohedron!
                #   i   = which( out$face$idx[ ,1] == idx[1]   &   out$face$idx[ ,2] == idx[2] )    this is slower
                i   = lookup2D[ idx[1], idx[2] ]        # this is faster
                
                if( i==0  ||  !(i %in% cluster) )
                    {
                    # something is wrong
                    print( generatorvec )
                    print( i )
                    print( cluster )
                    print( idx )
                    log.string( FATAL, "no match for idx=%d,%d.  Try reducing tol2=%g to a smaller value", 
                                        idx[1], idx[2], tol2 )                    
                    return(NULL)
                    }
                                     
                    
                out$face$center[i, ]    = facecenter3D[kk, ]
                
                lookupface[kk]  = i
                }
                
            out$lookupface[[k]] = lookupface
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
    betamin = min( out$face$beta  )
    if( betamin <= 0 )
        {
        log.string( FATAL, "Internal Error.  min(beta)=%g <= 0.", betamin )
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
    
    
#   a list with
#       area    total surface area = sum of the area of all the parallelograms
#       volume  total volume  = sum of the volume of all the pyramids, with vertex at the center of the zonohedron

metrics.zonohedron <- function( x )
    {
    out = list()
    
    out$area    = 2*sum( x$face$area )
    out$volume  = (2/3) * sum( x$face$area * x$face$beta )
    
    return( out )
    }
    
    
#   x           a zonohedron object
#   base        a numeric vector of length 3, the basepoint of all the rays
#               base must be in the interior of x,
#               or if x is non-negative, base can also be the black or white point on the boundary(x)
#   direction   an Nx3 matrix with non-zero directions in the rows
#
#   value   a dataframe N rows, and these columns
#           base        given basepoint of the ray (all the same)
#           direction   given direction of the ray
#           faceidx     of the parallelogram face where ray exits the zonohedron. The row number of x$face
#           tmax        ray parameter of intersection with face
#           boundary    the intersection with face
#           alpha       2 coordinates of boundary point in the parallelogram coords.  both in [0,1]
#           clusteridx  index of compound face cluster containg boundary intersection, or 0 if a simple face
#           faces       1 if face faceidx is a simple face in no cluster, and the number of p-grams in the compound face otherwise
#           timetrace   the time it took to compute intersection and associated stuff

raytrace.zonohedron <- function( x, base, direction )
    {
    #cat( "raytrace()", 'base=', base, 'direction=', direction, '\n' )
    
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
    faceidx     = rep(NA_integer_,n)    
    sign        = rep(NA_integer_,n)        
    boundary    = matrix(NA_real_,n,3)
    alpha       = matrix(NA_real_,n,2)
    faces       = rep(NA_integer_,n)
    clusteridx  = rep(NA_integer_,n)
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
            bound   = abs(numerator) < rtol * sum(abs(x$center))    #; print( sum(is.na(bound)) )
            
            bound[ is.na(bound) ] = FALSE
            #numerator[bound]  = 0

            tvec[bound]  = ifelse( 0 < hv[bound], 0, Inf )    # 0 here will later make the corresponding boundary = NA 
            }

        tvec[ ! is.finite(tvec) ]   = Inf
        
        j   = which.min( tvec )     # this ignores Infs
        
        tmax[k]     = tvec[j]       # tmax[k] is not negative, and might be 0 if base is black or white

        # idx = x$face$idx[j, ] 
        
        if( tmax[k] <= 0 )  next    # failed to intersect properly
        

        #   cat( "--------------", v, "---------------------------\n" )
        
        optcentered     = gcentered  +  tmax[k] * v             #;  print( optcentered )         
        boundary[k, ]   = optcentered + x$center
        
        if( is.null(x$clusteridx) )
            cidx = 0
        else
            cidx =  x$clusteridx[j]    
        
        if( is.na(cidx) )   next    # should not happen
        
        clusteridx[k]   = cidx
        
        theSign = sign(hv[j])     
        
        if( cidx == 0 )
            {
            #   ray intersects a simple p-face, a singleton cluster of normals
            faces[k]    = 1
            
            faceidx[k]  = j
        
            facecenter      = theSign * x$face$center[j, ]          #;  print( facecenter ) 
            
            edges           = t( x$Wcond[ x$face$idx[j, ], ] )         # 2 edges of the parallelogram, as 3x2 matrix
            M               = cbind( edges, x$face$normal[j, ] )    #; print( M )    # M is 3x3
            y               = base::solve( M, optcentered - facecenter )    #; print(y)
            
            #   test for inside parallelogram, with a tolerance
            if( all(abs(y[1:2]) <= 0.5 + 5.e-7 ) )
                {
                #   found it
                #   idx[k, ]    = x$face$idx[j, ]   # override above
                sign[k]     = theSign                
                alpha[k, ]  = pmin( pmax( y[1:2] + 0.5, 0), 1 )     # translate from [-0.5,0.5] to [0,1] and clamp
                }
            }
        else
            {
            #   ray intersects a compound face, a non-trivial cluster of normals
            
            cluster = x$cluster[[ cidx ]]  
            
            faces[k]    = length(cluster)   # length(cluster) is always > 1
            sign[k]     = theSign
                
            #   cat( sprintf( "Searching in cluster %d, with %d subfaces.  Not working yet.\n", cidx, length(cluster) ) )

            # get the zonogon3 and trace that instead
            
            #   project optcentered to 2D
            zono3   = x$zonogon3[[ cidx ]]    #; plot( zono3 )
            
            frame3x2    = x$frame3x2[[ cidx ]]    #   frame3x2fun( x$face$normal[j, ] )

            center3D    = theSign * x$center3D[cidx, ]
            opt2D   = (optcentered - center3D) %*% frame3x2  +  zono3$center  #; print(opt2D)

            res = raytrace2( zono3$zonohedron, c(opt2D,0), c(0,0,1) )
            # print( res )
            # print( res$source )
            
            faceidx[k]  = x$lookupface[[cidx]][ res$faceidx ]
            
            alpha[k, ]  = res$alpha
            
            #cat( "faceidx1=", j, "   normal1=", x$face$normal[j, ], '\n' )           
            #cat( "faceidx2=", faceidx[k],  "    normal2=", x$face$normal[faceidx[k], ], '\n' )
            #cat( "faceidx=", faceidx[k], "   cidx=", cidx, "  theSign=", theSign, "  frame3x2=", frame3x2, '\n' )

            #return(NULL)
            
            }

        timetrace[k]    = gettime() - time_start
        }
    
    out = data.frame( row.names=1:n )
    
    out$base        = matrix( base, n, 3, byrow=TRUE )  # replicate base to all rows
    out$direction   = direction
    out$faceidx     = faceidx    
    out$sign        = sign
    out$tmax        = tmax
    out$boundary    = boundary
    out$alpha       = alpha   
    out$clusteridx  = clusteridx
    out$faces       = faces    
    out$timetrace   = timetrace

    #if( blackwhite )    tmax[ tmax==0 ] = NA_real_      # so boundary becomes NA too !
    #out$boundary     = out$base  +  tmax * direction   # tmax is replicated to 3 columns

    cnames  = colnames(base)
    if( is.null(cnames) )   cnames = colnames(direction)    
        
    colnames(out$boundary)   = cnames
    
    return( out )
    }
    
    
    
#   in this one, base is not required to be inside the zonohedron
#   the returned intersection is with the boundary, and is the first one along the ray
#   this one assumes no compound faces, and does *NOT* handle clusters.
#
#   value   a dataframe N rows, and these columns
#           base        given basepoint of the ray (all the same)
#           direction   given directions of the ray; usually only one of them.
#           tint        ray parameter of intersection with face
#           faceidx     of the parallelogram face where ray first intersects the zonohedron. The row number of x$face
#           idx         pair of face indexes  (redundant)
#           boundary    the intersection with face
#           alpha       2 coordinates of boundary in the parallelogram coords.  both in [0,1]
#           source      point in source cube that maps to boundary point.  It always inverts boundary.
#           delta       difference between L(source) and boundary


raytrace2.zonohedron <- function( x, base, direction )
    {
    # cat( "raytrace2()", 'base=', base, 'direction=', direction, '\n' )
        
    
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
    #blackwhite  = ifelse( x$nonnegative, all(gcentered == x$center) || all(gcentered == -x$center), FALSE )  #; print( blackwhite )

    hg  = as.numeric( x$face$normal  %*%  gcentered )      #; print( str(hg) )
    
       
    n       = nrow(direction)   
    faces   = nrow( x$face$normal )
        
    tint        = rep(NA_real_,n)           # t at intersection
    faceidx     = rep(NA_integer_,n)
    idx         = matrix(NA_integer_,n,2)
    boundary    = matrix(NA_real_,n,3)
    alpha       = matrix(NA_real_,n,2)
    source      = matrix(NA_real_, n, nrow(x$Wcond) )
    
    for( k in 1:n )
        {
        v   = direction[k, ]
        
        if( any( ! is.finite(v) ) )   next
        
        if( sum(v*v) == 0 ) next    # 0-vector
        
        hv      = x$face$normal  %*%  v
                
        signvec = sign( hv )
        
        tmin    =   (-signvec * x$face$beta  -  hg) / hv
        tmax    =   ( signvec * x$face$beta  -  hg) / hv
        
        #mask0   = (hv == 0)        
        #tmin[mask0] = -Inf
        #tmax[mask0] =  Inf
    
        jmax    = which.min( tmax )     # this ignores NAs, but not +/-Inf
        tmax    = tmax[jmax]
        
        if( tmax <= 0 )  next    # failed to intersect the ray
        
        
        jmin    = which.max( tmin )     # this ignores NAs, but not +/-Inf   
        tmin    = tmin[jmin]
        
        #   the next check tests whether the ray intersects the zonohedron.
        #   If the ray is parallel to a slab, and outside it, then either tmax=-Inf or tmin=+Inf.
        if( tmax < tmin )   next
        
        if( 0 < tmin )
            {
            #   basepoint is outside zonohedron
            tint[k] = tmin
            j       = jmin
            theSign = -signvec[j]            
            }
        else
            {
            #   basepoint is inside zonohedron
            tint[k] = tmax
            j       = jmax
            theSign = signvec[j]
            }

        faceidx[k]      = j            
        idx[k, ]        = x$face$idx[j, ] 

        theNormal       = theSign * x$face$normal[j, ]
        
        optcentered     = gcentered  +  tint[k] * v             #;  print( optcentered )         
        boundary[k, ]   = optcentered + x$center                

        facecenter      = theSign * x$face$center[j, ]          #;  print( facecenter ) 
        
        edges           = t( x$Wcond[ x$face$idx[j, ], ] )         # 2 edges of the parallelogram, as 3x2 matrix
        M               = cbind( edges, x$face$normal[j, ] )    #; print( M )    # M is 3x3
        y               = base::solve( M, optcentered - facecenter )    #; print(y)
        
        #   test for inside parallelogram, with a tolerance
        ok  = all(abs(y[1:2]) <= 0.5 + 5.e-7 ) 
        if( ! ok )  next    # something went wrong

        #   sign[k]     = theSign                
        # translate from [-0.5,0.5] to [0,1] and clamp
        alpha[k, ]  = pmin( pmax( y[1:2] + 0.5, 0), 1 )     
        
        #   assign the source point in the cube [corresponds to the reflectance spectrum]
        source[ k, ]    = (sign(x$Wcond %*% theNormal) + 1) / 2

        #   override the coordinates in the parallelogram
        source[ k, idx[k, ]  ]   = alpha[k, ]
        }
    
    out = data.frame( row.names=1:n )
    
    out$base        = matrix( base, n, 3, byrow=TRUE )  # replicate base to all rows
    out$direction   = direction
    out$tint        = tint    
    out$faceidx     = faceidx        
    out$idx         = idx    
    out$boundary    = boundary
    out$alpha       = alpha
    out$source      = source
    out$delta       = source %*% x$Wcond  -  out$boundary       # compare mapped point to boundary
    
    if( TRUE )
        {
        delta   = max( abs(out$delta) ) #; print(delta)
        if( 1.e-10 < delta )
            {
            # cat( out$delta, '\n' )
            log.string( ERROR, "Failed verification, delta=%g.", delta )
            return(NULL)
            }
        }
    
    cnames  = colnames(base)
    if( is.null(cnames) )   cnames = colnames(direction)    
        
    colnames(out$boundary)   = cnames
    
    #   print( out )
    
    return( out )    
    }
    
    

#   given a point on the boundary, return a point in the unit n-cube that maps to it, under W  
#       x           the zonogon
#       boundary    Mx3 matrix with points on the boundary of x
#                   Such points typically come in 1 of 2 ways:
#                       1) as computed by raytrace().  In this case the next variable - data - is available.
#                       2) for a cyclic zonohedron, as a 2-transition combination of the generators.  The variable - data - is NOT available.
#
#       data    optional data.frame as returned from raytrace(), which is used to speed up the inversion.
#               There are M rows, and columns used are:
#                   faceidx     index of the face in x$face
#                   sign        sign of the face = +-1
#                   alpha       coordinate along the face, in [0,1]
#                   clusteridx  index of clustered compound face, or 0 in case the faces is simple  (redundant)
#       tol     tolerance for verification, or NULL to skip verification.  Set to NULL for RELEASE.
#
#   returns a data.frame with M rows and these columns:
#       boundary    the original given matrix
#       clusteridx  index of clustered compound face, or 0 in case the faces is simple
#       source      an mxn matrix, where m=nrow(data)  n=number of rows in x$W 
#                   each row of the matrix is a point in the n-cube
#                   that maps to the boundary point on the zonohedron
#       delta       L^1-difference between the mapped point and the given boundary point.  
#                   Only present when verification is done.

invertboundary.zonohedron <- function( x, boundary, data=NULL, tol=NULL )      #1.e-9 )
    {
    #cat( "invertboundary()", 'boundary=', boundary,   '\n' )
            
    ok  = is.numeric(boundary)  &&  is.matrix(boundary)  &&  0<nrow(boundary)  &&  ncol(boundary)==3
    if( ! ok )
        {
        log.string( ERROR, "argument boundary is invalid." )
        return(NULL)
        }
        

    if( is.null(data) )
        {
        #   compute data from x and boundary
        direction   = boundary - matrix( x$center, nrow(boundary), 3, byrow=TRUE )
        
        data    = raytrace( x, x$center, direction )
        if( is.null(data) ) return(NULL)
        #print(data)   
        
        #   check tmax
        delta   = max( abs(data$tmax-1) )
        eps = 5.e-10
        if( eps < delta )
            {
            log.string( WARN, "boundary delta = %g > %g.", delta, eps )
            }

        #return(NULL)
        }
        
    #   get the number of boundary points
    m   = nrow(boundary)      
    
    
    ok  = is.data.frame(data)  &&  nrow(data)==m
    ok  = ok   &&   !is.null(data$faceidx)  &&  !is.null(data$sign)  &&  !is.null(data$alpha)
    if( ! ok )
        {
        log.string( ERROR, "argument data is invalid." )
        return(NULL)
        }
    
    n   = nrow(x$W)
    
    #   faceindex   = apply( data$idx, 1, function(pair) {  which( x$face$idx[ ,1] == pair[1]  &  x$face$idx[ ,2] == pair[2] ) } )
    
    source  = matrix( NA_real_, m, n )
    
    for( i in 1:m )
        {
        if( data$clusteridx[i] == 0 )
            {
            #   a simple p-face; there are 2 generators
            k   = data$faceidx[i]
                    
            normalk             = x$face$normal[k, ]
            
            functional          = as.numeric( x$Wcond %*% normalk )       #; print(functional)
            pair                = as.integer( x$face$idx[k, ] )
            functional[pair]    = 0    # exactly
            pcube               = 0.5 * ( data$sign[i] * sign(functional) + 1 )      #; print(face_center) # in the n-cube         
            pcube[pair]         = data$alpha[i, ]                   # overwrite 0.5 and 0.5
            }
        else
            {
            #   a compound face, so there are more than 2 generators involved
            cidx    = data$clusteridx[i]
            
            k   = data$faceidx[i]
                    
            normalk     = x$face$normal[k, ]        #;             print( str(normalk) )
            
            frame3x2    = x$frame3x2[[ cidx ]]  #; frame3x2fun( normalk )            
            
            functional          = as.numeric( x$Wcond %*% normalk )       #; print(functional)
            tuple               = x$lookupgen[[ cidx ]]  #; print( tuple )
            functional[tuple]   = 0    # exactly            
            
            zono3   = x$zonogon3[[ cidx ]]    #; plot( zono3 )
                                    
            center3D    = data$sign[i] * x$center3D[cidx, ]
            
            # cat( "faceidx=",  k, "   cidx=", cidx, "  theSign=", data$sign[i], "   normalk=", normalk,  "  frame3x2=", frame3x2, '\n' )
                                    
            
            bcentered   = boundary[i, ]  -  x$center
            
            opt2D   = (bcentered - center3D) %*% frame3x2  +  zono3$center  #; print(opt2D)

            pmat    = matrix( NA_real_, 2, n )
            #   trans   = integer(2)

            for( j in 1:2 )
                {
                if( j == 2 )
                    opt2D   = 2*zono3$center - opt2D
                    
                res = raytrace2( zono3$zonohedron, c(opt2D,0), c(0,0,1) )
                # print( res )
                # print( res$source )            
                

                pmat[j, ]         = 0.5 * ( data$sign[i] * sign(functional) + 1 )
                pmat[j,tuple]     = res$source
                
                if( j == 2 )
                    pmat[j,tuple] = 1 - pmat[j,tuple]
                    
                #   trans[j]    = counttransitions( pmat[j, ] )
                }
                
            #print( trans )
            
            #   choose the option with fewer transitions
            if( counttransitions(pmat[1, ]) < counttransitions(pmat[2, ]) )
                pcube   = pmat[1, ]
            else
                pcube   = pmat[2, ]
            }
            
        #   pcube is in the "condensed" cube for the condensed generators
        #   now expand to the original generators
        source[ i, ]    = expandcoeffs( x, pcube )
        }
        
    rnames  = rownames(boundary)
    if( is.null(rnames) )   rnames = 1:m
    
    colnames(source)    = rownames( x$W )
    #   rownames(source)    = as.character( data$idx )    
    
    out = data.frame( row.names=rnames )
    out$clusteridx  = data$clusteridx
    out$boundary    = boundary
    out$source      = source

    if( is.numeric(tol) && length(tol)==1 && is.finite(tol) && 0<=tol )
        {
        #   map from boundary of n-cube to the boundary of zonohedron
        bpoint  = source %*% x$W
        delta   = abs(bpoint - data$boundary)   #; print( delta )
        
        #   attr( source, 'delta' ) = delta
        
        out$delta   = rowSums(delta)
        
        if( tol <= max(out$delta,na.rm=TRUE) )
            log.string( WARN, "Verification failed. max(delta)=%g > %g.", max(delta), tol )
        }
           

        
    return( out )
    }
    
    
#   x       a zonohedron, or any list with W, Wcond, group, groupidx
#   pcube   point in the condensed cube, thought of as a reflectance spectrum
#
#   returns a point out in (possibly) bigger cube, so that
#       out %*% x$W  =  pcube %*% x$Wcond   (up to roundoff)
#
#   If there are non-trivial groups, there are many ways to do this,
#   but try to minimize the number of transitions.

expandcoeffs <- function( x, pcube, tol=NULL )    
    {
    if( length(pcube) != nrow(x$Wcond) )
        {
        log.string( FATAL, "mismatch %d != %d", length(pcube), nrow(x$Wcond) )
        return(NULL)
        }
   
    n   = nrow(x$W)
    if( nrow(x$Wcond) == n )    return( pcube )     # no condensation
        
   
    knext   = c( 2:length(pcube), 1 )
    kprev   = c( length(pcube), 1:(length(pcube)-1) )
        
    out = numeric( n )
    
    
    #   examine the non-trivial groups
    #   the easy way is to assign the weight of the condensed generator to all members of the group.
    #   but this way would not minimize the number of transitions, so we must work harder.
    for( k in 1:length(pcube) )
        {
        group   = x$group[[k]]
        alpha   = pcube[k]        
        
        if( length(group)==1  ||  alpha==0  ||  alpha==1 )
            {
            #   trivial cases
            out[ group ]   = alpha
            next
            }
       
        #    nontrivial splitting
        # cat( "group", k, '\n' )
        
        w       = x$Wcond[k, ]
       
        betavec  = (x$W[group, ] %*% w) / sum(w*w)
        # print( betavec )
        
        m       = length(group)        
        csum    = c(0,cumsum(betavec))  #;     print( csum )
        upper   = csum[ 2:(m+1) ]
        
        alphasplit  = numeric(m)
        
        #   choose output form to minimize the number of transitions
        if( pcube[kprev[k]]  >  pcube[knext[k]] )
            {
            #   output form is 1111...X00000...
            alphasplit[ upper < alpha ]  = 1            
            i   = findInterval( alpha, csum )            
            alphasplit[i]   = (alpha - csum[i]) / betavec[i]
            }
        else
            {
            #   output form is 0000...X11111...            
            alphasplit[ upper < 1-alpha ]  = 1            
            i   = findInterval( 1-alpha, csum )            
            alphasplit[i]   = (1-alpha - csum[i]) / betavec[i]
            alphasplit  = 1 - alphasplit
            }
            
        #   print( alphasplit )
            
        out[ group ]    = alphasplit
        }
        
    #    now examine the 0-generators.  
    #   They can be assigned any coefficient, but assign to minimize # of transitions.
    idx     = which( is.na(x$groupidx) )
    inext   = c( 2:n, 1 )
    iprev   = c( n, 1:(n-1) )
    for( i in idx )
        {
        if( out[iprev[i]]==1  ||  out[inext[i]]==1 )
            out[i]  = 1
        }
        
    if( ! is.null(tol) )
        {
        #   verify the results
        delta   = out %*% x$W  -  pcube %*% x$Wcond
        mad = max( abs(delta) ) ; cat( "mad=", mad, '\n' )
        if( tol < mad )
            {
            log.string( WARN, "%g < %g = max(abs(delta))", tol, mad )
            }
        }

   return( out )
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
    
    n   = nrow( x$Wcond )

    #   compute functional in R^n
    functional  = as.numeric( x$Wcond %*% normal )   #; print( functional ) 
    
    #   find vertex of centered zonohedron where functional is maximized
    vertex_max  = 0.5 * sign( functional )      #; print( vertex_max )  # this a vertex of the n-cube, translated by -1/2
    #   vertex  = crossprod( vertex_max, x$Wcond )
    betamax = sum( functional * vertex_max )    # the maximum of <x,normal> for x in the zonotope
    betamin = -betamax  # by symmetry
    
    # print( betamax )
    
    
    #   find points on boundary of centered zonohedron where <x,normal> is maximized and minimized
    bp_max  = as.numeric( crossprod( x$Wcond, vertex_max ) )
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
    frame3x2    = frame3x2fun(normal)   #; base::svd( normal, nu=3 )$u[ , 2:3 ]   #; print(frame3x2)
    
    
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
            out[[k]]    = list( beta=beta[k],  section= matrix( vertex_max %*% x$Wcond  + x$center, 1, 3 ) )
            next
            }
        if( abs(beta_k - betamin) < tol )
            {
            #   special case - only one point of intersection
            out[[k]]    = list( beta=beta[k],  section= matrix( -vertex_max %*% x$Wcond  + x$center, 1, 3 ) )       
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
            W2x3    = x$Wcond[ x$face$idx[imod, ], ]    #;print( W2x3 )                        
            
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
        
   
    
    print( x$face[k, ] )
    
    w   = x$Wcond[ x$face$idx[k, ], ]
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
    functional  = as.numeric( x$Wcond %*% normalk )
    names(functional)    = 1:length(functional)
    func    = functional
    func[ c(i,j) ] = 0
    vertex  = 0.5*( sign(func) + 1 )
    functional  = rbind( functional, vertex )
    print(functional)    
    
    cat( "\n" )    
    mess = sprintf( "Transitions: %d\n", counttransitions(vertex) )
    cat( mess )
    
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
    generatorvec  = makeconsecutive( unique( as.integer(idx) ), nrow(x$Wcond) )
    mess    = sprintf( "Face %d(%d,%d) is in cluster %d, with %d faces and %d generators.\n",  
                        k, i, j, cidx, length(cluster), length(generatorvec) )
    cat( '\n' )                        
    cat( mess )                        
    dumpcluster( x, cidx )
    

    
    #cat( '\n' )    
    #df  = x$face[cluster, ]
    #print( df )
    
    
    #cat( "\nSpread of normal vectors:\n" )
    #print( as.numeric( diff( apply( df$normal, 2, range, na.rm=T ) ) ) )

    #cat( "\nSegments:", segmentvec, '\n' )
    
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
    generatorvec    = makeconsecutive( unique( as.integer(idx) ), nrow(x$Wcond) )    
    mess    = sprintf( "Cluster %d has %d faces and %d generators (condensed).\n",  
                        k, length(cluster), length(generatorvec) )
    cat( '\n' )                        
    cat( mess )
    
    cat( '\n' )    
    df  = x$face[cluster, ]
    print( df, max=maxrows*ncol(df) )
    
    
    cat( "\nSpread of normal vectors:\n" )
    print( as.numeric( diff( apply( df$normal, 2, range, na.rm=T ) ) ) )

    cat( "\nGenerators (condensed):", generatorvec, '\n' )        
    
    generators.org  = sapply( generatorvec, function(u) { paste( x$group[[u]], collapse='+' ) } )
    cat( "Generators  (original):", generators.org, '\n' )        
    
    Wsub    = x$Wcond[ generatorvec, ]
    rownames(Wsub)  = sprintf( "%d [%s]", generatorvec, generators.org )
    colnames(Wsub)  = 1:3
    print( Wsub )
    
    if( ! is.null(x$center3D) )
        {
        cat( "Center:", x$center3D[k, ], '\n' ) 
        }
        
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
    frame3x2    = x$frame3x2[[k]]        #;base::svd( normalk, nu=3 )$u[ , 2:3 ]   #; print(frame3x2)
    
    idx = x$face$idx[ cluster, ]
    segmentvec  = makeconsecutive( unique( as.integer(idx) ), nrow(x$Wcond) ) #; print(segmentvec)
    
    #   find center of the compound face in 3D
    functional  = as.numeric(x$Wcond %*% normalk)       #; print(functional)
    functional[ segmentvec ] = 0    # exactly
    face_center = 0.5 * sign(functional)            #; print(face_center) # in the unit cube
    centerface  = crossprod( face_center, x$Wcond )     #; print(centerface)  # in the centered zonohedron
        
        
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
#   type    'w' for wireframe, 'p' for points drawn at the center of each p-face, 'f' for filled faces
#   both    draw both symmetric halves,  (in black and red)
#   axes    c(1,2), or 1, or 2
#   cidx    draw only one specific cluster  (disabled)
#   faceidx draw only one face (disabled)
#   labels  label the faces

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

    if( 1 <= cidx  &&  cidx <= length(x$cluster) )
        cluster = x$cluster[[cidx]]
    else
        cluster = NULL

        

    rgl::points3d( 0, 0, 0, col='black', size=10, point_antialias=TRUE )
    rgl::points3d( white[1], white[2], white[3], col='white', size=10, point_antialias=TRUE )
        
        
        
    if( grepl( 'w', type ) )
        {
        rgl::wire3d( cube, lit=FALSE )    
        
        #   exact diagonal of box        
        rgl::lines3d( c(0,white[1]), c(0,white[2]), c(0,white[3]), col=c('black','white'), lwd=3, lit=FALSE )       
        
        for( tp in c('simple','compound') )
            {
            if( tp == 'simple' )
                {
                color   = 'black'                
                if( is.null(x$clusteridx) )
                    facedata    = x$face
                else
                    facedata    = x$face[ x$clusteridx == 0, ]
                }
            else
                {
                #   compound
                color   = 'green'
                if( is.null(x$clusteridx) )
                    facedata    = NULL
                else
                    facedata    = x$face[ 0 < x$clusteridx, ]
                }
                
            if( is.null(facedata) ) next
                
            faces   = nrow(facedata)
            
            if( faces == 0 )    next
            
            step    = 2 * length(axes)  # 2 * the number of segments per face

            mat     = matrix( 0, step*faces, 3 )

            # quad    = matrix( 0, 4, 3 )
            
            for( i in 1:faces )
                {
                center  = facedata$center[i, ] 

                edge    = x$W[ facedata$idx[i, ], ]
                
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
            
            rgl::segments3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col=color )    #, front=polymode, back=polymode, col=col, lit=FALSE )
            
            if( both )
                {
                #   draw opposite half
                xyz = offset - mat
                
                if( tp == 'simple' )  color = 'red'
                
                rgl::segments3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col=color )
                }            
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
        
        
    if( grepl( 'f', type ) )
        {
        #   draw filled quads
        facedata    = x$face 
        clusteridx  = x$clusteridx
        
        faces   = nrow(facedata)
        step    = 4   
        mat     = matrix( 0, step*faces, 3 )

        colvec  = character( step*faces )
        
        for( i in 1:faces )
            {
            center  = facedata$center[i, ] 

            edge    = x$W[ facedata$idx[i, ], ]
            
            k       = step*(i-1)
            
            mat[k+1, ] = center - 0.5 * edge[1, ] - 0.5*edge[2, ]            
            mat[k+2, ] = center - 0.5 * edge[1, ] + 0.5*edge[2, ]
            mat[k+3, ] = center + 0.5 * edge[1, ] + 0.5*edge[2, ]
            mat[k+4, ] = center + 0.5 * edge[1, ] - 0.5*edge[2, ]
            
            col = 'black'
            
            if( clusteridx[i] == 0 )
                #   a simple face, a parallelogram
                col = 'green'
            else
                {
                #   a compound face
                pfaces  = length(x$cluster[[ clusteridx[i] ]])

                if( pfaces == 3 )
                    col = 'blue'
                else if( pfaces == 6 ) 
                    col = 'yellow'                
                else if( 10 <= pfaces ) 
                    col = 'red'                                    
                }
                
                
            colvec[ (k+1):(k+4) ]   = col
            }
            
        offset  = matrix( x$center, step*faces, 3, byrow=TRUE )
        xyz = offset + mat
        
        rgl::quads3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col=colvec )    #, front=polymode, back=polymode, col=col, lit=FALSE )

        if( both )
            {
            #   draw opposite half
            xyz = offset - mat

            rgl::quads3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col=colvec )
            }                    
        }
        
        
        
    if( grepl( 'p', type ) )
        {
        #   draw first half in 'black'
        offset  = matrix( x$center, nrow(x$face), 3, byrow=TRUE )
            
        xyz = x$face$center + offset
        rgl::points3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col='black', size=6, point_antialias=TRUE )
        
        if( both )
            {
            #   draw 2nd half in 'red'
            xyz = -x$face$center + offset
            rgl::points3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col='red', size=6, point_antialias=TRUE )
            }
            
        if( ! is.null( x$center3D ) )
            {
            offset  = matrix( x$center, nrow(x$center3D), 3, byrow=TRUE )            
            xyz = x$center3D + offset            
            rgl::points3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col='yellow', size=10, point_antialias=TRUE )
            
            if( both )
                {
                xyz = -x$center3D + offset            
                rgl::points3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col='yellow', size=10, point_antialias=TRUE )
                }
            }
        }

    return( invisible(TRUE) )
    }
    
    
#   x       a zonohedron    
#   normal  projection is to the plane orthogonal to this vector
#   side    of the zonohedron to plot, + means the "top side" and - means the "bottom side"
#
#   Let normal[] point toward a light source at infinity.
#   Then the "top side" is the illuminated side, and the "bottom side" is the dark side.

plotprojection.zonohedron <- function( x, normal=c(0,0,1), side=+1 )
    {
    frame3x2    = frame3x2fun( normal )    
    if( is.null(frame3x2) ) return(FALSE)
    
    center2D    = x$face$center  %*%  frame3x2
    edge2D      = x$W  %*%  frame3x2
    
    center0     = x$center  %*%  frame3x2
    
    faces   = nrow(x$face)

    step    = 4  # the number of vertices to draw per face

    xy  = matrix( NA_real_, step*faces, 2 )  
    
    signvec = side * as.numeric( sign( x$face$normal %*% normal ) )

    for( i in 1:faces )
        {
        #dot     = sum( x$face$normal[i, ] * normal )
        #s       = sign( dot )
        
        if( is.na( signvec[i] ) )  next
        
        center  = signvec[i] * center2D[i, ]  +  center0
        
        edge    = edge2D[ x$face$idx[i, ], ]
        
        k       = step*(i-1)
        
        xy[k+1, ]   = center - 0.5 * edge[1, ] - 0.5*edge[2, ]
        xy[k+2, ]   = center - 0.5 * edge[1, ] + 0.5*edge[2, ]
        xy[k+3, ]   = center + 0.5 * edge[1, ] + 0.5*edge[2, ]
        xy[k+4, ]   = center + 0.5 * edge[1, ] - 0.5*edge[2, ]
        }
        
    
    xlim    = range( xy[ ,1], na.rm=TRUE )  #; print( xlim )
    ylim    = range( xy[ ,2], na.rm=TRUE )  #; print( ylim )
    
    plot( xlim, ylim, type='n', las=1, asp=1 )
    grid( lty=1 )
    abline( h=0, v=0 )


    for( i in 1:faces )
        {
        k       = step*(i-1)
        if( is.na( xy[k+1,1] ) )  next
        
        col = ifelse( 0 < signvec[i], 'pink', 'skyblue' )
        
        polygon( xy[ (k+1):(k+4), 1]  ,  xy[ (k+1):(k+4), 2] , col=col )
        }
        
        
    return( invisible(TRUE) )
    }
    
    
    
#   x       a zonohedron
#   normal  inward normal of halfspace that contains all the generators
#           if NULL, then one is computed automatically if possible
#
#   makes a 2D plot of the polygon formed by the generators    
plotpolygon.zonohedron <- function( x, normal=c(1,1,1) )
    {
    ncond   = nrow( x$Wcond )
    
    if( is.null(normal) )
        {
        if( ! requireNamespace( 'quadprog', quietly=TRUE ) )
            {
            log.string( ERROR, "Required package 'quadprog' could not be imported."  )
            return(NULL)
            }     
        
        res = try( quadprog::solve.QP( diag(3), numeric(3), t(x$Wcond), rep(1,ncond) ), silent=TRUE )
        #   print(res)
        
        if( inherits(res,"try-error") )    
            {
            log.string( ERROR, "zonohedron x is not salient." )
            return(FALSE)
            }
        
        normal  = res$solution
        normal  = normal / sqrt( sum(normal^2) )
        
        cat( "normal=", normal, '\n' )
        }    
    
    #   find denominator
    denom   = as.numeric( x$W %*% normal )
    if( ! all(0 < denom) )
        {
        log.string( ERROR, "normal=%g,%g,%g is invalid.", normal[1], normal[2], normal[3] )
        return(FALSE)
        }
        
    #   compute vertices in space
    vert        = x$W / denom   # denom is replicated to all columns

    
    #   rotate to 2D
    vert        = vert %*% frame3x2fun(normal)  ; print(vert)

    xlim    = range(vert[ ,1])  #; print(xlim)
    ylim    = range(vert[ ,2])  #; print(ylim)
    
    plot( xlim, ylim, type='n', xlab='x', ylab='y', asp=1, las=1 )
    grid( lty=1 )
    #   abline( h=0, v=0 )
    polygon( vert[ ,1], vert[ ,2] )
    points( vert[ ,1], vert[ ,2], pch=20 )
    
    return( invisible(TRUE) )
    }

    
    
print.zonohedron  <-  function( x, maxrows=50, ... )
    {
    cat( "original generators: ", nrow(x$W), '\n' )
    cat( "zero generators:     ", sum( is.na(x$groupidx) ), '\n' )    
    
    groupidx    = x$groupidx[ is.finite(x$groupidx) ]
    redundant   = sum( duplicated(groupidx) )
    nontriv = sum( 1 < sapply( x$group, length ) )
    mess    = sprintf( "redundant generators: %d  (generating %d lines)", redundant, nontriv )
    cat( mess, '\n' )    
    
    cat( "condensed generators:", nrow(x$Wcond), '\n' )
        
    mess    = sprintf( "total faces:          %d  [%d*(%d-1)/2 p-faces. only half the total, due to central symmetry]\n",
                        nrow(x$face), nrow(x$Wcond), nrow(x$Wcond) )
    cat( mess )                        

    cat( "simple faces:        ", sum(x$clusteridx==0), ' [distinct parallelograms, only half the total, due to central symmetry]\n' )
    
    
    #mask_degen  = is.na( x$face$beta )
    #count   = sum( mask_degen )
    #   cat( "degenerate faces:    ", count, ' [angle too small, or undefined]\n' )
   
    mets    = metrics( x )

    #area_degen  = 2*sum( x$face$area[ mask_degen ] ) 
    #mess    = sprintf( "total area:           %g  (degenerate face area: %g.    %g of total area).\n",
    #                    area_total, area_degen, area_degen/area_total )
    #cat( mess )  
    cat( "total area:          ", mets$area, '\n' )    
    
    cat( "total volume:        ", mets$volume, '\n' )    

    cat( "center:              ", x$center, "    [ white:", 2*x$center, ']\n' )
              
    cat( "object size:         ", object.size(x), "(bytes)\n" )
              
              
              
    if( length(x$cluster) == 0 )
        {
        cat( "No clusters (compound faces) were computed.\n" )
        return( invisible(TRUE) )        
        }

    cat( '\n' )
    cat( "non-trivial clusters:", length(x$cluster), '   [clustered by face normals]\n' )
    
    n       = length(x$cluster) 

    faces       = integer(n)
    generators  = integer(n)
    fgmatch     = logical(n)
    area        = numeric(n)
    normal      = matrix( NA_real_, n, 3 )
    
    for( k in 1:length(x$cluster) )
        {
        cluster     = sort( x$cluster[[k]] )
        faces[k]    = length(cluster)
        
        idx = x$face$idx[ cluster, ]
        generators[k] = length( unique( as.integer(idx) ) )
        
        fgmatch[k]  = (generators[k]*(generators[k]-1)/2  ==  faces[k])
                
        area[k]     = sum( x$face$area[cluster] )   
        
        normal[k, ] = x$face$normal[ cluster[1], ]
        
        #mess1 = sprintf( "#%d (faces=%d, generators=%d, area=%g, spread=%g)", 
        #            k, length(cluster), length(generators), area_cluster, x$spread[k] )
        #mess2 = sprintf( "%d(%d,%d)", cluster, idx[ ,1], idx[ ,2] )
        #cat( "    ", mess1, "  normal=", x$face$normal[ cluster[1], ], "   ", mess2, '\n' )
        }
        
    data    = data.frame( row.names=1:n )
    data$faces          = faces   
    data$generators     = generators  
    data$fgmatch        = fgmatch
    data$area           = area  
    data$normal         = normal  
    data$spread         = x$spread
    data$center3D       = x$center3D
    
    cat( '\n' )        
    print( data, max=maxrows*ncol(data) )
    
    area_compound =  2 * sum( data$area )
    cat( '\n' )    
    cat( "area of clustered (compound faces): ", area_compound, '     ', area_compound/mets$area, 'of total area\n' )        

    cat( '\n' )
    mess = sprintf( "maximum cluster normal spread: %g  (for cluster %d).\n", 
                        max(data$spread), which.max(data$spread) )
    cat( mess )
    
    return( invisible(TRUE) )
    }    
    
    
#   data    a data.frame with a row for each p-face, and with columns
#           idx         pair of generators
#           clusteridx  cluster to which p-face belongs, or 0 if a simple p-gram
dumpclusters <- function( data )
    {
    print( str(data) )
    
    count   = sum( is.na(data$clusteridx) )
    cat( "degenerate faces:    ", count, ' [angle too small, or undefined]\n' )
   
    #   find the cluster indexes
    cidx    = sort( unique(data$clusteridx) )
    if( cidx[1] == 0 )  cidx    = cidx[-1]
    
    m   = length(cidx)
    if( m == 0 )
        {
        cat( "There are no compound clusters.\n" )
        return(TRUE)
        }
        

    faces       = integer(m)
    generators  = integer(m)
    fgmatch     = logical(m)
    
    for( k in 1:m )
        {
        mask        = data$clusteridx == cidx[k]
        faces[k]    = sum( mask, na.rm=TRUE )
        
        idx             = data$idx[ mask, ]
        generators[k]   = length( unique( as.integer(idx) ) )
                
        fgmatch[k]  = (generators[k]*(generators[k]-1)/2  ==  faces[k])
        }
        
    df  = data.frame( row.names=cidx )
        
    df$faces        = faces   
    df$generators   = generators  
    df$fgmatch      = fgmatch
    print( df )
    
    out = all( fgmatch )
    if( ! out )
        {
        k   = which( ! fgmatch )[1]
        
        cat( "cluster ", k, "is bad!\n" )
        
        mask    = data$clusteridx == cidx[k]
        
        print( data[mask, ] )
        
        idx             = data$idx[ mask, ]        
        generatorvec    = sort( unique( as.integer(idx) ) )
        mat = t( utils::combn( generatorvec, 2 ) ) ; print( mat )
        
        mask    = logical( nrow(data) )
        for( j in 1:nrow(mat) )
            {
            i   = which( data$idx[ ,1]==mat[j,1]  &  data$idx[ ,2]==mat[j,2] )
            mask[i] = TRUE
            }
        print( data[mask, ] )
        }
    
        
    return( invisible(out) )
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
    
    
frame3x2fun <- function( normal )
    {
    ok  = is.numeric(normal)  &&  length(normal)==3  &&  0<sum(abs(normal))
    if( ! ok )
        {
        log.string( ERROR, "argument normal is invalid." )
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
            log.string( ERROR, "frame3x2fun() failed orthogonal test = %g!", max(abs(test))  )
            }
        
        test    = t(out) %*% out  -  diag(2)
        if( 1.e-14 < max(abs(test)) )
            {
            log.string( ERROR, "frame3x2fun() failed product test = %g!", max(abs(test))  )
            }
        }
        
    return(out)
    }
    
isaxis <- function( vec )
    {
    n   = length(vec)
    
    return( sum(vec==0)==n-1  &&  sum(vec==1)==1 )
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
        log.string( ERROR, "W is invalid; all rows are 0." )
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
                log.string( ERROR, "Generators are invalid. One generator is a negative multiple of another one." )
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
    
#   ivec    a vector of unique integers, always taken from the set {1,...,n}
#   n       the maximum possible value appearing in ivec
#
#   return ivec in consecutive order (possibly with wraparound) if possible
#   if not not possible, then return ivec unchanged
makeconsecutive  <-  function( ivec, n )    
    {
    m       = length( ivec )
    out     = sort( ivec )
    diffvec = diff( out )
    
    #   count the number of 1s
    count   = sum( diffvec==1 )
    
    if( count == m-1 )    return(out)  # out is in consecutive order
    
    if( count==m-2  &&  out[1]==1  &&  out[m]==n )
        {
        #   this is in consecutive order, if we allow wraparound
        k   = which( 1 < diffvec )
        out = c( out[(k+1):m], 1:k )
        return( out )
        }
    
    return( ivec )
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
    
    
    
#   x       a zonohedron
#   rays    # of random rays to trace, from the center    
#   tol     inversion tolerance
testinvertRT  <-  function( x, rays, tol=1.e-9 )
    {
    set.seed(0)
    
    direction   = matrix( rnorm(rays*3), rays, 3, byrow=TRUE )
    
    res = raytrace( x, x$center, direction  )
    if( is.null(res) )  return(NULL)
    #   print( res )
    
    compounds   = sum( 0 < res$clusteridx )
    mess    = sprintf( "compound face intersections: %d of %d.", compounds, rays )
    cat( mess, '\n' )
    
    out  = invertboundary( x, res$boundary, res, tol=tol )
    if( is.null(out) )   return(NULL)
    
    #   out = data.frame( row.names=1:rays )
    #out$vertex  = vertex
    #out$delta   = rowSums( abs(attr(vertex,'delta')) )
    
    return( invisible(out) )
    }
    
#   x       a zonohedron, which must be in cyclic order
#   count   # of random points to generate on the boundary
#   tol     inversion tolerance
testinvertCyclic  <-  function( x, count, tol=1.e-9 )
    {
    set.seed(0)
    
    n   = nrow( x$W )    
    
    source  = matrix( NA_real_, count, n )

    #   make random points on the boundary
    for( i in 1:count )
        source[i, ]   = random2trans(n) 
        
    #   print( source )
        
    boundary    = source %*% x$W

    res  = invertboundary( x, boundary, NULL, tol=tol )
    if( is.null(res) )   return(NULL)
    #   print( res )
    
    count   = sum( 0 < res$clusteridx )
    mess    = sprintf( "compound faces: %d.", count )
    cat( mess, '\n' )
            
    
    #   compare source and res$source
    delta   = apply( abs(source - res$source), 1, max )
    count   = sum( tol < delta )
    mess    = sprintf( "%d violations (delta > %g).    max(delta)=%g", count, tol, max(delta) )
    cat( mess, '\n' )
    

    return( invisible(delta) )
    }
        
    
    
########    a few classic zonohedra from https://www.ics.uci.edu/~eppstein/junkyard/ukraine/ukraine.html    ######
    
#   6 generators    
truncatedoctahedron <- function( perf=FALSE )
    {
    W   = matrix( c(1,1,0, 1,-1,0, 1,0,1, 1,0,-1, 0,1,1, 0,1,-1), 6, 3, byrow=TRUE )
    
    return( zonohedron(W,perf=perf) )
    }
    
#   9 generators        
truncatedcuboctahedron  <-  function(perf=FALSE)
    {
    s2  = sqrt(2)
    W   = matrix( c(1,1,0, 1,-1,0, 1,0,1, 1,0,-1, 0,1,1, 0,1,-1, s2,0,0, 0,s2,0, 0,0,s2),  9, 3, byrow=TRUE )
    
    return( zonohedron(W,perf=perf) )
    }
    
    
#   15 generators        
TruncatedIcosidodecahedron   <-  function(  perf=FALSE )  
    {
    GoldenRatio = (1 + sqrt(5))/2
        
    W   = c( 1,GoldenRatio,GoldenRatio-1,   1,-GoldenRatio,GoldenRatio-1,
            1,-GoldenRatio,1-GoldenRatio,   1,GoldenRatio,1-GoldenRatio,
            GoldenRatio,1-GoldenRatio,      1,GoldenRatio,1-GoldenRatio,-1,
            GoldenRatio,GoldenRatio-1,-1,   GoldenRatio,GoldenRatio-1,1,
            GoldenRatio-1,1,GoldenRatio,    GoldenRatio-1,-1,-GoldenRatio,
            GoldenRatio-1,1,-GoldenRatio,   GoldenRatio-1,-1,GoldenRatio,
            2,0,0, 0,2,0, 0,0,2 )

    W   = matrix( W, ncol=3, byrow=TRUE )   #; print(W)
    
    return( zonohedron(W,perf=perf) )
    }
    
    
#   21 generators        
truncatedsmallrhombicosidodecahedron <-  function(  perf=FALSE )
    {
    GoldenRatio = (1 + sqrt(5))/2
    
    W   = c( 1,0,-GoldenRatio,  1,0,GoldenRatio,
            0,-GoldenRatio,1,   0,GoldenRatio,1,
            -GoldenRatio,1,0,   GoldenRatio,1,0,
            1,GoldenRatio,GoldenRatio-1,    1,-GoldenRatio,GoldenRatio-1,
            1,-GoldenRatio,1-GoldenRatio,   1,GoldenRatio,1-GoldenRatio,
            GoldenRatio,1-GoldenRatio,1,    GoldenRatio,1-GoldenRatio,-1,
            GoldenRatio,GoldenRatio-1,-1,   GoldenRatio,GoldenRatio-1,1,
            GoldenRatio-1,1,GoldenRatio,    GoldenRatio-1,-1,-GoldenRatio,
            GoldenRatio-1,1,-GoldenRatio,   GoldenRatio-1,-1,GoldenRatio,
            2,0,0, 0,2,0, 0,0,2 )
            
    W   = matrix( W, ncol=3, byrow=TRUE )   #; print(W)
    
    return( zonohedron(W,perf=perf) )
    }
    
    
#--------       UseMethod() calls           --------------#    
    
inside <- function( x, g )
    {
    UseMethod("inside")
    }    
    
metrics <- function( x )    
    {
    UseMethod("metrics")
    }    
    
raytrace <- function( x, base, direction )
    {
    UseMethod("raytrace")
    }    
    
raytrace2 <- function( x, base, direction )
    {
    UseMethod("raytrace2")
    }    
        
section <- function( x, normal, beta )
    {    
    UseMethod("section")
    }    
    
invertboundary <- function( x, boundary, data=NULL, tol=NULL )     #1.e-9 )    
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
    
plotprojection  <- function( x, normal=c(0,0,1), side=+1 )    
    {
    UseMethod("plotprojection")
    }    

plotpolygon  <- function( x, normal=c(1,1,1)  )    
    {
    UseMethod("plotpolygon")
    }    
    

#   already in utils.R
#gettime <- function()
#    {
#    return( microbenchmark::get_nanotime() * 1.e-9 )
#    }
    