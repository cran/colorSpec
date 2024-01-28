
#   .quad   4x2 matrix defining simple quadrilateral.
#           clockwise or counter-clockwise. convex or concave. 
#           Non-degenerate, except that 2 consecutive vertices can be the same.
#   .point  inside the quadrilateral
#   .val    4xn matrix of vector function values at each vertex

#   returns:
#           n-vector value of the function at .point

interpQuad  <-  function( .quad, .point, .val )
    {
    if( length(.val) == 4 )
        dim(.val) = c(4,1)
    
    #   compute edges
    quad_next  = .quad[ c(2,3,4,1), ]
    
    edge    = quad_next - .quad

    e2  = rowSums( edge*edge )
    
    #   count the # of 0s
    zeros   = sum( e2 == 0 )    
    if( 2 <= zeros )
        {
        log_string( ERROR, ".quad is vertex-degenerate." )
        return(NULL)
        }

    #   assign RHS for barycentric calculation below
    b   = c( .point, 1 )  
        
    if( zeros == 1 )
        {
        #   quadrangle degenerates to a triangle - a special case
        log_string( WARN, ".quad degenerates to a triangle." )
        del = which( e2 == 0 )
        tri = .quad[ -del, ]    ; log_object( TRACE, tri )
        
        A   = rbind( t(tri), 1 )
        lambda  = solve(A,b)    ; log_object( TRACE, lambda )
        
        if( all(-1.e-6 < lambda ) )
            {
            #   got it !
            out = lambda %*% .val[ -del, ]
            return( as.numeric(out) )
            }        
        else
            {
            log_string( ERROR, "Point %g,%g is not inside the degenerate quadrangle.", .point[1], .point[2] )
            return(NULL)        
            }
        }
    
    #   now normalize
    #   edge    = edge / matrix( sqrt(e2), 4, 2 )   ; print( edge )
    
    edge_prev  = edge[ c(4,1,2,3), ]
    
    #   cross   = numeric(4)
    
    #for( k in 1:4 )
    #    cross[k]    = edge[k,1]*edge_prev[k,2]  -  edge[k,2]*edge_prev[k,1]
        
    cross   = edge[ ,1]*edge_prev[ ,2]  -  edge[ ,2]*edge_prev[ ,1]    #; print( cross )
    
    
    #   count the # of 0s
    zeros   = sum( cross == 0 )
    
    if( 2 <= zeros )
        {
        log_string( ERROR, ".quad is edge-degenerate." )
        return(NULL)
        }
    
    if( zeros == 1 )
        {
        #   split at this single "inline vertex"
        split   = which( cross == 0 )
        }
    else
        {
        sign.cross  = sign( cross )
        sum.cross   = sum( sign.cross )
        
        if( sum.cross == 0 )
            {
            log_string( ERROR, ".quad is not a simple polygon." )
            return(NULL)
            }
        
        if( abs(sum.cross) == 4 )
            {
            #   quadrilateral is convex
            #   split along the shortest diagonal
            dist2   = numeric(2)
            dist2[1]    = sum( (.quad[1, ] - .quad[3, ])^2 )
            dist2[2]    = sum( (.quad[2, ] - .quad[4, ])^2 )
            split       = which.min( dist2 )
            }
        else if( abs(sum.cross) == 2 )
            {
            #   quadrilateral must have a single concave vertex
            split   = which( sign.cross == -sum.cross/2 )
            }
        else
            {
            log_string( FATAL, "Internal error.  sum.cross=%g.", sum.cross )
            return(NULL)
            }
        }
        
    log_string( TRACE, "Split polygon at vertex=%d.", split )

    
    #   determine which triangle contains .point, using barycentric coords

    #   1st triangle
    del = ( (split-2) %% 4 ) + 1
    tri = .quad[ -del, ]    ; log_object( TRACE, tri )
    
    A   = rbind( t(tri), 1 )
    lambda  = solve(A,b)    ; log_object( TRACE, lambda )
    
    if( all(-1.e-6 < lambda ) )
        {
        #   got it !
        out = lambda %*% .val[ -del, ]
        return( as.numeric(out) )
        }
        
        
    
    #   2nd triangle
    del = ( split %% 4 ) + 1
    tri = .quad[ -del, ]    ; log_object( TRACE, tri )
    
    A   = rbind( t(tri), 1 )
    lambda  = solve(A,b)    ; log_object( TRACE, lambda )
    
    if( all(-1.e-6 < lambda ) )
        {
        #   got it !
        out = lambda %*% .val[ -del, ]
        return( as.numeric(out) )
        }
    
    log_string( ERROR, "Point %g,%g is not inside the quadrangle.", .point[1], .point[2] )
    
    return(NULL)
    }
    
    
#   .mat        input spectra, n x m
#   .wavemin    wavelength for the 1st row of .mat
#   .wavestep   defines the rest of the wavelengths for .mat.  Must be positive
#   .wavenew    the new wavelengths
#
#   returns interpolated matrix:  length(.wavenew) x m
interpSprague  <-  function( .mat, .wavemin, .wavestep, .wavenew )
    {
    #   time_start  = as.double( Sys.time() )
    
    n   = nrow(.mat)
    wavemax = .wavemin + (n-1)*.wavestep
    
    wav = (.wavenew - .wavemin) / .wavestep 
    i   = floor( wav )  #   ; print(i - 2)      # i is an index into .mat
    h   = wav - i
    
    #   make the 6 weights of the kernel
    h2  = h^2
    #h3  = h * h2
    #h4  = h2 * h2
    ch  = 1 - h
    hch = h*ch
    #   h5  = 5*h
    h25 = 25*h
    
    weight  = matrix( 0, length(.wavenew), 6 )
    
    weight[ ,1] = hch*ch^2 * (5*h + 2)/24
    weight[ ,2] = -hch * ((h25 - 39)*h2 + 16)/24
    weight[ ,3] = ch * ( (((h25 - 38)*h -3)*h + 12)*h + 12)/12
    weight[ ,4] = h * ( (((h25 - 62)*h + 33)*h + 8)*h + 8)/12
    weight[ ,5] = hch * ( ((h25 - 36)*h - 3)*h - 2)/24
    weight[ ,6] = -hch*h2 * (5*h - 7)/24
    
    #   return( rowSums(weight) )

    
    out = matrix( NaN, length(.wavenew), ncol(.mat) )
    
    #   extend by constant on the low end, to the left
    idx = which( .wavenew <= .wavemin )
    if( 0 < length(idx) )
        out[ idx, ] = matrix( .mat[1, ], length(idx), ncol(.mat), byrow=T )
    
    #   extend by constant on the high end, to the right
    idx = which( wavemax <= .wavenew )
    if( 0 < length(idx) )
        out[ idx, ] = matrix( .mat[n, ], length(idx), ncol(.mat), byrow=T )
    
    
    #   do the middle part, 3 available points on both sides
    #   this is the main part of out[,]
    idx = which( 1<=(i-1)  &  (i+4)<=n )    #; print(idx)
    for( k in idx )
        {
        #   print( i[k] )
        out[ k, ]   = weight[ k, ] %*% .mat[ (i[k]-1):(i[k]+4),  ]
        }
        
    #   make 7-row extrapolation matrix on the low side 
    extrap  = matrix( .mat[1, ,drop=F], 2, ncol(.mat), byrow=T )    # 2 rows
    extrap  = rbind( extrap, .mat[ 1:5, ,drop=F] )                  # 5 rows
    idx = which( 0<=i  &  i<=1 )    #; print(idx)
    for( k in idx )
        {
        #   print( i[k] )
        out[ k, ]   = weight[ k, ] %*% extrap[ (i[k]+1):(i[k]+6),  ]
        }    
        
    #   make 7-row extrapolation matrix on the high side 
    extrap  = matrix( .mat[n, ,drop=F], 2, ncol(.mat), byrow=T )    # 2 rows
    extrap  = rbind(  .mat[ (n-4):n, ,drop=F ], extrap )            # 5 rows
    idx = which( (n-3)<=i  &  i<=(n-2) )    #; print(idx)
    for( k in idx )
        {
        #   print( i[k] )
        out[ k, ]   = weight[ k, ] %*% extrap[ (i[k]-n+4):(i[k]-n+9),  ]
        }    
    
    #   print( c( "interpSprague().  elapsed: ", as.double(Sys.time()) - time_start) )
        
    if( any( is.nan(out) ) )
        {
        log_string( FATAL, "Failed to fill all the %dx%d entries !", nrow(out), ncol(out) )
        return(NULL)
        }
        
    return( invisible(out) )
    }
    
    
######################      testing functions   ############################    
    
testQuad <- function()    
    {
    point   = c(1,1)
    val     = c( 1, 2, 2, 3 )
    
    time_start  = as.double( Sys.time() )
    

    cat( "-----------------   vertices #2 and #3 are the same, drop vertex #2 and interpolate over triangle\n" )
    interpQuad( matrix( c(0,0, 2,1,  2,1,  2,5), 4, 2, byrow=T ), point, val )
    
    cat( "-----------------   concave polygon, split at vertex #2\n" )
    interpQuad( matrix( c(0,0, 2,1,  3,0,  2,5), 4, 2, byrow=T ), point, val )
    
    cat( "-----------------   convex polygon, split at vertex #2\n" )
    interpQuad( matrix( c(0,0, 1,3,  3,3,  2,0), 4, 2, byrow=T ), point, val )
        
    cat( "-----------------   inline polygon, at vertex #2\n" )
    interpQuad( matrix( c(0,0, 2,0,  3,0,  1,-2), 4, 2, byrow=T ), point, val )
            
    cat( "-----------------   non-simple polygon, a bowtie\n" )
    interpQuad( matrix( c(0,0, 1,2,  3,0,  3,2), 4, 2, byrow=T ), point, val )
                
    cat( "-----------------   degenerate polygon, on y-axis, fails\n" )
    interpQuad( matrix( c(0,0, 0,2,  0,3,  0,5), 4, 2, byrow=T ), point, val )
    
    cat( "-----------------   degenerate polygon, on diagonal line, fails\n" )
    interpQuad( matrix( c(0,0,  3,3,  2,2,  5,5), 4, 2, byrow=T ) , point, val )
        
    mess    = sprintf( "Elapsed Time: %g msec\n", 1000 * ( as.double(Sys.time() - time_start) ) )
    cat(mess)
    
    }
    

    