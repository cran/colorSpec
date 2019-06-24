

#   returns polyhedral realization of S^2, by squares

makeTransitions2  <-  function( n )
    {
    pairmat = as.matrix(  expand.grid(1:n,1:n) )
    
    vertex  = matrix( 0L, nrow(pairmat), n )
    
    for( i in 1:nrow(pairmat) )
        {
        pair    = pairmat[i, ]
        
        v   = vertex[i, ]
        
        v[ pair[1]:pair[2] ]    = 1L
        
        if( pair[2] < pair[1] )
            {
            #   fixup
            v   = 1L - v
            v[ pair ]   = 1L  #; print(pair) ; print(v)
            }
            
        vertex[ i, ]    = v
        }
    

    
    #   removed duplicated all 1s
    mask    = ! duplicated(vertex,MARGIN=1)
    
    pairmat = pairmat[ mask,  ]
    vertex  = vertex[ mask, ]
    
    #   prepend all 0s
    pairmat = rbind( c(0,0), pairmat )
    vertex  = rbind( integer(n), vertex )
    
    #   sort by number of 1s
    perm    = order( rowSums(vertex) )
    pairmat = pairmat[ perm,  ]
    vertex  = vertex[ perm, ]
    
    #   print( cbind( pairmat, vertex ) )   
    
    
    #   now make the edges

    #   the first n edges are a special case
    edgemat0    = matrix( c( rep(1,n), 2:(n+1) ), n, 2 )
    
    edgemat = matrix( NA_integer_, 2 * (nrow(pairmat)-n-2), 2 )
    
    k   = 0
    
    for( i in 2L:(nrow(pairmat)-(n+1)) )
        {
        pair    = pairmat[i, ]
        
        if( 1 )
        {
        #   generate new pair by decrementing the 1st one
        p1  = pair[1] - 1
        if( p1 == 0 )   p1 = n
        pairnew = c( p1, pair[2] )
        idx = which( pairmat[ ,1]==pairnew[1]  &  pairmat[ ,2]==pairnew[2] )
        if( length(idx) != 1 )  
            {
            cat( "Cannot find ", pairnew, '\n' )
            return(NULL)
            }
        k   = k+1
        edgemat[ k, ]   = c(i,idx)
        }
        

        #   generate new pair by incrementing the 2nd one
        p2  = pair[2] + 1
        if( p2 == n+1 )   p2 = 1
        pairnew = c( pair[1], p2 )
        idx = which( pairmat[ ,1]==pairnew[1]  &  pairmat[ ,2]==pairnew[2] )
        if( length(idx) != 1 )  
            {
            cat( "Cannot find ", pairnew, '\n' )
            return(NULL)
            }
        k   = k+1
        edgemat[ k, ]   = c(i,idx)
        }
        
    #   the last n edges are a special case
    m   = nrow(vertex)
    edgemat1    = matrix( c( (m-n):(m-1), rep(m,n) ), n, 2 )

    out = list()
    out$vertex  = vertex
    out$edge    = rbind( edgemat0, edgemat, edgemat1 )
    
    return( invisible(out) )
    }
    
    
#   W     nx2 matrix defining a linear map from the n-cube to R^2    
#
#   the 2-transition vertices of the n-cube, plus relevevant edges a squares, form S^2
#   the function plots the image of all these edges

plotImageEdges  <- function( W )    
    {
    n   = nrow( W )
    if( n < 3 ) return(FALSE)
    
    
    zono    = zonogon( W )
    if( is.null(zono) ) return(FALSE)
    
    plot( zono, interior=FALSE )
    
    #   get the vertices and edges
    trans2  = makeTransitions2( n )
    if( is.null(trans2) )   return(FALSE)
    
    #   map the vertices
    vert2   = trans2$vertex  %*%  W
    
    #xlim    = range( vert2[ ,1] )
    #ylim    = range( vert2[ ,2] )    
    
    #plot( xlim, ylim, type='n', las=1, asp=1 )
    #grid( lty=1 )
    #abline( h=0, v=0 )
    
    points( vert2[ ,1], vert2[ ,2], pch=19 )
    
    x0  = vert2[ trans2$edge[ ,1]  ,1]
    y0  = vert2[ trans2$edge[ ,1]  ,2]
    
    x1  = vert2[ trans2$edge[ ,2]  ,1]
    y1  = vert2[ trans2$edge[ ,2]  ,2]
    
    graphics::segments( x0, y0, x1, y1 )

    return( invisible(TRUE) )
    }
    
    
randomCauchy  <-  function( n, m=1, seed=NA )
    {
    if( ! is.na(seed) ) set.seed(seed)
    
    out = matrix( NA_real_, n, m )
    x   = seq( 1/(2*n), 1 - 1/(2*n), length.out=n )
    
    for( j in 1:m )
        {
        x0  = runif( 1, min=1/8, max=7/8 )
        amp = runif( 1, min=1, max=2 )
        gam  = runif( 1, min=0.125, 0.5 )
        
        
        out[ ,j]    = amp / (1 + ((x-x0)/gam)^2 )
        }
        
    if( m == 1 )    dim(out) = NULL
        
    return( out )
    }
    