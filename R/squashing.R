

#   evaluate log( (exp(t)-1)/t )  with series expansion near 0
#   tvec    should be real 
#   rad     threshold radius
Kfun <- function( tvec, rad=0.01 )
    {
    dim.saved   = dim(tvec)
    dim(tvec)   = NULL
    
    out     = tvec  # get the size and type right    
    
    mask    = (abs(tvec) < rad)
    if( any(mask) )
        {
        #   handle z's near 0 using power series
        s   = tvec[mask]
        s2  = s*s
        
        #   even powers
        coeff   = c( 1/24, -1/2880, 1/181440 )  
        
        K  = 0
        for( k in length(coeff):1 )
            {
            K  = coeff[k] + K*s2
            }
        out[mask]   = 0.5*s   +  K*s2  # first term is linear, but all even powers after that
        }
    
    out[ ! mask ]   = log( (exp(tvec[!mask]) - 1) / tvec[!mask] )
    
    dim(out)    = dim.saved
    
    return( out )
    }
    
    
    
#   evaluate (1/2)*coth(x/2) - 1/x + 1/2,  with series expansion near 0
#   zvec    can be real or complex, but emphasis on real
#   rad     threshold radius
Kpfun <- function( zvec, rad=0.01 )
    {
    dim.saved   = dim(zvec)
    dim(zvec)   = NULL
    
    mask    = (abs(zvec) < rad)
    
    out     = zvec  # get the size and type right
    if( any(mask) )
        {
        #   handle z's near 0 using power series
        s   = zvec[mask]
        s2  = s*s
        
        #   odd powers
        coeff   = c( 1/12, -1/720, 1/30240 ) #, -1/1209600 )  #, 1/47900160 )
        
        Kp  = 0
        for( k in length(coeff):1 )
            {
            Kp  = coeff[k] + s2*(Kp)
            }
        out[mask]   = s * Kp
        }
    
    out[ ! mask ]   = 1/(2*tanh(zvec[ ! mask ]/2)) - 1/zvec[ ! mask ]

    out = out + 1/2
    
    dim(out)    = dim.saved
    
    return( out )
    }
    
plotKpfun <- function()
    {
    par( mai=c(0.5,0.6,0.1,0.5) )
    
    plot( c(-30,30), c(0,1), type='n', xlab='', ylab='', lab=c(8,8,1), las=1, tcl=0, mgp=c(3, 0.25, 0) )    
    title( xlab='t', line=1.1 )
    title( ylab=expression( sigma(t) ), line=1.5 )
    
    grid( lty=1 )
    abline( h=c(0,1) , v=0 )

    x = seq( -35,35, by=0.1 )
    lines( x, Kpfun(x) )

    return(TRUE)
    }
    
#   evaluate 1/x^2 - (1/4)*cosech(x/2)^2 ,  with series expansion near 0
#   zvec    can be real or complex, but emphasis on real
#   rad     threshold radius
#   this function is even
Kppfun <- function( zvec, rad=0.01 )
    {
    dim.saved   = dim(zvec)
    dim(zvec)   = NULL
    
    mask    = (abs(zvec) < rad)
    
    out     = zvec  # get the size and type right
    if( any(mask) )
        {
        #   handle z's near 0 using power series
        s   = zvec[mask]
        s2  = s*s
        
        #   even powers
        coeff   = c( 1/12, -1/240, 1/6048 ) 
        
        Kpp  = 0
        for( k in length(coeff):1 )
            {
            Kpp  = coeff[k] + s2*(Kpp)
            }
        out[mask]   = Kpp
        }
    
    out[ ! mask ]   = 1/zvec[! mask]^2 - 0.25 /( sinh(zvec[ ! mask ]/2)^2 )
    
    dim(out)    = dim.saved
    
    return( out )
    }
    
    
#   this P(s) corresponds to rho(t) = t * 1_[0,1]    
Pprimefun <- function( zvec )
    {
    dim.saved   = dim(zvec)
    dim(zvec)   = NULL
    
    mask    = (zvec == 0)
    
    out     = zvec  # get the size and type right
    
    if( any(mask) )
        out[mask]     = -1/2
        
    out[ ! mask ] = - ( (zvec[!mask]-1)*exp(zvec[!mask]) + 1) # / zvec[!mask]^2
    
    dim(out)    = dim.saved
    
    return( out )
    }
    
    
#   wfun    function from interval to R^2    
#   rad     radii in the domain
#   res     # of points for the integration
#
#   We hope these curves are convex ovals
plotOvals <- function( wfun, rad=c(0.1, 1,2,3,4,5,10,20,40,80,200,400,100000), res=256 )
    {    
    #   make wmat with 2 columns
    xvec    = (1:res)/res - 1/(2*res)
    wmat    = wfun( xvec )    
    
    #   compute Gram matrix, and make it have det = 1
    G   = t(wmat) %*% wmat
    G   = G / sqrt(det(G))
    
    
    #   make points on the unit circle - a matrix with 2 columns
    theta   = (0:359) * pi / 180
    circle  = cbind( cos(theta), sin(theta) )
    
    if( 0 )
        {
        #   change circle to square
        for( i in 1:nrow(circle) )
            circle[i, ] = circle[i, ] / max( abs(circle[i, ]) )
        }
    
    if( 0 )
        {
        #   change circle to diamond
        for( i in 1:nrow(circle) )
            circle[i, ] = circle[i, ] / sum( abs(circle[i, ]) )
        }    
        
    if( 0 )
        #   change circle to ellipse
        circle  = circle %*% solve( G )
    
    oval    = circle    # just to get the size right
    
    for( k in length(rad):1 )
        {
        r   = rad[k]
        
        for( i in 1:nrow(circle) )
            {
            #   make linear combination of the 2 w's
            w   = r * ( wmat %*% circle[i, ] )
        
            #   apply soft clamp to [0,1]
            w   = Kpfun( w )
            
            #   apply wmat itself
            oval[i, ] = t(wmat) %*% w
            }
            
        oval    = oval / res
        
        if( k == length(rad) )
            {
            xlim    = range( oval[ ,1] ) / 1
            ylim    = range( oval[ ,2] ) / 1
            
            plot( xlim, ylim, type='n', las=1, xlab='w1', ylab='w2' )  #, asp=1 )
            grid( lty=1 )
            abline( h=0, v=0 )
            #   points( mean(xlim), mean(ylim) )            
            }
        
        polygon( oval[ ,1], oval[ ,2], border='black' )
        }

    return(TRUE)
    }


##------------------  gsl::hypergeo().   to be used later ?  -----------------------##    

if( 0 )
{
#   xvec        real vector, complex is OK if param==c(1,1)
#   param       parameters alpha and beta in the beta density

Kprime <- function( xvec, param=c(1,1) )
    {
    if( all( param==c(1,1) ) )
        #   uniform density
        return( Kpfun(xvec) )
        
    #   require( gsl )
    
    if( ! requireNamespace( 'gsl', quietly=TRUE ) )
        {
        log.string( ERROR, "Package '%s' could not be loaded.  It is required for gsl::hypergeo().",  'gsl' )
        return(NULL)
        }      
    
    #   use the 1F1 function - the Kummer M(a,b,s)
    a       = param[1]
    b       = sum(param)
    numer   = (a/b) * gsl::hypergeo( a+1, b+1, xvec )
    denom   = gsl::hypergeo( a, b, xvec )
    
    out     = numer / denom
    
    return( out )
    }
    
plotKprime <- function( param=c(1,1), xmax=100 )
    {
    x   = seq( -xmax, xmax, len=2^12+1 )

    y   = Kprime( x, param=param )
    
    plot( range(x), range(y,0,1), type='n', xlab='x', ylab='Kprime(x)', las=1 )
    grid( lty=1 )
    abline( h=c(0,1), lty=2 )
    points( 0, Kprime( 0, param=param ) )
    
    lines( x, y )
    
    
    return(TRUE)
    }    
}
