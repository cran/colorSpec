
########################        K(x) = log( (exp(x) - 1)/x  )       #######################################
    
#   power series for K(x) = log( (exp(x) - 1)/x  )
#   near 0.  When x=0 the result is 0.    
pseriesK <- function( x, n=10 )
    {    
    #   even powers starting with constant 1/24, computed in Maple
    coeff   = c( 1/24, -1/2880, 1/181440, -1/9676800, 1/479001600, -691/15692092416000, 1/1046139494400,  
                    -3617/170729965486080000,  43867/91963695909076992000, -174611/16057153253965824000000 )  

    x2  = x*x
        
    K  = 0
    for( k in n:1 )
        {
        K  = coeff[k] + K*x2
        }
    out = 0.5*x   +  K*x2  # first term is linear, but all even powers after that
    
    return( out )
    }
    
    
#   this version of K(x) divides into 3 pieces    
#       time required:  small  <   neg ~ pos
Kfun <- function( x, rad=1 )
    {
    dim.saved   = dim(x)
    dim(x)      = NULL
    
    out = numeric( length(x) )
        
    pos = rad < x
    if( any(pos) )
        {
        xpos        = x[pos]
        out[pos]    = xpos + base::log1p( -exp(-xpos) ) - log(xpos) 
        }
        
    neg = x < -rad
    if( any(neg) )
        {
        xneg        = x[neg]
        out[neg]    = base::log1p( -exp(xneg) ) - log(-xneg)
        }
    
    small   = !pos & !neg
    if( any(small) )
        out[small]  = pseriesK( x[small] )
    
    dim(out)    = dim.saved
    
    return( out )
    }
    
    
           
       
########################        K'(x) = ( coth(t/2) - 2/t + 1 ) / 2       #######################################
        
    
#   power series for K'(x) = ( coth(t/2) - 2/t + 1 ) / 2
#   near 0.  When x=0 the result is 0.    
pseriesKp <- function( x, n=10 )
    {    
    #   odd powers starting at x , computed in Maple
    coeff   = c( 1/12, -1/720, 1/30240, -1/1209600, 1/47900160, -691/1307674368000, 1/74724249600,  
                    -3617/10670622842880000,  43867/5109094217170944000, -174611/802857662698291200000 )  

    x2  = x*x
        
    Kp  = 0
    for( k in n:1 )
        {
        Kp  = coeff[k] + Kp*x2
        }
    out = 0.5  +  Kp*x   
    
    return( out )
    }
    
#   this version of K'(x) divides into 3 pieces    
#       time required:  small < neg < pos   (108 < 158 < 200)
Kpfun <- function( x, rad=1 )
    {
    dim.saved   = dim(x)
    dim(x)      = NULL
    
    out = numeric( length(x) )
        
    pos = rad < x
    if( any(pos) )
        {
        out[pos]    = ( 1/tanh(x[pos]/2) - 2/x[pos] + 1 ) / 2
        }
        
    neg = x < -rad
    if( any(neg) )
        {
        expx        = exp( x[neg] )
        out[neg]    = expx/(expx - 1) - 1/x[neg]    # alternate expression is more accurate for very negative x
        }
    
    small   = !pos & !neg
    if( any(small) )
        out[small]  = pseriesKp( x[small] )
    
    dim(out)    = dim.saved
    
    return( out )
    }
    

           
           
########################        K''(x) = 1/x^2 - (1/4)*cosech(x/2)^2    (an even function)     #######################################
        
pseriesKpp <- function( x, n=11 )
    {    
    #   even powers, starting with the constant 1/12, computed in Maple
    coeff   = c( 1/12, -1/240, 1/6048, -1/172800, 1/5322240, -691/118879488000, 1/5748019200,  
                    -3617/711374856192000,  43867/300534953951232000, -174611/42255666457804800000, 77683/671480954256752640000 )        

    x2  = x*x
        
    out = 0
    for( k in n:1 )
        {
        out = coeff[k] + out*x2
        }

    return( out )
    }
    
#   this version of K''(x) divides into 2 pieces, small and large  
#       time required:  
Kppfun <- function( x, rad=1 )
    {
    dim.saved   = dim(x)
    dim(x)      = NULL
    
    out = numeric( length(x) )
        
    large   = rad < abs(x)
    if( any(large) )
        {
        xlarge      = x[large]
        out[large]  = 1/xlarge^2  -  0.25 /( sinh(xlarge/2)^2 )
        #   out[large]  = 1/xlarge^2  -  exp(xlarge) / base::expm1(xlarge)^2
        #a   = 1 / xlarge
        #b   = 0.5 / sinh(xlarge/2)
        #out[large]  = (a + b)*(a - b)
        }
       
    small   = !large
    if( any(small) )
        out[small]  = pseriesKpp( x[small] )
    
    dim(out)    = dim.saved
    
    return( out )
    }
    