
#   invert.colorSpec()
#
#   x           a material responder, or a light responder
#   response    a matrix of responses, with a response in each row to be processed. ncol(response) == numSpectra(x)
#   method      'centroid' or 'Hawkyard' or 'TLSS'.  'Hawkyard' is only for material reflectances.
#   alpha       weighting coefficients for equalization ('centroid') or normalization ('Hawkyard').
#
#   return value:  a colorSpec object with spectra that have the given responses
invert.colorSpec <- function( x, response, method="centroid", alpha=1 )
    {   
    #   check x
    if( ! is.regular(x) )
        {
        log.string( ERROR, "responder x does not have regular wavelength step, which is necessary to speed and simplify calculations." )        
        return(NULL)
        }
        
    ok  = grepl( "^responsivity", type(x) )
    if( ! ok )
        {
        log.string( ERROR, "type(x) = '%s', but it must be 'responsivity.material' or 'responsivity.light'.", type(x) )        
        return(NULL)
        }
    
    
    #   check response
    spectra = numSpectra(x)     
    
    response    = prepareNxM( response, M=spectra )
    if( is.null(response) )
        return(NULL)
        
    #   check method
    methodvec   = c("centroid","Hawkyard","TLSS")
    
    idx = pmatch( tolower(method[1]), tolower(methodvec) )
    if( is.na(idx) )
        {
        log.string( ERROR, "method='%s' is invalid.", method )
        return(NULL)
        }
    method  = methodvec[idx]
        
        
    if( method=='centroid'  ||  method=='Hawkyard' )
        {
        #   check alpha
        if( length(alpha) == 1 )
            alpha = rep(alpha[1],spectra)
            
        ok  = is.null(alpha)  ||  ( is.numeric(alpha) && length(alpha)==spectra  &&  is.null(dim(alpha)) )
        if( ! ok )
            {
            log.string( ERROR, "alpha is invalid.  It must be a numeric vector with length = %d, but length(alpha)=%d.", spectra, length(alpha) )
            return(NULL)
            }        
        }
        
    if( method == 'centroid' )
        {
        if( type(x) == 'responsivity.material' )
            out = invertReflectanceCentroid( x, response, alpha=alpha )
        else if( type(x) == 'responsivity.light' )
            out = invertEnergyCentroid( x, response, alpha=alpha )
        }
    else if( method == 'Hawkyard' )
        {
        if( type(x) == 'responsivity.light' )
            {
            log.string( ERROR, "Method 'Hawkyard' is invalid when type(x) == 'responsivity.light'."  )
            return(NULL)
            }  
        out = invertReflectanceHawkyard( x, response, alpha=alpha )
        }
    else if( method == 'TLSS' )
        {
        out = invertTLSS( x, response )
        }
        
    if( is.null(out) )  return(NULL)
    
    count   = sum( is.na(out$estim.precis) )
    if( 0 < count )
        {
        log.string( WARN, "%d of %d responses could not be inverted.", count, nrow(response) )
        }        
        
    return( out )
    }
        
    

invertReflectanceCentroid <- function( x, response, alpha )
    {
    for( p in c('rootSolve') )     # ,'microbenchmark'
        {
        if( ! requireNamespace( p, quietly=TRUE ) )
            {
            log.string( ERROR, "Required package '%s' could not be imported.",  p )
            return(NULL)
            }           
        }

    spectra = numSpectra(x)         
            
        
    n   = nrow(response)  # number of responses to process
        
    if( is.null(rownames(response)) )
        namevec = sprintf( "estimate.%d", 1:n)
    else
        namevec = sprintf( "%s.estimate", rownames(response) )
        
    #   variable mat will hold the computed spectra
    mat = array( NA_real_, c(numWavelengths(x),n) )
    colnames(mat)   = namevec

    step    = step.wl(x)
    
    mat.responsivity            = step * as.matrix(x)  #; print( str(mat.responsivity) )
    
    mat.responsivity.equalized  = mat.responsivity
    ysum                        = 1
    
    if( ! is.null(alpha) )
        {
        ysum            = mat.responsivity %*% alpha
        if( any(ysum <= 0) )
            {
            log.string( ERROR, "The alpha-weighted sum of the %d responsivities must be all positive.", spectra )
            return(NULL)
            }

        for( j in 1:spectra )
            mat.responsivity.equalized[ ,j]    =   mat.responsivity[ ,j] / ysum
        }
        
    iters           = rep( NA_integer_, n )
    # precision       = rep( NA_real_, n )
    time            = rep( NA_real_, n )
    tau0            = array( NA_real_, c(n,spectra) )    
    estim.precis    = rep( NA_real_, n )

    #   myfun() is the function to be minimized
    myfun <- function( tau, y0 )
        {
        out = sum( Kfun( mat.responsivity.equalized %*% tau ) * ysum )
        
        return( out - sum( tau * y0 ) )
        }

    myGrad <- function( tau, y0 )
        {
        spec        = Kpfun(mat.responsivity.equalized %*% tau) #; print( str(spec) )
        dim(spec)   = NULL
        out         = spec %*%  mat.responsivity
        
        #cat( "-------------------\n" )
        #print( tau )
        #print( out  )
        #   cat( "myGrad()", y0, out-y0, '\n', file=stderr() )
      
        return( out - y0 )
        }
        
    myJac <- function( tau, y0 )
        {
        Kppterm = Kppfun(mat.responsivity.equalized %*% tau)
        dim(Kppterm) = NULL
        
        #out = matrix( 0, spectra, spectra )
        #for( j in 1:spectra )
        #    out[ ,j]    =  (Kppterm * mat.responsivity.equalized[ ,j]) %*% mat.responsivity
        
        mat = Kppterm * mat.responsivity.equalized  # this multiplies every column of mat.responsivity.equalized by the vector Kppterm
        
        out = crossprod( mat, mat.responsivity )
        
        # print(out)    this is always symmetric
        #s   = max(abs(out))
        #if( s < 1.e-8 )
        #    out = 10*(1.e-8/s) * out
        
        return(out)
        }
        
        
    atol    = 1.e-7     # 1.e-7
    rtol    = 1.e-9     # 1.e-9
        
    for( k in 1:n )
        {
        time[k] = gettime()         # microbenchmark::get_nanotime()  #as.double( Sys.time() )
        
        y0  = response[k, ]     # the k'th row of the input matrix

        if( all(y0 == 0) )
            {
            #  response = 0 is a special case
            time[k]         = gettime() - time[k]       
            estim.precis[k] = 0            
            mat[ ,k]        = 0  
            next
            }
        
        start0  = rep(0,spectra)  # c(12,16,4)
        
        if( FALSE )
            {
            #   COULD NOT GET THIS TO WORK WELL
            res = try( stats::optim( par=start0, fn=myfun, gr=myGrad, y0=y0, method="BFGS", control=list(abstol=atol,reltol=rtol,trace=3,REPORT=1) ) )
            
            cat( 'stats::optim()  res = ', str(res), '\n', file=stderr() )
                
            if( class(res) == "try-error"  || res$convergence != 0 )    
                {
                cat( 'stats::optim()  res = ', str(res), '\n', file=stderr() )
                next
                }        

            #   add vars to res so it looks like res came from multiroot()
            res$root            = res$par
            res$f.root          = myGrad( res$root, y0 )
            res$iter            = res$counts[1]
            res$estim.precis    = mean( abs(res$f.root) )
            }
        else
            {
            # res =  try( rootSolve::multiroot( myGrad, start0, rtol=rtol, atol=atol, ctol=0, verbose=F, y0=y0 ),  silent=F )
                    
            res =  try( rootSolve::multiroot( myGrad, start0, rtol=rtol, atol=atol, ctol=0, verbose=F, jacfunc=myJac, jactype='fullusr', y0=y0 ),  silent=F )
            
            if( class(res) == "try-error" )    
                {
                cat( 'rootSolve::multiroot()  res = ', str(res), '\n', file=stderr() )
                next
                }
            # print(res)
            }
            
        time[k]         = gettime() - time[k]                   
        iters[k]        = res$iter                  
        
        ok  = all( abs(res$f.root) < rtol*abs(res$root) + atol )
        if( ! ok )
            {
            cat( 'rootSolve::multiroot()  res = ', str(res), '\n', file=stderr() )               
            cat( 'res$root   = ', res$root, '\n', file=stderr() )            
            cat( 'res$f.root = ', res$f.root, '\n', file=stderr() )
            log.string( WARN, "Root-finding failed because ! all( abs(f.root) < rtol*abs(root) + atol ). estim.precis = %g.",  res$estim.precis )
            next 
            }  

        estim.precis[k] = res$estim.precis
        tau0[k, ]       = res$root
        
        #if( 1.e-6 < res$estim.precis )  # max(abs(res$f.root) ) )
        #   {
        #   mess = paste( sprintf( "%g", res$f.root ), collapse=' ' )
        #   log.string( WARN, "Root-solving error too large: %s, iters=%d, name='%s'\n", mess, res$iter, namevec[k] )
        #   next
        #   }

        mat[ ,k]    = Kpfun(mat.responsivity.equalized %*% res$root)
        }
        

    #   return(NULL)

        
    out = colorSpec( mat, wavelength=wavelength(x), quantity='reflectance', organization='df.row' )    

    extra   = data.frame( row.names=1:n )
    
    colnames(response)  = toupper( specnames(x) )
    extra$response      = response
    extra$time.msec     = 1000 * time
    extra$iters         = iters
    # extra$tau0        = tau0     nobody cares about tau0 but me
    extra$estim.precis  = estim.precis
    
    extradata(out)  = extra     # data.frame( response=response, time.msec=1000*time, iters=iters, tau0=tau0,  estim.precis=estim.precis )   #, precision=precision )

    return( out )
    }
    
        

#   response    matrix with response in each row
invertReflectanceHawkyard <- function( x, response, alpha=alpha )
    {   
    if( type(x) != 'responsivity.material' )
        {
        log.string( ERROR, "type(x) = '%s', but it must be 'responsivity.material'", type(x) )        
        return(NULL)    
        }
        
    for( p in c('microbenchmark') )     
        {
        if( ! requireNamespace( p, quietly=TRUE ) )
            {
            log.string( ERROR, "Required package '%s' could not be loaded.",  p )
            return(NULL)
            }           
        }        
    
    time_start  = gettime()         # microbenchmark::get_nanotime()       # as.double( Sys.time() )
    
    spectra = numSpectra(x)       
    
    step    = step.wl(x)    
        
    mat.responsivity            = step * as.matrix(x)    
    
    mat.responsivity.normalized = mat.responsivity

    if( ! is.null(alpha) )
        {
        ysum    = mat.responsivity %*% alpha
        if( any(ysum <= 0) )
            {
            log.string( ERROR, "The alpha-weighted sum of the %d responsivities must be all positive.", spectra )
            return(NULL)
            }

        for( j in 1:spectra )
            mat.responsivity.normalized[ ,j]    =   mat.responsivity[ ,j] / ysum
        }

    #   A is square matrix, with dimensions spectra x spectra
    A = t(mat.responsivity) %*% mat.responsivity.normalized
    
    #   check that coredata is full-rank
    singular    = svd( A, nu=0, nv=0 )$d   # ; print( singular )
    thresh  = spectra * singular[1] * 2^(-52)
    rank = sum( thresh < singular )
    if( rank <  spectra )
        {
        log.string( ERROR, "The responsivity matrix is rank-deficient (rank=%d < %d).",  rank, spectra )        
        return(NULL)
        }        
    
    A_inv   = solve(A)
    
    n   = nrow( response )    
    
    #   compute all estimated spectra in one line, and put them in mat
    mat = mat.responsivity.normalized  %*%  A_inv  %*% t(response)

    if( is.null(rownames(response)) )
        namevec = sprintf( "estimate.%d", 1:n)
    else
        namevec = sprintf( "%s.estimate", rownames(response) )
        
    colnames(mat)   = namevec
    

    feasible    = logical( n )
    for( j in 1:n )
        {
        spec    = mat[ ,j]
        
        feasible[j] = all( 0 <= spec  &  spec <= 1 )
        
        if( ! feasible[j] )
            #   fix it !
            mat[ ,j]    = pmin( pmax( spec, 0 ), 1 )
        }

    err = t(mat.responsivity) %*% mat - t(response)
    
    precision   = colMeans( abs(err) )
        
    time_elapsed    = gettime()  - time_start
    time.msec   = rep( 1000 * time_elapsed/n, n )
                
    out = colorSpec( mat, wavelength=wavelength(x), quantity='reflectance', organization='df.row' )    

    extra   = data.frame( row.names=1:n )
    
    colnames(response)  = toupper( specnames(x) )
    extra$response      = response
    extra$estim.precis  = precision
    extra$time.msec     = time.msec
    extra$clamped       = !feasible
    extradata(out)      = extra     # data.frame( time.msec=time.msec, clamped = !feasible, estim.precis=precision )

    return( out )
    }    
            
            
            
invertEnergyCentroid <- function( x, response, alpha )
    {   
    for( p in c('rootSolve') )     # ,'microbenchmark'
        {
        if( ! requireNamespace( p, quietly=TRUE ) )
            {
            log.string( ERROR, "Required package '%s' could not be imported.",  p )
            return(NULL)
            }           
        }

    #  check the responsivities
    mat.responsivity            = as.matrix(x)  #; print( str(mat.responsivity) )
    if( any( mat.responsivity < 0 ) )
        {
        log.string( ERROR, "Light responsivities are invalid because one or more values is negative." )
        return(NULL)
        }
        
    if( any( rowSums(mat.responsivity) <= 0 ) )
        {
        log.string( ERROR, "Light responsivities are invalid because their sum is not everywhere positive." )
        return(NULL)
        }
    
    spectra = numSpectra(x)     
    
    step    = step.wl(x)    
        
    mat.responsivity = step * mat.responsivity

    n   = nrow(response)  # number of responses to process
        
    if( is.null(rownames(response)) )
        namevec = sprintf( "estimate.%d", 1:n)
    else
        namevec = sprintf( "%s.estimate", rownames(response) )
        
    #   variable mat will hold the computed spectra
    mat = array( NA_real_, c(numWavelengths(x),n) )
    colnames(mat)   = namevec

    #   equalize or not ?
    mat.responsivity.equalized  = mat.responsivity

    if( ! is.null(alpha) )
        {
        ysum    = mat.responsivity %*% alpha
        if( any(ysum <= 0) )
            {
            log.string( ERROR, "The alpha-weighted sum of the %d responsivities must be everywhere positive.", spectra )
            return(NULL)
            }

        for( j in 1:spectra )
            mat.responsivity.equalized[ ,j]    =   mat.responsivity[ ,j] / ysum
        }
        
    iters           = rep( NA_integer_, n )
    # precision       = rep( NA_real_, n )
    time            = rep( NA_real_, n )
    tau0            = array( NA_real_, c(n,spectra) )    
    estim.precis    = rep( NA_real_, n )
    
    posfun  <- function(x)  {exp(x)}    #{ ifelse( 0<=x, 1+x, (1-x)^(-1) ) }
    
    dposfun <- function(x)  {exp(x)}    #{ ifelse( 0<=x, 1, (1-x)^(-2) ) }
    
    
    myGrad <- function( tau, y0 )
        {
        postau      = posfun(tau)   # now all positive
        
        lincomb     = mat.responsivity.equalized %*% postau # linear combination
                    
        spec        = 1/lincomb
        dim(spec)   = NULL
        out         = spec %*%  mat.responsivity
        
        #cat( "-------------------\n" )
        #print( tau )
        #print( out  )
        
        return( out - y0 )
        }
        
    myJac <- function( tau, y0 )
        {
        postau      = posfun(tau)   # now all positive
        
        dpostau     = dposfun(tau)   # now all positive
        
        # cat( 'tau = ', tau, '     exptau = ', exptau, '\n' )
        
        Kppterm = -(mat.responsivity.equalized %*% postau)^(-2)
        dim(Kppterm) = NULL
        
        # the next line multiplies every column of mat.responsivity.equalized by the vector Kppterm
        mat = Kppterm * mat.responsivity.equalized  
        
        #   the next line  multiplies every column of t(mat) by the vector exptau
        # mat = exptau * t(mat)
        # out = t( mat %*%  mat.responsivity )
        
        
        out = crossprod( mat.responsivity, mat )  *  matrix( dpostau, length(tau), length(tau), byrow=TRUE )   # not symmetric
        
        # print(out)   
        #s   = max(abs(out))
        #if( s < 1.e-8 )
        #    out = 10*(1.e-8/s) * out
        
        return(out)
        }
        
        
    rtol    = 1.e-9    
    atol    = 1.e-7
    
    for( k in 1:n )
        {
        time[k] = gettime()     # as.double( Sys.time() )
        
        y0  = response[k, ]     # the k'th row of the input matrix

        if( all(y0 == 0) )
            {
            #  response = 0 is a special case
            time[k]         = gettime() - time[k]       
            estim.precis[k] = 0            
            mat[ ,k]        = 0  
            next
            }        
        
        start0  = rep(0,spectra)  # c(12,16,4)
        #print( start0 )
         
        # res =  try( rootSolve::multiroot( myGrad, start0, rtol=rtol, atol=atol, ctol=0, verbose=F, y0=y0 ),  silent=F )

        res =  try( rootSolve::multiroot( myGrad, start0, rtol=rtol, atol=atol, ctol=0, verbose=F, jacfunc=myJac, jactype='fullusr', y0=y0 ),  silent=F )

        # cat( 'res = ', str(res), '\n', file=stderr() )
        
        if( class(res) == "try-error" )        
            {
            cat( 'res = ', str(res), '\n', file=stderr() )
            next
            }
           
        time[k]         = gettime() - time[k]               
        iters[k]        = res$iter

        
        ok  = all( abs(res$f.root) < rtol*abs(res$root) + atol )
        if( ! ok )
            {
            cat( 'res = ', str(res), '\n', file=stderr() )            
            cat( 'res$root   = ', res$root, '\n', file=stderr() )            
            cat( 'res$f.root = ', res$f.root, '\n', file=stderr() )
            log.string( INFO, "Root-finding failed because ! all( abs(f.root) < rtol*abs(root) + atol ). estim.precis = %g.",  res$estim.precis )
            next 
            }  
            
        estim.precis[k] = res$estim.precis   
        
        postau0         = posfun(res$root)  # now all positive                           
        tau0[k, ]       = postau0  
        
        #if( 1.e-6 < res$estim.precis )  # max(abs(res$f.root) ) )
        #   {
        #   mess = paste( sprintf( "%g", res$f.root ), collapse=' ' )
        #   mess = sprintf( "Root-solving error too large: %s, iters=%d, name='%s'\n", mess, res$iter, namevec[k] )
        #   cat( mess )
        #   next
        #   }
           
        mat[ ,k]    = 1 / (mat.responsivity.equalized %*% postau0)
        }
        

    quant   = ifelse( is.radiometric(x), 'energy', 'photons' )
    
    out = colorSpec( mat, wavelength=wavelength(x), quantity=quant, organization='df.row' )    

    #class(response) = "model.matrix"
    #class(tau0)     = "model.matrix"
    #extradata(out)  = data.frame( response=response, time.msec=1000*time, iters=iters, tau0=tau0,  estim.precis=estim.precis )   #, precision=precision )
    
    extra   = data.frame( row.names=1:n )
    
    colnames(response)  = toupper( specnames(x) )
    extra$response      = response
    extra$time.msec     = 1000 * time
    extra$iters         = iters
    # extra$tau0        = tau0     nobody cares about tau0 but me
    extra$estim.precis  = estim.precis
    
    extradata(out)  = extra     # data.frame( response=response, time.msec=1000*time, iters=iters, tau0=tau0,  estim.precis=estim.precis )   #, precision=precision )
    
    return( out )
    }
                
invertTLSS <- function( x, response )
    {  
    spectra = numSpectra(x)     # = length of the response vector = ncol(response)
    
    mat.responsivity    = as.matrix(x)
    
    step    = step.wl(x)            
    mat.responsivity = step * mat.responsivity

    n   = numWavelengths(x)
    
    m   = nrow(response)  # number of responses to process
    

    
    #   load matrix D
    D   = diag( rep(4,n) )
    D[ cbind( 2:n,1:(n-1) ) ]   = -2
    D[ cbind( 1:(n-1),2:n ) ]   = -2
    D[ cbind(c(1,n),c(1,n)) ]   = 2
    #   print( D )

    squash  = list()
    
    if( type(x) == 'responsivity.material' )
        {
        quantity    = 'reflectance'
        
        squash$sigma    <- function( z ) { (tanh(z) + 1) / 2 }
        squash$sigma1   <- function( z ) { 1/cosh(z)^2 / 2 }
        squash$sigma2   <- function( z ) { -1/cosh(z)^2 * tanh(z) }
            
        sigmainv <- function( rho ) { atanh( 2*rho - 1 ) }
        
        #   make mat0 of initial estimates 
        # mat0    = as.matrix( invertReflectanceCentroid( x, response, alpha=rep(1,spectra) ) )
        
        mat0    = as.matrix( invertReflectanceHawkyard( x, response, alpha=rep(1,spectra) ) )
        mat0    = pmin( pmax( mat0,0.01), 0.99 )
 
        mat0    = sigmainv( mat0 )        
        }
    else
        {
        #   type(x) == 'responsivity.light' )
        quantity    = ifelse( is.radiometric(x), 'energy', 'photons' )        
        
        squash$sigma    <- exp
        squash$sigma1   <- exp
        squash$sigma2   <- exp
            
        sigmainv <- log
        
        #   make mat0 of initial estimates 
        # mat0 =  array( 0, c(n,m) )   #   just use 0s

        whitelevel  = mean( colSums(mat.responsivity) )  #;  cat( "whitelevel=", whitelevel, '\n' )
        if( whitelevel <= 0 )
            {
            log.string( ERROR, "Cannot proceed because response to all 1s is non-positive. mean whitelevel=%g.", whitelevel )
            return(NULL)
            }
            
            
       mat0 = matrix( NA_real_, n, m )
       
       for( k in 1:m )
            {
            y0  = response[k, ]     # the k'th row of the input matrix
            
            meanresp    = mean( y0 )
            if( meanresp <= 0 )
                {
                log.string( WARN, "Cannot invert response %d, because mean response is non-positive.   meanresp=%g.", k, meanresp )
                next
                }            
                
            z0  = sigmainv( meanresp / whitelevel ) #; cat( "k=", k, "   z0=", z0, '\n' )
            mat0[ ,k]   = z0
            }
        }
        

    iters           = rep( NA_integer_, m )
    time            = rep( NA_real_, m )
    estim.precis    = rep( NA_real_, m )
    
    if( is.null(rownames(response)) )
        namevec = sprintf( "estimate.%d", 1:m )
    else
        namevec = sprintf( "%s.estimate", rownames(response) )
            
    
    #   variable mat will hold the computed spectra
    mat = array( NA_real_, c(n,m) )
    colnames(mat)   = namevec    
    
    for( k in 1:m )
        {
        time[k] = gettime()     # as.double( Sys.time() )
        
        y0  = response[k, ]     # the k'th row of the input matrix

        #   equality contraint function depends on y0
        #   heq <- function( z )    { as.numeric( trans(z)  %*%  mat.responsivity ) - y0 }        
        
        #   extract initial guess z0 from mat0
        z0  = mat0[ ,k]
        if( is.na(z0[1]) )
            #   y0 has already been determined to be bad, and warning issued, so skip it
            next           
            
        res = TLSS( z0, squash, t(mat.responsivity), y0 )
        
        # res = alabama::auglag( par=z0, fn=fun, gr=grad, heq=heq, heq.jac=heq.jac, control.outer=list(trace=FALSE) )
        #  res = alabama::constrOptim.nl( par=z0, fn=fun, gr=grad, heq=heq, heq.jac=heq.jac, control.outer=list(trace=FALSE) )
        #   print( str(res) )
        
        
        time[k]         = gettime() - time[k]        
        iters[k]        = ifelse( is.null(res$iterations), res$outer.iterations, res$iterations )
        
        if( ! is.null(res$convergence)  &&  res$convergence != 0 ) next    # convergence failed
        
        mat[ ,k]        = res$rho   # trans( res$par )
        
        #if( is.null(res$equal) )    res$equal = heq( res$par )
        # equal   = res$equal
            
        estim.precis[k] = mean( abs(res$equal) )
        
        #cat( "k=", k, "   iters=", iters[k], "  lambda=", res$lambda, "     counts=", res$counts, '\n' )
        }
        
    out = colorSpec( mat, wavelength=wavelength(x), quantity=quantity, organization='df.row' )    
    
    extra   = data.frame( row.names=1:m )
    
    colnames(response)  = toupper( specnames(x) )
    extra$response      = response
    extra$time.msec     = 1000 * time
    extra$iters         = iters
    extra$estim.precis  = estim.precis
    
    extradata(out)  = extra      
        
    return( out )
    }
    
    
    
#  from http://scottburns.us/wp-content/uploads/2019/05/LHTSS-text-file.txt    
#        translated from MATLAB/octave to R by Glenn Davis
#
#   This one function will do both LHTSS and LLSS, depending on the argument 'squash'
#
#   This is a non-linear minimization with non-linear constraints, solved by Lagrange multipliers.
#   The solution method is "full Newton"  (not quasi-Newton, e.g. augmented Lagrange).
#
#   The function works with the transformed variable z, and then rho = sigma(z).
#       for reflectance: z := atanh( 2*rho - 1 ), so rho = (tanh(z) + 1)/2.  LHTSS method
#                               rho is in (0,1)
#       for energy:      z := log(energy), so energy = exp(z)                LLSS method
#                                energy is in (0,Inf)
#
#   z0      the initial estimate for z, an n-vector.  n = # of wavelengths
#   squash  a list with 3 functions:
#               sigma   the squashing function itself, typically (tanh(z)+1)/2  or  exp(z)
#               sigma1  the 1st derivative of sigma
#               sigma2  the 2nd derivative of sigma
#   T       matrix defining the linear map from R^n to R^r,  an r x n matrix.  T is wide.
#   y0      the target response, an r-vector
#   ftol    solution tolerance
#
#   sigma() is called a "squashing function".
#
#
#
#   return value:   a list with items:
#                   convergence integer code, 0 means success
#                   zopt        optimal z, whether convergence succeeded or not
#                   rho         optimal computed reflectance or energy
#                   lambda      optimal lambda
#                   F1          the 1st part of the optimal F   n-vector, all should be small
#                   equal       the 2nd part of the optimal F   r-vector, all should be small
#                   iterations
#
#   a little bit of the notation was taken from:
#           OPTIMIZATION WITH CONSTRAINTS
#           2nd Edition, March 2004
#           K. Madsen, H.B. Nielsen, O. Tingleff

TLSS  <- function( z0, squash, T, y0, ftol=1.0e-8 )
    {
#   This is the Least Hyperbolic Tangent Slope Squared (LHTSS) algorithm for
#   generating a "reasonable" reflectance curve from a given sRGB color triplet.

#   Written by Scott Allen Burns, May 2019.
#   Licensed under a Creative Commons Attribution-ShareAlike 4.0 International
#   License (http://creativecommons.org/licenses/by-sa/4.0/).
#   For more information, see http://scottburns.us/reflectance-curves-from-srgb/            
#   [****  transcribed from MATLAB/octave to R by Glenn Davis  ****]

    n   = ncol(T)
    r   = nrow(T)
    
    #   D is the n x n  difference matrix for Jacobian
    #   having 4 on main diagonal and -2 on off diagonals,
    #   except first and last main diagonal are 2.
    D   = diag( rep(4,n) )
    D[ cbind( 2:n,1:(n-1) ) ]   = -2
    D[ cbind( 1:(n-1),2:n ) ]   = -2
    D[ cbind(c(1,n),c(1,n)) ]   = 2            
            
    #   initialize Newton's method
    z       = z0                # starting guess for z
    lambda  = numeric( r )      # starting Lagrange multiplier, all 0s
    maxit   = 50                # max number of iterations

    dslice  = cbind( 1:n, 1:n )
    
    #   Newton's method iteration
    convergence = 1     #  iteration limit reached
    
    for( count in 1:maxit )
        {
        d0 = squash$sigma(z)    # (tanh(z) + 1)/2
        d1 = squash$sigma1(z)   # 1/cosh(z)^2 / 2
        d2 = squash$sigma2(z)   # -1/cosh(z)^2 * tanh(z)
        
        Tlambda = as.numeric( crossprod(T,lambda) )     #   =   t(T) %*% lambda, but slightly faster
        
        F = rbind( D %*% z  +  d1*Tlambda, T %*% d0 - y0 )  # (n+r)x1 F vector

        if( all( abs(F) < ftol ) )
            {
            #   solution found !
            convergence = 0
            break
            }

        W           = D
        W[dslice]   = W[dslice]  +  d2*Tlambda
        
        Jct = d1 * t(T)     # == diag(d1) %*% t(T)  which is an R trick
        
        top = cbind( W, Jct )                   # n x (n+r) ; print( str(top) )
        bot = cbind( t(Jct), matrix(0,r,r) )    # r x (n+r) ; print( str(bot) )
        J   = rbind( top, bot  )                # (n+r)x(n+r) 
        
        delta   = try( base::solve( J, -F ) )   # solve Newton system of equations J*delta = -F.  This is the bottleneck.
        
        if( ! is.numeric(delta) )       # class(delta) == "try-error" )    
            {
            #   cat( 'base::solve()  res = ', str(delta), '\n', file=stderr() )
            convergence = 2  # J is bad,  too close to singular
            log.string( INFO, "Convergence failed because J is too close to singular." )
            break
            }
            
        if( ! all( is.finite(delta) ) )
            {
            #   cat( 'base::solve()  res = ', str(delta), '\n', file=stderr() )
            convergence = 3  # J is bad,  too close to singular
            log.string( INFO, "Convergence failed because %d of %d components of delta are not finite.", 
                                    sum(! is.finite(delta)), length(delta) )            
            break
            }

        z       = z + delta[ 1:n ]                  # update z
        lambda  = lambda + delta[ (n+1):(n+r) ]     # update lambda
        }
    
    out = list()
    out$convergence = convergence
    out$zopt        = z
    out$rho         = d0
    out$lambda      = lambda
    out$F1          = F[ 1:n ]
    #   out$F2          = F[ (n+1):(n+r) ]
    out$equal       = F[ (n+1):(n+r) ]
    out$iterations  = count

    return( out )
    }                


#--------       UseMethod() calls           --------------#            

invert <- function( x, response, method="centroid", alpha=1 )
    {
    UseMethod("invert")
    }

    
##############################  deadwood below  ##########################################    

        