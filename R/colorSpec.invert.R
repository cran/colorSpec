
#   invert.colorSpec()
#
#   x           a material responder, or a light responder
#   response    a matrix of responses, with a response in each row to be processed. ncol(response) == numSpectra(x)
#   method      'centroid' or 'Hawkyard'.  The latter is only for material reflectances.
#   alpha       weighting coefficients for equalization ('centroid') or normalization ('Hawkyard').
#
#   return value:  a colorSpec object with spectra that have the given responses
invert.colorSpec <- function( x, response, method="centroid", alpha=1 )
    {   
    if( ! is.regular(x) )
        {
        log.string( ERROR, "responder x does not have regular wavelength step, which is necessary to speed and simplify calculations." )        
        return(NULL)
        }
    
    spectra = numSpectra(x)     
    
    response    = prepareNxM( response, M=spectra )
    if( is.null(response) )
        return(NULL)
        
    if( length(alpha) == 1 )
        alpha = rep(alpha[1],spectra)
        
    ok  = is.null(alpha)  ||  ( is.numeric(alpha) && length(alpha)==spectra  &&  is.null(dim(alpha)) )
    if( ! ok )
        {
        log.string( ERROR, "alpha is invalid.  It must be a numeric vector with length = %d, but length(alpha)=%d.", spectra, length(alpha) )
        return(NULL)
        }        
              
    if( type(x) == 'responsivity.material' )
        return( invertReflectance( x, response, method=method, alpha=alpha ) )
    else if( type(x) == 'responsivity.light' )
        return( invertEnergyCentroid( x, response, alpha=alpha ) )
    
    log.string( ERROR, "type(x) = '%s', but it must be 'responsivity.material' or 'responsivity.light'.", type(x) )        
    
    return(NULL)
    }
        
invertReflectance <- function( x, response, method, alpha )
    {
    methodvec   = c("centroid","hawkyard")
    idx = pmatch( tolower(method), methodvec )
    if( is.na(idx) )
        {
        log.string( ERROR, "method='%s' invalid.", method )
        return(NULL)
        }
        
    if( idx == 1 )
        return( invertReflectanceCentroid( x, response, alpha=alpha ) )
    else
        return( invertReflectanceHawkyard( x, response, alpha=alpha ) )
    }
    

invertReflectanceCentroid <- function( x, response, alpha )
    {   
    spectra = numSpectra(x)         
    
    for( p in c('rootSolve') )     # ,'microbenchmark'
        {
        if( ! requireNamespace( p, quietly=TRUE ) )
            {
            log.string( ERROR, "Required package '%s' could not be imported.",  p )
            return(NULL)
            }           
        }

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
            log.string( ERROR, "Failed because ! all( abs(f.root) < rtol*abs(root) + atol ). estim.precis = %g.",  res$estim.precis )
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

    #class(response) = "model.matrix"
    #class(tau0)     = "model.matrix"
    
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
    
    #   compute all estimated spectra in one line, and put them in mat
    mat = mat.responsivity.normalized  %*%  A_inv  %*% t(response)

    if( is.null(rownames(response)) )
        namevec = sprintf( "estimate.%d", 1:n)
    else
        namevec = sprintf( "%s.estimate", rownames(response) )
        
    colnames(mat)   = namevec
    
    n   = nrow( response ) 
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
            log.string( ERROR, "Failed because ! all( abs(f.root) < rtol*abs(root) + atol ). estim.precis = %g.",  res$estim.precis )
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
                
            
            
            
    
#--------       UseMethod() calls           --------------#            

invert <- function( x, response, method="centroid", alpha=1 )
    {
    UseMethod("invert")
    }

    
##############################  deadwood below  ##########################################    

        