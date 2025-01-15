    
    
        
#   .obj    a colorSpec with responsivity, lambda in [lambda_min.lambda_max]
#   returns list with
#       omega.from.lambda()
#       lambda.from.omega()
#       responsivity.from.omega[[ ]]    a list with length = channels.  No longer needed, just take deriv=1 of the next one !
#       integral.from.omega[[ ]]        a list with length = channels
#   omega is in [0,1]

makeReparamFunctionList <- function( .obj, .name, .mode='equalize', .pnorm=2 )
    {
    #print( .name )
    #print( deparse(substitute(.obj) ) )
    
    if( ! is.colorSpec( .obj ) )
        {
        log_level( ERROR, "%s is not a valid colorSpec object", .name )
        return(NULL)
        }
        
    if( ! is.regular( .obj ) )
        {
        log_level( ERROR, "%s wavelengths are not regular.", .name )        
        return(NULL)
        }  
        
    coredata    = coredata( .obj, forcemat=TRUE )
    
    if( FALSE   &&  min(coredata) < 0 )
        {
        log_level( ERROR, "responsivity of %s is invalid; some values are negative.", .name )        
        return(NULL)
        }    

    #   check that coredata is full-rank
    singular    = svd( coredata, nu=0, nv=0 )$d
    
    #   log_level( DEBUG, "SVD time = %g sec", as.double(Sys.time()) - time_svd )
    
    thresh  = max( dim(coredata) ) * singular[1] * 2^(-52)
    rank = sum( thresh < singular )
    if( rank <  min(dim(coredata)) )
        {
        log_level( ERROR, "The responsivity matrix of %s is rank-deficient (rank=%d < %d).", 
                                .name, rank, min(dim(coredata)) )        
        return(NULL)
        }    
        
    n           = nrow(coredata)
    channels    = ncol(coredata)
    
    step.wl     = step.wl(.obj)
    
    lambda.center   = wavelength(.obj)
    lambda.min      = lambda.center[1]  
    lambda.max      = lambda.center[n]  
        
    #   wavelengths defining the bins
    lambda.break    = seq( lambda.min - step.wl/2, lambda.max + step.wl/2, len=(n+1) )        
        
    out = list()
           
    if( .mode == 'equalize' )
        {
        if( .pnorm == 2 )
            omega   = sqrt( rowSums(coredata*coredata) )
        else if( .pnorm == 1 )
            omega   = rowSums( abs(coredata) )
        else
            omega   = rowSums( abs(coredata)^.pnorm )^(1/.pnorm)
        
        if( any(omega == 0) )
            {
            log_level( ERROR, "Scanner '%s' has responsivity=0 at 1 or more wavelengths, which is invalid.", .name )
            return(NULL)
            }
            
        omega   = c( 0, cumsum(omega) ) # n+1 of these.
        omega   = omega / omega[n+1]    # n+1 values from 0 to 1.  Not regular.  #;     print( omega )
        
        out$omega.from.lambda   = splinefun( lambda.break, omega, method='monoH.FC' )      # hyman  monoH.FC  natural
        
        out$lambda.from.omega   = splinefun( omega, lambda.break, method='monoH.FC' )      # hyman  monoH.FC  natural
        }
    else if( .mode == 'linear' )
        {
        omega   = (0:n) / n  # n+1 values in regular steps
        
        out$omega.from.lambda   = function( lambda )    { (lambda - lambda.break[1]) / (lambda.break[n+1] - lambda.break[1]) }
        
        out$lambda.from.omega   = function( omega )     {  (1-omega)*lambda.break[1]  +  omega*lambda.break[n+1] }
        }
    else
        {
        log_level( ERROR, ".mode='%s' is invalid.", .mode )
        return(NULL)
        }
        
    out$lambda.break    = lambda.break
    out$omega           = omega
             

    out$integral.from.omega     = list()        
    #   out$responsivity.from.omega = list()    
    
    for( j in 1:channels )
        {
        integral    = c( 0, step.wl * cumsum( coredata[ , j ] ) )   #   ; print( range(integral) )
        
        out$integral.from.omega[[j]] <- splinefun( omega, integral, method='monoH.FC' )   # hyman  monoH.FC  natural
        
        #   omega.center    = out$omega.from.lambda( lambda.center )
        #   out$responsivity.from.omega[[j]]    = splinefun( omega.center, coredata[ ,j] )  # a tiny bit of extrapolation
            
        #   the next line does not work -- only 1 new function is created  -- I think maybe it is a bug
        #   out$responsivity.from.omega[[j]]    <- function( om ) { return( out$integral.from.omega[[j]]( om, deriv=1 ) ) }
        #   print( str(out$responsivity.from.omega[[j]]) )
        }        
        
    return( out )
    }
    
    
    
#   first the normal plot, and then the reparameterized one
    
plotReparam <- function( .obj, .mode='equalize', .full=F, .pnorm=1 )
    {
    par( mfcol= c(2,1) )
    
    #   fake the color by changing quantity to energy
    obj.copy  = .obj
    quantity(obj.copy) = 'energy'
    mat.rgb = product( obj.copy, colorSpec::BT.709.RGB, wavelength='auto' )  #; print( mat.rgb )
    mat.rgb = mat.rgb / max(mat.rgb)                        #; print( mat.rgb )    
    mat.rgb = DisplayRGBfromLinearRGB( mat.rgb )            #; print( mat.rgb )
    color_vec   = rgb( mat.rgb )      
    
    #   top plot is the normal plot of .obj, with extra spectrum sum
    ysum            = rowSums( as.matrix(.obj) )
    dim(ysum)       = c(length(ysum),1)
    colnames(ysum)  = 'sum'
       
    obj.plus    = bind( .obj, colorSpec( ysum, wavelength=wavelength(.obj), quantity=quantity(.obj) ) )
    
    color_vec   = c( color_vec, 'black' )
    
    plot( obj.plus, color=color_vec )

    
    #   bottom plot has the reparameterized responses
    theList = makeReparamFunctionList( .obj, deparse(substitute(.obj)), .mode, .pnorm )    #  ; print( str(theList) )
    
    if( is.null(theList) )  
        {
        par( mfcol= c(1,1) )        
        return(NULL)
        }
    
    spectra = numSpectra( .obj )

    omega   = theList$omega     # (0:200)/200
    
    ylim    = 0
    ysum    = 0
    for( j in 1:spectra )
        {
        yvec    = theList$integral.from.omega[[j]]( omega, deriv=1 )
        ysum    = ysum + yvec        
        ylim    = range( yvec, ylim ) 
        }
    ylim    = range( ysum, ylim ) 
    #   print( range(ysum) )
        
    plot.default( c(0,1), ylim, type='n', las=1, xlab='', ylab='', lab=c(10,8,7),
                             tcl=0, mgp=c(3, 0.25, 0) )
    title( xlab=expression(omega), line=1.5 )
    title( ylab='Responsivity', line=2 )                             
    grid( lty=1 )
    abline( h=0, v=0 )

    for( j in 1:spectra )
        {
        #   y   = theList$responsivity.from.omega[[j]]( omega ) ; print( range(y) )
        
        y   = theList$integral.from.omega[[j]]( omega, deriv=1 )    #; print( range(y) )
        
        lines( omega, y, col=color_vec[j] ) 
        # points( omega, y, col=color_vec[j], pch=20, cex=0.5 )     
        }
    lines( omega, ysum, col='black' )
    
        
    legend( "topright", c(specnames(.obj),'sum'), col=c(color_vec,'black'), bty='n', lwd=9 )
    
    
    if( .full )
        {
        #   add wavelength to omega plot
        wave    = theList$lambda.break
        plot.default( range(wave), c(0,1), type='n', las=1, xlab='', ylab='', lab=c(10,8,7), tcl=0, mgp=c(3, 0.25, 0) )
        title( ylab=expression(omega), line=2 )
        title( xlab='Wavelength (nm)', line=1.5 )      
        title( main="Equalized Wavelength Reparameterization" )
        grid( lty=1 )
        abline( h=0, v=0 )
        #print( str(wave) )
        #print( str(theList$omega) )
        lines( wave, theList$omega, col='red'  )
        
        
        #   plot integral.from.omega[[]]
        ylim    = 0
        for( j in 1:spectra )
            {
            ylim    = range( theList$integral.from.omega[[j]]( omega ), ylim ) 
            }
        
        plot.default( c(0,1), ylim, type='n', las=1, xlab='', ylab='', lab=c(10,8,7),
                                 tcl=0, mgp=c(3, 0.25, 0) )
        title( xlab=expression(omega), line=1.5 )
        title( ylab='Integral of Responsivity', line=2 )                             
        grid( lty=1 )
        abline( h=0, v=0 )
        
        for( j in 1:spectra )
            {
            y   = theList$integral.from.omega[[j]]( omega )     #; print( range(y) )
            
            lines( omega, y, col=color_vec[j] )   
            }
            
        legend( "topleft", specnames(.obj), col=color_vec, bty='n', lwd=9 )
        }
    
    
    par( mfcol=c(1,1) )
    
    return( invisible(TRUE) )
    }
    
    