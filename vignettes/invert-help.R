


#   .data   a list of colorSpec objects - with the same specnames() and same number of spectra
#   the first one is supposed to be the original (the "true") spectrum
#   and the rest are estimates.    length(.data) <= 4.
#
plotOriginalPlusEstimates  <-  function( .data, mfrow=c(2,3), ymax=1 )
    {
    plots_max   = prod(mfrow)    
    
    plots   = numSpectra( .data[[1]] )
    
    #   check that there is space for all plots
    if( plots_max < plots )
        {
        mess = sprintf( "First object has %d spectra, but there is only space for %d = %d*%d.",
                        plots, plots_max, mfrow[1], mfrow[2] )
        ca(mess)
        return(FALSE)
        }

    # check that all objects have the same number of spectra
    for( j in 1:length(.data) )
        {
        if( numSpectra(.data[[j]]) != plots )
            {
            mess    = sprintf( "# spectra = %d (!= %d).\n", numSpectra(.data[[j]]), plots )
            cat(mess)
            return(FALSE)
            }
        }
        
    if( length(ymax) == 1 ) ymax = rep( ymax, plots )
    
    par( mfrow=mfrow )
    par( mai=rep(0.2,4) )
    par( omi=c(0.2,0.1,0,0) )

    wave    = wavelength(.data[[1]]) 
    xlim    = range(wave)
    
    namevec = specnames(.data[[1]])
    #   colvec  = ifelse( grepl('white',namevec,ignore=T), 'black', NA ) ; print(colvec)
    
    
    if( organization(.data[[1]])=='df.row'  &&   ! is.null(.data[[1]]$SAMPLE_ID) )
        specnames(.data[[1]]) = sprintf( "%d:  %s", .data[[1]]$SAMPLE_ID, specnames(.data[[1]]) )
    
    ltype   = c( 'solid', '42', '11', '22' )    
    lwid    = c( 0.75, 1.5, 2, 1 )
    for( k in 1:plots )
        {
        #plot.default( xlim, c(0,1), type='n', las=1, xlab='', ylab='', lab=c(10,8,7), tcl=0, mgp=c(3, 0.25, 0) )
        #title( ylab=expression(omega), line=2 )
        #title( xlab='Wavelength (nm)', line=1.5 )      
        #title( main="Equalized Wavelength Reparameterization" )
        #grid( lty=1 )
        #abline( h=c(0,1) ) #, v=0 )
        #print( str(wave) )
        #print( str(theList$omega) )
        #lines( wave, as.matrix(.data$original)[ ,j] )
        
        if( grepl('white',namevec[k],ignore=T)  ||  type(.data[[1]])=='light' )
            color = 'black'
        else
            color = NULL
            
        if( is.na(ymax[k]) )
            {
            #   compute it from the data
            ymax[k] = 0
            
            for( j in 1:length(.data) )
                ymax[k] = max( ymax[k], as.matrix( subset(.data[[j]], k) ) )
                
            ymax[k] = 1.1 * ymax[k]    # add a little more
            }
            
            
        plot( subset(.data[[1]], k), color=color, main=F, ylim=c(0,ymax[k]), legend='topleft', lty=ltype[1], lwd=lwid[1] )

        if( 2 <= length(.data) )
            {
            for( j in 2:length(.data) )
                plot( subset(.data[[j]], k), color=color, add=T, legend=F, lty=ltype[j], lwd=lwid[j]  )
            }
        }
        
    title( sub="Wavelength (nm)", line=0, outer=T )
    
    par( mfrow=c(1,1) )
    
    return( invisible(TRUE) )
    }
    
    
    

#   .obj    a colorSpec with responsivity, lambda in [lambda_min.lambda_max]
#   returns list with
#       omega.from.lambda()
#       lambda.from.omega()
#       integral.from.omega[[ ]]        a list with length = channels
#   omega is in [0,1]

createReparamList <- function( .obj, .name, .alpha=1 )
    {
    #print( .name )
    #print( deparse(substitute(.obj) ) )
    
    if( ! is.colorSpec( .obj ) )
        {
        log.string( ERROR, "%s is not a valid colorSpec object", .name )
        return(NULL)
        }
        
    if( ! is.regular( .obj ) )
        {
        log.string( ERROR, "%s wavelengths are not regular.", .name )        
        return(NULL)
        }  
        
    coredata    = coredata( .obj, forcemat=TRUE )
    
    if( FALSE   &&  min(coredata) < 0 )
        {
        log.string( ERROR, "responsivity of %s is invalid; some values are negative.", .name )        
        return(NULL)
        }    

    #   check that coredata is full-rank
    singular    = svd( coredata, nu=0, nv=0 )$d
    
    #   log.string( DEBUG, "SVD time = %g sec", as.double(Sys.time()) - time_svd )
    
    thresh  = max( dim(coredata) ) * singular[1] * 2^(-52)
    rank = sum( thresh < singular )
    if( rank <  min(dim(coredata)) )
        {
        log.string( ERROR, "The responsivity matrix of %s is rank-deficient (rank=%d < %d).", 
                                .name, rank, min(dim(coredata)) )        
        return(NULL)
        }    
        
    n           = nrow(coredata)
    channels    = ncol(coredata)
    
    if( length(.alpha) == 1 )
        .alpha = rep( .alpha[1], channels )
    
    step.wl     = step.wl(.obj)
    
    lambda.center   = wavelength(.obj)
    lambda.min      = lambda.center[1]  
    lambda.max      = lambda.center[n]  
        
    #   wavelengths defining the bins
    lambda.break    = seq( lambda.min - step.wl/2, lambda.max + step.wl/2, len=(n+1) )        
        
    out = list()
           
    omega   = coredata %*% .alpha

    if( any(omega <= 0) )
        {
        log.string( ERROR, "Scanner '%s' has non-positive linear combination at some wavelengths, which is invalid.", .name )
        return(NULL)
        }
        
    omega   = c( 0, cumsum(omega) ) # n+1 of these.
    omega   = omega / omega[n+1]    # n+1 values from 0 to 1.  Not regular.  #;     print( omega )
    
    out$omega.from.lambda   = splinefun( lambda.break, omega, method='monoH.FC' )      # hyman  monoH.FC  natural
    
    out$lambda.from.omega   = splinefun( omega, lambda.break, method='monoH.FC' )      # hyman  monoH.FC  natural
           
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
        
    
compileText  <-  function( text )
    {
    return( eval( parse( text=text ) ) )    
    }    
    
    
#   .obj    responder
#   1)  the normal plot, plus the sum
#   2)  the wavelength -> omega function
#   3)  the reparameterized responses
    
plotReparam3 <- function( .obj, .alpha=1 )
    {        
    par( mfcol= c(1,3) )
    par( mai=c(0.4,0.4,0.3,0.1) )
    
    #   fake the color by changing quantity to 'energy'
    obj.copy  = .obj
    quantity(obj.copy) = 'energy'
    mat.rgb = product( obj.copy, colorSpec::BT.709.RGB, wave='auto' )  #; print( mat.rgb )
    mat.rgb = mat.rgb / max(mat.rgb)                        #; print( mat.rgb )    
    mat.rgb = DisplayRGBfromLinearRGB( mat.rgb )          #; print( mat.rgb )
    color_vec   = rgb( mat.rgb )      
    
    #   plot #1 is the normal plot of .obj, with extra spectrum sum
    ysum            = rowSums( as.matrix(.obj) )
    dim(ysum)       = c(length(ysum),1)
    colnames(ysum)  = 'sum'
       
    obj.plus    = bind( .obj, colorSpec(ysum,wave=wavelength(.obj), quant=quantity(.obj) ) )
    
    color_vec   = c( color_vec, 'black' )
    
    plot( obj.plus, main=FALSE, color=color_vec, ylab='Responsivity', legend=FALSE )
    
    mess = paste( sprintf( "bar(%s)", specnames(.obj) ), collapse=',' )
    mess = sprintf( "expression( %s, 'sum' )", mess )
    legend = compileText( mess )    # "expression( bar(x), bar(y), bar(z), 'sum' )" )
        
    legend( 'topright', legend, bty='n', col=color_vec, lty=1, lwd=10 )  #, seg.len=4 )
    
    
    mtext( "(a)", side=3, line=0.5, at=par("usr")[1], adj=1.3, cex=1 )
    
    
    theList = createReparamList( .obj, deparse(substitute(.obj)), .alpha=.alpha )    #  ; print( str(theList) )
        
    if( 0 )
    {
    #   plot #2    wavelength to omega plot

    wave    = theList$lambda.break
    plot.default( range(wave), c(0,1), type='n', las=1, xlab='', ylab='', lab=c(10,8,7), tcl=0, mgp=c(3, 0.25, 0) )
    title( ylab=expression(omega ~ '(reparameterized wavelength)'), line=2 )
    title( xlab='Wavelength (nm)', line=1.5 )      
    # title( main="Equalized Wavelength Reparameterization" )
    grid( lty=1, lwd=0.5 )
    abline( h=0, v=0 )
    #print( str(wave) )
    #print( str(theList$omega) )
    lines( wave, theList$omega, col='red'  )
    }
    
    if( 1 )
    {
    #   plot #2     omega to wavelength plot
    wave    = theList$lambda.break
    plot.default( c(0,1), range(wave),  type='n', las=1, xlab='', ylab='', lab=c(10,8,7), tcl=0, mgp=c(3, 0.25, 0) )
    title( xlab=expression(omega ~ '(reparameterized wavelength)'), line=1.5 )
    title( ylab='Wavelength (nm)', line=2 )      
    # title( main="Equalized Wavelength Reparameterization" )
    grid( lty=1, lwd=0.5 )
    abline( h=0, v=c(0,1) )
    lines( theList$omega, wave, col='red'  )
    legend( 'topleft', expression( lambda ~ '=' ~ varphi(omega) ), lty=1, bty='n', col='red', inset=c(0.1,0.1) )
    }
    
    mtext( "(b)", side=3, line=0.5, at=par("usr")[1], adj=1.3, cex=1 )
        
    
    #   plot #3 has the reparameterized responses

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
    title( xlab=expression(omega ~ '(reparameterized wavelength)' ), line=1.5 )
    title( ylab='Responsivity (reparameterized)', line=2 )                             
    grid( lty=1, lwd=0.5 )
    abline( h=0, v=0 )

    for( j in 1:spectra )
        {
        #   y   = theList$responsivity.from.omega[[j]]( omega ) ; print( range(y) )
        
        y   = theList$integral.from.omega[[j]]( omega, deriv=1 )    #; print( range(y) )
        
        lines( omega, y, col=color_vec[j] ) 
        # points( omega, y, col=color_vec[j], pch=20, cex=0.5 )     
        }
    lines( omega, ysum, col='black' )

    #if( specnames(.obj)[1] == 'x' )
    #    legend =  expression( widehat(x), widehat(y), widehat(z), 'sum' )
    #else
    #    legend =  expression( widehat(l), widehat(m), widehat(s), 'sum' )
            
    mess = paste( sprintf( "widehat(%s)", specnames(.obj) ), collapse=',' )
    mess = sprintf( "expression( %s, 'sum' )", mess )
    legend = compileText( mess )
                    
    legend( "topright", legend, col=c(color_vec,'black'), bty='n', lwd=10, inset=c(0,0.04) )
    
    
    mtext( "(c)", side=3, line=0.5, at=par("usr")[1], adj=1.3, cex=1 )
        

    par( mfcol=c(1,1) )
    
    return( invisible(TRUE) )
    }
    
        
#   RGB     linear RGB any sort of array
#           it is OK if val is a matrix, and then the return value is a matrix of the same shape
#   gamma   name of transfer function, supported are 'sRGB'
#           or numeric gamma of the display, so output is (linear)^(1/gamma)
#
#   return  first clips to [0,1], and then maps [0,1] to [0,1].
#           in case of ERROR it logs a message and returns the clipped values only
DisplayRGBfromLinearRGB <- function( RGB, gamma='sRGB' )
    {
    out = as.numeric(RGB)
    
    out = pmin( pmax( RGB, 0 ), 1 )
    
    if( is.character(gamma) )
        {
        if( tolower(gamma[1]) == 'srgb' )
            out = ifelse( out <= 0.0031308,    12.92 * out,  1.055 * out^(1/2.4) - 0.055 )
        #else
        #    log.string( ERROR, "gamma is invalid" )
        }
    else if( is.numeric(gamma) && 0 < gamma[1] )
        out = out ^ (1/gamma[1])
    else
        {
        #log.string( ERROR, "gamma is invalid" )
        }
    
    dim(out)        = dim(RGB)
    dimnames(out)   = dimnames(RGB) #; print( out )
    
    return( out )
    }    