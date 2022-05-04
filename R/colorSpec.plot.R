
#   x          colorSpec object
#   xlim        for plot.default
#   ylim        for plot.default
#   add        add to existing plot
#   legend     apply a legend, ignored if add is TRUE
#   color      'natural' means simulated true color 
#   subset     of YYYY given by name, a single regexp, or a vector of exact strings 
    
plot.colorSpec  <- function( x, color=NULL, subset=NULL, main=TRUE, legend=TRUE, CCT=FALSE, add=FALSE, ... )
    {
    m   = numSpectra(x)
    if( m == 0 )
        {
        log.string( ERROR, "x has no spectra !" )
        return(FALSE)
        }

    theName = deparse( substitute(x) )
        
    vararg  = list(...)
    


    wave    = wavelength(x)
  
    color_vec = color
                
    if( is.null(color)  )
        {
        #   compute color from the spectra themselves
        #   note that BT.709.RGB is a global object
        if( type(x) == "material" )
            {
            #   compute accurate color for sRGB display
            mat.rgb = product( colorSpec::D65.1nm, x, colorSpec::BT.709.RGB, wavelength='auto' )
            }
        else
            {
            #   fake it by changing quantity to energy
            x.energy    = x    #   resample( x, wavelength(BT.709.RGB) )    #   we know BT.709.RGB has regular wavelength steps
            quantity(x.energy) = 'energy'
            mat.rgb = product( x.energy, colorSpec::BT.709.RGB, wave='auto' )         #; print( mat.rgb )
            
            if( all( is.na(mat.rgb) ) )
                {
                log.string( ERROR, "Cannot determine colors for plot.  Is all data NA ?" )                
                return( invisible(FALSE) )                
                }            
                
            theMax  = max( mat.rgb, na.rm=TRUE ) 

            if( 0 < theMax )
                mat.rgb = mat.rgb / theMax              #   normalize so max RGB is 1; print( mat.rgb )
            }
            
        mat.rgb = DisplayRGBfromLinearRGB( mat.rgb )    #; print( mat.rgb )
        #   mat.rgb = spacesRGB::SignalRGBfromLinearRGB( mat.rgb, space='sRGB', which='scene' )$RGB    #; print( mat.rgb )
         
        color_vec   = myrgb( mat.rgb )                  #;   print( color_vec )
        }
    else if( color_vec[1] == 'auto' )
        {
        jc <- colorRampPalette( rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2")))
        color_vec   = jc(m)
        }
    else
        {
        color_vec = rep( color_vec, len=m )     #; print(color_vec)
        }
        
    #   print( color_vec )
    
    #   add colors as extradata, so we can subset and match everything up OK
    organization(x) = 'df.row' 
    
    rhs = data.frame( color=color_vec, stringsAsFactors=FALSE ) #; print( str(rhs) )
    #print( is.data.frame(rhs) )
    
    extradata(x)    = rhs

    if( is.null(subset) )
        {
        x.sub   = x     # no change, so avoid slow subsetting
        }  
    else
        {
        x.sub = subset( x, subset )
        if( is.null(x.sub) )
            {
            log.string( ERROR, "subset is invalid." )
            return(FALSE)
            }
            
        if( numSpectra(x.sub) == 0 )
            {
            log.string( ERROR, "subset is empty; nothing to plot." )
            return(FALSE)
            }        
        }
        
    #   plot just the subset
    mat = as.matrix( x.sub )        
    
    if( all( is.na(mat) ) )
        {
        log.string( ERROR, "Cannot plot because all data in subset to be plotted is NA." )                
        return( invisible(FALSE) )                
        }              
        
    m   = ncol(mat)
            
        
        
    if( ! add )
        {
        #   get arguments needed for a new plot
        
        xlim    = vararg$xlim
        if( is.null(xlim) ) xlim=c(NA,NA)
        
        ylim    = vararg$ylim
        if( is.null(ylim) ) ylim=c(NA,NA)
        
        ylab    = vararg$ylab
        if( is.null(ylab) ) ylab=NA
        
        log     = vararg$log
        if( is.null(log) )  log=''
    
    
        #   init for just the subset
        if( ! initPlot.colorSpec( x.sub, xlim, ylim, ylab, log ) ) return(FALSE)
        
        if( is.logical(main)  &&  main )
            {
            #date    = metadata( x, "date" )    #; print( date )
            #if( is.null(date) ) date = file.info(path)$mtime
            #.title = sprintf( "%s\n[%s]", path, format(date,usetz=T) )
            main  = metadata( x, "path" )
            if( is.null(main) )
                main = theName      # last resort
            else
                main = basename(main)   # much shorter
            }
            
        if( is.character(main) )
            {
            title( main=main, cex.main=1 )
            }
        }
        
    
    color_vec   =   x.sub$color
    
    type    = vararg$type
    if( is.null(type) ) type = 'l'
    
    lty = vararg$lty
    lwd = vararg$lwd
    
    if( type == 'step' )
        {
        #   use segments()
        xvec    = breakandstep(wave)$breakvec
        n       = length(wave)
        
        if( is.null(lwd) )  lwd = 1
        
        if( is.null(lty) )  lty = 1
        
        if( length(lty) == 1 )  lty = rep(lty,2)
        
        for( j in 1:m )
            {
            y   = mat[ ,j]
            
            #   draw n horizontal segments, using lwd[1] and lty[1]
            segments( xvec[1:n], y, xvec[2:(n+1)], y, col=color_vec[j], lwd=lwd[1], lty=lty[1], lend='butt' )
            
            if( 2 <= length(lwd)  &&  ! is.na(lwd[2]) )
                #   draw n-1 vertical segments, using lwd[2] and lty[2]
                segments( xvec[2:(n+1)], y[1:(n-1)], xvec[2:(n+1)], y[2:n], col=color_vec[j], lwd=lwd[2], lty=lty[2], lend='round' )
            }
            
        lty.legend  = lty[1]    # for possible use in legend() below.  All spectra will get the same lty.
        }
    else
        {
        #   use lines()
        lty = suppressWarnings( rep( lty, len=m ) )
        lwd = suppressWarnings( rep( lwd, len=m ) )
        
        pch = suppressWarnings( rep( vararg$pch, len=m ) )
        
        for( j in 1:m )
            {
            lines.default( wave, mat[ ,j], type=type, col=color_vec[j], lty=lty[j], lwd=lwd[j], pch=pch[j] )
            }
            
        lty.legend  = lty       # for possible use in legend() below
        }
        
    #   the legend
    if( is.logical(legend) && legend )  
        legend = 'topright'
    
    if( is.character(legend) )
        {
        location    = legend
        
        legend  = specnames( x.sub )
        
        legend  = sub( "[.](Energy|Power)$", '', legend, ignore.case=TRUE )    # cleanup, maybe remove later
        
        if( CCT  &&  type(x) == "light" )
            {
            legend  = paste( legend, sprintf( "    [CCT=%g K]", round( computeCCT(x) ) ) )     # degree symbol °
            }
        
        if( 1 < length(color_vec)  &&  length( unique(color_vec) ) == 1 )   #  1 < length( unique(vararg$lty) ) 
            lwd = 2     # thinner so we can distinguish the lty
        else
            lwd = 11    # thicker so we can see colors instead
        
        legend( location, legend, col=color_vec, bty='n', lty=lty.legend, lwd=lwd, seg.len=4 )        
        }

        
    return( invisible(TRUE) )
    }
    
    



#   .x      a colorSpec object
#   .ymax   on the plot
#   .ymin   on the plot, but only when    
initPlot.colorSpec  <- function( .x, .xlim=c(NA,NA), .ylim=c(NA,NA), .ylab=NA, .log='' )
    {    
    #par( omi=rep(0,4) )
    #par( mai=c(1,3/4,1/2,1/4) )
    
    specnames   = specnames(.x)
    quantity    = quantity(.x)
    
    if( ! is.expression(.ylab)  &&  is.na(.ylab) )
        {
        #   assign a good default
        .ylab = "AU"
        
        if( quantity == "energy" || quantity == "power" )
            .ylab = "Radiant Energy (AU)"
        else if( quantity == "photons" || quantity == "photons/sec" )
            .ylab = "Photon Count  (AU)"
        else if( grepl(  "^(energy|power)->electrical$", quantity ) )
            .ylab = "Electrical Response / Radiant Energy"
        else if( grepl(  "^(energy|power)->neural$", quantity )  )
            .ylab = "Neural Response / Radiant Energy"
        else if( grepl(  "^(energy|power)->action$", quantity ) )
            .ylab = "Action Response / Radiant Energy"
        else if( quantity == "photons->electrical" )
            .ylab = "Electrical Response / Photon Count (QE)"
        else if( quantity == "photons->neural" )
            .ylab = "Neural Response / Photon Count"
        else if( quantity == "photons->action" )
            .ylab = "Action Response / Photon Count"
        else if( quantity == "reflectance" ||  quantity == "transmittance"  ||  quantity == "absorbance" )
            .ylab = paste( toupper( substr(quantity,1,1) ), substr( quantity, 2, nchar(quantity) ), sep='' )
        else if( quantity == "material->electrical" )
            .ylab = "Electrical Response / Material reflectance or transmittance"
        else if( quantity == "material->neural" )
            .ylab = "Neural Response / Material reflectance or transmittance"
        else if( quantity == "material->action" )
            .ylab = "Action Response / Material reflectance or transmittance"
        }
        
    data    = coredata( .x, forcemat=T )
    
    n   = ncol(data)
    
    wave    = wavelength(.x)
        
    xmin    = .xlim[1]
    xmax    = .xlim[2]
        
    if( is.na(xmin) )
        xmin = min( wave, na.rm=T )           
        
    if( is.na(xmax) )
        xmax = max( wave, na.rm=T )        

    ymin    = .ylim[1]
    ymax    = .ylim[2]
    
    if( is.na(ymax) )
        {
        ymax = max( data, na.rm=T )
        
        if( is.na(ymax) )
            {
            log.string( ERROR, "Cannot determine ymax limit for plot.  Is all data NA ?" )
            return( invisible(FALSE) )
            }
        }
    if( is.na(ymin) )
        {
        if( .log == 'y' )
            ymin = max( min(data), 1.e-4, na.rm=T )
        else
            ymin = min( data, 0, na.rm=T )
            
        if( is.na(ymin) )
            {
            log.string( ERROR, "Cannot determine ymin limit for plot.  Is all data NA ?" )
            return( invisible(FALSE) )
            }
        }
        
    marginx.axislabels   = 0.25
    
    plot.default( c(xmin,xmax), c(ymin,ymax), las=1, xlab='', ylab='', type='n', lab=c(10,8,7), log=.log,  tcl=0, mgp=c(3, marginx.axislabels, 0) )
    
    title( xlab="Wavelength (nm)", line=1.25 + marginx.axislabels )
        
    line_count = lines_for_ylab( as.character( axTicks(2) ) ) + 0.5   #; print( line_count )
    
    title( ylab=.ylab, line=line_count )  # or 2.5 or 2.0
    grid( lty=1, lwd=0.5, equilogs=F )
    abline( h=0 )    
    
    return( invisible(TRUE) )
    }
    
    
lines_for_ylab <- function( str )
    {
    #cat( "------------------\n")
    
    #   str = as.character( axTicks(2) )    #; print(str)
    width.str   = max(strwidth(str,units="inches"))
    height.str  = max(strheight(str,units="inches"))   #   should all have the same height !
    #cat( width.str, " ", height.str, '\n' )
    
    #usr = par('usr')    #; print(usr)
    #width.plot  = abs( usr[2] - usr[1] )
    #height.plot = abs( usr[4] - usr[3] )
    
    #width.plot  = par('fin')[1]
    #height.plot = par('fin')[2]
    
    # height.str is now in y-units, convert to x-units
    # height.str  = (width.plot/height.plot) * height.str ; print( height.str )
    
    fudge   = 0.70    
    line_count  = fudge * width.str / height.str   #; print(line_count)

    return( line_count )
    }
     

#   mat             an Nx3 numeric matrix, possibly with NAs
#   maxColorValue   same as rgb()
#   returns a character vector, with NAs where any input row is NA

myrgb   <- function( mat, maxColorValue=1 ) 
    {
    n   = nrow( mat )
    out = rep( NA_character_, n )
    
    #   this is the fastest way I found to test a whole row, even though we do not need the sums
    finite  = is.finite( .rowSums( mat, n, 3 ) )
    
    out[finite] = rgb( mat[finite, ,drop=F], maxColorValue=maxColorValue )
    
    return(out)
    }
    
    
    
    
    
#   obj         mx3 matrix of linear RGBs, with rownames and colnames assigned
#               or a data.frame with a matrix column RGB, optional columns LEFT,TOP,WIDTH,HEIGHT
#   normalize   scale the RGBs so the maximum is 1
#   gamma       of the target display, can be "sRGB" or a positive number (e.g. 2.2)
#   background  color of the background
#   shape       'full', 'half', etc.
plotPatchesRGB  <-  function( obj, normalize=FALSE, gamma='sRGB', background='gray50', labels=TRUE, shape='full', add=FALSE )
    {
    if( is.matrix(obj) )
        {
        if( ncol(obj) != 3 )
            {
            log.string( ERROR, "'%s' has %d columns, but it must be 3.", deparse(substitute(obj)), ncol(obj) )
            return(FALSE)
            }    
            
        #   test the column names
        diag    = ( diag( rep(T,3) ) == 1 )
        
        mat1    =   multiPatternMatch( c("^R","^G","^B"), colnames(obj), .ignore=T )  #; print(mat1)
        mat2    =   multiPatternMatch( c("R$","G$","B$"), colnames(obj), .ignore=T )  #; print(mat2)
        
        ok  = all( mat1==diag ) || all( mat2==diag )
        #   ok  = all( toupper(colnames(obj))  ==  c('R','G','B') )
        
        if( ! ok )
            {
            log.string( ERROR, "'%s' is not RGB.", deparse(substitute(.spec)) )
            return(FALSE)
            }    
                    
                    
        #   convert the matrix to a data.frame with that matrix as the only column
        
        #obj   = as.data.frame.model.matrix( obj )
        #colnames(obj) = 'RGB'
        
        df  = data.frame( row.names=1:nrow(obj) )
        df$RGB  = obj
        obj     = df
        }
        
    if( ! is.data.frame(obj) )
        {
        log.string( ERROR, "data is invalid; neither matrix nor data.frame" )
        return(FALSE)
        }

    #   look for column RGB
    if( is.null(obj$RGB)  ||  is.null(dim(obj$RGB))  ||  ncol(obj$RGB) != 3 )
        {
        log.string( ERROR, "data is invalid; there is no column RGB with 3 columns" )
        return(FALSE)
        }    

    
    n = nrow(obj)
        
    
    # look for columns LEFT,TOP,WIDTH,HEIGHT
    ok  = all( c("LEFT","TOP","WIDTH","HEIGHT") %in% colnames(obj) )
    
    if( ok )
        {
        xlim = range( obj$LEFT, obj$LEFT + obj$WIDTH )
        ylim = range( obj$TOP, obj$TOP + obj$HEIGHT )
        ylim = ylim[2:1]
        
        aspect  = 1        
        }
    else
        {
        #   make vertically stacked patches on the left, with lots of room for labels on the right
        #  add extra columns LEFT,TOP,WIDTH,HEIGHT        
        obj   = cbind( obj, LEFT=0, TOP=0:(n-1), WIDTH=1, HEIGHT=1 )
        
        xlim = c( 0, 3 )
        ylim = c( n, 0 )
        aspect  = NA
        
        #print( obj )
        #return(FALSE)
        }
        
    if( normalize )
        {
        obj$RGB = obj$RGB / max(obj$RGB)
        #   names(obj)    = theNames
        }

    if( gamma == 'sRGB' )
        obj$RGB = DisplayRGBfromLinearRGB( obj$RGB, gamma )
    else if( is.numeric(gamma)  &&  0 < gamma )
        obj$RGB = DisplayRGBfromLinearRGB( obj$RGB, gamma )
    else
        {
        log.string( ERROR, "gamma is invalid" )
        return(FALSE)
        }    
        
    colvec  = rgb( obj$RGB )      #;    print( colvec )
    

    
    #par( omi=rep(0,4) )
    #par( mai=c( 0.25, 0.25, 0.25, 0.25) )   

    if( ! add )
        {
        #   check background
        if( is.numeric(background) )
            {
            if( length(background) == 1 )   background = rep( background, 3 )
            
            if( length(background) != 3 )
                {
                log.string( ERROR, "background is invalid, because length(background)==%d", length(background) )
                return(FALSE)
                }
                
            dim(background) = c(1,3)            
            background  = rgb( DisplayRGBfromLinearRGB( background, gamma ) )
            }
            
        if( ! is.character( background ) )
            {
            log.string( ERROR, "background is invalid" )
            return(FALSE)
            }        

        bg.prev = par( bg=background )
        
        par( mgp=c(0, 0.5, 0) )
        
        plot.new()    

        plot.window( xlim, ylim, asp=aspect )    
        }
        
        
    #   print( dimnames(obj) )
    
    rectangle   = shape %in% c('full','left','right','bottom','top','half')
    
    #   triangle    = shape %in% c("bottomright", "bottomleft", "topleft", "topright")
    #   triangle    = match( shape, c("topleft", "topright",  "bottomleft","bottomright") , nomatch=0 )
 
    triangle    = match( shape, c("bottomright", "bottomleft", "topright", "topleft") , nomatch=0 )
 
 
    #   set the full size
    left    = obj$LEFT
    right   = left + obj$WIDTH
    top     = obj$TOP
    bottom  = top + obj$HEIGHT
    
    if( rectangle )
        {
        #   optionally move 1 or all sides
        if( shape == 'left' )
            right   = obj$LEFT + 0.5 * obj$WIDTH
        else if( shape == 'right' )
            left    = obj$LEFT + 0.5 * obj$WIDTH
        else if( shape == 'top' )
            bottom  = obj$TOP  + 0.5 * obj$HEIGHT
        else if( shape == 'bottom' )
            top     = obj$TOP  + 0.5 * obj$HEIGHT
        else if( shape == 'half' )
            {
            x   = 0.5*(left + right)
            y   = 0.5*(top + bottom)
            
            left    = 0.5*(left + x)
            right   = 0.5*(right + x)
            top     = 0.5*(top + y)
            bottom  = 0.5*(bottom + y)
            }
            
        rect( left, bottom, right, top, col=colvec, border=NA )
        }
    else if( 0 < triangle )
        {
        for( i in 1:n )
            {
            #   assign all 4 vertices            
            xy  = c( left[i],top[i],  right[i],top[i],  left[i],bottom[i],  right[i],bottom[i] )
            dim(xy) = c(2,4)
            
            #   now drop one column
            xy  = xy[ ,-triangle]

            polygon( xy[1, ], xy[2, ], col=colvec[i], border=NA )
            }
        }
    else
        {
        log.string( ERROR, "shape='%s' unknown.", shape )
        }
        
    if( labels )
        {
        text( right + 0.1*obj$WIDTH, top + 0.5*obj$HEIGHT, rownames(obj), adj=c(0,0.5) )
        }
            
    if( ! add )
        {            
        par( bg=bg.prev )   # restore previous background       
        }
        
        
    return( invisible(T) )
    }

        
#--------       UseMethod() calls           --------------#                    
        
if( 0 )
{        
plot <- function( .x,  .xlim=c(NA,NA), .ylim=c(NA,NA), .add=FALSE, .legend=TRUE, .title=T, 
                        .color='natural', .subset=NA, .ylab=NA, .lty=1, .log='', .CCT=FALSE, .points=FALSE  )
    {
    UseMethod("plot")
    }
}    
