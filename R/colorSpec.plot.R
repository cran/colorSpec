
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
    
    xlim    = vararg$xlim
    if( is.null(xlim) ) xlim=c(NA,NA)
    
    ylim    = vararg$ylim
    if( is.null(ylim) ) ylim=c(NA,NA)
    
    ylab    = vararg$ylab
    if( is.null(ylab) ) ylab=NA
    
    log     = vararg$log
    if( is.null(log) )  log=''
    
    pch     = vararg$pch
    
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
            #   fake it by changing quantity to power
            x.power  = x    #   resample( x, wavelength(BT.709.RGB) )    #   we know BT.709.RGB has regular wavelength steps
            quantity(x.power) = 'power'
            mat.rgb = product( x.power, colorSpec::BT.709.RGB, wave='auto' )         #; print( mat.rgb )
            mat.rgb = mat.rgb / max(mat.rgb)            #   normalize so max RGB is 1; print( mat.rgb )
            }
            
        mat.rgb = DisplayRGBfromLinearRGB( mat.rgb )  #; print( mat.rgb )
         
        color_vec   = rgb( mat.rgb )                #;   print( color_vec )
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
        
    if( ! add )
        {
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
        
    #   plot just the subset
    mat = as.matrix( x.sub )
    
    color_vec   =   x.sub$color
        
    lty = suppressWarnings( rep( vararg$lty, len=ncol(mat) ) )
    lwd = suppressWarnings( rep( vararg$lwd, len=ncol(mat) ) )

    for( j in 1:ncol(mat) )
        {
        lines.default( wave, mat[ ,j], col=color_vec[j], lty=lty[j], lwd=lwd[j] )      # lwd=0.5
        
        if( is.numeric(pch) )
            points.default( wave, mat[ ,j], col=color_vec[j], pch=pch )
        }
        
        
    #   the legend
    if( is.logical(legend) && legend )  
        legend = 'topright'
    
    if( is.character(legend) )
        {
        location    = legend
        
        legend  = specnames( x.sub )
        
        legend  = sub( "[.]Power$", '', legend, ignore.case=TRUE )    # cleanup, maybe remove later
        
        if( CCT  &&  type(x) == "light" )
            {
            legend  = paste( legend, sprintf( "    [CCT=%d K\u00B0]", round( computeCCT(x) ) ) )     # degree symbol °
            }
        
        if( 1 < length(color_vec)  &&  length( unique(color_vec) ) == 1 )   #  1 < length( unique(vararg$lty) ) 
            lwd = 2     # thinner so we can distinguish the lty
        else
            lwd = 11    # thicker so we can see colors instead
        
        legend( location, legend, col=color_vec, bty='n', lty=vararg$lty, lwd=lwd, seg.len=4 )        
        }

        
    return( invisible(T) )
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
    
    if( is.na(.ylab) )
        {
        .ylab = "AU"
        
        if( quantity == "power" )
            .ylab = "Radiant Power (AU)"
        else if( quantity == "photons" || quantity == "photons/sec" )
            .ylab = "Photon Flux  (photons per second) (AU)"
        else if( quantity == "power->electrical" )
            .ylab = "Electrical Response / Radiant Power"
        else if( quantity == "power->neural" )
            .ylab = "Neural Response / Radiant Power"
        else if( quantity == "power->action" )
            .ylab = "Action Response / Radiant Power"
        else if( quantity == "photons->electrical" )
            .ylab = "Electrical Response / Photon Flux  (QE)"
        else if( quantity == "photons->neural" )
            .ylab = "Neural Response / Photon Flux"
        else if( quantity == "photons->action" )
            .ylab = "Action Response / Photon Flux"
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
        ymax = max( data, na.rm=T )
        
    if( is.na(ymin) )
        {
        if( .log == 'y' )
            ymin = max( min(data), 1.e-4, na.rm=T )
        else
            ymin = min( data, 0, na.rm=T )
        }
        
    margin.axislabels   = 0.25
    
    plot.default( c(xmin,xmax), c(ymin,ymax), las=1, xlab='', ylab='', type='n', lab=c(10,8,7), log=.log,  tcl=0, mgp=c(3, margin.axislabels, 0) )
    

    title( xlab="Wavelength (nm)", line=1.25 + margin.axislabels )
    
    line_count = lines_for_ylab() + margin.axislabels   #; print( line_count )
    
    title( ylab=.ylab, line=line_count )  # or 2.5 or 2.0
    grid( lty=1, equilogs=F )
    abline( h=0 )    
    
    return( invisible(TRUE) )
    }
    
    
#   no argument - use the current plot context    
lines_for_ylab <- function()
    {
    #cat( "------------------\n")
    
    str = as.character( axTicks(2) )    #; print(str)
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
    line_count  = fudge * width.str / height.str #; print(line_count)

    return( line_count )
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
                    
        obj   = as.data.frame.model.matrix( obj )
        colnames(obj) = 'RGB'
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
            xy  = c( left[i],top[i],  right[i],top[i],  right[i],bottom[i],  left[i],bottom[i] )
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
