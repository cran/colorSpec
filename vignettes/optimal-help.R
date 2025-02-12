

#   reproduce Fig. 3(3.7) page 182, and similar ones, in Wyszecki & Stiles


#   seclist     list of sections
#   Ylevel      beta, the plane constant
#   obj         the light responder
#   white       white XYZ, of the illuminant

plotSections <- function( seclist, Ylevel, obj, white, col='black', add=FALSE )
    {
    if( ! add ) plotChromaticityDiagram( obj, white )
    
    lwd = 1   # ifelse( add, 0.5, 1 )
    
    for( k in 1:length(seclist) )
        {
        section = seclist[[k]]$section
        
        denom   = rowSums( section )
        x   = section[ ,1] / denom
        y   = section[ ,2] / denom
        
        polygon( x, y, border=col, lwd=lwd )
        
        if( ! add )
            {
            idx     = which.min( x + y )
            gray    = Ylevel[k]
            if( gray < 1 )
                # display as a percentage
                gray = 100 * gray

            text( x[idx], y[idx], sprintf( "%g", gray ), adj=c(-0.25,0), cex=0.6 )
            }
        }

    if( ! add ) plotWavelengthPoints( obj )    

    return( invisible(TRUE) )
    }

    
plotChromaticityDiagram  <-  function( .xyz=xyz1931.1nm, white )
    {
    coredata    = coredata( .xyz )
    denom       = rowSums( coredata )
    x   = coredata[ ,1] / denom
    y   = coredata[ ,2] / denom    
    
    xylab   =  tolower( specnames(.xyz) )   #tolower( substr(specnames(.xyz),1,1) )
    
    plot.default( range(x), range(y), type='n', las=1, xlab='', ylab='', asp=1, lab=c(10,8,7), tcl=0, mgp=c(3, 0.25, 0)  )
    title( xlab=xylab[1], line=1.5 )
    title( ylab=xylab[2], line=2 )    
    grid( lty=1 )
    abline( h=0, v=0 )    
    
    polygon( x, y, col='white' )
    
    #   put a black dot at white point
    # denom   = sum( white )
    xy      = white[1:2] / sum(white)
    points( xy[1], xy[2], pch=20 )
    
    legend( 'topright', c("Schrödinger","optimal"), col=c("black","red"), bty='n', lty=1, lwd=5 )  #, seg.len=4 )
        
    
    return( TRUE )
    }

    
#   plot nice points, pch=21, at multiples of 20 nm
    
plotWavelengthPoints <- function( .obj )
    {
    wave    = wavelength(.obj)
    idx = which( wave %% 20 == 0 )
    
    if( length(idx) <= 1 )  return(FALSE)

    coredata    = coredata( .obj )
    denom       = rowSums( coredata )
    x   = coredata[ ,1] / denom
    y   = coredata[ ,2] / denom    
           
    xy      = cbind( x[idx], y[idx] )   #; print( xy )
    n       = nrow(xy)
    dist    = xy[ 2:n, ] - xy[1:(n-1), ]    #; print( dist )
    dist    = sqrt( rowSums(dist*dist) )    #; print( dist )
    idx     = idx[ 0.75*strheight("560") < dist ]
    
    points( x[idx], y[idx], pch=21, bg='white', cex=0.8 )

    for( j in idx )
        {    
        tangent = c( x[j+1], y[j+1] ) - c( x[j], y[j] )
        normal  = c( -tangent[2], tangent[1] )
        normal  = normal / sqrt( sum(normal*normal) )
        adj     = -normal + c(0.4,0.5)
        text( x[j], y[j], as.character(wave[j]), adj=adj, cex=0.75 )
        }
        
    return(T)
    }
