
#   x           a colorSpec object, with type "responsivity.material".  The number of spectra must be 3.
#   gray        vector of numbers in (0,1) that define neutral gray on segment from black to white
#               this neutral gray point is the base of the probe ray
#   direction   a non-zero 3-vector defining the ray, or a matrix with 3 columns defining multiple rays
#   tol         convergence tolerance. Iterations continue until root-finding error < .tol.  
#   aux         return auxiliary performance data
#
#   return value
#   data.frame with a row for each direction and these columns:
#       gray        input
#       direction   input
#       s           position along the ray that intersects boundary
#       optimal     the optimal color on the boundary, where the ray intersects boundary
#       lambda      lambda.1 and lambda.2 of the ideal material producing the optimal color
#                       lambda.1 < lambda.2 => bandpass
#                       lambda.1 > lambda.2 => bandstop
#       do          delta and omega, the Logvinenko parameters - analogous to latitude and longitude
#
#   and if aux is TRUE, these auxiliary columns:
#       time_grid   time to find initial estimate point on boundary, in seconds
#       iters       # of interations taken to find XYZ
#       btracks     total # of backtracks in "damped" Newton's method
#       time_newt   time spent in Newton iterations, in seconds
#       error       root-finding error, in coordinates of the optimal color
#
#   in case of ERROR returns NULL

probeOptimalColors.colorSpec <- function( x, gray, direction, tol=1.e-6, aux=FALSE )
    {
    time_start  = as.double( Sys.time() )
    
    theName = deparse(substitute(x))
    
    spectra = numSpectra(x)
    
    if( spectra != 3 )
        {
        log.string( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )        
        return(NULL)
        }        

    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )        
        return(NULL)
        }    
                
        
    funlist = makeFunctionList( x, theName, 'equalize' )
    
    if( is.null(funlist) ) return(NULL)


    step.wl = step.wl(x)

    range.wl    =  range(wavelength(x))

    mask    =  0<gray & gray<1 
    if( ! all(mask) )
        {
        log.string( ERROR, "gray = %g is invalid.", gray[ which(!mask)[1] ] )
        return(NULL)
        }    
    
    if( length(direction) == 3 )
        direction = matrix( direction, 1, 3 )
        
    ok  = is.matrix(direction) &&  ncol(direction)==3
    if( ! ok )
        {
        log.string( ERROR, "direction is invalid." )             
        log.object( ERROR, direction, addname=F )
        return(NULL)
        }       
        

    gray    = rep( gray, each=nrow(direction) )
        
    rays    = length(gray)
    
    direction  = matrix( t(direction), rays, 3, byrow=T )
        
    colnames(direction) = toupper( specnames(x) )
    class(direction) = "model.matrix"
    

    wavelengths = numWavelengths(x)
          
        
    #   log.object( TRACE, sfun, type='str' )
    
    grid.boundary   = computeGridOfOptimals( funlist$integral.from.omega )
    
    log.object( TRACE, grid.boundary, type='str' )
        
    XYZ.white   = numeric( 3 )
    for( k in 1:3 )
        XYZ.white[k]    = funlist$integral.from.omega[[k]]( 1 )
    
    #   XYZ.white           = mat.cum[ wavelengths, ]     #; log.object( TRACE, XYZ.white ) ;  log.object( TRACE, XYZ.white / XYZ.white[2] )
    
    #base                = matrix( gray * XYZ.white, 1, 3 )
    #colnames(base)      = colnames(direction)
    
    #basemat             = matrix( base, directions, 3, byrow=T )
    #colnames(basemat)   = colnames(direction)      
    #class(basemat)      = "model.matrix"

    XYZmat              = matrix( as.numeric(NA), rays, 3 )
    colnames(XYZmat)    = colnames(direction)      #; print( colnames(XYZmat) )
    class(XYZmat)       = "model.matrix"
    
    lambda  = matrix( as.numeric(NA), rays, 2 )   # in nanometers
    class(lambda)       = "model.matrix"
    
    dol = matrix( as.numeric(NA), rays, 3 )
    colnames(dol)       = c( "delta", "omega", "lambda" )     # reparameterized width and center. Both in [0,1].  plus lambda equivalent of omega
    class(dol)          = "model.matrix"
    
    
    
    out = data.frame( gray=gray, direction=direction, s=as.numeric(NA), optimal=XYZmat, lambda=lambda, dol=dol )
    
    if( aux )
        df.aux  = data.frame( time_grid=rep(as.numeric(NA),nrow(out)), iters=as.integer(NA),  btracks=as.integer(NA), time_newt=as.numeric(NA), error=as.numeric(NA) )   #;   print( str(out) )
    

    #XYZ = rep( as.numeric(NA), 3 )  # matrix( as.numeric(NA), 1, 3 )    
    #dim(XYZ)    = NULL
    #dim(base)   = NULL
    
    #   jac2x2  = matrix( c(1,-1/2,1,1/2), 2, 2, byrow=T )
    
    time_elapsed    = as.double( Sys.time() ) - time_start
    
    log.string( INFO, "Preprocessing time: %g sec.", time_elapsed )
    
    time_start  = as.double( Sys.time() )
    
    quad.start  = c(1,1)
    
    for( i in 1:rays )
        {
        base    = out$gray[i] * XYZ.white
    
        vec     = out$direction[i, ]
        
        log.string( INFO, "##########  gray=%g           base=%g,%g,%g      direction=%g,%g,%g   ##################", 
                                out$gray[i],    base[1], base[2], base[3],     vec[1], vec[2], vec[3] )

        if( all( vec == 0 ) )
            {
            log.string( WARN, "Direction %d is 0.  No ray tracing is possible.", i )
            next
            }
        
        rot = t( svd(vec,3,3)$u )   ; log.object( TRACE, rot )
        
        iters_max   = 10

        ok          = FALSE        

        time0   = as.double( Sys.time() )
        
        res = gridSearchBoundary( grid.boundary, base, rot, quad.start )    #; log.object( TRACE, cw_min )

        #cw0     = c(1/2,1/2)
        #radius  = c(16,16)
        #step    = c(1/32,1/32)
                
        #   cw_min     = gridSearchCW( cw0, radius, step, sfun, base, rot )   ; print( cw_min )  obsolete
            
        #   cw_min  = c( 0.330883, 0.545647 )   # very accurate override    
        #   cw_min  = c( 0.3238551, 0.5558677 ) #  test override

        time1               = as.double( Sys.time() )
        
        if( aux )
            df.aux$time_grid[i]    =  time1 -  time0        
        
        if( is.null(res)  )   next
        
        
        cw_min      = res$cw.root
        quad.start  = res$quad          # to be used for next direction
        
        time0   = time1
        
        #   compute data at cw_min
        data  = dataFromCW( cw_min, funlist$integral.from.omega, base, rot ) #; print( data$XYZ )

        ok          = FALSE
        iters_max   = 20
        btracks     = 0
        
        for( iter in 1:iters_max )
            {
            log.string( DEBUG, "----------------  Newton iteration %d  ----------", iter )
            log.object( DEBUG, data )

            if( data$uv.norm2 < tol )   
                {
                #   all done !
                ok  = TRUE
                break
                }
                

            #   solve for delta
            uv  = data$test[2:3] 

            delta     = try( solve( data$jac2x2, -uv ) )  #; log.object( WARN, delta )
            if( class(delta) == "try-error" )
                {
                log.string( WARN, "Resorting to pseudo-inverse..." )

                log.object( TRACE, data$jac2x2 )
                
                requireNamespace( 'MASS', quietly=TRUE )
                
                delta = MASS::ginv(data$jac2x2) %*% (-uv)
                
                if( all(delta==0) )
                    {
                    log.string( WARN, "Pseudo-inverse failed." )
                    log.string( ERROR, "delta=0.  uv=%g %g", uv[1], uv[2] )
                    break
                    }
                }
                
            log.string( TRACE, "delta = %g %g", delta[1], delta[2] )
            
            #   clamp delta so 0.001 < w < 0.999
            if( delta[2] != 0 )
                {
                for( lim in c(0.001,0.999) )
                    {
                    s   = (lim - cw_min[2]) / delta[2]
                    if( 0<s  &&  s<1 ) 
                        {
                        log.string( TRACE, "clamping delta so 0.001<w && w<0.999" )        
                        delta = s * delta
                        log.string( TRACE, "delta = %g %g", delta[1], delta[2] )
                        }
                    }
                }
            
            delta.norm2 = sqrt( sum(delta*delta) )
            if( 0.1 < delta.norm2 )
                {
                log.string( TRACE, "clamping delta to length 0.1" )
                delta = 0.1 * delta / delta.norm2
                log.string( TRACE, "delta = %g %g", delta[1], delta[2] )         
                }
                
            #   using damping technique to find new point that reduces objective
            reduction   = FALSE
            
            for( j in 0:9 )
                {
                cw      = data$cw + (2^-j)*delta
                cw[1]   = cw[1] %% 1   
                
                data.next   = dataFromCW( cw, funlist$integral.from.omega, base, rot ) #; print( data$XYZ )
                
                if( data.next$test[1] <= 0 )
                    {
                    log.string( INFO, "For direction=%d,  wandered into wrong halfspace. test[1]=%g.  j=%d. Try cutting step in half.", 
                                        i, data.next$test[1], j )
                    btracks = btracks + 1
                    next
                    }
                    
                if( data.next$uv.norm2 < data$uv.norm2  )
                    {
                    #   jump to new point
                    if( 0 < j )                    
                        {
                        log.string( TRACE, "Backtrack Successful at lambda=%g.  %g < %g", 
                                            2^-j, data.next$uv.norm2, data$uv.norm2  )                  
                        }
                    data        = data.next                    
                    reduction   = TRUE
                    break
                    }
                    
                #   backtrack
                if( j == 0 )
                    {
                    log.string( TRACE, "Backtrack necessary.    %g < %g", data$uv.norm2, data.next$uv.norm2 )
                    }

                btracks = btracks + 1                
                }
                
            if( ! reduction )
                {
                log.string( ERROR, "Backtrack failed. Cannot reduce objective using damping technique." )
                break
                }
            }   
            
        time1   = as.double( Sys.time() )
        
        if( ok )
            {
            out$optimal[i, ]    = data$XYZ
            out$lambda[i, ]     = funlist$lambda.from.omega( data$omega )
            out$dol[i, ]        = c( deltaFromWidth( data$cw[2] ), data$cw[1], funlist$lambda.from.omega( data$cw[1] ) )
            difference          = data$XYZ - base
            out$s[i]            = sum(difference*vec) / sum(vec*vec)
            }
            
        if( aux )
            {
            df.aux$time_newt[i] = time1  -  time0    
            df.aux$iters[i]     = iter      
            df.aux$btracks[i]   = btracks
            df.aux$error[i]     = data$uv.norm2               
            }
        }
        
    time_elapsed    = as.double( Sys.time() ) - time_start
    log.string( INFO, "Processed %d rays in %g sec  (%g sec per rays)",
                        rays, time_elapsed, time_elapsed/rays ) 

    failures    = sum( is.na( out$s ) )
    if( 0 < failures )
        log.string( WARN, "There were %d failures out of %d rays.\n", failures, nrow(out) )
    
    if( aux )
        out = cbind( out, df.aux )
    
    return(out)
    }
    
twoTransitionGrid.colorSpec <- function( .obj, .size=c(33,33) )
    {
    theName = deparse(substitute(.obj))

    if( type(.obj) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(.obj) )        
        return(NULL)
        }        
    
    theList = makeFunctionList( .obj, theName )
    
    if( is.null(theList) ) return(NULL)
        
    out = computeGridOfOptimals( theList$integral.from.omega, .size )
    
    return( out )
    }
        
    
makeSplineFunctions <- function( .obj, .name )
    {
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
        
    step.wl = step.wl(.obj)
        
    coredata    = coredata( .obj )
    if( min(coredata) < 0 )
        {
        log.string( ERROR, "responsivity of %s is invalid; some values are negative.", .name )        
        return(NULL)
        }    
    
    spectra = ncol( coredata )
       
    #   time_svd    = as.double( Sys.time() )
    
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
        
    
    wavelengths = nrow(coredata)
    
    #   make mat.cum and 3 spline interpolators
    mat.cum = matrix( 0, wavelengths, spectra )
    sfun    = list()
    for( j in 1:spectra )
        {
        mat.cum[ ,j]    = step.wl * cumsum( coredata[ , j ] )
        
        sfun[[j]]       = splinefun( (0:wavelengths)/wavelengths, c(0,mat.cum[ ,j]), method='monoH.FC' )      # hyman  monoH.FC  natural
        }
        
    return( sfun )
    }
        
            
    
    
#   .grid.boundary  list returned from computeGridOfOptimals()
#   .radius     half-size of grid, both directions
#   .step       step size of grid, both directions    
#   .base       of the search ray
#   .rotation   orthogonal frame at the .base, with direction of ray in 1st row
#   .sfun       spline function integral of XYZ responsivity    

#   returns    cw.root
#       

gridSearchBoundary.old <- function( .grid.boundary, .base, .rotation )
    {
    c.seq   = .grid.boundary$c.seq
    w.seq   = .grid.boundary$w.seq
    
    uv   = array( as.numeric(NA), dim=c(length(c.seq), length(w.seq), 2) )
    
    for( x in 1:length(c.seq) )
        {
        for( y in 1:length(w.seq) )
            {
            test    = .rotation %*% ( .grid.boundary$XYZ[x,y, ] - .base )

            if( test[1] < 0 )  next    #in the wrong half-plane

            uv[x,y, ]  = test[2:3]
            }
        }
    
    cw.mat  = matrix( 0, 4, 2 )
    uv.mat  = matrix( 0, 4, 2 )
    for( x in 1:(nrow(uv)-1) )
        {
        for( y in 1:(ncol(uv)-1) )
            {
            uv.mat[1, ] = uv[x,y, ]
            uv.mat[2, ] = uv[x+1,y, ]
            uv.mat[3, ] = uv[x+1,y+1, ]
            uv.mat[4, ] = uv[x,y+1, ]
            
            if( any(is.na(uv.mat)) )  next   

            #   polygon( uv.mat[1, ], uv.mat[2, ] ) ; next            

            w   = windingNumber( c(0,0), uv.mat[ ,1], uv.mat[ ,2] )
            
            if( w == 0 )    next
            
            log.string( TRACE, "Found inside:  x=%d  y=%d.  Winding Number = %d.", x, y, w )  
            log.object( TRACE, uv.mat )                            
            
            cw.mat[1, ] = c( c.seq[x], w.seq[y] )
            cw.mat[2, ] = c( c.seq[x+1], w.seq[y] )    
            cw.mat[3, ] = c( c.seq[x+1], w.seq[y+1] ) 
            cw.mat[4, ] = c( c.seq[x], w.seq[y+1] ) 

            cw.root = interpQuad( uv.mat, c(0,0), cw.mat ) # ; print( cw.root )
            
            return( cw.root )
            }
        }
        
    log.string( ERROR, "Found no quadrilateral that covers 0." )
    
    return( NULL )
    }
        

#   returns a list with
#       cw.root     (center,width) root
#       quad        (i,j) of quadrilateral where root was found        
        
gridSearchBoundary <- function( .grid.boundary, .base, .rotation, .start )
    {
    c.seq   = .grid.boundary$c.seq
    w.seq   = .grid.boundary$w.seq
    
    #   in this uv array
    #       NA      means that the value has not been computed yet
    #       Inf     means it has been computed, but boundary point is in the wrong half-space
    uv   = array( as.numeric(NA), dim=c(length(c.seq), length(w.seq), 2) )
    
    if( FALSE )
    {
    for( i in 1:length(c.seq) )
        {
        for( j in 1:length(w.seq) )
            {
            test    = .rotation %*% ( .grid.boundary$XYZ[i,j, ] - .base )

            if( test[1] < 0 )  next    #in the wrong half-space

            uv[i,j, ]  = test[2:3]
            }
        }
    }   
    
    cw.mat  = matrix( 0, 4, 2 )
    uv.mat  = matrix( 0, 4, 2 )
    
    #time0   = as.double( Sys.time() )
    
    mat.search  = ringOrder( c(length(c.seq),length(w.seq)) - 1L, .start )   # about 1 msec  
    
    #cat( sprintf( "ringOrder() elapsed = %g sec\n",  as.double( Sys.time() ) - time0 ) )

    i4  = c(0,1,1,0)
    j4  = c(0,0,1,1)
 
    for( k in 1:nrow(mat.search) )
        {
        i   = mat.search[k,1]
        j   = mat.search[k,2]
        
        #uv.mat[1, ] = uv[i,j, ]
        #uv.mat[2, ] = uv[i+1,j, ]
        #uv.mat[3, ] = uv[i+1,j+1, ]
        #uv.mat[4, ] = uv[i,j+1, ]
        
        for( m in 1:4 )
            {
            ip  = i + i4[m]
            jp  = j + j4[m]
            
            if( is.na( uv[ip,jp,1] ) )
                {
                # ip,jp has not been visited yet
                test    = .rotation %*% ( .grid.boundary$XYZ[ip,jp, ] - .base )

                if( 0 <= test[1] )
                    #   in the right half-space
                    uv[ip,jp, ] = test[2:3]
                else
                    #   in the wrong half-space
                    uv[ip,jp, ] = Inf
                }
                
            uv.mat[m, ] = uv[ip,jp, ]               
            }

        if( ! all( is.finite(uv.mat) ) )    next    # quadrilateral is invalid
        
        
        #   polygon( uv.mat[1, ], uv.mat[2, ] ) ; next            

        w   = windingNumber( c(0,0), uv.mat[ ,1], uv.mat[ ,2] )
        
        if( w == 0 )    next
        
        log.string( TRACE, "Found inside:  i=%d  j=%d.  Winding Number = %d.", i, j, w )  
        log.object( TRACE, uv.mat )                            
        
        cw.mat[1, ] = c( c.seq[i], w.seq[j] )
        cw.mat[2, ] = c( c.seq[i+1], w.seq[j] )    
        cw.mat[3, ] = c( c.seq[i+1], w.seq[j+1] ) 
        cw.mat[4, ] = c( c.seq[i], w.seq[j+1] ) 
        
        cw.root = interpQuad( uv.mat, c(0,0), cw.mat ) # ; print( cw.root )
        
        if( is.null(cw.root) )  return(NULL)    # degenerate quad
        
        out = list()        
        out$cw.root = cw.root
        out$quad    = c(i,j)
        
        return( out )
        }
        
        
    log.string( ERROR, "Found no quadrilateral that covers 0." )
    
    return( NULL )
    }

    
    
#   .n              number in center and width directions
#   .funintegral    integral.from.omega list
#
#   returns list with
#       c.seq       from 0 to 1, with length .n[1]
#       w.seq       from 0 to 1, with length .n[2]
#       XYZ[ , ,3]  computed XYZ values    
computeGridOfOptimals <- function( .funintegral, .n=c(33,33) )
    {
    if( length(.n) != 2 )   .n = rep( .n[1], 2 )
    
    c.seq   = seq( 0, 1, len=.n[1] )
    w.seq   = seq( 0, 1, len=.n[2] )
    
    #w.seq[ .n[2] ]      = 0.98
    #w.seq[ .n[2]+1 ]    = 1
    #print( w.seq )
    
    channels    = length( .funintegral )

    XYZ.white = numeric(channels)
            
    for( k in 1:channels )
        XYZ.white[k]   = .funintegral[[k]]( 1 )     #  .mat.cum[ nrow(.mat.cum), ]
    
    XYZ = array( 0, dim=c(length(c.seq), length(w.seq), channels) )
    
    
    for( y in 1:length(w.seq) )
        {
        width   = w.seq[y]
        
        if( width == 0 )
            #   special case XYZ=0, bandwidth=0  =>  all black
            next
        
        for( x in 1:length(c.seq) )
            {
            center  = c.seq[x]

            if( width < 1 )
                {
                #   usual case
                omega   = omegaFromCW( center, width )

                for( k in 1:channels )
                    {
                    r1  =   evaluateD0D1( omega[1], .funintegral[[k]] )    # .mat.cum[ ,k] )
                    r2  =   evaluateD0D1( omega[2], .funintegral[[k]] )    # .mat.cum[ ,k] )      
                    if( omega[1] < omega[2] )
                        #   bandpass                    
                        XYZ[x,y,k]  = r2$v - r1$v
                    else
                        #   bandstop
                        XYZ[x,y,k]  = r2$v - r1$v + XYZ.white[k]
                    }
                }
            else
                {
                #   special case, bandwidth=1  => all white
                XYZ[x,y, ]  =  XYZ.white
                }
            }
        }
        
    out = list( c.seq=c.seq, w.seq=w.seq, XYZ=XYZ )
        
    return( out )
    }
    
    
    
    
#   .cw         center in [0,1] but taken as in the circle
#               width in [0,1]
#   .sfun       spline function integral of xyz responsivity    
#   .base       base of ray
#   .rotation   3x3 rotation matrix
dataFromCW <- function( .cw, .sfun, .base, .rotation  )
    {
    omega   = omegaFromCW( .cw[1], .cw[2] )

    XYZ = numeric(3)
    XYZ.white = XYZ
            
    for( k in 1:3 )
        XYZ.white[k]   = .sfun[[k]]( 1 )     #  .mat.cum[ nrow(.mat.cum), ]

    jac3x2  = matrix( 0, 3, 2 )
        
    if( omega[1] < omega[2] )
        {
        #   bandpass
        for( k in 1:3 )
            {
            r1  =   evaluateD0D1( omega[1], .sfun[[k]] )    # .mat.cum[ ,k] )
            r2  =   evaluateD0D1( omega[2], .sfun[[k]] )    # .mat.cum[ ,k] )                        
            XYZ[k]  = r2$v - r1$v
            jac3x2[k,1]  = -r1$dv
            jac3x2[k,2]  =  r2$dv                    
            }
        }
    else
        {
        #   bandstop
        for( k in 1:3 )   
            {
            r1  =   evaluateD0D1( omega[1], .sfun[[k]] )    # .mat.cum[ ,k] )
            r2  =   evaluateD0D1( omega[2], .sfun[[k]] )    # .mat.cum[ ,k] )   
            XYZ[k]  = r2$v +  XYZ.white[k] - r1$v
            jac3x2[k,1]  = -r1$dv
            jac3x2[k,2]  =  r2$dv                    
            }
        }
        
    #   print( jac3x2 )        
        
    a   = 0.25 * pi * sin( pi * .cw[2])
    
    jac2x2  = matrix( c(1,-a, 1,a), 2, 2, byrow=T )
    
    difference  = XYZ - .base
    test        = .rotation %*% difference
     
    uv  =   test[2:3]
    uv.norm2    = sqrt( sum(uv*uv) )        
        

    out = list()
    out$cw      = as.numeric(.cw)
    out$omega   = omega
    out$XYZ     = XYZ
    #   out$jac3x2  = jac3x2 %*%  jac2x2
    
    out$test        = as.numeric(test)
    out$uv          = uv
    out$uv.norm2    = uv.norm2        
    out$jac2x2      = .rotation[2:3, ] %*% jac3x2 %*%  jac2x2 #  ; print(jac2x2) ; print( solve(jac2x2) )

    return(out)
    }
    
        


#   .center     in [0,1] but taken as in the circle
#   .width      in (0,1)

omegaFromCW <- function( .center, .width )
    {
    a   = 0.5 * deltaFromWidth( .width )
    
    omega1  = .center - a   # ; 0.5*delta
    omega2  = .center + a   # ; 0.5*delta
      
    out     = c( omega1 %% 1, omega2 %% 1 ) #; print(out)
    
    if( out[1] == out[2] )
        {
        log.string( ERROR, "Endpoints are equal = %g\n", out[1] )
        return(NULL)
        }

    return( out )
    }

#   deltaFromWidth()
#   a special non-linear function from [0,1] to [0,1] which compresses at the endpoints
deltaFromWidth <- function( .width )
    {
    return( 0.5 * (1 - cos(pi*.width)) )
    }

    
    
#   evaluate spline function on [0,1]    

evaluateD0D1 <- function( .x, .sfun )
    {    
    valid   = 0<=.x  &&  .x<=1
    if( ! valid )   return(NULL)
    
    v   = .sfun( .x )
    dv  = .sfun( .x, deriv=1 )
    
    return( list( v=v, dv=dv ) )
    }
        
    

    
gridPlotF <- function( .gray, .direction, .x.seq=NULL, .y.seq=NULL, .illum=colorSpec::C.5nm, .xyz=colorSpec::xyz1931.1nm )
    {
    if( ! identical( .illum$Wavelength, .xyz$Wavelength ) )
        {
        return(NULL)
        }
        
    wavelengths = length(.illum$Wavelength)
    
    #   make mat.cum and 3 spline interpolators
    mat.cum = matrix( 0, wavelengths, 3 )
    sfun    = list()
    for( j in 1:3 )
        {
        mat.cum[ ,j]    = cumsum( .illum$Power * .xyz[ , 1+j ] )
        
        sfun[[j]]   =   splinefun( (0:wavelengths)/wavelengths, c(0,mat.cum[ ,j]), method='hyman' )      # hyman  monoH.FC  natural
        }
        
    XYZ.white   = mat.cum[ wavelengths, ]   ;  log.object( TRACE, XYZ.white ) ;  log.object( TRACE, XYZ.white / XYZ.white[2] ) 
    
    base            = .gray * XYZ.white
    colnames(base)  = colnames(.direction)        

    rot = t( svd(.direction,3,3)$u )   ; log.object( TRACE, rot )
    
    n   = 32
    
    if( is.null(.x.seq) )
        .x.seq  = (0:n)/n

    if( is.null(.y.seq) )
        .y.seq  = (1:(n-1))/n
    
    uv   = array( as.numeric(NA), dim=c(length(.x.seq), length(.y.seq), 2) )
    
    uv.det   = array( as.numeric(NA), dim=dim(uv)[1:2] )
    
    count.pos   = 0
    count.neg   = 0
    count.0     = 0
    for( x in 1:nrow(uv) )
        {
        for( y in 1:ncol(uv) )
            {
            data    = dataFromCW( c(.x.seq[x], .y.seq[y]), sfun, base, rot )   #; print(data)
            
            if( data$test[1] < 0 )  next    #in the wrong half-plane

            uv[x,y, ]  = data$uv
            
            d   = det( data$jac2x2 )
            uv.det[x,y]  = d
                
            if( d < 0 )
                count.neg = count.neg + 1
            else if( 0 < d )
                count.pos = count.pos + 1             
            else
                {
                count.0 = count.0 + 1
                #print( data$jac2x2 )
                #stop( "jac2x2 is singular !" )
                }
            }
        }
        
    log.object( TRACE, uv, type="str" )
    log.string( TRACE, "count.neg=%d  count.pos=%d   count.0=%d", count.neg, count.pos, count.0 )
    
    
    uvdist   = abs( uv[ , ,1] )  +  abs( uv[ , ,2] )  
    idx     = arrayInd( which.min(uvdist), dim(uvdist) ) #; print(idx)
    
    cw.min  = c( .x.seq[idx[1]], .y.seq[idx[2]] )
    log.object( TRACE, cw.min )
    
    
    plot( c(-30,30), c(-30,30), type='n', asp=1, las=1, xlab='u', ylab='v' )
    abline( h=0, v=0 )
    
    #   plot the cols of uv
    i3 =   (idx[1]-1):(idx[1]+1)   
    j3 =   (idx[2]-1):(idx[2]+1)    
    
    for( j in j3 )      # 1:ncol(uv) )
        {
        #   lines( uv[ , j, 1], uv[ , j,2] )
        
        lines( uv[ i3, j, 1], uv[ i3, j,2], lwd=0.1 )
        }        
            
    #   plot the rows of uv
    mask.col    = logical( ncol(uv.det) )
    mask.col[ j3 ] = TRUE

    for( i in i3 )      #1:nrow(uv) )
        {        
        lines( uv[i, j3, 1], uv[ i, j3, 2], lwd=0.1  )
        # points( uv[ i, j3, 1], uv[ i, j3,2], pch=20, col='red'  )
        
        mask    = uv.det[i, ] < 0
        points( uv[i,mask ,1], uv[i,mask ,2], pch=20, col='cyan' )        
        
        mask    = uv.det[i, ] == 0
        points( uv[i,mask&mask.col,1], uv[i,mask&mask.col,2], pch=20, col='black', cex=10 )        
        
        mask    = 0 < uv.det[i, ]
        points( uv[i,mask&mask.col,1], uv[i,mask&mask.col,2], pch=20, col='red' )        
        }
        
    points( uv[ idx[1], idx[2], 1], uv[ idx[1], idx[2], 2],   col='red' )
    
    
    #   find the quadrilateral that covers 0
    cw.mat  = matrix( c(0,0, 1,0,  1,1, 0,1), 4, 2, byrow=T )
    
    uv.mat  = matrix( 0, 4, 2 )
    for( x in 1:(nrow(uv)-1) )
        {
        for( y in 1:(ncol(uv)-1) )
            {
            #bad = is.na( uv[x,y, ] )  |  is.na( uv[x+1,y, ] )  |  is.na( uv[x,y+1, ] )  |  is.na( uv[x+1,y+1, ] )
            #if( any(bad) )  next
            
            uv.mat[1, ] = uv[x,y, ]
            uv.mat[2, ] = uv[x+1,y, ]
            uv.mat[3, ] = uv[x+1,y+1, ]
            uv.mat[4, ] = uv[x,y+1, ]
            
            if( any(is.na(uv.mat)) )  next   

            #   polygon( uv.mat[1, ], uv.mat[2, ] ) ; next            

            if( windingNumber( c(0,0), uv.mat[ ,1], uv.mat[ ,2] ) != 0 )
                {
                log.string( TRACE, "Found inside:  x=%d  y=%d.  Winding Number != 0.", x, y )
                log.object( TRACE, uv.mat )                           
                polygon( uv.mat[, 1], uv.mat[ ,2], border='red'  )
                points( uv.mat[, 1], uv.mat[ ,2] )       
                
                cw.mat[1, ] = c( .x.seq[x], .y.seq[y] )
                cw.mat[2, ] = c( .x.seq[x+1], .y.seq[y] )    
                cw.mat[3, ] = c( .x.seq[x+1], .y.seq[y+1] ) 
                cw.mat[4, ] = c( .x.seq[x], .y.seq[y+1] ) 

                cw.root = interpQuad(  uv.mat, c(0,0), cw.mat )
                log.object( TRACE, cw.root )
                
                #B = projectiveMatrixGeneral( uv.mat, cw.mat )
                #print( B )
                #print( B %*% rbind( uv.mat, 1 ) ) 
                #cw.root = B[1:2,3] / B[3,3] ; print( cw.root )                
                #stop()
                }
            
            }
        }
    
    return(TRUE)
    }    
    
    
windingNumber <- function( iPoint, iX, iY )
    {
    n = length(iX)  #;    assert( length(iY)==n )
    
    #   look for intersections of edges with ray from iPoint
    
    x = iPoint[1]
    y = iPoint[2]
    
    winding = 0
    
    
    for( i in 1:n )
        {
        i_next = (i %% n) + 1

        y_diff = iY[i] - y
        
        y_diff_next = iY[i_next] - y
        
        if( 0 <= y_diff * y_diff_next ) next
        
        #  find intersection of edge and horizontal line through iPoint
        t = (y - iY[i]) / (iY[i_next] - iY[i])
        
        x_int = iX[i] + t * (iX[i_next] - iX[i])
        
        if( x < x_int )
            winding = winding + sign(y_diff_next)
        }
    
    winding
    }
    
    

testJacobian2x2 <- function( .gray, .direction, .cw0, .delta, .illum=colorSpec::C.5nm, .xyz=colorSpec::xyz1931.1nm )
    {
    if( ! identical( .illum$Wavelength, .xyz$Wavelength ) )
        {
        return(NULL)
        }
        
    wavelengths = length(.illum$Wavelength)
    
    #   make mat.cum and 3 spline interpolators
    mat.cum = matrix( 0, wavelengths, 3 )
    sfun    = list()
    for( j in 1:3 )
        {
        mat.cum[ ,j]    = cumsum( .illum$Power * .xyz[ , 1+j ] )
        
        sfun[[j]]   =   splinefun( (0:wavelengths)/wavelengths, c(0,mat.cum[ ,j]), method='hyman' )      # hyman  monoH.FC  natural
        }
        
    XYZ.white   = mat.cum[ wavelengths, ]   ;  print( XYZ.white ) #;  print( XYZ.white / XYZ.white[2] ) 
    
    base            = .gray * XYZ.white    

    rotation    = t( svd(.direction,3,3)$u )   ; print( rotation )
            
    data0   = dataFromCW( .cw0, sfun, base, rotation  )
    
    cat( "\nAnalytical Jacobian2x2:\n" )
    print( data0$jac2x2 )

    #   now compute numerical estimate
    jac2x2  = matrix(0,2,2)
    
    data0.c = dataFromCW( .cw0 + c(.delta,0), sfun, base, rotation  )
    jac2x2[ ,1] = (data0.c$uv - data0$uv) / .delta

    data0.w = dataFromCW( .cw0 + c(0,.delta), sfun, base, rotation  )
    jac2x2[ ,2] = (data0.w$uv - data0$uv) / .delta

    cat( "\nNumerical Jacobian2x2:\n" )
    print( jac2x2 )
    
    return(TRUE)
    }
    
    
    
    
    
testIntegral <- function( .spectrum=c(1,2,3,4) )
    {
    x   = (0:100)/100
    v   = x
    
    csum    = cumsum(.spectrum)
    
    for( k in 1:length(x) )
        {
        res = evaluateD0D1( x[k], csum )
        
        v[k]    = res$v
        }
    
    plot( x, v, type='l' )
    }
               
               

evalAtCW <- function( .gray, .direction, .cw, .illum=colorSpec::C.5nm, .xyz=colorSpec::xyz1931.1nm ) 
    {
    wavelengths = length(.illum$Wavelength)
    
    #   make mat.cum and 3 spline interpolators
    mat.cum = matrix( 0, wavelengths, 3 )
    sfun    = list()
    for( j in 1:3 )
        {
        mat.cum[ ,j]    = cumsum( .illum$Power * .xyz[ , 1+j ] )
        
        sfun[[j]]   =   splinefun( (0:wavelengths)/wavelengths, c(0,mat.cum[ ,j]), method='hyman' )      # hyman  monoH.FC  natural
        }
        
    print( str(sfun) )
        
    XYZ.white   = mat.cum[ wavelengths, ]   ;  print( XYZ.white ) ;  print( XYZ.white / XYZ.white[2] ) 
    
    base            = .gray * XYZ.white
    
    rot = t( svd(.direction,3,3)$u )   ; print( rot )
        
    out    = dataFromCW( .cw, sfun, base, rot )
    
    return( out )
    }
    
     
    
    
#   .data   as returned from twoTransitionGrid.colorSpec()

plotParallels <- function( .data, .index=NULL )
    {
    dim         = dim( .data$XYZ )
    stopifnot( length(dim)==3  &&  dim[3]==3 )
    
    parallels   = dim[2]
    
    #   at every meridian the top parallel should be exactly the same
    XYZ.white   = .data$XYZ[ 1, parallels, ]    ; print( XYZ.white )
    
    rot = svd(XYZ.white,3,3)$u  ; print( rot )
    
    if( is.null(.index) )
        .index  = as.integer( parallels/2 ) # equator
    
    equator = .data$XYZ[ , parallels/2,  ,drop=T] ;  print( str(equator) )
    
    uv  = equator %*% rot[ , 2:3 ]  
    
    xlim    = range( uv[ ,1] )
    ylim    = range( uv[ ,2] )
    
    plot( xlim, ylim, type='n', las=1, xlab='u', ylab='v', asp=1 )
    abline( h=0, v=0 )
    
    for( k in 1:length(.index) )
        {
        print( .data$w.seq[ .index[k] ] )
        
        print( .data$XYZ[ , .index[k], ] )
        
        uv  = .data$XYZ[ , .index[k], ] %*% rot[ , 2:3 ]
        
        #   print( uv )
        lines(  uv[ ,1],  uv[ ,2] )
        points(  uv[ ,1],  uv[ ,2] )        
        }
    
    return( TRUE )
    }
    

plotOptimals3D.colorSpec <- function( x, size=c(33,33) )
    {
    theName = deparse( substitute(x) )
    
    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(FALSE)
        }        
    
    if( numSpectra(x) != 3 )
        {
        log.string( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )
        return(FALSE)
        }     
    
    if( ! requireNamespace( 'rgl', quietly=TRUE ) ) 
        {
        log.string( ERROR, "Package 'rgl' is required.  Please install it." )        
        return(FALSE)
        }
        
    #   attachNamespace( 'rgl' )
    
    grid    = twoTransitionGrid.colorSpec( x, size )    
    
    if( is.null(grid) ) return(FALSE)
    
    dim         = dim( grid$XYZ )
    
    ok  = length(dim)==3  &&  dim[3]==3
    if( ! ok )
        {
        log.string( ERROR, "Cannot compute the 2-Transition grid." )
        return(FALSE)
        }

    meridians   = dim[1]
    parallels   = dim[2]
    
    #   at every meridian the top parallel should be exactly the same
    XYZ.white   = grid$XYZ[ 1, parallels, ]    #; print( XYZ.white )
    center  = 0.5 * XYZ.white
    
    #   start 3D drawing
    rgl::bg3d("gray50")
    rgl::light3d()
    
    
    cube    = rgl::scale3d( rgl::cube3d(col="white"), center[1], center[2], center[3] )
    cube    = rgl::translate3d( cube, center[1], center[2], center[3]  )
    rgl::wire3d( cube, lit=FALSE )    
    
    #   exact diagonal of box        
    rgl::lines3d( c(0,XYZ.white[1]), c(0,XYZ.white[2]), c(0,XYZ.white[3]), col=c('black','white'), lwd=3, lit=FALSE )    
    
    label   = toupper( specnames(x) )
    cex     = 1.5
    rgl::text3d( center[1], 0, 0, label[1], col='white', cex=cex )
    rgl::text3d( 0, center[2], 0, label[2], col='white', cex=cex )
    rgl::text3d( 0, 0, center[3], label[3], col='white', cex=cex )
    
    radius  = 0.01 * sqrt( sum(XYZ.white^2) )
    
    ball    = rgl::scale3d( rgl::icosahedron3d(), radius, radius, radius )
    rgl::shade3d( ball, lit=T, col='black' )    
    rgl::shade3d( rgl::translate3d( ball, XYZ.white[1], XYZ.white[2], XYZ.white[3] ), lit=T, col='white'  )
    rgl::wire3d( rgl::translate3d( ball, center[1], center[2], center[3] ), lit=F, col='black'  )
    
    #   plot the meridians
    for( i in 1:(meridians-1) )
        {
        rgl::lines3d( grid$XYZ[i, ,1],  grid$XYZ[i, ,2],  grid$XYZ[i, ,3] )
        }
        
    #   plot the parallels
    for( j in 2:(parallels-1) )
        {
        rgl::lines3d( grid$XYZ[ ,j,1],  grid$XYZ[ ,j,2],  grid$XYZ[ ,j,3] )
        }
    
    return( invisible(TRUE) )
    }
    
computeADL.colorSpec <- function( x, response )
    {
    theName = deparse( substitute(x) )
    
    if( type(x) != "responsivity.material" )
        {
        log.string( ERROR, "type(%s) = '%s' != 'responsivity.material'", theName, type(x) )
        return(NULL)
        }        
    
    if( numSpectra(x) != 3 )
        {
        log.string( ERROR, "numSpectra(%s) = %d != 3", theName, numSpectra(x) )
        return(NULL)
        }     
        
    if( ! is.numeric(response) )
        {
        log.string( ERROR, "response must be numeric, but typeof(response) = '%s'.", typeof(response) )
        return(NULL)
        }     
    
    if( is.null(dim(response))  &&  length(response) == 3 )
        dim(response)   = c(1,3)
        
    ok  = length(dim(response))==2  &&  ncol(response)==3
    if( ! ok )
        {
        log.string( ERROR, "response is not a matrix with 3 columns." )
        return(NULL)
        }     
    
            
    response.white  = step.wl(x) * colSums( as.matrix(x) )
    
    gray    = 0.5 * response.white  #;    print(gray)
    
    direction   = response - matrix( gray, nrow(response), 3, byrow=T )
    
    data    = probeOptimalColors( x, 0.5, direction )
    
    if( is.null(data) ) return(NULL)
    
    #print( data )
    
    ADL = cbind( 1/data$s, data$dol[ ,1], data$dol[ ,3] )
    rownames(ADL)   = NULL
    colnames(ADL)   = c('alpha','delta','lambda')
    class(ADL)      = "model.matrix"    
    
    colnames(response)  = toupper( specnames(x) )
    class(response)     = "model.matrix"

    out = data.frame( response=response, ADL=ADL, omega=data$dol[ ,2],
                        lambda=data$lambda, row.names=rownames(response) )  # as.data.frame.model.matrix
    
    return( out )
    }
    
    
        
    
    
        
#   .obj    a colorSpec with responsivity, lambda in [lambda_min.lambda_max]
#   returns list with
#       omega.from.lambda()
#       lambda.from.omega()
#       responsivity.from.omega[[ ]]    a list with length = channels
#       integral.from.omega[[ ]]        a list with length = channels
#   omega is in [0,1]

makeFunctionList <- function( .obj, .name, .mode='equalize' )
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
    
    step.wl     = step.wl(.obj)
    
    lambda.center   = wavelength(.obj)
    lambda.min      = lambda.center[1]  
    lambda.max      = lambda.center[n]  
        
    #   wavelengths defining the bins
    lambda.break    = seq( lambda.min - step.wl/2, lambda.max + step.wl/2, len=(n+1) )        
        
    out = list()
           
    if( .mode == 'equalize' )
        {
        omega   = sqrt( rowSums(coredata*coredata) )
        
        if( any(omega == 0) )
            {
            log.string( ERROR, "Scanner '%s' has 0 responsivity=0 at 1 or more wavelengths.", .name )
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
        log.string( ERROR, ".mode='%s' is invalid.", .mode )
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
    
plotReparam <- function( .obj, .mode='equalize', .full=F )
    {
    par( mfcol= c(2,1) )
    
    #   top plot is the normal plot of .obj
    plot( .obj )
    
    
    #   bottom plot has the reparameterized responses
    theList = makeFunctionList( .obj, deparse(substitute(.obj)), .mode )    #  ; print( str(theList) )
    
    if( is.null(theList) )  
        {
        par( mfcol= c(1,1) )        
        return(NULL)
        }
    
    spectra = numSpectra( .obj )

    omega   = theList$omega     # (0:200)/200
    
    ylim    = 0
    for( j in 1:spectra )
        {
        ylim    = range( theList$integral.from.omega[[j]]( omega, deriv=1 ), ylim ) 
        }
        
    plot.default( c(0,1), ylim, type='n', las=1, xlab='', ylab='', lab=c(10,8,7),
                             tcl=0, mgp=c(3, 0.25, 0) )
    title( xlab=expression(omega), line=1.5 )
    title( ylab='Responsivity', line=2 )                             
    grid( lty=1 )
    abline( h=0, v=0 )
    
    
    #   fake the color by changing quantity to power
    obj.copy  = .obj
    quantity(obj.copy) = 'power'
    mat.rgb = product( obj.copy, colorSpec::BT.709.RGB, wave='auto' )  #; print( mat.rgb )
    mat.rgb = mat.rgb / max(mat.rgb)                        #; print( mat.rgb )    
    mat.rgb = DisplayRGBfromLinearRGB( mat.rgb )            #; print( mat.rgb )
    color_vec   = rgb( mat.rgb )  
    
    for( j in 1:spectra )
        {
        #   y   = theList$responsivity.from.omega[[j]]( omega ) ; print( range(y) )
        
        y   = theList$integral.from.omega[[j]]( omega, deriv=1 )    #; print( range(y) )
        
        lines( omega, y, col=color_vec[j] ) 
        #   points( omega, y, col=color_vec[j] )     
        }
        
    legend( "topright", specnames(.obj), col=color_vec, bty='n', lwd=9 )
    
    
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
    
    
    
#--------       UseMethod() calls           --------------#            
              
probeOptimalColors <- function( x, gray, direction, tol=1.e-6, aux=FALSE )
    {
    UseMethod("probeOptimalColors")
    }
        
plotOptimals3D <- function( x, size=c(33,33) )
    {
    UseMethod("plotOptimals3D")
    }
                           
computeADL <- function( x, response )
    {
    UseMethod("computeADL")
    }
 
        
        
############################################        deadwood below  ###########################################        

#   .cw0        center of grid search
#   .radius     half-size of grid, both directions
#   .step       step size of grid, both directions    
#   .base       of the search ray
#   .rotation   orthogonal frame at the .base, with direction of ray in 1st row
#   .sfun       spline function integral of xyz responsivity    

#   returns    cw_min
#       

gridSearchCW <- function( .cw0, .radius, .step, .sfun, .base, .rotation )
    {
    c_seq   = seq( .cw0[1] - .radius[1]*.step[1],  .cw0[1] + .radius[1]*.step[1], by=.step[1] )
    w_seq   = seq( .cw0[2] - .radius[2]*.step[2],  .cw0[2] + .radius[2]*.step[2], by=.step[2] )
    
    #   only allow valid w's
    w_seq   = w_seq[ 0 < w_seq  &  w_seq < 1 ]
    
    uv   = array( as.numeric(NA), dim=c(length(c_seq), length(w_seq), 2) )
    
    for( x in 1:length(c_seq) )
        {
        for( y in 1:length(w_seq) )
            {
            data    = dataFromCW( c(c_seq[x], w_seq[y]), .sfun, .base, .rotation )

            if( data$test[1] < 0 )  next    #in the wrong half-plane

            uv[x,y, ]  = data$uv
            }
        }
    
    cw.mat  = matrix( 0, 4, 2 )
    uv.mat  = matrix( 0, 4, 2 )
    for( x in 1:(nrow(uv)-1) )
        {
        for( y in 1:(ncol(uv)-1) )
            {
            uv.mat[1, ] = uv[x,y, ]
            uv.mat[2, ] = uv[x+1,y, ]
            uv.mat[3, ] = uv[x+1,y+1, ]
            uv.mat[4, ] = uv[x,y+1, ]
            
            if( any(is.na(uv.mat)) )  next   

            #   polygon( uv.mat[1, ], uv.mat[2, ] ) ; next            

            w   = windingNumber( c(0,0), uv.mat[ ,1], uv.mat[ ,2] )
            
            if( w == 0 )    next
            
            mess = sprintf( "Found inside:  x=%d  y=%d.  Winding Number = %d\n", x, y, w ) ; cat(mess)
            print( uv.mat )                            
            
            cw.mat[1, ] = c( c_seq[x], w_seq[y] )
            cw.mat[2, ] = c( c_seq[x+1], w_seq[y] )    
            cw.mat[3, ] = c( c_seq[x+1], w_seq[y+1] ) 
            cw.mat[4, ] = c( c_seq[x], w_seq[y+1] ) 

            cw.root = interpQuad(  uv.mat, c(0,0), cw.mat ) ; print( cw.root )
            
            return( cw.root )
            #B = projectiveMatrixGeneral( uv.mat, cw.mat )
            #print( B )
            #print( B %*% rbind( uv.mat, 1 ) ) 
            #cw.root = B[1:2,3] / B[3,3] ; print( cw.root )                
            #stop()
            }
        }
    
    return( NULL )
    }
    
    
                   