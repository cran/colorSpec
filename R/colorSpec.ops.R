

#   x           object of class colorSpec
#   wavelength  the *new* wavelengths
#   method      'auto',  'sprague', 'spline', 'loess', or 'linear'
#   span        smoothing factor passed to loess().  0 means no smoothing
#
#   returns     a  colorSpec object with the same organization

resample.colorSpec <-  function( x, wavelength, method='auto', span=0.02 )
    {
    #   partial matching
    table   = c('auto','sprague','spline','loess','linear')
    k   = pmatch( tolower(method), table )   
    if( is.na(k) )
        {
        log.string( ERROR, "method='%s' does not uniquely match any of '%s'.",
                            method, paste(table,collapse=',' ) )
        return(NULL)
        }
         
    method  = table[k]
    
    if( method=='loess'  &&  span<=0 )
        {
        log.string( WARN, "method='%s' with span=%g is invalid; changed to method='auto'.",
                            method, span )
        method = 'auto'
        }
        
    if( method == 'auto' )
        #   use CIE recommendation
        method  = ifelse( is.regular(x), 'sprague', 'spline' )
    
    wave    = wavelength(x)
        
    if( method!='loess' &&  identical(wavelength,wave) )
        {
        #   no change !   nothing to do !
        #   log.string( TRACE, "new wavelengths are identical to current wavelengths, and no smoothing. So nothing to do." )
        return(x)
        }
        
    if( method=='sprague' &&  ! is.regular(x) )
        {
        log.string( ERROR, "Sprague interpolation cannot be used with irregular wavelengths." )
        return(NULL)
        }
                
    if( ! isStrictlyIncreasingSequence(wavelength) )
        {
        log.string( ERROR, "New wavelength sequence is not strictly increasing." )
        return(NULL)
        }
        
    #   count the number of extrapolated
    extrap  = (wavelength < wave[1])  |  (wave[length(wave)] < wavelength)
    if( 0 < sum(extrap) )
        {
        log.string( TRACE, "For object %s, wavelength extrapolation occurred at %d points (outside [%g,%g] nm).",
                            deparse(substitute(x)), sum(extrap), wave[1], wave[length(wave)] )
        }
        
    mat.in  = as.matrix( x )
    
    theType = type(x)
    
    if( method == 'sprague' )
        {
        #   convolution weights are complicated to compute
        #   so best to apply to all columns in mat.in at the same time

        mat.out = interpSprague( mat.in, wave[1], step.wl(x), wavelength )
        
        if( is.null(mat.out) )  return(NULL)
        }
    else
        {
        time.start  = as.double( Sys.time() )
        
        mat.out = matrix( NA_real_, length(wavelength), ncol(mat.in) )
    
        for( k in 1:ncol(mat.out) )
            {
            y = mat.in[ , k]
            
            if( is.na(y[1]) )   next
            
            mat.out[ ,k]    = resampleXY( wave, y, wavelength, method, span )
            }
            
        #   print( c( "resampleXY().  elapsed: ", as.double(Sys.time()) - time.start) )        
        }
        
    #   check for negative undershoot in mat.out[ , ]
    for( k in 1:ncol(mat.out) )
        {
        if( is.na( mat.out[1,k]) )  next    #   resampleXY() failed
            
        y = mat.in[ , k]        
        
        if( all( -1.e-8 * max(y) <= y ) )
            #   clip to 0
            mat.out[ ,k]   = pmax( mat.out[ ,k], 0 )

        #if( theType == 'material' )
        #    #   prevent overshoot
        #    data.out[ ,k]   = pmin( data.out[ ,k], 1 )            
        }
        
        
    colnames( mat.out )   = specnames( x )        
        
    org = organization(x) 

    if( org == "df.row" )
        {
        #   preserve the initial columns, and change the final one
        #   an easy overwrite in this case
        mat.out             = t( mat.out )   
        class( mat.out )    = "model.matrix"
        out = x
        out[[ ncol(out) ]]  = mat.out
        
        wavelength(out)     = wavelength
        }
    else
        {
        #   create a new object
        out = colorSpec( mat.out, wavelength, quantity(x), org )
        }
        
    attr( out, "metadata" )  = attr( x, "metadata" )
    metadata(out)   = list( resampled=TRUE )
    
    if( method=='loess'  &&  0 < span )
        metadata(out)   = list( span=span )
    
    return( out )
    }
    
    
#   .x      vector of current wavelengths, length n
#   .y      current values,  n x q matrix
#   .xnew   new wavelengths, length n'
#   .span   smoothing parameter
#
#   value   new matrix,  n' x q

resampleXY <- function( .x, .y, .xnew, .method, .span )
    {
    #   if( any( is.na(.y) ) )  return( rep(NA_real_,length(.xnew)) )
    
    if( .method == 'spline' )
        {
        #   use simple spline
        out = spline( .x, .y, xout=.xnew, method="natural" )$y
        
        if( .x[1] == .x[2]  &&  .y[1] == .y[2] )
            #   duplicated knot => custom extrapolation - a constant
            out[ .xnew < .x[1] ] = .y[1]
            
        n   = length(.x)
        
        if( .x[n-1] == .x[n]  &&  .y[n-1] == .y[n] )
            #   duplicated knot => custom extrapolation - a constant
            out[ .x[n] < .xnew ] = .y[n]
        }
    else if( .method == 'loess' )
        {
        #   use loess smoother:
        #df = data.frame( X=.x, Y=.y )
        xy.lo = try( loess( .y ~ .x, span=.span ) )     #, family="symmetric" )
        
        if( class(xy.lo) == "try-error" )        
            {
            log.string( WARN, "loess smoothing with span=%g failed !  Probably span is too small.  Returning all NA.", .span )
            return( rep(NA_real_,length(.xnew)) )       #return( resampleXY( .x, .y, .xnew, .span=0 ) )
            }
        
        out = predict( xy.lo, .xnew )
        }
    else if( .method == 'linear' )
        {
        out = approx( .x, .y, .xnew )$y
        }
        
    return( out )
    }
        
                
        
#   .list       a list of colorSpec objects, with names        
#   value       TRUE iff all objects have the same quantity and wavelength, and distinct specnames

areSpectraBindable <- function( .list )
    {
    if( ! is.list( .list ) )    return(FALSE)
    
    n   = length(.list)
    
    if( n == 0 )    return(TRUE)
    
    for( k in 1:n )
        {
        if( ! is.colorSpec( .list[[k]] ) )   
            {
            log.string( ERROR, "The list of %d spectra are not bindable, because '%s' is not a valid colorSpec object.", 
                                    n, names(.list)[k] )
            return(FALSE)
            }
        }

    if( n == 1 )    return(TRUE)
    
    #   identical quantities
    qvec    = sapply( .list, quantity.colorSpec )
    
    if( 2 <= length(unique(qvec)) ) 
        {
        log.string( ERROR, "The list of %d spectra are not bindable, because they do not have the same quantity.", n )
        return(FALSE)
        }
        
        
    #   identical wavelengths
    wave    = wavelength( .list[[1]] )
    
    for( k in 2:n )
        {
        if( ! identical( wavelength( .list[[k]] ), wave ) ) 
            {
            log.string( ERROR, "The list of %d spectra are not bindable, because they do not have the same wavelengths.", n )
            return(FALSE)
            }
        }
        
    #   duplicated specnames
    namevec = unlist( sapply( .list, specnames.colorSpec ) )
    if( any(duplicated(namevec)) ) 
        {
        log.string( ERROR, "The list of %d spectra are not bindable, because the specnames are not all distinct.", n )
        return(FALSE)
        }
        
    return(TRUE)
    }
        
        
    
#   ...     colorSpec objects with identical wavelength and quantity
#
#   returns new colorSpec with spectra combined
#           organization is the most complex from the inputs
#           metadata is taken from the 1st spectrum
bind.colorSpec  <-  function( ... ) #, .removedups=FALSE )
    {
    theList =  list( ... ) 
    
    n   = length(theList)
    if( n == 0 )
        {
        log.string( ERROR, "No arguments." )
        return(NULL)
        }
        
    log.string( TRACE, "Found %d objects in '...'", n )
            
    theNames = as.character( substitute(list(...)) )    # ; print( theNames )
    if( length(theNames) == n+1 )
        theNames    = theNames[ 2:(n+1) ]
    else
        {
        log.string( WARN, "length(theNames) = %d != %d.  Using fake names.", length(theNames), n+1 )
        theNames    = sprintf( "Name%d", 1:n )
        }
        
    names(theList)  = theNames  #; print( str(theList) )

    return( bindSpectra( theList ) )
    }
    
    
#   .list   of colorSpec objects, with names    
bindSpectra <- function( .list )
    {
    if( ! areSpectraBindable(.list) ) return(NULL)    
     
    n   = length( .list )
    
    if( n == 1 )    return( .list[[1]] )  # nothing to do !   

    #   find the output organization
    orgvec  = sapply( .list, organization.colorSpec )
    
    for( org in c('df.row','df.col','matrix','vector') )
        {
        if( any( orgvec == org ) )  break
        }
    if( org == "vector" )   
        org = "matrix"  # the bind must have more than 1 spectrum in it !
        

    mat = coredata( .list[[1]], forcemat=T )
    
    specnames   = specnames( .list[[1]] )
    extradata   = extradata( .list[[1]] )

    for( k in 2:n )
        {
        mat = cbind( mat, coredata( .list[[k]] ) )
        
        specnames   = c( specnames, specnames( .list[[k]] ) )
        
        if( org == 'df.row' )
            extradata   = rbind.super( extradata, extradata( .list[[k]] ) )
        }
    
    colnames(mat)   = specnames
    
    out     = colorSpec( mat, wavelength( .list[[1]] ), quantity( .list[[1]] ), organization=org )
        
    if( org == 'df.row' )        
        extradata(out)  = extradata
        
    metadata(out)   = metadata( .list[[1]] )
    
    return( out )
    }
    
    
#   subset.colorSpec()
#
#   x       a colorSpec object
#   subset  a set of integer indexes, no duplicates
#           a logical mask with length(subset) = #(spectra in x)
#           a regular expression, matching the specnames, case insensitive
#           
#   returns a subset of the spectra in x
subset.colorSpec  <-  function( x, subset, ... )
    {   
    spectra = numSpectra(x)
    
    if( spectra == 0 )  return(x)
    
    if( is.logical(subset) )
        {
        if( length(subset) != spectra )
            {
            log.object( ERROR, subset )
            log.string( ERROR, "subset is logical, and length(subset) = %d != %d spectra.", length(subset), spectra )
            return( NULL )
            }    
            
        subset  = which( subset )
        }    
    else if( is.character(subset) )
        {
        # interpret subset as a regular expression        
        subset    = which( grepl( subset, specnames(x), ignore.case=T ) )
        }   
        
    if( is.numeric(subset) )
        {
        subset = as.integer(subset)
        
        if( anyDuplicated(subset) )        
            {
            log.object( ERROR, subset )
            log.string( ERROR, "subset indexes are invalid.  Duplicates are not allowed."  )
            return( NULL )
            }
            
        ok      = all( 1L <= subset & subset <= spectra )
        if( ! ok )
            {
            log.object( ERROR, subset )
            log.string( ERROR, "subset indexes are invalid.  One or more are outside the interval [%d,%d]", 1, spectra )
            return( NULL )
            }
        }
    else
        {
        log.object( ERROR, subset )
        log.string( ERROR, "subset argument is invalid." ) 
        return(NULL)        
        }

    org = organization(x)   
    
    if( org == "vector" )
        {
        #   length(subset) must be 0 or 1
        if( length(subset) == 1 )  return(x)   # no change
        org = "matrix"          # force empty matrix return
        }
    
    mat = as.matrix( x )
        
    mat = mat[ , subset, drop=F ]  # the actual subsetting happens here   print( colnames(mat) )

    #   colnames(mat)   = specnames(x)[subset]  previous line does this
    
    out = colorSpec( mat, wavelength(x), quantity(x), org )

    extradata(out)  = extradata(x)[ subset, , drop=F]   # the actual subsetting happens here, and is in synch with the one in mat[,]
    metadata(out)   = metadata(x)
            
    metadata(out)   = list( subsetted=TRUE )
    
    return( out )
    }
            
        

    
#   x   colorSpec object
#   returns mean of all spectra with 'vector' organization
mean.colorSpec <-  function( x, ...  )
    {
    spectra     = numSpectra( x )

    if( spectra <= 1 )  return( x )
    
    wavelength  = wavelength( x )
        
    mat = coredata( x, forcemat=T )
    
    vec = rowMeans( mat )


    out = colorSpec( vec, wavelength, quantity(x), 'vector' )
    
    specnames(out)  = sprintf( "mean.%s", deparse(substitute(x) ) )
        
    metadata(out)   = metadata(x)
    metadata(out)   = list( samples=sprintf( "mean of %d spectra", spectra ) )
    
    return( out )
    }
    
    
#   x       colorSpec object
#   norm    desired norm
#   returns colorSpec, with all spectra scaled to have norm 1

normalize.colorSpec  <-  function( x, norm='L1'  )
    {
    coremat = coredata( x, forcemat=T )
    
    if( is.character(norm) )
        {
        step.wl = step.wl(x)
        
        if( grepl( '1', norm ) )
            normvec = step.wl * colSums( abs(coremat) )
        else if( grepl( '2', norm ) )
            normvec = step.wl * sqrt( colSums( coremat^2 ) )
        else if( grepl( 'inf', norm, ignore.case=T ) )
            normvec = base::apply( abs(coremat), 2, max )       #  fun <- function( y ) {  y / max(abs(y)) }
        else
            {
            log.string( ERROR, "norm = '%s' is invalid.", norm )
            return( x )
            }
        }
    else if( is.numeric(norm) )
        {
        #   interpret as a wavelength
        wavelength  = wavelength(x)
        
        i   = which( wavelength == norm )
        if( length(i) == 0 )
            {
            log.string( ERROR, "norm = %g nm is an invalid wavelength.", norm )
            return( x )
            }        
        normvec = coremat[i, ]
        }
    else
        {
        log.string( ERROR, "norm = '%s' is invalid.", as.character(norm) )
        return( x )
        }        
        
    #   prevent division by 0
    normvec[ normvec==0 ]   = 1
        
    out = multiply( x, 1/normvec )
    
    return( out )
    }
    
    

#   linearize.colorSpec()
#
#   force any spectrum to be ready for colorimetry
#   At this time, the only conversion is absorbance to transmittance
#
linearize.colorSpec <- function( x )
    {
    quantity    = quantity( x )
    
    if( quantity == 'absorbance' )
        {
        #   log.string( TRACE, "Converting '%s' from 'absorbance' to 'transmittance'.", deparse(substitute(x)) )
        myfun <- function( y )  { 10^(-y) }
 
        quantity    = "transmittance"  
        }
    else
        {
        return( x )  #  no change needed
        }
        
    #   apply the function to all spectra
    out = applyspec.colorSpec( x, myfun )   
    
    #   and change the quantity
    quantity(out)   = quantity
    
    metadata(out)   = metadata(x)
        
    return( out )
    }
    
    
#   x   a colorSpec object with N spectra
#   interval   vector with 2 wavelength values - giving the blending interval for lo and hi
#   adj         adjustment parameter in [0,1]
#   returns:    a colorSpec object with 2*N spectra:  y.lo, y.hi, ...
#               where .y = y.lo + y.hi
chop.colorSpec <- function( x, interval, adj=0.5 )    
    {    
    specnames   = specnames(x) 
    
    spectra     = length(specnames)
    
    if( spectra == 0 )
        {
        log.string( ERROR, "'%s' has 0 spectra !", deparse(substitute(x)) )
        return(x)
        }
        
        
    #theList =  list(...)     
    #n   = length(theList)

    wave    = wavelength(x)
        
    i1  = which.min( abs(interval[1]-wave) )
    i2  = which.min( abs(interval[2]-wave) )
    
    if( i2 - i1 < 2 )
        {
        log.string( ERROR, ".interval endpoints %g and %g are too close (or swapped).", 
                        interval[1], interval[2] )
        return(NULL)
        }    
        
    out = matrix( 0, length(wave), 2 * spectra )

    core    = coredata(x)
    
    for( j in 1:spectra )
        {
        mat =   splitSpectrum( core[ ,j], c(i1,i2), adj )   
        if( is.null(mat) )  return(NULL)
        
        pair    = c(2*j-1,2*j)
        out[  , pair ]  = mat
        colnames( out )[ pair ] = c( sprintf("%s.lo",specnames[j]), sprintf("%s.hi",specnames[j]) )
        }
    
    out = colorSpec( out, wave, quantity(x) )
    
    return( out )
    }
    
    
    
    
    
    
#--------       UseMethod() calls           --------------#            

        
#   x   a colorSpec object with M spectra
#   s   a scalar
#       an M-vector
#       an MxM matrix

multiply.colorSpec   <-  function( x, s )
    {
    if( ! is.numeric(s) )
        {
        log.string( ERROR, "s is not numeric. type(s)='%s'", typeof(s) )
        return(x)
        }
        
    #   print( length(s) )
    
    spectra = numSpectra(x) 
    
    if( length(dim(s)) == 2 )
        ok  = nrow(s) == spectra
    else
        ok  = length(s)==1  ||  length(s)==spectra

    if( ! ok )
        {
        log.string( ERROR, "Size of s is invalid for %d spectra.", spectra )
        return(x)
        }
    

    
    org = organization(x)
    
    if( length(s) == 1 )
        {
        #   simple scalar multiplication
        out = x        
        
        if( org == 'vector'  ||  org == 'matrix' )
            out = s * out
        else if( org == 'df.col' )
            out[ 2:ncol(out) ]  = s * out[ 2:ncol(out) ]
        else if( org == 'df.row' )        
            out[ ncol(out) ]    = s * out[ ncol(out) ]
            
        return( out )
        }
        

    if( length(dim(s)) == 2 )
        mat = s
    else
        {
        mat = diag(s)
        colnames(mat)   = specnames(x)
        }

    if( ncol(mat) == spectra )
        {
        #   mat is square, so avoid unpacking and repacking when possible
        out = x
    
        if( org == 'matrix' )
            {
            out = out %*% mat
            colnames(out)   = specnames(x)
            out = colorSpec( out, wavelength(x), quantity(x), "matrix" )
            }
        else if( org == 'df.col' )
            out[ 2:ncol(out) ]  = as.matrix.data.frame( out[ 2:ncol(out) ] ) %*% mat
        else if( org == 'df.row' )        
            out[[ ncol(out) ]]    = t(mat) %*% out[[ ncol(out) ]]
        }
    else
        {
        #   mat is not square, so the number of spectra in the output is different
        #   must unpack and repack
        core    = coredata( x, forcemat=T ) %*% mat
        
        out     = colorSpec( core, wavelength(x), quantity(x), org )
        }
        
    return(  out  )
    }    

    

    
applyspec.colorSpec <- function( x, FUN, ... )    
    {
    mat = as.matrix( x )    #; print( str(mat) )

    out = base::apply( mat, 2, FUN, ... )     #;  print( str(mat) )

    if( nrow(out) != nrow(mat) )
        {
        log.string( ERROR, "Function FUN mapped %d-vector to a %d-vector.", nrow(mat), nrow(out) )
        return(NULL)
        }
    
    out = colorSpec( out, wavelength(x), quantity(x), organization(x) )

    if( organization(out) == "df.row" )
        extradata(out)  = extradata(x)    
        
    for( a in c('metadata','sequence','calibration') )
        attr(out,a) = attr(x,a)    
    
    return(out)
    }
    

    
convolvewith.colorSpec <- function( x, coeff )    
    {
    if( is.character(coeff) )
        {
        if( coeff == "SS3" )
            coeff = c(-1,14,-1)/12
        else if( coeff == "SS5" )
            coeff = c(1,-12,120,-12,1)/98
        else
            {
            log.string( ERROR, "Unknown coeff='%s'.", coeff )
            return(NULL)
            }
        }

    k   = length(coeff) 
    ok  = is.numeric(coeff)   &&  (k %% 2L == 1L)
    
    if( ! ok )
        {
        log.string( ERROR, "coeff is not a numeric vector of odd length.  length=%d.", length(coef) )
        return(NULL)
        }
        
    if( k == 1L )
        {
        #   not likely
        return( multiply(x,coeff) )
        }

    #   out = applyspec( x, stats::filter, filter=coeff,  method='convolution', sides=2 )  this works
    

    half    = as.integer( k/2 )

    myfun   <- function( y )    
        { 
        y   = stats::filter( y, filter=coeff,  method='convolution', sides=2 )
        
        #   use constant extrapolation at endpoints        
        y[1:half] = y[half+1]
        
        n   = length(y)
        
        y[ (n-half+1):n ] = y[ n-half ]
        
        return(y) 
        }
    
    out = applyspec( x, myfun )
    
    return(out)
    }
    
    
        
#--------       UseMethod() calls           --------------#            
        
        
resample <- function(  x, wavelength, method='auto', span=0.02 )
    {
    UseMethod("resample")
    }
    
bind <- function( ... )
    {
    UseMethod("bind")
    }    
    
chop <- function( x, interval, adj=0.5  )    
    {       
    UseMethod("chop")
    } 
   
multiply <- function( x, s )
    {
    UseMethod("multiply")
    }    
    
applyspec <- function( x, FUN, ... )        
    {
    UseMethod("applyspec")
    }    
            
normalize <- function( x, norm='L1'  )            
    {
    UseMethod("normalize")
    }        
    
linearize <- function( x )            
    {
    UseMethod("linearize")
    }            
 
convolvewith <- function( x, coeff )          
    {
    UseMethod("convolvewith")
    }             

 
        