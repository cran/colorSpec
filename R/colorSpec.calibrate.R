

#
#   x           a colorSpec object 
#               with type  'responsivity.light'  or  'responsivity.material'
#               with M spectra (output channels)
#
#   stimulus    a single spectrum input, of type 'light' or 'material' respectively
#
#   response    a vector with length M.  Must be all positive.
#
#   method      'scaling', 'Bradford', or 'Von Kries', or 'MCAT02', or an MxM matrix
#
#   return:     a new colorSpec object, with same size etc.
#               always = multiply( x, mat )    where mat is an MxM matrix
#               In case of ERROR, it returns x unchanged !
#
#   if an argument value is NULL, the function attempts to find an appropriate default
#
#   another possible name for "calibrate()" is "whitebalance()"

calibrate.colorSpec <- function( x, stimulus=NULL, response=NULL, method=NULL )      
    {
    theType = type(x)
    
    ok  = theType=='responsivity.light'  ||  theType=='responsivity.material'
    if( ! ok )
        {
        log.string( ERROR, "type(x) = '%s' is invalid.", theType )
        return(x)
        }
        
    m       = numSpectra(x)
    wave    = wavelength(x)
    
    is.neural   = grepl( 'neural$', quantity(x) )

    is.XYZ      = is.neural  &&  m==3  &&  all( pmatch( c('x','y','z'), tolower(specnames(x))) == 1:3 )
    
    #   change NULLs to appropriate defaults
    if( is.null(stimulus) )
        {
        if( theType == 'responsivity.light' )
            {
            log.string( TRACE, "Set stimulus to Illuminant E." )
            stimulus = illuminantE( 1, wave )
            }
        else
            {
            log.string( TRACE, "Set stimulus to the Perfect Reflecting Diffuser." )            
            stimulus = neutralMaterial( 1, wave )   # the perfect-reflecting-diffuser
            }
        }
    
    if( is.null(response) )
        {
        if( is.neural )        
            {
            log.string( ERROR, "Since quantity(x)='%s', an explicit response is required.", quantity(x) )
            return(x)
            }

        response    = rep( 1, m )   # all 1s is conventional
        
        names(response) = toupper( specnames(x) )
        
        log.string( TRACE, "Set desired response to all 1s." )              
        }
    
    if( is.null(method) )
        {
        if( is.XYZ )
            method  = "Bradford"
        else
            method  = "scaling"
            
        log.string( TRACE, "Set method to '%s'.", method )        
        }
        
    #   change characters to numbers
    Ma  = method
    
    if( is.character(Ma) )
        {
        if( Ma == "scaling" )
            Ma = diag(m)      
        else if( m == 3 )
            {
            if( Ma == "Bradford" )
                Ma = matrix( c(0.8951,0.2664,-0.1614,  -0.7502,1.7135,0.0367,  0.0389,-0.0685,1.0296), 3, 3, byrow=T )
            else if( Ma == "Von Kries" )
                Ma = matrix( c(0.40024,0.7076,-0.08081,  -0.2263,1.16532,0.0457,  0,0,0.91822), 3, 3, byrow=T )
            else if( Ma == "MCAT02" )
                Ma = matrix( c( 0.7328, 0.4296, -0.1624,  -0.7036, 1.6975, 0.0061, 0.0030, 0.0136, 0.9834 ), 3, 3, byrow=T )
                
            if(  ! is.character(Ma)  &&  is.character(method)  &&  ! is.XYZ  )
                log.string( WARN, "method='%s' is not really appropriate for non-XYZ object x.", method )
            }
            
        if( is.character(Ma) )
            {
            log.string( ERROR, "method='%s' unknown, for M=%d.", Ma, m )
            return(x)
            }                 
        }
        
        
    #   check validity of stimulus
    if( ! is.colorSpec( stimulus ) )
        {
        log.string( ERROR, "stimulus is not a valid colorSpec." )
        return(x)
        }                    
    
    if( ! identical( wave, wavelength(stimulus) ) )
        {
        log.string( ERROR, "x and stimulus do not have the same wavelengths." )
        return(x)
        }           

    ok1 = theType == 'responsivity.light'    && type(stimulus) == 'light'
    ok2 = theType == 'responsivity.material' && type(stimulus) == 'material' 
    
    if( ! ok1  &&  ! ok2 )
        {
        log.string( ERROR, "type(stimulus) = '%s' is invalid for x.", type(stimulus) )
        return(x)
        }     
        
    if( numSpectra( stimulus ) != 1 )
        {
        log.string( ERROR, "numSpectra(stimulus) = %d is invalid; it must be 1.", numSpectra(stimulus) )
        return(x)
        }

    #   check validity of response
    if( length(response) == 1 ) response = rep( response, m )
    
    ok  = is.numeric(response) && length(response)==m  &&  all(0 < response)
    if( ! ok )
        {
        log.string( ERROR, "response is invalid for x. It must be numeric with length %d, and all positive.", m )
        return(x)
        }     
        
    #   check validity of Ma, it must be an MxM matrix
    ok  = is.numeric(Ma)  &&  length(dim(Ma)==2)  &&  all( dim(Ma) == c(m,m) )
    if( ! ok )
        {
        log.string( ERROR, "adaption matrix is invalid for x. It must be a %dx%d matrix.", m, m )
        return(x)
        }     
        
    #   ensure that both are radiometric        
    x           = radiometric( x, warn=TRUE )
    stimulus    = radiometric( stimulus, warn=TRUE )
        
    #   compute response from the original object x 
    response.src    = product( stimulus, x )
    
    if( is.null(response.src) ) return(x)   # some ERROR unforseen
    
    if( any(response.src <= 0) )
        {
        log.object( ERROR, response.src )
        log.string( "Cannot continue, because 1 or more x response values is non-positive." )
        return(x)
        }
        
    if( m == 1 )
        {
        #   special case
        mat = as.numeric( response/response.src )
        out = multiply( x, mat )        
        }
    else
        {
        mat = makeMappingMatrix( Ma, response.src, response ) 
        if( is.null(mat) )  return(x)
        
        out = multiply( x, t(mat) )        
        }
    
    if( ! is.character(method) )    method=as.character(NA)
    

    #   add useful data to the attribute list.  This will be printed in summary().
    attr( out, "calibrate" )  = list( method=method, Ma=Ma, response.before=response.src, response.after=response, gain=mat )
    
    return( out )
    }
    
    
#   Ma     the adaption matrix
#   white  the reference white    
#   returns a matrix that maps white to (1,1,1,...)
unitMappingMatrix  <-  function( Ma, white )
    {
    lms = Ma %*% as.numeric(white)
    
    if( any( lms <= 0 ) )
        {
        log.object( lms )
        log.string( ERROR, "One component of Ma*white is <= 0" )
        return(NULL)
        }
      
    #   this special test is not necessary, above code takes care of it
    #if( length(lms) == 1 )
    #    #   special case, all are scalars
    #    return( Ma / lms )
        
    return( diag( as.double(1/lms) ) %*% Ma )
    }
    
    
#   Ma  the adaption matrix
#   w1  the source white   
#   w2  the destination white     
#   returns a matrix that maps w1 to w2
makeMappingMatrix  <- function( Ma, w1, w2 )
    {
    crm_src     = unitMappingMatrix( Ma, w1 )
    
    crm_dest    = unitMappingMatrix( Ma, w2 )
    
    if( is.null(crm_src)  ||  is.null(crm_dest) )   return(NULL)
    
    #   print( lms_dest )
    
    return( solve(crm_dest) %*% crm_src )
    }
    
    
    
#--------       UseMethod() calls           --------------#                    
        
calibrate <- function(  x,  stimulus=NULL, response=NULL, method=NULL ) 
    {
    UseMethod("calibrate")
    }    
    
    
    
        