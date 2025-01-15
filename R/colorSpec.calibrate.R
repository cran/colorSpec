

#
#   x           a colorSpec object
#               with type  'responsivity.light'  or  'responsivity.material'
#               with M spectra (output channels)
#
#   stimulus    a single spectrum input, of type 'light' or 'material' respectively
#
#   response    a vector with length M.  Must be all positive, or exactly one positive and the remainder NA.
#               When the latter, the method in
#                   ASTM E308-01 Standard Practice for Computing the Colors of Objects by Using the CIE System.
#                               sec. 7.1.2 is used  (not sec 7.1.1) 
#                   CIE Technical Report 15:2004 3rd Edition, section 7.1
#
#   method      'scaling', 'Bradford', 'Von Kries', 'MCAT02', or 'Bianco+Schettini' or an MxM matrix
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
        log_level( ERROR, "type(x) = '%s' is invalid.", theType )
        return(x)
        }

    m       = numSpectra(x)
    wave    = wavelength(x)

    is.neural   = grepl( 'neural$', quantity(x) )

    is.XYZ      = is.neural  &&  m==3  &&  all( pmatch( c('x','y','z'), tolower(specnames(x)), nomatch=0L  ) == 1:3 )


    #   assign stimulus, if necessary
    if( is.null(stimulus) )
        {
        if( theType == 'responsivity.light' )
            {
            log_level( TRACE, "Set stimulus to Illuminant E." )
            stimulus = illuminantE( 1, wave )
            }
        else
            {
            log_level( TRACE, "Set stimulus to the Perfect Reflecting Diffuser." )
            stimulus = neutralMaterial( 1, wave )
            }
        }

    #   check validity of stimulus
    if( ! is.colorSpec( stimulus ) )
        {
        log_level( ERROR, "stimulus is not a valid colorSpec." )
        return(x)
        }

    if( ! identical( wave, wavelength(stimulus) ) )
        {
        log_level( ERROR, "x and stimulus do not have the same wavelengths." )
        return(x)
        }

    ok1 = theType == 'responsivity.light'    && type(stimulus) == 'light'
    ok2 = theType == 'responsivity.material' && type(stimulus) == 'material'

    if( ! ok1  &&  ! ok2 )
        {
        log_level( ERROR, "type(stimulus) = '%s' is invalid for x.", type(stimulus) )
        return(x)
        }

    if( numSpectra( stimulus ) != 1 )
        {
        log_level( ERROR, "numSpectra(stimulus) = %d is invalid; it must be 1.", numSpectra(stimulus) )
        return(x)
        }




    #   assign response, if necessary

    if( is.null(response) )
        {
        if( is.neural )
            {
            log_level( ERROR, "Since quantity(x)='%s', an explicit response is required.", quantity(x) )
            return(x)
            }

        response    = rep( 1, m )   # all 1s is conventional

        names(response) = toupper( specnames(x) )

        log_level( TRACE, "Set desired response to all 1s." )
        }

    if( length(response) == 1 ) response = rep( response, m )


    #   check validity of response

    ok  = is.numeric(response)  &&  length(response)==m
    if( ! ok )
        {
        log_level( ERROR, "response is invalid for x. It must be numeric with length %d.", m )
        return(x)
        }

    idxfinite   = which( is.finite(response) )

    ok = length(idxfinite) %in% c(1,m)  &&  all(0 < response[idxfinite])
    if( ! ok )
        {
        log_level( ERROR, "response is invalid for x. response must have either 1 or %d non-NA components, which are positive.", m )
        return(x)
        }


    #   assign method, if necessary

    if( is.null(method) )
        {
        if( is.XYZ  &&  length(idxfinite)==m )
            method  = "Bradford"
        else
            method  = "scaling"

        log_level( TRACE, "Set method to '%s'.", method )
        }

    Ma  = method

    if( is.character(method) )
        {
        #   convert character name to a 3x3 matrix using global list p.Ma, which is lazy-loaded from sysdata.rda
        full    = names(p.Ma)

        idx     = pmatch( tolower(method), tolower(full) )
        if( is.na(idx) )
            {
            log_level( ERROR, "method='%s' unknown, for M=%d.", method, m )
            return(x)
            }

        method = full[idx]

        if( method == "scaling" )
            Ma = diag(m)
        else if( m == 3 )
            {
            Ma = p.Ma[[ idx ]]

            if( ! is.XYZ  )
                log_level( WARN, "method='%s' is not really appropriate for non-XYZ responder x.", method )
            }
        else
            {
            log_level( ERROR, "method='%s' invalid, for M=%d.", Ma, m )
            return(x)
            }
        }

    #   check validity of Ma
    ok  = is.numeric(Ma)  &&  length(dim(Ma)==2)  &&  all( dim(Ma) == c(m,m) )
    if( ! ok )
        {
        log_level( ERROR, "adaptation matrix is invalid for x. It must be a %dx%d matrix.", m, m )
        return(x)
        }

    #   for validity of Ma, there is some interaction with response
    if( length(idxfinite) < m  &&  ! is.identity(Ma) )
        {
        log_level( ERROR, "adaptation method is invalid for x. It must be 'scaling'." )
        return(x)
        }


    #   force both colorSpec objects to be radiometric
    x           = radiometric( x, warn=TRUE )
    stimulus    = radiometric( stimulus, warn=TRUE )

    #   compute response from the original object x
    response.src    = product( stimulus, x )

    if( is.null(response.src) ) return(x)   # some ERROR unforseen

    if( any(response.src <= 0) )
        {
        log_object( ERROR, response.src )
        log_level( ERROR, "Cannot continue, because 1 or more x response values are non-positive." )
        return(x)
        }

    #   compute gain matrix 'gmat' and out
    if( m == 1 )
        {
        #   special case, gmat is just a scalar
        gmat = as.numeric( response/response.src )
        out = multiply( x, gmat )
        }
    else if( length(idxfinite) == m )
        {
        gmat = makeMappingMatrix( Ma, response.src, response )
        if( is.null(gmat) )  return(x)

        out = multiply( x, t(gmat) )
        }
    else
        {
        #   we know that length(idxfinite) == 1, the special CIE and ASTM calibration
        g   = response[idxfinite] / response.src[idxfinite]

        gmat = g * diag(m)

        out = multiply( x, gmat )
        }

    if( ! is.character(method) )    method = as.character(NA)

    #   add useful data to the attribute list.  This will be printed in summary().
    attr( out, "calibrate" )  = list( method=method, Ma=Ma, response.before=response.src, response.after=response, gain=gmat )

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
        log_object( lms )
        log_level( ERROR, "One component of Ma*white is <= 0" )
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



