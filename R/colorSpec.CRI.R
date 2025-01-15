
#   x      a colorSpec object with a single spectrum and type 'light'.  Already verified.
#
#   if full is FALSE, just returns the CRI in interval (-Inf,100], or NA
#   if full is TRUE, returns a list of intermediate data, and also the CRI
#
#   requires private data frame TCSforCRI, which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
#
#   in case of error, returns NULL

computeCRIsingle   <- function( x, full, CCT=NULL, adapt=TRUE, tol=5.4e-3 )
    {
    if( is.null(CCT) )
        {
        if( ! requireNamespace( "spacesXYZ", quietly=TRUE ) )
            {
            log_level( WARN, "Package 'spacesXYZ' needs to be installed to compute CCT." )
            return( NULL )
            }

        CCT = computeCCT( x )
        }

    ok  = length(CCT)==1  &&  is.numeric(CCT)  &&  is.finite(CCT)  &&  0<CCT

    if( ! ok )
        {
        log_level( WARN, "CCT is invalid.  It must be a single positive number." )
        return(NULL)
        }

    wave    = wavelength(x)

    #   find the reference illuminant
    locus   = ifelse( CCT < 5000, 'planckian', 'daylight' )

    #log_level( DEBUG, "For CCT=%g, using %s radiator as reference.", CCT, locus )

    if( locus == 'planckian' )
        #   Planckian radiator
        illum.ref   = planckSpectra( CCT, wavelength=wave )
    else
        #   daylight radiator
        illum.ref   = daylightSpectra( CCT, wavelength=wave, roundMs=TRUE )

    if( is.null(illum.ref) )    return(NULL)    # error already issued


    XYZ.illum.test  = product( x, colorSpec::xyz1931.1nm, wavelength=wave )

    XYZ.illum.ref   = product( illum.ref, colorSpec::xyz1931.1nm, wavelength=wave )

    #   scale the reference illuminant so the Y's are the same
    illum.ref       = multiply( illum.ref, XYZ.illum.test[2]/XYZ.illum.ref[2] )
    XYZ.illum.ref   = product( illum.ref, colorSpec::xyz1931.1nm, wavelength=wave )

    #   convert both illuminants to uv
    xy.illum.test = xyY_from_XYZ(XYZ.illum.test)[ ,1:2, drop=F]
    uv.illum.test = uv_from_xy( xy.illum.test )     #; cat( "uv (test)", uv.illum.test, '\n' )

    #print( xyY_from_XYZ(XYZ.illum.test)    )
    #print( xyY_from_XYZ(XYZ.illum.test)[ ,1:2, drop=F]  )

    xy.illum.ref  = xyY_from_XYZ(XYZ.illum.ref)[ ,1:2, drop=F]
    uv.illum.ref  = uv_from_xy( xy.illum.ref )      #; cat( "uv (ideal)", uv.illum.ref, '\n' )

    Delta_uv  = sqrt( sum( (uv.illum.test - uv.illum.ref)^2 ) )      #; print(Delta_uv)

    if( tol < Delta_uv )
        {
        log_level( WARN, "For spectrum '%s', the distance from uv.test to uv.ref (on the %s locus) = %g > %g.  It is too large.",
                                specnames(x), locus, Delta_uv, tol )
        return( NULL )
        }

    #   compute UVW.ref for the first 8 test samples
    #   TCS8    = subset( colorSpec::TCSforCRI, 1:8 )

    XYZ.ref = product( illum.ref, TCSforCRI, colorSpec::xyz1931.1nm, wavelength=wave )
    UVW.ref = UVW_from_XYZ( XYZ.ref, XYZ.illum.ref )

    #   now compute UVW.test for the first 8 test samples

    XYZ.test    = product( x, TCSforCRI, colorSpec::xyz1931.1nm, wavelength=wave )

    #   prepend an 'R' to the test sample number
    rnames  = paste0( 'R', 1:nrow(TCSforCRI) )
    
    if( adapt )
        {
        #   compute uv for each test sample
        xyY = xyY_from_XYZ( XYZ.test )
        uv.before   = uv_from_xy( xyY[ ,1:2] )

        #   adapt uv from test to references, change is small
        uv  = chromaticityAdaption.uv( uv.before, uv.illum.test, uv.illum.ref )
        uvY = cbind( uv, XYZ.test[ ,2] )

        #class(uv.before)    = 'model.matrix'
        #class(uv)           = 'model.matrix'
        #table3  = data.frame( before=uv.before, after=uv, difference=uv-uv.before )

        table3              = data.frame( row.names=rnames )
        table3$before       = uv.before
        table3$after        = uv
        table3$difference   = uv - uv.before

        uvY.0   =  c( uv.illum.ref, XYZ.illum.test[2] )    # note use of uv ref illum here

        #   from uvY to UVW
        UVW.test = UVW_from_uvY( uvY, uvY.0 )
        }
    else
        {
        UVW.test = UVW_from_XYZ( XYZ.test, XYZ.illum.test )         # note use of test illum here
        table3  = NULL
        }

    #   now for the DeltaE's
    DeltaE  = sqrt( rowSums( (UVW.test- UVW.ref)^2 ) )  # a vector of length 14
    CRI     = 100 - 4.6*DeltaE                          # a vector of length 14

    if( ! full )
        {
        #   only a number is requested, so we can exit now
        #   the final CRI is the mean of the 1st 8 test samples only !!
        return( mean( CRI[1:8] ) )
        }

    #print( rbind(DeltaE,CRI) )

    tab1    = rbind( XYZ.illum.test, XYZ.illum.ref )
    tab2    = rbind( xy.illum.test, xy.illum.ref )
    tab3    = rbind( uv.illum.test, uv.illum.ref )

    table1  = cbind( tab1, tab2, tab3 )

    #class(XYZ.ref)  = 'model.matrix'
    #class(XYZ.test) = 'model.matrix'
    table2          = data.frame( row.names=rnames )
    table2$referen  = XYZ.ref
    table2$test     = XYZ.test


    table4          = data.frame( row.names=rnames )
    table4$referen  = UVW.ref
    table4$test     = UVW.test
    table4$DeltaE   = DeltaE
    table4$CRI      = CRI

    #   the final CRI is the mean of the 1st 8 test samples only !!
    out = list( CCT=CCT,
                illum.ref=illum.ref,
                table1=table1,
                Delta_uv=Delta_uv,
                table2=table2,     # data.frame( referen=XYZ.ref,  test=XYZ.test),
                table3=table3,
                table4=table4,
                CRI=mean( CRI[1:8] ) )      # data.frame( referen=UVW.ref,  test=UVW.test, DeltaE=DeltaE, CRI=CRI ) )

    return( out )
    }


#   x      a colorSpec object with type 'light'
#
#   returns a list of intermediate data, followed by CRI in interval (-Inf,100], or NA

#   requires private data frame TCSforCRI, which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
#
#   in case of error, returns NULL

computeCRIdata.colorSpec   <- function( x, CCT=NULL, adapt=TRUE, tol=5.4e-3 )
    {
    if( numSpectra( x ) != 1 )
        {
        log_level( WARN, "Object '%s' has %d spectra, but must have exactly 1.",
                    deparse(substitute(x)), numSpectra( x ) )
        return( NULL )
        }

    if( type(x) != 'light' )
        {
        log_level( WARN, "The type of of '%s' is '%s', but it must be 'light'.",
                    deparse(substitute(x)), type(x) )
        return( NULL )
        }

    return( computeCRIsingle( x, full=TRUE, CCT=CCT, adapt=adapt, tol=tol ) )
    }


#   x      a colorSpec object with type 'light'
#
#   returns a vector of CRIs in interval (-Inf,100], or NA

computeCRI.colorSpec <- function( x, CCT=NULL, adapt=TRUE, tol=5.4e-3, attach=FALSE )
    {
    n   = numSpectra( x )

    out = rep( NA_real_, n )

    if( n == 0 )
        {
        log_level( WARN, "Object '%s' has 0 spectra; it is empty.  Returning 0-length vector.", deparse(substitute(x)) )
        return( out )
        }

    names(out)  = specnames(x)


    if( type(x) != 'light' )
        {
        log_level( WARN, "The type of of '%s' is '%s', but it must be 'light'.",
                                deparse(substitute(x)), type(x) )
        return( out )
        }
        
    if( ! is.null(CCT) )
        {
        if( length(CCT)==1 )    CCT = rep( CCT, n )
        
        ok  = length(CCT)==n  &&  is.numeric(CCT)
        if( ! ok )
            {
            log_level( WARN, "CCT is invalid." )
            return( out )
            }
        }

    if( attach && n==1 )
        {
        #   a special case
        dat = computeCRIsingle( x, full=TRUE, CCT=CCT, adapt=adapt, tol=tol )    # dat is a big list

        if( ! is.null(dat) )
            {
            out[1]  = dat$CRI
            attr(out,'data')    = dat
            }

        return( out )
        }

    for( k in 1:n )
        {
        if( is.null(CCT) )
            CCTk = NULL
        else
            CCTk = CCT[k]
        
        CRI = computeCRIsingle( subset(x,k), full=FALSE, CCT=CCTk, adapt=adapt, tol=tol )     # CRI is just a number

        if( is.null(CRI) )  next

        out[k]  = CRI
        }

    return( out )
    }





#   .uv         Mx2 matrix of test uv's
#   .uv.test    uv of test illuminant
#   .uv.ref     uv of reference illuminant
#
#   value       Mx2 matrix of corrected test uv's
chromaticityAdaption.uv <- function( .uv, .uv.test, .uv.ref )
    {
    cd.test = cd_from_uv( .uv.test )
    cd.ref  = cd_from_uv( .uv.ref )

    cd  = cd_from_uv( .uv )     #; print( cd )

    #   modify the cd's by simple scaling, so that test maps to ref
    cd[ ,1] = cd[ ,1] * cd.ref[1] / cd.test[1]
    cd[ ,2] = cd[ ,2] * cd.ref[2] / cd.test[2]

    # and now back to uv
    #   this matrix is the inverse of the one in cd_from_uv(), up to a constant (and roundoff)
    mat = matrix( c(0.404,-4,10.872,  0,0,5.52,  1.481,-1,16.518), 3, 3 )

    uvw = cbind( cd, 1 ) %*% mat

    uv  = uvw[ , 1:2, drop=F] / uvw[ ,3]    # divides each column appropriately

    colnames(uv)    = c( 'u', 'v' )

    return( uv )
    }


#   .uv         Mx2 matrix of test uv's
cd_from_uv  <- function( .uv )
    {
    #cd  = .uv %*% matrix( c(-1,-10, -1.481, 1.708), 2, 2 )

    #cd  = cd + matrix( c(4, 0.404), nrow(.uv), 2, byrow=T )

    #cd  = cd / matrix( .uv[ ,2], nrow(.uv), 2 )

    #c   = (4 - .uv[ ,1] - 10*.uv[ ,2]) / .uv[ ,2]

    #d   = (0.404 - 1.481*.uv[ ,1] + 1.708*.uv[ ,2]) / .uv[ ,2]

    if( length(.uv) == 2 )  dim(.uv) = c(1,2)

    mat = matrix( c(-1,-10,4,  -1.481,1.708,0.404, 0,1,0), 3,3 )

    cde = cbind(.uv,1) %*% mat

    cd  = cde[ , 1:2, drop=F] / cde[ ,3]    # divides each column appropriately

    return( cd )
    }





#--------       UseMethod() calls           --------------#

computeCRIdata <- function( x, CCT=NULL, adapt=TRUE, tol=5.4e-3 )
    {
    UseMethod("computeCRIdata")
    }

computeCRI <- function( x, CCT=NULL, adapt=TRUE, tol=5.4e-3, attach=FALSE  )
    {
    UseMethod("computeCRI")
    }
