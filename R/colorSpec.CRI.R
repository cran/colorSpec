
    
        
#   x      a colorSpec object with type 'light'   
#
#   returns CRI in interval (-Inf,100], or NA 

#   requires private data frame TCSforCRI, which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()

computeCRI.colorSpec   <- function( x, adapt=TRUE, attach=FALSE, tol=5.4e-3  )
    {
    out = NA_real_
    
    if( numSpectra( x ) != 1 )
        {
        log.string( ERROR, "Object '%s' has %d spectra, but must have exactly 1.",
                    deparse(substitute(x)), numSpectra( x ) )
        return( out )
        }

    names(out)  = specnames(x)
    
    if( type(x) != 'light' )
        {
        log.string( WARN, "The type of of '%s' is '%s', but it must be 'light'.",
                    deparse(substitute(x)), type(x) )
        return( out )
        }
    
    CCT = computeCCT(x)     #;  print(CCT)
    
    if( is.na(CCT) )    return(out) # error message already issued
    
    wave    = wavelength(x)
    
    #   find the reference illuminant
    locus   = ifelse( CCT < 5000, 'planckian', 'daylight' )
    
    #log.string( DEBUG, "For CCT=%g, using %s radiator as reference.", CCT, locus )
    
    if( locus == 'planckian' )
        #   Planckian radiator
        illum.ref   = planckSpectra( CCT, wavelength=wave )
    else
        #   daylight radiator
        illum.ref   = daylightSpectra( CCT, wavelength=wave, roundMs=TRUE )
    
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
    
    dc  = sqrt( sum( (uv.illum.test - uv.illum.ref)^2 ) )      #; print(dc)
    
    if( tol < dc )
        {
        #   print(dc)
        log.string( ERROR, "The distance from uv.test to uv.ref (on the %s locus) = %g > %g.  It is too large.",
                                locus, dc, tol )
        return( out )
        }
        
    #   compute UVW.ref for the first 8 test samples
    #   TCS8    = subset( colorSpec::TCSforCRI, 1:8 )
        
    XYZ.ref = product( illum.ref, TCSforCRI, colorSpec::xyz1931.1nm, wavelength=wave )
    UVW.ref = UVW_from_XYZ( XYZ.ref, XYZ.illum.ref )    
        
    #   now compute UVW.test for the first 8 test samples

    XYZ.test    = product( x, TCSforCRI, colorSpec::xyz1931.1nm, wavelength=wave )
    
    if( adapt )  
        {
        #   compute uv for each test sample
        xyY = xyY_from_XYZ( XYZ.test )
        uv.before   = uv_from_xy( xyY[ ,1:2] )

        #   adapt uv from test to references, change is small
        uv  = chromaticityAdaption.uv( uv.before, uv.illum.test, uv.illum.ref )
        uvY = cbind( uv, XYZ.test[ ,2] )
        
        class(uv.before)    = 'model.matrix'        
        class(uv)           = 'model.matrix'
        table4  = data.frame( before=uv.before, after=uv, difference=uv-uv.before )
        
        uvY.0   =  c( uv.illum.ref, XYZ.illum.test[2] )    # note use of uv ref illum here

        #   from uvY to UVW
        UVW.test = UVW_from_uvY( uvY, uvY.0 )        
        }
    else
        {        
        UVW.test = UVW_from_XYZ( XYZ.test, XYZ.illum.test )         # note use of test illum here
        table4  = NULL
        }
        
    #   now for the DeltaE's
    DeltaE  = sqrt( rowSums( (UVW.test- UVW.ref)^2 ) )  # a vector of length 14
    CRI     = 100 - 4.6*DeltaE                          # a vector of length 14
    
    #print( rbind(DeltaE,CRI) )
    
    #   the output is the mean of the 1st 8 test samples only !!
    out = mean( CRI[1:8] )        
    
    names(out)  = specnames(x)
    
    if( attach )  
        {
        table1  = rbind( XYZ.illum.test, XYZ.illum.ref )
        table2  = rbind( xy.illum.test, xy.illum.ref )        
        table3  = rbind( uv.illum.test, uv.illum.ref )

        table1  = cbind( table1, table2, table3 )
        
        class(XYZ.ref)  = 'model.matrix'
        class(XYZ.test) = 'model.matrix'
                
        class(UVW.ref)  = 'model.matrix'
        class(UVW.test) = 'model.matrix'
        
        attr(out,'data')    = list( CCT=CCT,
                                    table1=table1,
                                    #   table2=table2,                                    
                                    table2=data.frame( referen=XYZ.ref,  test=XYZ.test),
                                    table3=table4,
                                    table4=data.frame( referen=UVW.ref,  test=UVW.test, DeltaE=DeltaE, CRI=CRI ) )
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
        
        
computeCRI <- function( x, adapt=TRUE, attach=FALSE, tol=5.4e-3    )
    {
    UseMethod("computeCRI")
    }
    