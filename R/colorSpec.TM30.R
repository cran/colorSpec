



#   These colors are in Annex B of the IES TM-30 spec.  page 22

p.bar_color =  c( 163,  92,  96,
                  204, 118,  94,
                  204, 129,  69,
                  216, 172,  98,
                  172, 153,  89,
                  145, 158,  93,
                  102, 139,  94,
                   97, 178, 144,
                  123, 186, 166,
                   41, 122, 126,
                   85, 120, 141,
                  112, 138, 178,
                  152, 140, 170,
                  115,  88, 119,
                  143, 102, 130,
                  186, 122, 142 )

p.bar_color = grDevices::rgb( matrix( p.bar_color, nrow=16, ncol=3, byrow=TRUE ), max=255 )


p.arrow_color =  c( 230,  40,  40,
                    231,  75,  75,
                    251, 129,  46,
                    255, 181,  41,
                    203, 202,  70,
                    126, 185,  76,
                     65, 192, 109,
                      0, 156, 124,
                     22, 188, 176,
                      0, 164, 191,
                      0, 133, 195,
                     59,  98, 170,
                     69, 104, 174,
                    106,  78, 133,
                    157, 105, 161,
                    167,  79, 129 )

p.arrow_color = grDevices::rgb( matrix( p.arrow_color, nrow=16, ncol=3, byrow=TRUE ), max=255 )



#   computeTM30.colorSpec() calcuations depend on the colorSpec object:
#       CESforTM30
#   which is stored in sysdata.rda and lazily loaded at startup


#   x           colorSpec object of type 'light', with 1 test spectra
#   reference   colorSpec object of type 'light' with 1 reference spectra.
#               if NULL, then the reference spectra is computed from the CCT of x
#   digits      for rounding the output
#   isotherms   isotherms for the CCT
#   locus       locus for the CCT
#
#   returns     an M-vector with SSI for the corresponding pair of test and reference spectra

computeTM30.colorSpec <-function( x, reference=NULL )
    {
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )
        {
        log_level( WARN, "Required package 'spacesXYZ' could not be loaded; returning NULL."  )
        return(NULL)
        }

    if( type(x) != 'light' )
        {
        log_level( ERROR, "type(x)='%s', but it must be 'light'.", type(x) )
        return(NULL)
        }

    if( numSpectra(x) != 1 )
        {
        log_level( ERROR, "Argument 'x' has %d spectra, but must have exactly 1", numSpectra(x) )
        return(NULL)
        }

    wave    = wavelength( CESforTM30 )      # 380:780       # later get this from CES object


    #   make sure x is radiometric
    x   = radiometric(x)

    #   resample to the desired wavelengths
    x   = resample( x, wave, extrapolation=0 )

    #   now scale so Y=100
    XYZ = product( x, colorSpec::xyz1964.1nm, wavelength=wave )
    x   = multiply( x, 100/XYZ[2] )

    #   the next function depends on package spacesXYZ
    CCT = computeCCT( x, isotherms='native', locus='precision', strict=FALSE )

    if( is.na(CCT) )    return(NULL)


    if( is.null(reference) )
        {
        #   make planck spectrum and scale so Y=100
        planck      = planckSpectra( CCT, wavelength=wave, normalize=FALSE )
        XYZ         = product( planck, colorSpec::xyz1964.1nm, wavelength=wave )
        planck      = multiply( planck, 100/XYZ[2] )
        specnames(planck)   = sprintf( "Planck reference CCT=%g", CCT )

        #   make daylight spectrum and scale so Y=100
        daylight    = daylightSpectra( CCT, wavelength=wave )
        XYZ         = product( daylight, colorSpec::xyz1964.1nm, wavelength=wave )
        daylight    = multiply( daylight, 100/XYZ[2] )
        specnames(daylight)   = sprintf( "daylight reference CCT=%g", CCT )

        if( CCT <= 4000 )
            reference   = planck
        else if( 5000 <= CCT )
            reference   = daylight
        else
            {
            #   compute weighted average of planck and daylight
            rho = (5000 - CCT) / (5000 - 4000)

            mat = rho * as.matrix(planck)  +  (1 - rho) * as.matrix(daylight)

            reference   = colorSpec( mat, wavelength=wave, quantity=quantity(planck), specnames="blended reference" )
            }
        }
    else
        {
        ok  = is.colorSpec(reference)  &&  type(reference)=='light'
        if( ! ok )
            {
            log_level( ERROR, "Argument 'reference' is not a colorSpec object, or else type(reference)!='light'." )
            return(NULL)
            }

        if( numSpectra(reference) != 1 )
            {
            log_level( ERROR, "numSpectra(reference)==%d, but it must be 1.", numSpectra(reference) )
            return(NULL)
            }

        #   make sure reference is radiometric
        reference   = radiometric( reference )

        reference   = resample( reference, wave,  extrapolation=0 )

        #   now scale so Y=100
        XYZ         = product( reference, colorSpec::xyz1964.1nm, wavelength=wave )
        reference    = multiply( reference, 100/XYZ[2] )
        }


    out = list()
    
    path    = metadata( x )$path
    if( TRUE )      # is.null(path) )
        name    = specnames(x)      # x has only one spectrum
    else
        name    = basename(path)
    
    out$name    = name

    out$CCT = CCT
    out$Duv = attr(CCT,'Duv')

    out$spectra_test_and_ref    = bind( x, reference )      # this works, because wavelength and quantity are the same

    #   for (u',v') use xyz1931.1nm
    XYZ = product( out$spectra_test_and_ref, colorSpec::xyz1931.1nm, wavelength=wave )

    out$uv_test_and_ref = spacesXYZ::uvfromXYZ( XYZ, space=1976 )


    #   for the CIECAM02 Jab calculations, use xyz1964.1nm

    #   test data
    XYZ_test        = product( x, CESforTM30, colorSpec::xyz1964.1nm, wavelength=wave )

    XYZ_test_illum  = product( x, colorSpec::xyz1964.1nm, wavelength=wave )            #;  cat( "XYZ_test_illum = ", XYZ_test_illum, '\n' )

    CAM_test        = spacesXYZ::CIECAM02fromXYZ( XYZ_test, XYZ_test_illum, discount=TRUE )

    if( is.null(CAM_test) ) return(NULL)



    #   reference data
    XYZ_ref         = product( reference, CESforTM30, colorSpec::xyz1964.1nm, wavelength=wave )

    XYZ_ref_illum   = product( reference, colorSpec::xyz1964.1nm, wavelength=wave )    #;  cat( "XYZ_ref_illum = ", XYZ_ref_illum, '\n' )

    CAM_ref         = spacesXYZ::CIECAM02fromXYZ( XYZ_ref, XYZ_ref_illum, discount=TRUE )

    if( is.null(CAM_ref) ) return(NULL)



    DeltaE  = spacesXYZ::DeltaE_CAM02( CAM_test, CAM_ref )

    #   compute Sample Color Fidelity (SCF)
    SCF = color_fidelity( DeltaE )        # 10 * log( exp( (100 - 6.73*DeltaE) /10 )  +  1 )

    #   make data.frame for the samples
    sample_data = data.frame( row.names = 1:length(SCF) )

    sample_data$SCF         = SCF
    sample_data$Jab_test    = CAM_test[ , c("Jp","abp") ]
    sample_data$Jab_ref     = CAM_ref[ , c("Jp","abp") ]
    sample_data$DeltaE      = DeltaE
    sample_data$hue_bin     = as.integer( floor(CAM_ref$h/22.5) + 1 )

    out$sample_data = sample_data

    #   Global Color Fidelity Index R_f
    out$R_f = color_fidelity( mean( DeltaE ) )      # 10 * log( exp( (100 - 6.73*mean(DeltaE)) / 10 )  +  1 )


    sample_count    = rep( NA_integer_, 16 )
    abp_test        = matrix( NA_real_, 16, 2 )
    abp_ref         = matrix( NA_real_, 16, 2 )
    xy_ref          = matrix( NA_real_, 16, 2 )
    lcf             = rep( NA_integer_, 16 )        # local_color_fidelity


    colnames(abp_test)  = c('ap','bp')
    colnames(abp_ref)  = c('ap','bp')


    hue_bin = sample_data$hue_bin

    for( i in 1:16 )
        {
        #   find the row indexes for this bin
        idx = which( hue_bin==i )

        sample_count[i] = length(idx)

        if( sample_count[i] == 0 ) next     # should never really happen

        mat_sub = sample_data$Jab_test[ idx, 'abp' ]

        abp_test[i, ]   = colMeans(mat_sub)

        mat_sub = sample_data$Jab_ref[ idx, 'abp' ]

        abp_ref[i, ]    = colMeans(mat_sub)

        hj      = atan2( mat_sub[ ,2], mat_sub[ ,1] )   # note the order:  2  then  1
        hj_mean = mean(hj)

        xy_ref[i, ] = c( cos(hj_mean), sin(hj_mean) )

        lcf[i]  = color_fidelity( mean( DeltaE[idx] ) )     # 10 * log( exp( (100 - 6.73*mean( DeltaE[idx] )) / 10 )  +  1 )
        }

    mask_empty  = (sample_count==0)

    if( any( mask_empty ) )
        {
        log_level( WARN, "%d of the %d hue bins are empty. The TM-30 measures for '%s' should be used with caution.",
                            sum(mask_empty), length(sample_count), specnames(x) )
        idx = which(mask_empty)

        log_level( WARN, paste( "    The empty bins are: ", paste(idx,collapse=', ') ) )
        }




    #   compute center of each of the 16 hue arcs
    theta               = seq( 22.5/2, by=22.5, len=16 ) * pi/180
    xy_center           = cbind( cos(theta), sin(theta) )
    colnames(xy_center) = c('x','y')


    bin_data = data.frame( row.names = 1:16 )

    abp_diff        = (abp_test - abp_ref) / sqrt( rowSums(abp_ref^2) )     # (60) and (61)

    bin_data$chroma_shift   = rowSums( abp_diff  *  xy_center )                                         # (62)
    bin_data$hue_shift      = rowSums( cbind( -abp_diff[ ,1], abp_diff[ ,2] )  *  xy_center[ ,2:1] )    # (63)
    bin_data$color_fidelity = lcf

    bin_data$sample_count   = sample_count
    bin_data$abp_test       = abp_test
    bin_data$abp_ref        = abp_ref
    bin_data$abp_diff       = abp_diff
    bin_data$xy_center      = xy_center

    # bin_data$xy_ref         = xy_ref

    out$bin_data = bin_data



    #   Gamut Index R_g = ratio of the area of 2 polygons
    poly_test   = abp_test[ ! mask_empty, ]
    poly_ref    = abp_ref[ ! mask_empty, ]

    out$R_g = 100 * areaOfPolygon( poly_test ) / areaOfPolygon( poly_ref )

    class( out )    = c( "TM30", class(out) )

    return(out)
    }


#   returns a vector the same length as delta
#   0 maps to 100.0005..
color_fidelity  <- function( delta )
    {
    10 * log( exp( (100 - 6.73*delta) /10 )  +  1 )
    }
    
    
#   data    as returned from computeTM30()
#
#   see Table E-2 Recommended Specifications Criteria

computePVF  <- function( data )
    {
    out = rep('-',3)
    
    names(out)  = c('P','V','F')
    
    R_f     = data$R_f
    R_g     = data$R_g
    
    R_cs1   = data$bin_data$chroma_shift[1]     # chroma shift for hue #1
    R_f1    = data$bin_data$color_fidelity[1]   # color_fidelity for hue #1
    
    #   Preference, in ascending order
    if( 78<=R_f  &&  95<=R_g  &&  -0.01<=R_cs1  &&  R_cs1<=0.15 )
        out['P']    = '1'
    else if( 75<=R_f  &&  92<=R_g  &&  -0.07<=R_cs1  &&  R_cs1<=0.19 )
        out['P']    = '2'
    else if( 70<=R_f  &&  89<=R_g  &&  -0.12<=R_cs1  &&  R_cs1<=0.23 )
        out['P']    = '3'        
        
    #   Vividness, in ascending order
    if( 118<=R_g  &&  0.15<=R_cs1 )
        out['V']    = '1'
    else if( 110<=R_g  &&  0.06<=R_cs1 )
        out['V']    = '2'
    else if( 100<=R_g  &&  0<=R_cs1 )
        out['V']    = '3'
        
    #   Fidelity, in ascending order
    if( 95 <= R_f )
        out['F']    = '1'
    else if( 90<=R_f  &&  90<=R_f1 )
        out['F']    = '2'
    else if( 85<=R_f  &&  85<=R_f1 )
        out['F']    = '3'
        
    return( out )
    }
    
    

#   x      TM30 object

print.TM30  <-  function( x, ... )
    {    
    cat( 'name:             "', x$name, '"\n', sep='' )

    cat( "CCT:              ", round(x$CCT), ' K\n', sep='' )
    cat( "Duv:              ", round(x$Duv,4), '\n', sep='' )

    cat( "R_f:              ", round(x$R_f), '\n', sep='' )
    cat( "R_g:              ", round(x$R_g), '\n', sep='' )

    cat( "R_cs,h1:          ", round(100*x$bin_data$chroma_shift[1]), '%\n', sep='' )
    cat( "R_hs,h1:          ", round(x$bin_data$hue_shift[1],2), '\n', sep='' )    
    cat( "R_f,h1:           ", round(x$bin_data$color_fidelity[1]), '\n', sep='' )

    PVF     = computePVF( x )
    label   = paste( names(PVF), PVF, sep='', collapse=' ' )
    cat( "Design Intent:    ", label, '\n', sep='' )

    return( invisible(TRUE) )
    }


plot.TM30   <- function( x, omi=c(0.75,0.2,1.25,0.25), ... )
    {
    figs    = c( 0, 0.5, 0.54,    1,
                 0, 0.5, 0.31, 0.54,
                0.5,  1, 0.77,    1,
                0.5,  1, 0.54, 0.77,
                0.5,  1, 0.31, 0.54,
                  0,  1,    0, 0.31 )

    figs    = matrix( figs, nrow=6, ncol=4, byrow=TRUE )

    graphics::split.screen( figs )

    omisaved    = par( omi=omi )
    #   cat( "omi = ", par('omi'), '\n' )
    
    title( main="ANSI/IES TM-30-20 Color Rendition Report", line=2.8, cex.main=1.4, outer=TRUE )
    
    
    mess    = sprintf( "Test Source: %s    Report Time: %s", x$name, format(Sys.time(), usetz=T) )
    title( main=mess, line=1.2, cex.main=0.9, outer=TRUE )
    
    info    = library( help='colorSpec' )      
    info    = format( info )
    mask    = grepl( "^(Version|Built):", info )     
    info    = gsub( "[ ]+", ' ', info[mask] )
    info    = paste( info, collapse='.  ' )
    mess    = sprintf( "report created by R package colorSpec. %s", info )
    title( sub=mess, line=1, cex.sub=0.8, outer=TRUE )
    
    graphics::screen(1)
    plot_CVG(x)

    graphics::screen(2)
    plot_spectra(x)
    #par( mar=c(3,3,0,0) )
    #spectra_test_and_ref = x$spectra_test_and_ref
    #specnames(spectra_test_and_ref) = c("Test","Reference")
    #plot( spectra_test_and_ref, color=c('red','black'), main=FALSE, legend="topright" )

    graphics::screen(3)
    plot_chroma_shift(x)

    graphics::screen(4)
    plot_hue_shift(x)

    graphics::screen(5)
    plot_color_fidelity(x)

    graphics::screen(6)
    plot_SCF(x)

    graphics::close.screen( all=TRUE )
    
    par( omi=omisaved )
    
    return( invisible(TRUE) )
    }


#
#   plot_CVG()
#
#   depends on nativeRaster object:
#       backgroundTM30
#   stored in sysdata.rda

plot_CVG    <- function( data )
    {
    if( FALSE )
        {
        if( ! requireNamespace( 'png', quietly=TRUE ) )
            {
            log_level( WARN, "Required package 'png' could not be loaded; returning FALSE."  )
            return(FALSE)
            }        

        path    = system.file( 'extdata/targets/backgroundTM-30.png', package='colorSpec' )

        if( nchar(path) == 0 )
            {
            log_level( WARN, "Cannot locate required file '%s'.", "backgroundTM-30.png" )
            return(FALSE)
            }

        #   TODO:  must put a try() around this one
        img     = png::readPNG( path, native=TRUE, info=FALSE )
        }


    par( mar=c(0,0,0,0) )

    xlim    = c(-1.5,1.5)
    ylim    = c(-1.5,1.5)

    plot( xlim, ylim, type='n', xlab='', ylab='', axes=F, asp=1 )

    graphics::rasterImage( backgroundTM30, xlim[1], ylim[1], xlim[2], ylim[2], interpolate=TRUE )


    #   the bin boundaries
    theta   = (pi/180) * 22.5 * (0:15)
    unit    = cbind( cos(theta), sin(theta) )
    start   = 0.1 * unit
    end     = 1.5 * unit
    segments( start[ ,1], start[ ,2], end[ ,1], end[ ,2], col="gray50", lty=2, lwd=1 )
    
    #   the little '+' at center
    rad = 0.03
    segments( c(-rad,0), c(0,-rad), c(rad,0), c(0,rad), col="gray50", lwd=1 )

    #   the bin labels
    xy_center   = data$bin_data$xy_center    # arc midpoints
    rad     = 1.4
    text( rad*xy_center[ ,1], rad*xy_center[ ,2], 1:16, col="gray50", adj=c(0.5,0.5) )


    #   the nested circles
    x   = numeric(5)
    graphics::symbols( x, x, circles=c(0.8,0.9,1,1.1,1.2) , fg=c('white','white','black','white','white'), inches=FALSE, add=TRUE )


    #   the test curve, solid red
    abp_test    = xy_center + data$bin_data$abp_diff
    abp         = rbind( abp_test, abp_test[1, ] )
    avec = stats::spline( 1:17, abp[ ,1], n=4*(17-1)+1, method='periodic' )$y
    bvec = stats::spline( 1:17, abp[ ,2], n=4*(17-1)+1, method='periodic' )$y
    lines( avec, bvec, col='red', lwd=3 )

    #   the arrows, from xy_center to abp_test
    graphics::arrows( xy_center[ ,1], xy_center[ ,2], abp_test[ ,1], abp_test[ ,2], col=p.arrow_color, length=0.07, lwd=1.5 )
    
    #   inner and outer circle labels, in bin 12
    rad = c( 0.75, 1.25 )
    text( rad*xy_center[12,1],  rad*xy_center[12,2], c("-20%", "+20%"), col='white', cex=0.6, adj=c(0.5,0.5) )

    #   R_f in left top 
    text( xlim[1], ylim[2], round(data$R_f), col='black', cex=1.25, adj=c(-0.25,1.5) )
    text( xlim[1], ylim[2], expression( italic(R)[f] ), col='black', cex=1.1, adj=c(-0.25,2.7) )
   
    #   R_g in right top 
    text( xlim[2], ylim[2], round(data$R_g), col='black', cex=1.25, adj=c(1.25,1.5) )
    text( xlim[2], ylim[2], expression( italic(R)[g] ), col='black', cex=1.1, adj=c(1.25,2.5) )
    
    #   local fidelity for hue #1 in left bottom
    R_fh1   = data$bin_data$color_fidelity[1]
    text( xlim[1], ylim[1], round(R_fh1), col='black', cex=1.25, adj=c(-1,-1) )
    text( xlim[1], ylim[1], expression( italic(R)['f,h1'] ), col='black', cex=1.1, adj=c(-0.25,-1.9) )
   
    #   local chroma shift for hue #1 in right bottom
    R_csh1  = data$bin_data$chroma_shift[1]
    label   = sprintf( "%g%%", round(100*R_csh1) )
    text( xlim[2], ylim[1], label, col='black', cex=1.25, adj=c(1.25,-1) )
    text( xlim[2], ylim[1], expression( italic(R)['cs,h1'] ), col='black', cex=1.1, adj=c(1.25,-1.9) )
   

   
    #   Duv above left-top corner
    text( xlim[1], ylim[2], expression( italic(D)[uv] ~ ':'), col='black', cex=0.7, adj=c(-0.25,-0.2), xpd=NA )
    label = sprintf( "%.4f", data$Duv )
    text( xlim[1]+0.3, ylim[2], label, col='black', cex=0.7, adj=c(0,-0.5), xpd=NA )
    
    #   CCT above left-top corner
    text( xlim[1], ylim[2], expression( CCT~':'), col='black', cex=0.7, adj=c(-0.25,-2.2), xpd=NA )
    label = sprintf( "%g K", round(data$CCT) )
    text( xlim[1]+0.31, ylim[2], label, col='black', cex=0.7, adj=c(0,-2.4), xpd=NA )
    
    #   PVF above right-top corner
    PVF = computePVF( data )
    label   = paste( names(PVF), PVF, sep='', collapse=' ' )
    text( xlim[2], ylim[2], label, col='black', cex=1, adj=c(1,-0.5), xpd=NA )

    return(TRUE)
    }
    
    
plot_spectra <- function( data )    
    {
    par( mar=c(1.5,2,0,0) )

    wave    = wavelength( data$spectra_test_and_ref )
    xlim    = range(wave)
    
    mat         = as.matrix( data$spectra_test_and_ref )
    spec_test   = mat[ ,1]
    spec_ref    = mat[ ,2]
    ylim        = c( 0, 1.15*max(mat) )  # add a little for the legend
    
    plot( xlim, ylim, type='n', xlab='', ylab='', axes=F )
    
    xvec    = seq( xlim[1], xlim[2], by=50 )
    # axis( side=1, at=xvec, col="gray28" )
    
    graphics::clip( xlim[1], xlim[2], ylim[1], ylim[2] ) 
    
    abline( v=xvec, col='gray93' )
    grid( nx=NA, ny=NULL, lty=1, col='gray93' )
    
    lines( wave, spec_test, col='red', lwd=2 )
    lines( wave, spec_ref, col='black', lwd=2 )

    # draw custom border
    rect( xlim[1], ylim[1], xlim[2], ylim[2], xpd=T, lwd=2, border='black' )
    
    #   x values
    text( xvec, 0, as.character(xvec), col='gray28', adj=c(0.5,2.1), cex=0.60, xpd=NA )
    
    #   tick marks
    eps = 0.02 * (ylim[2] - ylim[1])
    segments( xvec, 0, xvec, -eps, col='black', xpd=NA )
    
    #   axis labels
    title( xlab="Wavelength (nm)", line=0.4, xpd=NA, col.lab="gray28", cex.lab=0.8 )
    title( ylab="Radiant Power", line=0, xpd=NA, col.lab="gray28", cex.lab=0.8 )
    
    #   legend
    width   = xlim[2] - xlim[1]
    x0  = xlim[1] + 0.08*width
    x1  = x0 + 0.08*width
    
    height  = ylim[2] - ylim[1]
    y0  = ylim[2] - 0.06*height
    
    lines( c(x0,x1), c(y0,y0), col='red', lwd=2 )
    text( x1, y0, "Test", col='gray28', adj=c(-0.1,0.5), cex=0.75 )
    
    x0  = x0 + 0.2*width
    x1  = x0 + 0.08*width
    
    lines( c(x0,x1), c(y0,y0), col='black', lwd=2 )
    text( x1, y0, "Reference", col='gray28', adj=c(-0.1,0.5), cex=0.75 )
    
    
    return(TRUE)
    }


plot_color_fidelity <- function( data )
    {
    par( mar=c(1.5,2.5,0,0) )

    height  = data$bin_data$color_fidelity

    ylim    = c(0,100)

    mp = graphics::barplot( height, ylim=ylim, plot=FALSE )

    xlim = c( mp[1] - 0.7, mp[ length(mp) ] + 0.7 )

    #   set up the figure limits
    plot( xlim, ylim, type='n', xlab='', ylab='', axes=F )

    #   y values
    yvec    = 20 * (0:5)
    text( xlim[1], yvec, as.character(yvec), adj=c(1.2,0.5), xpd=NA, cex=0.75, col="gray28" )

    graphics::clip( xlim[1], xlim[2], ylim[1], ylim[2] )    # necessary for abline to clip to correctly
    abline( h=yvec, col="gray93", xpd=FALSE )

    #   segments( xlim[1], yvec, xlim[2], yvec, col="gray93" )   # clips correctly without the call to clip()

    #   put the bars *over* the abline segments
    graphics::barplot( height, ylim=ylim, col=p.bar_color, border=NA, las=1, tck=0, axes=F, add=TRUE )

    # draw custom border
    rect( xlim[1], ylim[1], xlim[2], ylim[2], xpd=T, lwd=2, border='black' )

    #   x values
    text( mp, 0, as.character( 1:length(mp) ), col='gray28', adj=c(0.5,1.4), cex=0.60, xpd=NA )
    text( mp, height, round(height), col='gray28', adj=c(0.5,-0.2), cex=0.75 )

    #   axis labels
    title( xlab="Hue-Angle Bin (j)", line=0, xpd=NA, col.lab="gray28", cex.lab=0.8 )
    title( ylab="Local Color Fidelity", line=0.5, xpd=NA, col.lab="gray28", cex.lab=0.8 )

    return(TRUE)
    }


plot_hue_shift <- function( data )
    {
    par( mar=c(1.5,2.5,0,0) )

    height  = data$bin_data$hue_shift

    ylim    = c(-0.4,0.4)

    mp = graphics::barplot( height, ylim=ylim, plot=FALSE )

    xlim = c( mp[1] - 0.7, mp[ length(mp) ] + 0.7 )

    #   set up the figure limits
    plot( xlim, ylim, type='n', xlab='', ylab='', axes=F )

    #   y values
    yvec    = (-4:4) / 10

    text( xlim[1], yvec, as.character(yvec), adj=c(1.2,0.5), xpd=NA, cex=0.75, col="gray28" )

    segments( xlim[1], yvec, xlim[2], yvec, col="gray93" )  # clips correctly

    # draw y=0 in the middle
    lines( xlim, c(0,0), col="gray28" )

    #   put the bars *over* the segments

    graphics::clip( xlim[1], xlim[2], ylim[1], ylim[2] )    # necessary for the bars

    if( TRUE )
        {
        graphics::barplot( height, ylim=ylim, xpd=FALSE, col=p.bar_color, border=NA, las=1, tck=0, axes=F, add=TRUE )
        }
    else
        {
        left    = mp - 0.5
        right   = mp + 0.5
        bottom  = rep(0,length(height))

        rect( left, bottom, right, height, border=NA, col=p.bar_color, xpd=FALSE )
        }




    # draw custom border
    rect( xlim[1], ylim[1], xlim[2], ylim[2], xpd=TRUE, lwd=2, border='black' )

    #   x values
    text( mp, ylim[1], as.character( 1:length(mp) ), col='gray28', adj=c(0.5,1.4), cex=0.60, xpd=NA )

    graphics::clip( xlim[1], xlim[2], ylim[1], ylim[2] )    # necessary for the bar labels

    #   label the bars with positive height
    mask    = (0 <= height )
    if( any(mask) )
        {
        y   = height[mask]
        text( mp[mask], y, sprintf("%.2f",y), col='gray28', adj=c(0.5,-0.2), cex=0.4, xpd=FALSE )
        }

    #   label the bars with negative height
    mask    = ! mask
    if( any(mask) )
        {
        y   = height[mask]
        text( mp[mask], y, sprintf("%.2f",y), col='gray28', adj=c(0.5,1.2), cex=0.4, xpd=FALSE )
        }

    #   axis labels
    title( xlab="Hue-Angle Bin (j)", line=0,  col.lab="gray28", cex.lab=0.8 )
    title( ylab="Local Hue Shift", line=0.8,  col.lab="gray28", cex.lab=0.8 )

    return(TRUE)
    }


plot_chroma_shift <- function( data )
    {
    par( mar=c(1.5,2.5,0,0) )

    height  = data$bin_data$chroma_shift

    ylim    = c(-0.4,0.4)

    mp = graphics::barplot( height, ylim=ylim, plot=FALSE )

    xlim = c( mp[1] - 0.7, mp[ length(mp) ] + 0.7 )

    #   set up the figure limits
    plot( xlim, ylim, type='n', xlab='', ylab='', axes=F )

    #   y values
    yvec    = (-4:4) / 10

    text( xlim[1], yvec, sprintf("%g%%",100*yvec), adj=c(1.2,0.5), xpd=NA, cex=0.6, col="gray28" )

    segments( xlim[1], yvec, xlim[2], yvec, col="gray93" )  # clips correctly

    # draw y=0 in the middle
    lines( xlim, c(0,0), col="gray28" )

    #   put the bars *over* the segments

    graphics::clip( xlim[1], xlim[2], ylim[1], ylim[2] )    # necessary for the bars

    if( TRUE )
        {
        graphics::barplot( height, ylim=ylim, xpd=FALSE, col=p.bar_color, border=NA, las=1, tck=0, axes=F, add=TRUE )
        }
    else
        {
        left    = mp - 0.5
        right   = mp + 0.5
        bottom  = rep(0,length(height))

        rect( left, bottom, right, height, border=NA, col=p.bar_color, xpd=FALSE )
        }


    # draw custom border
    rect( xlim[1], ylim[1], xlim[2], ylim[2], xpd=TRUE, lwd=2, border='black' )

    #   x values
    text( mp, ylim[1], as.character( 1:length(mp) ), col='gray28', adj=c(0.5,1.4), cex=0.60, xpd=NA )

    graphics::clip( xlim[1], xlim[2], ylim[1], ylim[2] )    # necessary for the bar labels

    #   label the bars with positive height
    mask    = (0 <= height )
    if( any(mask) )
        {
        y   = height[mask]
        text( mp[mask], y, sprintf("%g%%",round(100*y)), col='gray28', adj=c(0.5,-0.2), cex=0.4, xpd=FALSE )
        }

    #   label the bars with negative height
    mask    = ! mask
    if( any(mask) )
        {
        y   = height[mask]
        text( mp[mask], y, sprintf("%g%%",round(100*y)), col='gray28', adj=c(0.5,1.2), cex=0.4, xpd=FALSE )
        }

    #   axis labels
    title( xlab="Hue-Angle Bin (j)",  line=0,  col.lab="gray28", cex.lab=0.8 )
    title( ylab="Local Chroma Shift", line=1.1,  col.lab="gray28", cex.lab=0.8 )

    return(TRUE)
    }



plot_SCF <- function( data )
    {
    if( ! requireNamespace( 'spacesRGB', quietly=TRUE ) )
        {
        log_level( WARN, "Required package 'spacesRGB' could not be loaded; returning FALSE."  )
        return(FALSE)
        }

    par( mar=c(1.5,2.5,0,0) )

    #   first compute all 99 colors, for the spectra in colorSpec::CESforTM30
    wave    = wavelength( CESforTM30 )      # 380:780

    #   the test spectrum is spectrum #2
    spectrum_ref    = subset( data$spectra_test_and_ref, 2 )

    #   in the next line CESforTM30 is private data, stored in sysdata.rda
    #   and D50.5nm and BT.709.RGB is public data, stored in colorSpec.rda
    RGB = product( colorSpec::D50.5nm, CESforTM30, colorSpec::BT.709.RGB, wavelength=wave )   #;  print( str(RGB) )

    RGB = spacesRGB::SignalRGBfromLinearRGB( RGB )$RGB

    colvec  = grDevices::rgb( RGB, max=1 )    #; print( str(colvec) )

    #   do not put any space between these bars.
    space   = 0

    height  = data$sample_data$SCF

    ylim    = c(0,100)

    mp = graphics::barplot( height, ylim=ylim, space=space, plot=FALSE )

    xlim = c( mp[1] - 0.7, mp[ length(mp) ] + 0.7 )

    #   set up the figure limits
    plot( xlim, ylim, type='n', xlab='', ylab='', axes=F )

    #   y values
    yvec    = 10 * (0:10)
    text( xlim[1], yvec, as.character(yvec), adj=c(1.2,0.5), xpd=NA, cex=0.75, col="gray28" )

    graphics::clip( xlim[1], xlim[2], ylim[1], ylim[2] )    # necessary for abline to clip to correctly
    abline( h=yvec, col="gray93", xpd=FALSE )

    #   segments( xlim[1], yvec, xlim[2], yvec, col="gray93" )   # clips correctly without the call to clip()

    #   put the bars *over* the abline segments
    graphics::barplot( height, ylim=ylim, space=space, col=colvec, border=NA, las=1, tck=0, axes=F, add=TRUE )

    # draw custom border
    rect( xlim[1], ylim[1], xlim[2], ylim[2], xpd=T, lwd=2, border='black' )

    #   x values, only multiples of 5, but also add 1 and 99
    idx = c( 1L, seq( 5L, length(mp), by=5L ), length(mp) )
    text( mp[idx], ylim[1], as.character(idx), col='gray28', adj=c(0.5,1.5), cex=0.60, xpd=NA )
    # text( mp, height, round(height), col='gray28', adj=c(0.5,-0.2), cex=0.75 )

    #   axis labels
    title( xlab="Color Evaluation Sample", line=0, xpd=NA, col.lab="gray28", cex.lab=0.8 )
    title( ylab="Sample Color Fidelity", line=0, xpd=NA, col.lab="gray28", cex.lab=0.8 )

    return(TRUE)
    }


#--------       UseMethod() calls           --------------#

computeTM30 <- function( x, reference=NULL )
    {
    UseMethod("computeTM30")
    }
