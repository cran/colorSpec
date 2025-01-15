#   Specification   S-2018-001
#   Spectral Similarity Index (SSI)
#   The Academy of Motion Picture Arts and Sciences
#   Science and Technology Council
#   Solid State Lighting (SSL) Project
#
#   The Science and Technology Council wishes to acknowledge the following key contributors to the creation
#   of the SSI metric and the drafting and review of this document and the procedure and calculations described.
#   Jack Holm       Tom Maier       Paul Debevec    Chloe LeGendre
#   Joshua Pines    Jonathan Erland George Joblove  Scott Dyer
#   Blake Sloan     Joe di Gennaro  Dan Sherlock    Alex Forsythe



g.Trap30x301 = NULL     # matrix for trapezoidal integration.  Created one time in a session, on demand.




#   x           colorSpec object of type 'light', with M test spectra
#   reference   colorSpec object of type 'light' with 1 or M reference spectra.
#               if NULL, then M reference spectra are computed from the CCT of x
#   digits      for rounding the output
#   isotherms   isotherms for the CCT
#   locus       locus for the CCT
#
#   returns     an M-vector with SSI for the corresponding pair of test and reference spectra

computeSSI.colorSpec <-function( x, reference=NULL, digits=0, isotherms='mccamy', locus='robertson' )
    {
    if( type(x) != 'light' )
        {
        log_level( ERROR, "type(x)='%s', but it must be 'light'.", type(x) )
        return(NULL)
        }
        
    if( is.null(g.Trap30x301) )
        {
        base::unlockBinding( "g.Trap30x301", asNamespace('colorSpec') )          
        g.Trap30x301 <<- makeTrap30x301()     # only called one time, on demand        
        }
        
        
    lambda1     = 375:675        
        
    m   = numSpectra(x)
    if( m == 0 )    return( numeric(0) )
    
    numref = 0    
    
    if( is.null(reference) )
        {
        #   from Specification   S-2018-001.   These are the same values as in Wyszecki & Stiles  (1982)
        c   = 2.99792458e8  # m/s, the speed of light
        h   = 6.626176e-34  # J*s, the Planck constant
        k   = 1.380662e-23  # J/K, the Boltzmann constant
                
        c2  = h*c/k         # m*K,  we may need c2 later in planckSpectra()
        }
    else
        {
        ok  = is.colorSpec(reference)  &&  type(reference)=='light'
        if( ! ok )
            {
            log_level( ERROR, "Argument 'reference' is not a colorSpec object, or else type(reference)!='light'." )
            return(NULL)
            }
            
        numref  = numSpectra(reference)
        ok      = numref %in% c(1,m)
        if( ! ok )
            {
            log_level( ERROR, "numSpectra(reference)==%d, but it must be %d or %d.", numref, 1, m )
            return(NULL)
            }
            
        #   make sure reference is radiometric
        reference   = radiometric(reference)
        
        #   Section 4.1.1 Interpolation to 1 nm Increments
        reference   = resample( reference, lambda1, method='linear', extrapolation=0 )
            
        #   4.2 Integrate Spectra in 10 nm Intervals            
        ref30xNUMREF    = g.Trap30x301 %*% as.matrix(reference)       
        }
        
    #   make sure x is radiometric
    x   = radiometric(x)
    
    #   Section 4.1.1 Interpolation to 1 nm Increments
    x   = resample( x, lambda1, method='linear', extrapolation=0 )
        
    #   4.2 Integrate Spectra in 10 nm Intervals
    test30xM    = g.Trap30x301 %*% as.matrix(x)       
        
    if( is.null(reference) )
        {
        #   compute CCT for every spectrum in x, which we will need later
        CCTvec  = computeCCT( x, isotherms=isotherms, locus=locus )
        if( is.null(CCTvec) )   return(NULL)
        
        CCTvec  = round( CCTvec )
        }
        
    #   ready to iterate over each test spectrum
    out     = rep(NA_real_,m)
    refname = character(m)
    for( j in 1:m )
        {
        test30  = test30xM[ ,j]

        #   4.3 Area Normalize the Spectra
        sumtest30 = sum(test30)        
        if( is.na(sumtest30) || sumtest30 <= 0 )
            {
            log_level( WARN, "spectrum '%s' has invalid sum=%g. SSI is NA.", specnames(x)[j], sumtest30 )
            next
            }
        test30    = test30 / sumtest30
        
        #   find appropriate reference
        if( numref == 0 )
            {
            #   compute a reference from CCT
            CCT = CCTvec[j]
            
            if( is.na(CCT)  )
                {
                log_level( WARN, "Cannot compute CCT for spectrum '%s'. SSI is NA.", specnames(x)[j] )
                next
                }
            
            if( CCT < 4000 )
                {
                ref = planckSpectra( CCT, wavelength=lambda1, c2=c2 )
                refname[j]  = sprintf("P%g",CCT)
                }
            else if( CCT <= 25000 )
                {            
                ref = daylightSpectra( CCT, wavelength=lambda1 )
                
                # if( (CCT %% 100) == 0 ) CCT = CCT / 100     # D40, D55, D65, etc.
                
                refname[j]  = sprintf("D%g",CCT)
                }                
            else
                {
                log_level( WARN, "spectrum '%s' has CCT = %g > 25000. SSI is NA.", specnames(x)[j], CCT )
                next
                }
                
            ref30   = g.Trap30x301 %*% as.matrix(ref)   # a vector of length 30  
            }        
        else
            {
            k           = min( j, numref )
            ref30       = ref30xNUMREF[ ,k]
            refname[j]  = specnames(reference)[k]
            }
            
        #   4.3 Area Normalize 
        sumref30 = sum(ref30)        
        if( is.na(sumref30) || sumref30 <= 0 )
            {
            log_level( WARN, "reference spectrum %d has invalid sum=%g. SSI is NA.", j, sumref30 )
            next
            }
        ref30 = ref30 / sumref30
            
        out[j]  = computeSSI30( test30, ref30 )
        }
        
    out = round( out, digits=digits )   #  to comply with official SSI spec, set digits=0

    names(out)  = sprintf( "%s_SSI[%s]", specnames(x), refname )
    
    return(out)
    }
    
    
#   test30  30-vector of test samples (380 to 670 by 10nm).     Steps 4.3 and earlier
#   ref30   30-vector of ref  samples (380 to 670 by 10nm).     Steps 4.3 and earlier
computeSSI30 <- function( test30, ref30 )
    {
    #   4.4 Create Normalized Difference Vector
    d30 = test30 - ref30
    
    #   4.5 Create Relative Difference Vector
    dr30 = d30 / (ref30 + 1/length(ref30))
    
    #   4.6 Weighted Relative Difference Vector
    w    = c( c(12,22,32,40,44)/45, rep(1,23), c(11,3)/15 )     # about 3 microseconds
    dw30 = w * dr30
    
    #   4.7 Smooth the Weighted Relative Difference Vector
    s   = stats::filter( c(0,dw30,0), c(0.22,0.56,0.22), method="convolution", sides=2 )[2:31]   # filter() puts NAs on both ends, so subset [2:31]
    
    #   4.8 Calculate SSI Value (to full precision and do not round)
    SSI = 100 - 32*sqrt(sum(s^2))
    
    return(SSI)
    }
    

    
#   this function is designed to be called only 1 time per session
#   it takes about 170 usec    
makeTrap30x301 <- function()
    {
    lambda1     = seq(375,675,by=1)     # fine
    lambda10    = seq(380,670,by=10)    # coarse
    
    out = matrix( 0, length(lambda10), length(lambda1) )
    
    trap    = c( 0.5, rep(1,9), 0.5 )    # length is 1+9+1 = 11
    
    for( i in 1:nrow(out) )
        {
        out[ i, (10*(i-1)+1):(10*i+1) ] = trap
        }
        
    return( out )
    }
    
#--------       UseMethod() calls           --------------#            
        
computeSSI <-function( x, reference=NULL, digits=0, isotherms='mccamy', locus='robertson' )
    {
    UseMethod("computeSSI")
    }
    