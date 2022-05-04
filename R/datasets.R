  

saveDatasets  <- function( .path="../data/colorSpec.rda" )
    {
    savevec = character(0)
    

    ##-------------------------     illuminant power functions     -----------------------------##
    #   illuminant A
    path    =  "../inst/extdata/illuminants/A.1nm.txt"
    A.1nm = readSpectra( path )
    organization( A.1nm ) = 'vector'
    A.1nm   = A.1nm / 100
    specnames(A.1nm) = "A"
    savevec = c( savevec, "A.1nm" )
    
    #   illuminant B
    path    =  "../inst/extdata/illuminants/B.txt"
    B.5nm = readSpectra( path )
    organization( B.5nm ) = 'vector'
    B.5nm   = B.5nm / 100    
    specnames(B.5nm) = "B"
    savevec = c( savevec, "B.5nm" )
    
    #   illuminant C
    path    =  "../inst/extdata/illuminants/C.txt"
    C.5nm = readSpectra( path )
    organization( C.5nm ) = 'vector'
    C.5nm   = C.5nm / 100        
    specnames(C.5nm) = "C"
    savevec = c( savevec, "C.5nm" )
    
    
    #   D50 with Standard UV
    #D50.5nm = readSpectra( "../inst/extdata/illuminants/D50_1.0.sp" ) #; print( str(D50) )
    #organization( D50.5nm ) = 'vector'
    #   D50.5nm = D50.5nm / 100       #   we know that D50_1.0.sp is normalized to 100, but we prefer 1
    #specnames(D50.5nm) = "D50.Power"
    #savevec = c( savevec, "D50.5nm" )
    
    #   D50 with Standard UV
    correction  = 14388 / 14380     # note 5, page 69 in CIE:15:2004
    D50.10nm    = daylightSpectra( correction*5000, wavelength=seq(300,830,by=10), components=daylight1964, roundMs=TRUE )
    #wave        = seq(300,830,by=5)    
    #vec.5nm     = approx( wavelength(D50.10nm), coredata(D50.10nm), wave )$y
    D50.5nm     = resample( D50.10nm, seq(300,830,by=5) , method='linear')
    #   vec.5nm     = round( vec.5nm, 5 )
    D50.5nm     = round( D50.5nm, 5 )
    #   D50.5nm     = colorSpec( vec.5nm, wave, 'energy', 'vector' )
    specnames(D50.5nm)  = 'D50'
    desc                = "computed from the 1964 daylight components: S0, S1, and S2"   
    metadata(D50.5nm)   = list( description=desc )    
    savevec     = c( savevec, "D50.5nm" )

    #   illuminant D65
    path  = "../inst/extdata/illuminants/D65.1nm.txt" 
    D65.1nm = readSpectra( path )
    organization( D65.1nm ) = 'vector'
    D65.1nm   = D65.1nm / 100        
    specnames(D65.1nm) = "D65"
    savevec = c( savevec, "D65.1nm" )

    path  = "../inst/extdata/illuminants/D65.5nm.txt" 
    D65.5nm = readSpectra( path )
    organization( D65.5nm ) = 'vector'
    D65.5nm   = D65.5nm / 100            
    specnames(D65.5nm) = "D65"
    savevec = c( savevec, "D65.5nm" )

    
    
    
    ##---------------       series D daylight characteristic spectra      ---------------##
    path  = "../inst/extdata/illuminants/daylight1964.txt"
    daylight1964 = readSpectraXYY( path )
    organization(daylight1964)  = mostEfficientOrganization(daylight1964)    
    savevec = c( savevec, "daylight1964" )
    
    #  smoothed version of series D daylight from 2013
    path  = "../inst/extdata/illuminants/daylight2013.txt"
    daylight2013 = readSpectraXYY( path )
    organization(daylight2013)  = mostEfficientOrganization(daylight2013)        
    savevec = c( savevec, "daylight2013" )
    
    ##---------------     series F Fluorescent illuminants          ----------------------------##    
    path  = "../inst/extdata/illuminants/Fs.5nm.txt"
    Fs.5nm = readSpectraXYY( path )  
    organization( Fs.5nm ) = mostEfficientOrganization( Fs.5nm )
    savevec = c( savevec, "Fs.5nm" )
    
    
    ##--------------    Solar Irradiance and Atmospheric Transmittance     -----------------##
    path    = "../inst/extdata/illuminants/ASTMG173.txt"
    solar.irradiance    = readSpectra( path )
    #   solar.irradiance    = subset( solar.irradiance, c(1,3) )    # throw away the global which includes skylight
    solar.irradiance    = resample( solar.irradiance, 280:1000 )
    specnames(solar.irradiance) = c( "AirMass.0", "GlobalTilt", "AirMass.1.5" )
    #   print( summary(solar.irradiance) )
    organization(solar.irradiance)  = mostEfficientOrganization(solar.irradiance)
    savevec = c( savevec, "solar.irradiance" )
    
    atmosphere2003  = solar.irradiance[ ,3] / solar.irradiance[ ,1]
    atmosphere2003  = colorSpec( atmosphere2003, wavelength(solar.irradiance), quantity='transmittance', specnames="AirMass.1.5" )
    desc    = "from ASTM G163-03"
    desc    = c( desc, "Standard Tables for Reference Solar Spectral Irradiances: Direct Normal and Hemispherical on 37\u00B0 Tilted Surface" )     # degree symbol:  °
    desc    = c( desc, "transmittance of the atmosphere through Air Mass 1.5 is the quotient AM1.5/AM0" )    
    metadata(atmosphere2003)    = list( header=desc )
    #   print( summary(solar.transmittance) )
    organization(atmosphere2003)  = mostEfficientOrganization(atmosphere2003)      
    savevec = c( savevec, "atmosphere2003" )
    
    
    ##------------------        photons       --------------##
    path    = "../inst/extdata/sources/F96T12-GR8D.txt"
    F96T12  = readSpectra( path )
    specnames(F96T12) = "F96T12"
    #   print( summary(F96T12) )
    organization(F96T12)  = mostEfficientOrganization(F96T12)
    savevec = c( savevec, "F96T12" )
    
    
    ##-------------------------     human eyes     -----------------------------##
    #   the CMFs of 1931 - 2-degree
    path        = "../inst/extdata/eyes/ciexyz31_1.csv"
    xyz1931.1nm = readSpectraXYY( path )
    organization(xyz1931.1nm)  = mostEfficientOrganization(xyz1931.1nm)        
    savevec = c( savevec, "xyz1931.1nm" )

    path        = "../inst/extdata/eyes/xyz1931.5nm.txt"
    xyz1931.5nm = readSpectraXYY( path )
    organization(xyz1931.5nm)  = mostEfficientOrganization(xyz1931.5nm)        
    savevec = c( savevec, "xyz1931.5nm" )

    #   the CMFs of 1964 - 10-degrees
    path        = "../inst/extdata/eyes/ciexyz64_1.csv"
    xyz1964.1nm = readSpectraXYY( path )
    organization(xyz1964.1nm)  = mostEfficientOrganization(xyz1964.1nm)        
    savevec = c( savevec, "xyz1964.1nm" )

    path        = "../inst/extdata/eyes/xyz1964.5nm.txt"
    xyz1964.5nm = readSpectraXYY( path )
    organization(xyz1964.5nm)  = mostEfficientOrganization(xyz1964.5nm)        
    savevec = c( savevec, "xyz1964.5nm" )

    #   Judd and Vos cone fundamentals
    lms1971.5nm = readSpectraXYY( "../inst/extdata/eyes/lms1971.txt" )
    organization(lms1971.5nm)  = mostEfficientOrganization(lms1971.5nm)        
    savevec = c( savevec, "lms1971.5nm" )
    
    #   Stockman and Sharpe cone fundamentals
    lms2000.1nm = readSpectraXYY( "../inst/extdata/eyes/lms2000.1nm.csv" )
    organization(lms2000.1nm)  = mostEfficientOrganization(lms2000.1nm)        
    savevec = c( savevec, "lms2000.1nm" )
    
    #   luminsivity functions
    #   the union of all the wavelength vectors is 360:830, which is the wavelength of photopic1924.
    #   set extrap=0 to zero-pad the smaller ones    
    pathvec = c( "../inst/extdata/eyes/photopic1924.1nm.csv", "../inst/extdata/eyes/scotopic1951.1nm.csv",
                "../inst/extdata/eyes/photopic1978.1nm.csv", "../inst/extdata/eyes/photopic2008.1nm.csv" )
    luminsivity.1nm = readSpectra( pathvec, wave=360:830, extrap=0 )
    specnames(luminsivity.1nm)  = sub( "^([A-z0-9]+).*", '\\1', specnames(luminsivity.1nm) )
    #   print( specnames(luminsivity.1nm) )
    savevec = c( savevec, "luminsivity.1nm" )
    
    
    ##-------------------------     animal eyes     -----------------------------##    
    
    #   Higher Passerines bird vision - 4-channel
    HigherPasserines    = readSpectraXYY( "../inst/extdata/eyes/BirdEyes.txt" )
    HigherPasserines    = subset( HigherPasserines, 1:4 )
    HigherPasserines    = multiply( HigherPasserines, 0.01 )
    specnames(HigherPasserines) = c("UV","Short","Medium","Long")
    quantity(HigherPasserines)  = "photons->neural"
    organization(HigherPasserines)  = mostEfficientOrganization(HigherPasserines)         
    savevec = c( savevec, "HigherPasserines" )



    ##---------------       cameras     ----------------------------##        
    
    Flea2.RGB   = readSpectra( "../inst/extdata/cameras/Flea2-spectral.txt", seq(360,800,by=10) )
    organization(Flea2.RGB)  = mostEfficientOrganization(Flea2.RGB)       
    savevec = c( savevec, "Flea2.RGB" )
    
    # ideal BT.709 camera with negative lobes impossible to actually build.  D65 maps to RGB=(1,1,1)
    P = matrix( c(0.64,0.33,NA,  0.3,0.6,NA, 0.15,0.06,NA ), 3, 3, byrow=T )
    rownames(P) = c('R','G','B')    
    BT.709.RGB  = ptransform( xyz1931.1nm, P, D65.1nm )  
    quantity(BT.709.RGB)    = "energy->electrical"    
    desc    = "This is a theoretical RGB camera"
    desc    = c( desc, "They acquire RGB components for display using BT.709 primaries, which are the same as sRGB primaries." )
    desc    = c( desc, "This theoretical camera satisfies the Maxwell-Ives condition, but has negative lobes." )    
    desc    = c( desc, "Compare with Figure 26.5 on page 302 of:" )
    desc    = c( desc, "Poynton, Charles" )
    desc    = c( desc, "Digital Video and HD - Algorithms and Interfaces." )
    desc    = c( desc, "Second Edition. 2012." )
    metadata(BT.709.RGB,add=TRUE)   = list( header=desc )    
    metadata(BT.709.RGB,add=TRUE)   = list( path=NULL )         # erase path   
    organization(BT.709.RGB)  = mostEfficientOrganization(BT.709.RGB)         
    savevec = c( savevec, "BT.709.RGB" )
    
    # ideal Adobe.RGB camera with negative lobes impossible to actually build.  D65 maps to RGB=(1,1,1)
    # it's the same as BT.709.RGB, except for the matrix P
    P = matrix( c(0.64,0.33,NA,  0.21,0.71,NA, 0.15,0.06,NA ), 3, 3, byrow=T )
    rownames(P) = c('R','G','B')    
    Adobe.RGB  = ptransform( xyz1931.1nm, P, D65.1nm )  
    quantity(Adobe.RGB)    = "energy->electrical"        
    desc    = "This is a theoretical RGB camera."
    desc    = c( desc, "It acquires RGB components for display using Adobe RGB primaries." )
    desc    = c( desc, "This theoretical camera satisfies the Maxwell-Ives condition, but has negative lobes." )
    metadata(Adobe.RGB,add=TRUE)    = list( header=desc )
    metadata(Adobe.RGB,add=TRUE)    = list( path=NULL )         # erase path      
    organization(Adobe.RGB)  = mostEfficientOrganization(Adobe.RGB)         
    savevec = c( savevec, "Adobe.RGB" )
    
    
    # ideal ACES.RGB camera that encompasses all possible colors, without negative lobes.  D60.ACES maps to RGB=(1,1,1)
    P = matrix( c(0.73470,0.26530,NA,  0,1,NA,  0.00010,-0.07700,NA ), 3, 3, byrow=T )
    rownames(P) = c('R','G','B') 
    white = c(0.32168,0.33767)
    white = c( white, 1-sum(white) ) / white[2]           # D60.ACES
    ACES.RGB  = ptransform( xyz1931.1nm, P, white ) 
    ACES.RGB  = calibrate( ACES.RGB, illuminantE(1,wavelength=wavelength(ACES.RGB)), 1, method='scaling' )    
    quantity(ACES.RGB)    = "energy->electrical"        
    desc    = "This is a theoretical RGB camera."
    desc    = c( desc, "It acquires RGB components for display using ACES RGB primaries.")
    desc    = c( desc, "S-2008-001. Academy Color Encoding Specification (ACES)  Annex C." )
    desc    = c( desc, "This theoretical camera satisfies the Maxwell-Ives condition, and is everywhere non-negative." )
    metadata(ACES.RGB,add=TRUE)    = list( header=desc )
    metadata(ACES.RGB,add=TRUE)    = list( path=NULL )         # erase path      
    organization(ACES.RGB)  = mostEfficientOrganization(ACES.RGB)         
    savevec = c( savevec, "ACES.RGB" )

    
    ##---------------       materials     ----------------------------##    
    
    #ColorChecker = readSpectra( "../inst/extdata/targets/CC_Avg20_spectrum_XYY.txt" )
    #savevec = c( savevec, "ColorChecker" )
    
    Hoya = readSpectra( "../inst/extdata/objects/Hoya.txt" )
    organization(Hoya)  = mostEfficientOrganization(Hoya)         
    savevec = c( savevec, "Hoya" )
        
        
    ##------------------   material responders  (scanners)  ----------##
    
    scanner.ACES    = readSpectraXYY( "../inst/extdata/scanners/SMPTE-ST-2065-2.txt" )
    #scanner.ACES    = normalize( scanner.ACES, "L1" )   #  normalize so perfect-reflecting-diffuser response is RGB=(1,1,1)
    #metadata(scanner.ACES)  = list( normalized=TRUE )    
    scanner.ACES    = calibrate( scanner.ACES, stimulus=neutralMaterial(1,wavelength(scanner.ACES)), response=1, method='scaling' ) # so perfect-reflecting-diffuser response is RGB=(1,1,1)
    organization(scanner.ACES)  = mostEfficientOrganization(scanner.ACES)      
    #   summary( scanner.ACES )
    savevec = c( savevec, "scanner.ACES" )
        

    ##  finally ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE
    
    return( invisible(TRUE) )
    }
    
    
#   an advantage of the private data in "sysdata.rda" is that these
#   do not have to be documented, and therefore exposed    
#
#   p.Ma                    3x3 adaptation matrices


savePrivateDatasets  <- function( .path="sysdata.rda" )
    {
    savevec = character(0)
        
    if( FALSE )
    {
    #   these have been moved to package 'spacesXYZ'
    ##---------------       CCT table    ------------------------##
    path    = "../inst/extdata/illuminants/dataCCT.txt"
    dataCCT = read.table( path, sep='\t', header=T, stringsAsFactors=F )
    attr(dataCCT,"description") = readComments( path )
    savevec = c( savevec, "dataCCT" )
    }
    
    ##---------------       illuminants table    ------------------------##
    path    = "../inst/extdata/illuminants/illuminants.txt"
    dataIlluminants = read.table( path, sep='\t', header=T, stringsAsFactors=F )
    attr(dataIlluminants,"description") = readComments( path )
    savevec = c( savevec, "dataIlluminants" )

    
    #---------------       spectra for CRI     ------------------------##
    path    = "../inst/extdata/targets/TCSforCRI.txt"
    TCSforCRI = readSpectra( path )
    savevec = c( savevec, "TCSforCRI" )

    #---------------       lens absorbance dependence on age     -------##
    path    = "../inst/extdata/eyes/LensAbsorbance1987.txt"
    LensAbsorbance1987 = readSpectra( path )
    organization(LensAbsorbance1987)  = mostEfficientOrganization(LensAbsorbance1987)        
    savevec = c( savevec, "LensAbsorbance1987" )
    
    #   list of adaptation matrices
    p.Ma    = list()
    p.Ma[[ "Bradford" ]]            = matrix( c(0.8951,0.2664,-0.1614,  -0.7502,1.7135,0.0367,  0.0389,-0.0685,1.0296), 3, 3, byrow=T )
    p.Ma[[ "VonKries" ]]            = matrix( c(0.40024,0.7076,-0.08081,  -0.2263,1.16532,0.0457,  0,0,0.91822), 3, 3, byrow=T )
    p.Ma[[ "MCAT02" ]]              = matrix( c( 0.7328, 0.4296, -0.1624,  -0.7036, 1.6975, 0.0061, 0.0030, 0.0136, 0.9834 ), 3, 3, byrow=T )
    p.Ma[[ "Bianco+Schettini" ]]    = matrix( c( 0.8752, 0.2787, -0.1539,  -0.8904, 1.8709, 0.0195, -0.0061, 0.0162, 0.9899 ), 3, 3, byrow=T )
    p.Ma[[ "scaling" ]]             = diag(3)
    
    for( k in 1:length(p.Ma) )
        {
        rownames( p.Ma[[k]] )   = c('L','M','S')
        colnames( p.Ma[[k]] )   = c('X','Y','Z')
        }
        
    savevec = c( savevec, "p.Ma" )    

    ##  finally ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE
    
    return( invisible(TRUE) )
    }    
    
    
pingDatasets  <- function( .path="../data/colorSpec.rda", .verbose=FALSE )    
    {
    theName     = load(.path)
    print( theName )
    
    if( 0 < length(theName)  &&  .verbose )
        {
        for( k in 1:length(theName ) )
            {
            obj = get( theName[k] )
            cat( '\n', theName[k], '\n' )
            print( str(obj) )
            }
        }
        
    return( invisible(T) )
    }
    
    
loadDatasets  <- function( .path="../data/colorSpec.rda" )
    {
    load( .path, parent.frame(n=2) )
    }    
    

refreshDatasets  <- function( .path="../data/colorSpec.rda" )
    {
    saveDatasets( .path )
    load( .path, parent.frame(n=2) )
    }
  

readComments <- function( .path )
    {
    line    = readLines( .path, n=1024 )
    
    out = line[ grepl( "^[ \t]*#", line ) ]
    
    if( length(out) == 0 )  out = NULL
    
    return( out )
    }
    
mostEfficientOrganization  <- function( obj )
    {
    if( organization(obj) == 'df.row' ) 
        {
        if( ncol(obj) <= 1 )
            #   no extradata, only 1 spectrum
            return( 'vector' )
        else
            #   must preserve the extradata
            return( 'df.row' )
        }
        
    if( numSpectra(obj) == 1 )
        return( 'vector' )
    else
        return( 'matrix' )
    }