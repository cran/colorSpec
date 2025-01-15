


#   XYZ     XYZ values adapted to D65
#           can be a simple 3-vector, or a matrix with 3 columns
#
#   returns a matrix with 3 columns, or NULL in case of ERROR

sRGBfromXYZ  <-  function( XYZ )
    {
    if( length(XYZ) == 3 )
        XYZ = matrix( XYZ, 1, 3 )
    
    if( length(dim(XYZ))!=2  ||  ncol(XYZ)!=3 )
        {
        log_level( ERROR, "XYZ is invalid" )
        return(NULL)
        }
        
    #   load the matrix.  Note that we really create the transpose of published matrices.  byrow=FALSE
    
    #   matrix with only 4-digit precision, taken from  http://en.wikipedia.org/wiki/SRGB    
    #   mat = matrix(  c( 3.2406, -1.5372, -0.4986, -0.9689, 1.8758, 0.0415, 0.0557, -0.2040, 1.0570), 3, 3, byrow=FALSE )
    
    #   matrix with high precision, taken from   http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html    
    mat = matrix( c( 3.2404542, -1.5371385, -0.4985314,
                    -0.9692660,  1.8760108,  0.0415560,
                     0.0556434, -0.2040259,  1.0572252  ), 3, 3, byrow=FALSE )
        
    out =  XYZ  %*%  mat

    dim(out)    = dim(XYZ)
    
    rownames(out)   = rownames(XYZ)
    colnames(out)   = c('R','G','B')
    
    return(out)
    }    
        
        

#   XYZ     XYZ values adapted to D65
#           can be a simple 3-vector, or a matrix with 3 columns
#   space   RGB space specified with character string
#   returns a matrix with 3 columns, or NULL in case of ERROR
#
#   the right way to do this - a precomputed 3D array with names
RGBfromXYZ  <-  function( XYZ, space )
    {
    if( length(XYZ) == 3 )
        XYZ = matrix( XYZ, 1, 3 )
    
    if( length(dim(XYZ))!=2  ||  ncol(XYZ)!=3 )
        {
        log_level( ERROR, "XYZ is invalid" )
        return(NULL)
        }

    if( ! is.character(space) )
        {
        log_level( ERROR, "typeof(space) = '%s', but must be character.", typeof(space) )
        return(NULL)
        }

    space   = tolower( gsub( '[ ]+', '', space[1] ) )
    
    if( space == 'srgb' )
        {
        mat = matrix( c( 3.2404542, -1.5371385, -0.4985314,
                        -0.9692660,  1.8760108,  0.0415560,
                         0.0556434, -0.2040259,  1.0572252  ), 3, 3, byrow=TRUE )
        }
    else if( space == 'adobergb' )
        {
        #   matrix with high precision, taken from   http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html    
        mat = matrix(  c(    2.0413690, -0.5649464, -0.3446944,
                            -0.9692660,  1.8760108,  0.0415560,
                             0.0134474, -0.1183897,  1.0154096 ), 3, 3, byrow=TRUE )
        }
    else
        {
        log_level( ERROR, "RGB space '%s' unknown", space )
        return(NULL)
        }
                
    out =  XYZ  %*%  t(mat)

    dim(out)    = dim(XYZ)
    rownames(out)   = rownames(XYZ)
    colnames(out)   = c('R','G','B')
    
    return(out)
    }    
     
    


#   RGB     linear RGB any sort of array
#           it is OK if val is a matrix, and then the return value is a matrix of the same shape
#   gamma   name of transfer function, supported are 'sRGB'
#           or numeric gamma of the display, so output is (linear)^(1/gamma)
#
#   return  first clips to [0,1], and then maps [0,1] to [0,1].
#           in case of ERROR it logs a message and returns the clipped values only
DisplayRGBfromLinearRGB <- function( RGB, gamma='sRGB' )
    {
    out = as.numeric(RGB)
    
    out = pmin( pmax( RGB, 0 ), 1 )
    
    if( is.character(gamma) )
        {
        if( tolower(gamma[1]) == 'srgb' )
            out = ifelse( out <= 0.0031308,    12.92 * out,  1.055 * out^(1/2.4) - 0.055 )
        else
            log_level( ERROR, "gamma is invalid" )
        }
    else if( is.numeric(gamma) && 0 < gamma[1] )
        out = out ^ (1/gamma[1])
    else
        log_level( ERROR, "gamma is invalid" )
        
    
    dim(out)        = dim(RGB)
    dimnames(out)   = dimnames(RGB) #; print( out )
    
    return( out )
    }
    
    
    
xyY_from_XYZ <- function( XYZ )
    {
    out = XYZ
    if( is.null(dim(out))  &&  length(out) == 3 )
        {
        #   make into a matrix
        dim(out) = c(1,3)
        }    
        
    if( ncol(out) != 3 )    return(NULL)
    
    sums    = rowSums(out)
    
    out = cbind( out[ ,1]/sums, out[ ,2]/sums, out[ ,2] )
    
    dim(out)        = dim(XYZ)
    rownames(out)   = rownames(XYZ)
    colnames(out) = c( 'x', 'y', 'Y' )
        
    return( out )
    }
    
    
XYZ_from_xyY    <- function( x, y, Y )
    {
    out = cbind( x, y, 1 - x - y )
    
    out = out * (Y / y )
    
    colnames( out ) = c( 'X', 'Y', 'Z' )
    
    if( length(out) == 3 )  dim(out) = NULL
    
    return( out )
    }
    
    
    
#   requires private data frame dataIlluminants, which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
officialXYZ <- function( name )
    {
    XYZ             = as.matrix( subset( dataIlluminants, select=c('X','Y','Z') ) )
    rownames(XYZ)   = dataIlluminants$Name
    
    out = matrix( NA_real_, length(name), 3 )
    rownames(out)   = name
    colnames(out)   = c('X','Y','Z')
    
    common  = intersect( name, rownames(XYZ) )

    out[ common, ] = XYZ[ common, ]
    
    return(out)
    }
    

    
        
    
#   .xy     from 1931    a 2-vector or a matrix with 2 columns
uv_from_xy <- function( .xy, .year=1960 )
    {
    if( .year == 1931 )
        {
        #   no change !
        return(.xy)    
        }
    else if( .year == 1960 )
        {
        coeff = 6
        cnames  = c( "u", "v" )
        }
    else if( .year == 1976 )
        {
        coeff = 9
        cnames  = c( "u'", "v'" )
        }
    else
        {
        log_level( ERROR, 'year %g is invalid\n', .year )
        return(NULL)
        }   
        
    if( is.null(dim(.xy))  &&  length(.xy) == 2 )  
        dim(.xy) = c(1,2)   
    
    out = cbind( 4 * .xy[ ,1], coeff * .xy[ ,2] ) /  ( -2 * .xy[ ,1]  +  12 * .xy[ ,2] + 3 )
    

    #   dim(out) = dim(.xy)
    #   print( dim(out) )

    colnames(out)   = cnames
    rownames(out)   = rownames(.xy)

            
    #   print( out )
    
    return( out )
    }
    
    
#   XYZ a 3-vector or matrix with 3 columns    
UVW_from_XYZ <- function( XYZ, XYZ.white )
    {
    xy          = xyY_from_XYZ( XYZ.white )[1:2]
    uv          = uv_from_xy( as.double(xy) )
    uvY.white   = c( uv, XYZ.white[2] )
    
    if( length(XYZ) == 3 )  dim(XYZ) = c(1,3)
    
    xy  = xyY_from_XYZ( XYZ )[ ,1:2]
    uv  = uv_from_xy( xy )
    
    uvY = cbind( uv, XYZ[ ,2] )
    
    UVW = UVW_from_uvY( uvY, uvY.white )

    return( UVW )
    }    
    
    
#   uvY a 3-vector or matrix with 3 columns    
UVW_from_uvY <- function( uvY, uvY.white )
    {
    UVW = uvY   #  get the dimension right    
    
    if( length(UVW) == 3 )  dim(UVW) = c(1,3)
    
    Y.rel   = uvY[ ,3] / uvY.white[3]   # relative Y
    
    UVW[ ,3]    = 25 * (100*Y.rel)^(1/3) - 17
    
    UVW[ ,1]    = 13*UVW[ ,3]*( uvY[ ,1] - uvY.white[1] )
    UVW[ ,2]    = 13*UVW[ ,3]*( uvY[ ,2] - uvY.white[2] )
    
    colnames(UVW)   = c('U','V','W')
    
    return( UVW )
    }    

    