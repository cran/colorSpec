
#   ptransform.colorSpec()
#
#   x           a colorSpec responder (light or material) with M channels

#   primary     an MxM matrix.  The i'th row of primary is a point in the response vector space of x.
#               It is OK for each row to have a single value NA; the missing value is then
#               set so the row sum is 1.
#               The rownames() are the names of the new primaries, and must be present.

#   white       a vector of length M. 
#               white can also be a colorSpec object suitable as input to x;
#               in this case the vector white is set to product( white, x, wavelength='auto' )

#   For primary and white to be valid, the rows of primary must be a basis of R^m,
#   and the coordinates of white w.r.t. this basis must all be non-zero.

#   return value
#               a colorSpec responder with M spectra, each of which is a linear combination of the spectra in x.
#               This implies that
#                   out = multiply( x, mat )    for an appropriate matrix mat.
#               The matrix mat is chosen so that:
#                   1)  the i'th row of primary maps to a non-zero multiple of the elementary vector
#                           e_i = ( 0, 0, ..., 1, 0, ..., 0, 0 )       where the 1 is in the i'th coordinate
#                   2) the vector white maps to the M-vector of all 1s.
#
#               The specnames() of the returned object are set to tolower( rownames(primary) ).

ptransform.colorSpec  <-  function( x, primary, white, digits=Inf   )
    {
    #   check x
    valid   = is.colorSpec(x) &&  grep( "^responsivity", type(x) )
    if( ! valid )
        {
        log_level( ERROR, "x is not a colorSpec responder; type(x)='%s'.", type(x) )
        return(NULL)
        }
    
    m    = numSpectra(x)
    
    valid   = is.numeric(primary) && !is.null(dim(primary)) &&  !is.null(rownames(primary))  &&  all( dim(primary)==c(m,m) ) 
    if( ! valid )
        {
        log_level( ERROR, "primary is invalid. It must be a %dx%d matrix with rownames.", m, m )
        return(NULL)
        }

        
    for( i in 1:m )
        {
        row = primary[ i, ]
        idx = which( is.na(row) )
        if( length(idx) == 0 )  next
        
        if( min(2,m) <= length(idx) )
            {
            log_level( ERROR, "The %th row of primary has %d NAs, which is too many.", i, length(idx) )
            return(NULL)
            }
        
        row[idx]    = 1 - sum(row,na.rm=T)
        primary[ i, ]   = row
        }
    
    if( is.colorSpec(white) )
        {
        w   = product( white, x, wavelength='auto' )
        if( is.null(w) )    return(NULL)
        white   = w
        }
        
    valid   = is.numeric(white)  &&  (length(white)==m)  
    if( ! valid )
        {
        log_level( ERROR, "white is invalid. It must be a numeric vector of length %d.", m )
        return(NULL)
        }
        
    white   = as.numeric(white)     #   strip dim() attribute if present
    
    if( ! is.null(digits)  &&  is.finite(digits)  )
        {
        white.sum   = sum(white)
        if( white.sum == 0 )
            {
            log_level( ERROR, "sum(white) = %g, which is invalid.", white.sum )
            return(NULL)
            }
        
        white.chrom = roundAffine( white/sum(white.sum), digits )
        if( is.null(white.chrom) )  return(NULL)

        #   project white onto ray defined by white.chrom
        white   =   (sum(white*white.chrom) / sum(white.chrom*white.chrom)) * white.chrom
        
        #cat( white.chrom, '\n' )        
        #cat( white, '\n' )
        }

    B   = projectiveMatrix( t(primary), white )
    if( is.null(B) )
        {
        log_level( ERROR, "primary and/or white are degenerate." )
        return(NULL)
        }    
    
    A   = t( solve(B) )
    
    out = multiply( x, A )
    specnames(out) = tolower( rownames(primary) )
    
    attr( out, 'ptransform' )    = list( primary=primary, white=white, A=A )
    
    return( out )
    }
    
    
#   projectiveMatrix()
#    
#   .matrix     invertible matrix, for example a 3x3 matrix with columns the tristimulus coordinates of RGB primaries
#   .unit       non-zero vector.  For example the tristimulus coordinates of white.
#
#   return      square matrix  B, so that   
#               B = matrix  %*%  diag(a)  <=>   each column of B is a multiple of the corresponding column in .matrix
#               B %*% 1  =  .unit.      (1 is the vector of all 1s)
#
#   so for colors, B maps RGB to XYZ
#
#    Another way to write these properties:
#        B %*% I = matrix     up to multiples of the columns
#        B %*% 1  =  .unit
#   So I and 1 are the *standard* projective basis,
#   and .matrix and .unit are a different one

projectiveMatrix  <-  function( .matrix, .unit )
    {
    a   = try( solve( .matrix, .unit ), silent=TRUE )
    
    if( ! is.numeric(a) ) return(NULL)
    
    ran = range( abs(a) )   #; print(ran)
    
    if( ran[1] <= 1.e-6 * ran[2] ) return(NULL)
    
    return( .matrix  %*%  diag(a) )
    }
    
    
    
    
#--------       UseMethod() calls           --------------#                    
        
ptransform <- function( x,  primary, white, digits=Inf   )
    {
    UseMethod("ptransform")
    }    
        
    
