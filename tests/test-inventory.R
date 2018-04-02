


library( colorSpec )

inventory <- function()
    {
    pos = 2
    
    theNames    = ls(pos)
    
    theList = sapply( theNames, function(x) get(x, pos=pos )  )

    n = length(theList)    
    
    if( n == 0 ) { return(NULL) }    
    
    class1  <- function( x )    { return( class(x)[1] ) }
    
    classvec    = sapply(theList,class1)
    mask        = classvec == 'colorSpec' 
    theNames    = theNames[ mask ]
    theList     = theList[ mask ]
    
    out     = data.frame( name=theNames,     
                                size=sapply( theList, object.size ), 
                                N=sapply( theList, numWavelengths ),                                
                                M=sapply( theList, numSpectra ),
                                P=sapply( theList, function( x )  { numWavelengths(x)*(1+numSpectra(x)) } ),
                                bpp=sapply( theList, function( x )  { object.size(x)/(numWavelengths(x)*(1+numSpectra(x))) } ),
                                organization=sapply( theList, organization ),
                                type=sapply( theList, type ),                                           
                                quantity=sapply( theList, quantity ),         
                                calibrate=sapply( theList, function( x )  { ! is.null( attr(x,"calibrate") ) } ),                                
                                stringsAsFactors=F )
                                
    #out     = subset( out, class == 'colorSpec' )
    
    #out     = cbind( out, organization=sapply( out$name, function(x) { organization(get(x,pos=pos)) } ) )
    #out     = cbind( out, quantity=sapply( out$name, function(x) { quantity(get(x,pos=pos)) } ) )

    rownames(out)   = 1:nrow(out)
    
    cat( sprintf( "Total size: %d\n", sum(out$size) ) )
    
    return( out )
    }
    
options( width=144 )    
print( inventory(), width=144 )
