

#   the same as logger::layout_logging, except put fn() between timestamp and the msg    
layout_mine <- structure(
    function(level, msg, namespace="colorSpec",
                                    .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame())
        {
        # cat( "obj_addr()=", obj_addr( .topcall[[1L]] ), '\n' )
        # cat( "deparse1 =", deparse1( .topcall[[1L]] ), '\n' )
        
        fn  = deparse1( .topcall[[1L]] )
        
        paste0( attr(level, 'level'), ' [', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '] ', namespace, "::", fn, '(). ', msg )
        },
    generator = quote(layout_mine())
)



#   the same as logger::appender_console(), except maybe stop on ERROR or FATAL
appender_mine <- structure(
    function(lines)
        {
        cat(lines, file=stderr(), sep = '\n' )
        
        # cat( "stop =", .Options$colorSpec.stoponerror, '\n', file=stderr() )
        
        if( any( grepl("^FATAL",lines ) )  )
            {
            base::stop( "Stopping, because log level is FATAL.", call.=FALSE )
            }

        #   test for STOP
        if( .Options$colorSpec.stoponerror  &&  any( grepl("^ERR",lines ) )  )
            {
            base::stop( "Stopping, because log level is ERROR, and option colorSpec.stoponerror is TRUE.", call.=FALSE )
            }
        },
    generator = quote(appender_mine())
    )




#   similar to logger::formatter_logging() but handles non-character things with print()
#   this did not work because when multiple lines are returned, they are converted to a single string
formatter_mine <- function(...,
                              .logcall = sys.call(),
                              .topcall = sys.call(-1),
                              .topenv = parent.frame() ) {
  params <- list(...)
  .logcall <- substitute(.logcall)

  if (is.character(params[[1]])) {
    return(do.call(sprintf, params, envir = .topenv))
  }

  return( capture.output( print( params[[1]] ) ) )
}
attr(formatter_mine, "generator") <- quote(formatter_mine())
