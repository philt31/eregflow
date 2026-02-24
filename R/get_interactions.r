
#' @noRd
reorder.names.in.complexes <- function( complexes )
{
    complexes.split <- strsplit( complexes, ":" )
    complexes.sort <- lapply( complexes.split, sort )
    complexes.join <- lapply( complexes.sort, paste0, collapse = ":" )
    unlist( complexes.join )
}


#' Get interactions
#' 
#' It gets interactions from a database, by reading them from an R file 
#' imported from the database. Interactions are selected and formatted with 
#' the standard interaction names as given by 
#' [get.standard.names.in.interactions()]. A CSV table is generated as 
#' reference for the selected interactions. A cache is also generated, for 
#' optional reuse in later executions. 
#' 
#' @param database.name Character string with the name of the database to use. 
#'     Currently supported databases are given by [get.supported.databases()]. 
#' @param database.file.path Character string with the directory and name of an 
#'     R object with interactions imported from the database. This file will 
#'     usually have been generated with the corresponding function 
#'     `import.database.*()`. 
#' @param table.file.path Character string with the directory and name of a CSV 
#'     table to save processed interactions, with both standard and original 
#'     format. 
#' @param cache.file.path Character string with the directory and name of an R 
#'     data file to save and reuse processed interactions. 
#' @param cache.ignore Boolean indicating whether to force execution and 
#'     ignore a possible previous cache. 
#' 
#' @return A data frame with the processed interactions in standard format. 
#' 
#' @export

get.interactions <- function( database.name, database.file.path, 
    table.file.path, cache.file.path, cache.ignore )
{
    stopifnot( database.name %in% get.supported.databases() )
    
    if ( ! file.exists( cache.file.path ) || cache.ignore )
    {
        get.interactions.fname <- sprintf( "get.interactions.%s", 
            database.name )
        
        interactions.data <- do.call( get.interactions.fname, 
            list( database.file.path ) )
            
        # remove auto-interactions
        
        interactions.data <- interactions.data[ 
            interactions.data$source != interactions.data$target, 
        ]
        
        # remove duplicated entries
        
        interactions.edge <- sprintf( "%s|%s", 
            interactions.data$source, 
            interactions.data$target )
        
        interactions.edge.duplicated <- unique( interactions.edge[ 
            duplicated( interactions.edge ) ] )
        
        interactions.edge.remove.idx <- lapply( interactions.edge.duplicated, 
            function( ied ) {
                edge.idx <- which( interactions.edge == ied )
                edge.idx.keep.idx <- which.max( interactions.data$evidence[ 
                    edge.idx ] )
                edge.idx[ - edge.idx.keep.idx ]
            } )
        
        interactions.edge.remove.idx <- unlist( interactions.edge.remove.idx )
        
        interactions.data <- interactions.data[ 
            - interactions.edge.remove.idx, ]
        
        # reorder names in complexes
        
        interactions.data$source <- reorder.names.in.complexes( 
            interactions.data$source )
        
        interactions.data$target <- reorder.names.in.complexes( 
            interactions.data$target )
        
        # save table and cache of processed interactions
        
        interactions.data <- interactions.data[ 
            order( interactions.data$source, interactions.data$target ), ]
        
        write.csv( interactions.data, table.file.path, row.names = FALSE )
        
        interactions.data.std <- interactions.data[ , 
            get.standard.names.in.interactions() ]
        
        save( interactions.data.std, file = cache.file.path )
    }
    else
    {
        object.loaded <- load( cache.file.path )
        stopifnot( object.loaded == "interactions.data.std" )
    }
    
    interactions.data.std
}

