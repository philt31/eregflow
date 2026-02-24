
erfa.get.path.gene.new <- function( rg.subgraph, gene.from, gene.to.net, 
    path.gene, path.gene.current )
{

    if ( is.null( path.gene ) )
        return( NULL ) 
    
    gene.connected <- path.gene[ 
        apply( distances( rg.subgraph, path.gene, gene.from, mode = "in" ), 1, 
            function( dist.row ) any( dist.row < Inf ) ) & 
        apply( distances( rg.subgraph, path.gene, gene.to.net, mode = "out" ), 
            1, function( dist.row ) any( dist.row < Inf ) )
    ]
    
    if ( length( gene.connected ) == 0 )
        gene.connected <- NULL
    
    return( gene.connected )
}


erfa.calculate.regulation <- function( rg.pathway, regulation.step, gene.from, 
    gene.to.net )
{
    rg.pathway.step <- as_adjacency_matrix( rg.pathway )
    
    edge.end <- ends( rg.pathway, E( rg.pathway ) )
    
    for ( edge.idx in 1 : nrow( edge.end ) )
    {
        edge.v1 <- edge.end[ edge.idx, 1 ]
        edge.v2 <- edge.end[ edge.idx, 2 ]
        
        rg.pathway.step[ edge.v1, edge.v2 ] <- 
            regulation.step[ edge.v1, edge.v2 ]
    }
    
    rg.pathway.flow <- solve( Diagonal( nrow( rg.pathway.step ) ) - 
        rg.pathway.step )
    
    rg.pathway.flow[ match( gene.from, rownames( rg.pathway.step ) ), 
        match( gene.to.net, rownames( rg.pathway.step ) ) ]
}


#' Plot pathway representations obtain from ERFA model
#' 
#' The function plots partial pathway representations between genes or set of genes, in 
#' increasing amounts of genes included and flow explained. 
#' 
#' The call is ignored if `gene.from` or `gene.to`, after removing overlapping 
#' values from `gene.to`, is `NULL` or has no elements. 
#' 
#' @param rg An [igraph] graph object with the interaction graph, usually 
#'     obtained [build.regulation.graph()]. 
#' @param regulation.step A matrix with one-step regulation flow values, 
#'     usually obtained with [calculate.regulation.matrices()]. 
#' @param regulation.flow A matrix with all-step regulation flow values, 
#'     usually obtained with [calculate.regulation.matrices()]. 
#' @param gene.from Character vector with gene names as source of regulation. 
#' @param gene.to Character vector with gene names as target of regulation. 
#'     Elements also in `gene.from` are ignored. 
#' @param plot.dir Character string with the directory to save the plot. 
#' @param plot.file.name.base Character string with a prefix for the plot file 
#'     names. If `NULL`, default values are generated depending on the number 
#'     of genes in `gene.from` and `gene.to`. 
#' @param plot.title.base Character string with a prefix for the plot titles. 
#'     If `NULL`, default values are generated depending on the number of genes 
#'     in `gene.from` and `gene.to`. 
#' @param pathway.regulation.approx Float with the desired approximation of 
#'     pathway coverage.
#' @param pathway.n.max Integer with the maximum number of pathways to plot. 
#' @param pathway.gene.max Integer with the maximum number of genes to be 
#'     plotted in any pathway.
#' @param alpha a vector containing relative perturbation magnitudes for 
#' source of perturbation. If NULL, alpha is set to a vector with
#' components are equal to 1 for gene.from 
#' genes and 0 for other genes.
#'
#' @param plot.moa Should pathways be saved as plots to PDF files?
#' 
#' @return `NULL`, invisibly.
#' 
#' @examples
#' import.database.omnipath("./database_omnipath.rda", TRUE)
#' interacs <- get.interactions( "omnipath", 
#'                               "./database_omnipath.rda", "./table_interacs.csv",
#'                               "./cache_interacs.rda", FALSE )
#' 
#' reg.graph <- build.regulation.graph( interacs, NULL, "./cache_reg_graph.rda", 
#'                                      FALSE )
#' 
#' reg.matrix <- calculate.regulation.matrices( reg.graph, 
#'                                            "./cache_reg_matrix.rda", FALSE )
#' 
#' reg.step <- reg.matrix$step
#' reg.flow <- reg.matrix$flow
#' path <- system.file("extdata", "sample_data.RData", package = "regflow")
#' load(path)
#' 
#' # Compute ERFA model
#' regul.matrices.refined <- erfa( reg.graph, reg.step,  stimulationTrain, train,
#'                                 plot.dir=".", n.iteration=30,
#'                                 stepsize=5, 
#'                                 stimulationTest = stimulationTest,
#'                                 expressionTest = test, fractionValid = 0.5)
#' 
#' # Compute MoAs for one drug perturbation

#' completeset <- cbind(train,test)

#' drug_to_test <- "mitoxantrone"
#' id <- which(colnames(completeset) == drug_to_test)
#' 
#' largestLFC=order(abs(completeset[,id]),decreasing = T)[1:10]
#' to.genes=rownames(completeset)[largestLFC]
#' 
#' alpha <- as.matrix(regul.matrices.refined$alpha[[id]])
#' rownames(alpha) <- rownames(reg.step)
#' 
#' from.genes<- rownames(stimulationTest)[stimulationTest[,drug_to_test]!=0]
#' 
#' erfa.discovery(reg.graph,regul.matrices.refined$step, regul.matrices.refined$flow,
#'                from.genes, to.genes,
#'                pathway.n.max = 10, pathway.gene.max = 1000, 
#'                plot.dir="./pathways", alpha=alpha, plot.moa = T)
#'  
#' @export

erfa.discovery <- function( rg, regulation.step, regulation.flow, gene.from, 
    gene.to, plot.dir, plot.file.name.base = NULL, plot.title.base = NULL, 
    pathway.regulation.approx = 0.05, pathway.n.max = 10, 
    pathway.gene.max = 20, alpha = NULL, plot.moa = F )
{
    if(is.null(alpha)){
      alpha <- Matrix(0, nrow=nrow(regulation.step), ncol=1, 
             dimnames = list(rownames(regulation.step), NULL))
      alpha[gene.from,]<-1
    }
    gene.from <- correct.genes(gene.from, regulation.flow)
    gene.to <- correct.genes(gene.to, regulation.flow)
    
    gene.to.ex <- setdiff( gene.to, gene.from )
    
    gene.from.n <- length( gene.from )
    gene.to.n <- length( gene.to.ex )
    
    if ( gene.from.n == 0 || gene.to.n == 0 )
        return( invisible() )
    
    if ( gene.from.n > 1 )
        regulation.flow.from <-  as.vector( (t(alpha)%*%regulation.flow) )
    else
        regulation.flow.from <-  alpha[gene.from,1]*regulation.flow[ gene.from,]
    
    names(regulation.flow.from)<- rownames(regulation.flow)
    

    if ( gene.to.n > 1 )
        regulation.flow.to <- rowSums( as.matrix( 
            regulation.flow[ , gene.to.ex ] ) )
    else
        regulation.flow.to <- regulation.flow[ , gene.to.ex ]
    
    regulation.flow.auto <- Matrix::diag( regulation.flow )
    
    regulation.intermediate <- abs( regulation.flow.from * regulation.flow.to / 
        regulation.flow.auto )
    
    gene.to.net <- gene.to.ex[ regulation.intermediate[ gene.to.ex ] != 0 ]
    gene.to.net.n <- length( gene.to.net )
    
    gene.from.to.net <- union( gene.from, gene.to.net )
    
    regulation.intermediate <- regulation.intermediate[ !( 
        regulation.intermediate == 0 | is.infinite( regulation.intermediate ) 
    ) ]
    
    if ( length( regulation.intermediate ) == 0 )
        return( invisible() )
    
    regulation.intermediate.gene <- names( sort( 
        regulation.intermediate, decreasing = TRUE ) )
    
    regulation.intermediate.gene <- regulation.intermediate.gene[ 
        ! regulation.intermediate.gene %in% gene.from.to.net ]
    
    if ( ! dir.exists( plot.dir ) )
        dir.create( plot.dir, recursive = TRUE )
    
    if ( ! is.null( plot.file.name.base ) )
        plot.pw.file.name.base <- plot.file.name.base
    else if ( gene.from.n == 1 && gene.to.n == 1 )
        plot.pw.file.name.base <- 
            sprintf( "pathway_%s_2_%s", correct.file.name(gene.from), correct.file.name(gene.to.ex) )
    else if ( gene.from.n == 1 )
        plot.pw.file.name.base <- sprintf( "pathway_%s_2_gene_set", correct.file.name(gene.from) )
    else if ( gene.to.n == 1 )
        plot.pw.file.name.base <- sprintf( "pathway_gene_set_2_%s", correct.file.name(gene.to.ex) )
    else
        plot.pw.file.name.base <- "pathway_gene_set_2_gene_set"
    
    if ( ! is.null( plot.title.base ) )
        plot.pw.title.base <- plot.title.base
    else if ( gene.from.n == 1 && gene.to.n == 1 )
        plot.pw.title.base <- 
            sprintf( "Pathway from %s to %s", correct.file.name(gene.from), correct.file.name(gene.to.ex) )
    else if ( gene.from.n == 1 )
        plot.pw.title.base <- sprintf( "Pathway from %s to Gene Set", correct.file.name(gene.from) )
    else if ( gene.to.n == 1 )
        plot.pw.title.base <- sprintf( "Pathway from Gene Set to %s", correct.file.name(gene.to.ex) )
    else
        plot.pw.title.base <- "Pathway between Gene Sets"
    
    plot.pw.title.base <- paste(strwrap(plot.pw.title.base, width = 30), collapse = "\n")
    
    regulation.rg <- sum( (t(alpha)%*%regulation.flow)[,gene.to.net] )
    gene.next.max <- min( pathway.gene.max, 
        length( regulation.intermediate.gene ) )
    
    path.gene <- NULL
    path.gene.current <- NULL
    rg.subgraph <- induced_subgraph( rg, gene.from.to.net )
    rg.pathway <- induced_subgraph( rg, path.gene.current )
    gene.next.idx <- 1
    pathway.next.idx <- 1
    regulation.fraction <- 0
    plotting.pathway <- FALSE
    
    while ( pathway.next.idx <= pathway.n.max && 
            abs( regulation.fraction - 1 ) >= pathway.regulation.approx )
    {
        while( TRUE )
        {
            if ( plotting.pathway )
            {
                plotting.pathway <- FALSE
            }
            else
            {
                path.gene.new <- erfa.get.path.gene.new( rg.subgraph, gene.from, 
                    gene.to.net, path.gene, path.gene.current )
                
                if ( ! identical( path.gene.new, path.gene.current ) )
                {
                    path.gene.current <- path.gene.new
                    
                    rg.pathway <- induced_subgraph( rg, 
                        c( gene.from.to.net, path.gene.current ) )
                    
                    regulation.rg.pathway <- sum( 
                      matrix(alpha[gene.from,],ncol=length(gene.from), nrow=1)%*%
                        erfa.calculate.regulation( rg.pathway, regulation.step, 
                            gene.from, gene.to.net ) 
                          
                    )
                    regulation.fraction.new <- regulation.rg.pathway / 
                        regulation.rg
                    
                    if ( regulation.fraction.new != regulation.fraction ) {
                        regulation.fraction <- regulation.fraction.new
                        break
                    }
                }
            }
            
            if ( gene.next.idx > gene.next.max )
                return( invisible() )
            
            path.gene <- c( path.gene, 
                regulation.intermediate.gene[ gene.next.idx ] )
            rg.subgraph <- induced_subgraph( rg, c( gene.from.to.net, 
                path.gene ) )
            gene.next.idx <- gene.next.idx + 1
        }
        
        plotting.pathway <- TRUE
        
        rg.pathway.vertex.size <- sapply( V( rg.pathway )$name, function( rgv ) 
        {
            vertex.size <- abs(regulation.flow.from)[ rgv ] /( 
                mean( abs(regulation.flow.from)[ gene.from ] )+1e-6)
            
            10 * vertex.size^0.15
        } )
        
        rg.pathway.vertex.color <- ifelse( 
            colMeans( as.matrix( regulation.flow[ gene.from, 
                V( rg.pathway )$name, drop = FALSE ] ) ) > 0, 
            "blue3", "red3" )
        
        rg.pathway.label.font <- ifelse( 
            V( rg.pathway )$name %in% gene.from.to.net, 2, 1 )
        
        rg.pathway.edge.width <- sapply( E( rg.pathway ), function( erp ) {
            edge.vertex <- ends( rg.pathway, erp )[ 1, ]
            edge.width <- abs( regulation.step[ edge.vertex[ 1 ], 
                edge.vertex[ 2 ] ] )
            10 * edge.width^0.5 
        } )
        names( rg.pathway.edge.width ) <- NULL
        
        rg.pathway.edge.color <- ifelse( E( rg.pathway )$sign == 1, "green4", 
            "orange2" )
        
        rg.pathway.edge.lty <- sapply( E( rg.pathway ), function( erp ) {
            edge.vertex <- ends( rg.pathway, erp )[ 1, ]
            edge.sign <- E( rg.pathway )$sign[ erp ]
            ifelse( 
                colMeans( as.matrix( regulation.flow[ gene.from, 
                    edge.vertex[ 1 ], drop = FALSE ] ) ) * 
                colMeans( as.matrix( regulation.flow[ gene.from, 
                    edge.vertex[ 2 ], drop = FALSE ] ) ) * 
                edge.sign > 0, 
                1, 2 )
        } )
        names( rg.pathway.edge.lty ) <- NULL
        
        rg.pathway.layout.edge.weight <- apply( 
            ends( rg.pathway, E( rg.pathway ) ), 1, function( ev ) 
                1 + abs( regulation.step[ ev[ 1 ], ev[ 2 ] ] )
        )
        
        edge.weight.boost.bool <- apply(
            ends( rg.pathway, E( rg.pathway ) ), 1, function( ev )
                ev[ 1 ] %in% gene.from | ev[ 2 ] %in% gene.to.ex
        )
        
        rg.pathway.layout.edge.weight[ edge.weight.boost.bool ] <-
            rg.pathway.layout.edge.weight[ edge.weight.boost.bool ] +
            mean( rg.pathway.layout.edge.weight )
        
        rg.pathway.layout <- layout_with_fr( rg.pathway, 
             weights = 1 / rg.pathway.layout.edge.weight )
        
        plot.pw.file.name <- sprintf( "%s_%02d.pdf", plot.pw.file.name.base, 
            pathway.next.idx )
        
        plot.pw.title <- sprintf( "%s (#%02d, %dv, %.1f%%)", 
            plot.pw.title.base, pathway.next.idx, length( V( rg.pathway ) ), 
            100 * regulation.fraction )
        
        if(plot.moa){
          pdf(file = file.path(plot.dir, plot.pw.file.name), width = 18, 
              height = 18)
          par(mar = c(1, 2, 6, 2))
          plot(rg.pathway, layout = rg.pathway.layout, vertex.size = rg.pathway.vertex.size, 
               vertex.color = rg.pathway.vertex.color, vertex.label.cex = 2, 
               vertex.label.font = rg.pathway.label.font, vertex.label.color = "black", 
               vertex.label.dist = 2, edge.width = rg.pathway.edge.width, 
               edge.color = rg.pathway.edge.color, edge.lty = rg.pathway.edge.lty, 
               edge.arrow.size = 2)
          title(plot.pw.title, cex.main = 2, font.main = 1)
          dev.off()
        }
        

        graph.pw.file.name <- sprintf( "%s/%s_%02d.RData",plot.dir, plot.pw.file.name.base, 
                                      pathway.next.idx )
        
        rg.pathway$fraction <-100 * regulation.fraction
        
        cat("Saving file...\r\n")
        
        save(rg.pathway, file = graph.pw.file.name)
        
        pathway.next.idx <- pathway.next.idx + 1
    }
    
    invisible()
}

correct.genes <- function (gene, reg.matrix){
  genes <- colnames(reg.matrix)
  
  gene <- sapply(gene, function(g) {
    if (g %in% genes ){
      return(g)
    } else {
      pattern <- paste("\\b", g, "\\b", sep = "")
      in_complex <- grepl(pattern, genes)
      
      if(sum(in_complex) == 0) stop ("One or more genes are not present in the knowledge graph")
      
      return(genes[in_complex][1])
    }
    
  })
  
  return(gene)
}

correct.file.name <- function(gene.name){
  if(grepl(":", gene.name)){
    return(paste(strsplit(gene.name, ":")[[1]][1], "complex", sep="_"))
  } else{
    return(gene.name)
  }
}



