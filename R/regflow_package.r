#' regflow: A package for carrying out Regulatory Flow Analysis
#'
#' The `regflow` package provides utilities to build regulation graphs, 
#' obtain regulation matrices, perform calculations about regulation flow 
#' and pathways, and generate figures with flow distributions and 
#' representations of pathways. 
#' 
#' @importFrom grDevices dev.off pdf 
#' 
#' @importFrom ggplot2 .data aes coord_fixed element_text geom_abline 
#'     geom_hline geom_point geom_vline ggplot ggsave labs lims theme 
#'     theme_classic
#' 
#' @importFrom ggrepel geom_text_repel
#' 
#' @importFrom graphics abline hist mtext par text title 
#' 
#' @importFrom igraph E E<- V add_edges as_adjacency_matrix components degree 
#' @importFrom igraph delete_vertices distances ends graph_from_data_frame 
#' @importFrom igraph induced_subgraph is_connected layout_with_fr 
#' 
#' @importFrom Matrix Diagonal Matrix solve t
#' 
#' @importFrom moments moment
#' 
#' @importFrom RSpectra eigs 
#' 
#' @importFrom stats density runif 
#' 
#' @importFrom utils read.csv read.delim write.csv 

"_PACKAGE"
