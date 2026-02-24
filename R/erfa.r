library(ggplot2)
#' Estimation of ERFA model from the observed expression data
#' 
#' ERFA method uses Projected Gradient Descent algorithm with backtracking 
#' to make predicted regulatory flow values closer to the observed log-fold changes. 
#' The method returns updated step and flow (one-step and all-step) matrices.
#' Optimization is based on minimizing distance between the log fold changes and 
#' the predicted regulatory flow.
#' 
#' @param rg An [igraph] graph object with the interaction graph, usually 
#'     obtained [build.regulation.graph()]. 
#' @param regulation.step Named numeric matrix with one-step regulation values, 
#'     usually obtained with [calculate.regulation.matrices()]. 
#' @param stimulationTrain Named numeric matrix, defining stimulation values for one 
#'     or more genes over multiple experiments in the training data: genes in rows, experiments in column. 
#'     A value in the matrix should be equal to 0 if the gene is not stimulated, and non-zero
#'     if the gene is stimulated.
#' @param stimulationTest Named numeric matrix, defining stimulation values for one 
#'     or more genes over multiple experiments in the test data: genes in rows, experiments in column. 
#'     A value in the matrix should be equal to 0 if the gene is not stimulated, and non-zero
#'     if the gene is stimulated.     
#' @param expressionTrain Named numeric matrix, providing log-foldchanges of 
#'     expression values for training data:
#'     genes in rows, experiments in columns. 
#' @param expressionTest Named numeric matrix, providing log-foldchanges of 
#'     expression values for test data:
#'     genes in rows, experiments in columns.  
#' @param fractionValid Should training data be split into training and validation? 
#' The number specifies a proportion of split.
#' @param stepsize Numeric value that specifies the step size in the 
#'     Gradient Descent algorithm
#' @param n.iteration Integer with the number of optimization iterations to be 
#'     run. 
#' @param plot.dir Character string with the directory to save plots to. If 
#'     `NULL`, no plotting is carried out. 
#' @param plot.title Character string with the title for all the plots. 
#' 
#' @param plotIDtrain index of perturbations in the training data
#' that shall be used
#' to plot dependence of log-foldchanges versus predicted regulatory flow. 
#' 
#' @param plotIDtest index of perturbations in the test data
#' that shall be used
#' to plot dependence of log-foldchanges versus predicted regulatory flow. 
#' 
#' @param plot.title Character string with the title for all the plots. 
#' @param no_cores Number of available threads for parallellization
#' 
#' 
#' @return A list with two named numeric matrices: 
#' * `step`: Matrix with one-step regulation, refined. 
#' * `flow`: Matrix with all-step regulation, refined. 
#' 
#' @importFrom magrittr `%>%`
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

erfa <- function( rg, regulation.step, stimulationTrain, 
                                           expressionTrain, stepsize=0.1, n.iteration=20, plot.dir=NULL, 
                                           plot.title="Plot", stimulationTest=NULL, expressionTest=NULL, fractionValid=0.5, regularisation=0, 
                                           plotIDtrain=NULL, plotIDtest=NULL, no_cores=8 )
{
  library(foreach)
  library(doParallel)

  cl <- makeCluster(no_cores,outfile="")
  registerDoParallel(cl)
  
  reg.step <- regulation.step
  treg.step <- t(reg.step)
  
  nTr <- ncol(expressionTrain)
  nTe <- 0
  if(!is.null(stimulationTest)){
    expression.profile <- cbind(expressionTrain, expressionTest)
    stimulation <- cbind(stimulationTrain, stimulationTest)
    nTe <- ncol(stimulationTest)
  }
  
  ##  verify argument invariants
  
  stopifnot( nrow( reg.step ) == ncol( reg.step ) && 
               nrow( reg.step ) == length( V( rg ) ) )
  
  stopifnot( rownames( reg.step ) == colnames( reg.step ) & 
               rownames( reg.step ) == sort( V( rg )$names ) )
  
  graph.gene <- rownames( reg.step )
  graph.n <- length( graph.gene )
  
  ##  change max eigenvalue
  
  reg.step.lm.eigenval <- eigs( treg.step, 1 )$values
  
  reg.step <- 0.99 / Mod( reg.step.lm.eigenval ) * reg.step
  
  
  ##  identify genes to fit
  
  profile.gene <- rownames( expression.profile )
  
  stimulationNames <- rownames(stimulation)
  stimulation.genes <- apply(stimulation, MARGIN = 2, FUN=function(x) stimulationNames[x!=0], simplify = F)
  stopifnot( unlist(stimulation.genes) %in% graph.gene )
  
  graph.reg.edge <- ends( rg, E(rg))
  graph.reg.factor <- graph.reg.edge[ , 1 ]
  graph.reg.gene <- graph.reg.edge[ , 2 ]

  stimulation.Matrix <- Matrix( 0, graph.n, ncol(stimulation) )
  rownames( stimulation.Matrix ) <- graph.gene
  stimulatedGenes <- unique(unlist(stimulation.genes))
  stimulation.Matrix[ stimulatedGenes, ] <- stimulation[stimulatedGenes,]
  

  fit.genes <- lapply( stimulation.genes, FUN = function(stimulation.gene) 
    intersect(  profile.gene, graph.gene ))
  

  x <- stimulation.Matrix
  response.vector <- list()
  yp <- list()
  where_expressed <- list()
  where_not_expressed <- list()
  where_perturbed <- list()
  W <- list()
  Dm <- list()

  Sn <- treg.step

  alpha <- list()
  obj0 <- list()
  objv <- list()
  objT0 <- list()
  objV0 <- list()
  objT <- list()
  objV <- list()
  Mask <- (Sn!=0)
  

  Omatrix <- Matrix( 0,  graph.n, 1 )
  Gmatrix <- Matrix(data=0, nrow=nrow(Sn), ncol=ncol(Sn))

  for (i in 1:ncol(expression.profile)){

    response.vector[[i]] <- Omatrix
    rownames( response.vector[[i]] ) <- graph.gene  
    response.vector[[i]][fit.genes[[i]],] <- expression.profile[ fit.genes[[i]], i ]
    
    
    yp[[i]] <- response.vector[[i]]
    
  
    where_not_expressed[[i]]  <-  (response.vector[[i]]==0)
    where_expressed[[i]] <- as.vector(response.vector[[i]]!=0)
    where_perturbed[[i]] <- as.vector(stimulation.Matrix[,i]!=0)
    
    W[[i]] <- Diagonal(graph.n, x=where_expressed[[i]])
    
  
    Dm[[i]] <- Matrix::Diagonal(graph.n, as.vector(stimulation.Matrix[,i]))
    De <- Matrix::Diagonal(graph.n, as.vector(yp[[i]]))
  

    Z <- Matrix::solve(Diagonal( graph.n )-Sn, Dm[[i]])
  
    Zprime <- Z[where_expressed[[i]], where_perturbed[[i]], drop=F]
    
    ZZ <- t(Zprime)%*%Zprime
    
    if(Matrix::rcond(ZZ)>1e-15) {
      regul=0
    }  else {
      regul= Matrix::Diagonal(ncol(Zprime), 1e-6)
    } 
    
    alphac <- Matrix::solve(ZZ+ regul, t(Zprime)%*%yp[[i]][where_expressed[[i]]])
    
    
    alpha[[i]] <- Matrix( 0,  graph.n, 1 )
    alpha[[i]][where_perturbed[[i]],1] <- alphac
    

    if(i<=nTr) {
      obj0[[i]] <- Matrix::mean(t(yp[[i]]-Z%*%alpha[[i]])%*%W[[i]]%*%(yp[[i]]-Z%*%alpha[[i]]))
    } else if (i<= floor(nTr+fractionValid*nTe)){
      objV0[[i-nTr]] <- Matrix::mean(t(yp[[i]]-Z%*%alpha[[i]])%*%W[[i]]%*%(yp[[i]]-Z%*%alpha[[i]]))
    } else{
      objT0[[i-floor(nTr+fractionValid*nTe)]] <- Matrix::mean(t(yp[[i]]-Z%*%alpha[[i]])%*%W[[i]]%*%(yp[[i]]-Z%*%alpha[[i]]))
    }
    
  }
  
  obj_value0 <- mean(unlist(obj0))
  obj_valueT0 <- mean(unlist(objT0))
  obj_valueV0 <- mean(unlist(objV0))

  cat("Initial objective:", obj_value0, "\r\n")
  if(fractionValid>0) cat("Initial objective Valid:", obj_valueV0, "\r\n")
  cat("Initial objective Test:", obj_valueT0, "\r\n")
  TrainE <- c(obj_value0)
  ValidE <- c(obj_valueV0)
  TestE <- c(obj_valueT0)
  
  
  beta <- stepsize/4

  
  
  for ( j in 1 : n.iteration )
  {

    Sn0 <- Sn
    alpha0 <- alpha
    Dsn <- Diagonal( graph.n )-Sn

    Grad <- Gmatrix

    res <-
    foreach(i=seq(1,nTr),.combine = "+",.packages=c("foreach", "Matrix")) %dopar%{
      G1 <- Grad
      zi <- Matrix::solve(Dsn,alpha[[i]]*x[,i])
      gradI  <-  Matrix::solve(t(Dsn),W[[i]]%*%(yp[[i]]-zi))%*%t(zi)#t(R)%*%W[[i]]%*%(yp[[i]]-R%*%(alpha[[i]]*x[,i]))%*%t(alpha[[i]]*x[,i])%*%t(R)
      G1[Mask] <- G1[Mask]-gradI[Mask]
      G1
    }
    Grad[Mask] <- res[Mask]

    Grad <- (1/nTr*Grad)+regularisation*Sn
    Grad <- Grad/Matrix::norm(Grad, type="F")
    
    beta <- beta*4
    obj_value <- Inf


    while((beta >1e-6) && (obj_value>obj_value0)) {
      Sn1 <- Sn-beta*Grad
      Sn1 <- Sn1*(sign(Sn1) == sign(treg.step))


      for (i in 1:ncol(expression.profile)){
        Z <- Matrix::solve(Diagonal( graph.n )-Sn1, Dm[[i]]) 
        Zprime <- Z[where_expressed[[i]], where_perturbed[[i]], drop=F]
        ZZ <- t(Zprime)%*%Zprime
        
        if(Matrix::rcond(ZZ)>1e-15) {
          regul=0
        }  else {
          regul= Matrix::Diagonal(ncol(Zprime), 1e-6)
        } 
        
        alphac <- Matrix::solve(ZZ+regul, t(Zprime)%*%yp[[i]][where_expressed[[i]]])
        alpha[[i]] <- Omatrix
        alpha[[i]][where_perturbed[[i]],1] <- alphac

        if(i<=nTr) {
          objv[[i]] <- Matrix::mean(t(yp[[i]]-Z%*%alpha[[i]])%*%W[[i]]%*%(yp[[i]]-Z%*%alpha[[i]]))
        } else if (i<= floor(nTr+fractionValid*nTe)){
          objV[[i-nTr]] <- Matrix::mean(t(yp[[i]]-Z%*%alpha[[i]])%*%W[[i]]%*%(yp[[i]]-Z%*%alpha[[i]]))
        } else{
          objT[[i-floor(nTr+fractionValid*nTe)]] <- Matrix::mean(t(yp[[i]]-Z%*%alpha[[i]])%*%W[[i]]%*%(yp[[i]]-Z%*%alpha[[i]]))
        }
      }
      obj_value <- mean(unlist(objv))
      obj_valueV <- mean(unlist(objV))
      obj_valueT <- mean(unlist(objT))

      beta <- beta/2
      cat("Backtracking:", obj_value, "\r\n")
      
    }
    
    Sn <- Sn1

    if (obj_value>=obj_value0) break;
    
    obj_value0 <- obj_value
    cat("Objective function value iteration", j, ": ", obj_value, "\r\n")
    if(fractionValid>0)  cat("Objective function value Valid iteration", j, ": ", obj_valueV, "\r\n")
    cat("Objective function value Test iteration", j, ": ", obj_valueT, "\r\n")
    TrainE <- c(TrainE, obj_value)
    ValidE <- c(ValidE, obj_valueV)
    TestE <- c(TestE, obj_valueT)

  }
  reg.step <- t(Sn0)
  reg.flow <- Matrix::solve(Diagonal( graph.n )-reg.step)
  
  plot_i <- c()
  if(!is.null(plotIDtrain)){
    plot_i <- c(plot_i, plotIDtrain)
  }
  
  if(!is.null(plotIDtest)){
    plot_i <- c(plot_i, nTr+plotIDtest)
  }
  
  dir.create(plot.dir, showWarnings = FALSE)
  
  for (i in plot_i){

    PredResponseOriginal <- Matrix::solve(Diagonal( graph.n )-treg.step, x[,i, drop=F])[where_expressed[[i]]]
    PredResponseAdjusted <- Matrix::solve(Diagonal( graph.n )-Sn0, alpha0[[i]]*x[,i, drop=F])[where_expressed[[i]]]

    df <- data.frame(DiffExpression=response.vector[[i]][where_expressed[[i]]],
                  PredResponseOriginal=PredResponseOriginal, PredResponseAdjusted=PredResponseAdjusted,
                  GeneNames=rownames(response.vector[[i]])[where_expressed[[i]]])

    p0  <-  df%>%plotly::plot_ly(x=~DiffExpression, y=~PredResponseOriginal, hovertext=~GeneNames)%>%plotly::add_markers()
    p1  <-  df%>%plotly::plot_ly(x=~DiffExpression, y=~PredResponseAdjusted, hovertext=~GeneNames)%>%plotly::add_markers()

    if(! is.null( plot.dir )){
      htmlwidgets::saveWidget(plotly::as_widget(p0), paste(plot.dir, "/", plot.title,"_", i, "_", "BeforeAdjustment.html", sep="" ))
      htmlwidgets::saveWidget(plotly::as_widget(p1), paste(plot.dir, "/", plot.title,"_", i, "_", "AfterAdjustment.html", sep=""))
    }

  }
  
  df <- data.frame(TrainE=TrainE, ValidE=ValidE, TestE=TestE)
  df$iteration <- 1:nrow(df)
  df1 <-  df%>%tidyr::pivot_longer(TrainE:TestE, names_to="Type", values_to = "Error")
  itp <- ggplot(df1, aes(x=iteration, y=Error, color=Type))+geom_line()

  if(! is.null( plot.dir )){
    ggsave(paste(plot.dir, "/", "convergence_plot.pdf", sep=""), plot = itp, width = 8, height = 5)
  }
  
  stopCluster(cl)
  
  list( step = reg.step, flow = reg.flow , alpha=alpha0, errors=df1)
}
