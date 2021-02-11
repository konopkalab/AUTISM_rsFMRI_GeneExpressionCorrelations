FastCor <- function(exp,sme,method="pearson",alternative="two.sided",cores=1,override=FALSE,...){
  suppressPackageStartupMessages({library(parallel)})
  mat = as.matrix(exp)
  nr = nrow(mat)
  if(nr > 1000){
    if(override){
      show("Wait for it!")
    } else {
      stop("Set override to TRUE")
    }
  }
  if(is.null(rownames(mat))){
    nm = as.character(1:nr)
  } else { 
    nm = rownames(mat)
  }
  corrAndP = mclapply(1:nr,function(i){
      res = cor.test(mat[i,],as.numeric(sme),method=method,alternative=alternative)
      cbind(res$estimate,res$p.value)
    },mc.cores=cores,...)
  Rho = unlist(sapply(corrAndP,function(i){i[,1]}))
  Pval  = unlist(sapply(corrAndP,function(i){i[,2]}))
  results = cbind(Rho,Pval)
  rownames(results) <- rownames(exp)
  return(results)
}

# counts = gene x sample matrix
# meta = sample x covariate matrix
# threshold = number of PCA to consider (e.g. 5)
# inter = interaction between covariates (e.g. Age:Sex)
VarExp <- function(counts, meta, threshold, inter){
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(optimx))
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)

  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}

  pred.list <- colnames(meta)
  meta <- droplevels(meta)

  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}

  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit,control = lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}

# Bar plot for data visualization
plotVarExp <- function(pvca.res, title){
  suppressPackageStartupMessages(library(ggplot2))
  plot.dat <- data.frame(eff=names(pvca.res), prop=pvca.res)
  p <- ggplot2::ggplot(plot.dat, aes(x=eff, y=prop))
  p <- p + ggplot2::ggtitle(title)
  p <- p + ggplot2::geom_bar(stat="identity", fill="steelblue", colour="steelblue")
  p <- p + ggplot2::geom_text(aes(label=round(prop,3), y=prop+0.04), size=4)
  p <- p + ggplot2::scale_x_discrete(limits=names(pvca.res))
  p <- p + ggplot2::scale_y_continuous(limits = c(0,1))
  p <- p + ggplot2::labs(x= "Effects", y= "WAPV")
  p <- p + ggplot2::theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
}


# Clean genomic data function
cleaningP = function(y, mod, svaobj,  P=ncol(mod)) {
        X=cbind(mod,svaobj$sv)
        Hat=solve(t(X)%*%X)%*%t(X)
        beta=(Hat%*%t(y))
        cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
        return(cleany)
}

# HCP function
hcp <- function(F,Y,k,lambda,lambda2,lambda3,iter=100,tol = 1e-6){
      #### Standardize/Mean center inputs ####
      Yn = scale(Y)
      Fn = scale(F)
      #### Set default values for outputs and initialise B ####
      U = matrix(0,dim(F)[2],k)
      Z = matrix(0,dim(F)[1],k)
      B = matrix(runif(k*dim(Y)[2]),k,dim(Y)[2])
      #### Error checking ####
      if (dim(F)[1] != dim(Y)[1])
        stop('number of rows in F and Y must agree')
  
      if (k<1 | lambda<tol | lambda2<tol | lambda3<tol)
        stop('lambda, lambda2, lambda3 must be positive and/or k must be an integer');
      #### Coefficients esitmation ####
      o = matrix(0,1,iter)
      for (ii in 1:iter){
        o[ii] = norm(Yn-Z %*% B,type='E') + norm(Z-Fn %*% U, type='E')*lambda + norm(B,type='E')*lambda2 + norm(U,type='E')*lambda3
        Z = ((Yn %*% t(B)) + lambda * (Fn %*% U)) %*% solve((B %*% t(B)) + lambda*diag(dim(B)[1]))
      B = solve(((t(Z) %*% Z) + lambda2 * diag(dim(Z)[2])),(t(Z) %*% Yn))
        U = solve(((t(Fn) %*% Fn) * lambda + lambda3 * diag(dim(U)[1])),(lambda* t(Fn) %*% Z))
      if (ii > 1)
            if ((abs(o[ii]-o[ii-1])/o[ii]) < tol)
        break 
      }
      if (ii>=iter)
        writeLines("\nPre-mature convergence: Consider increasing the number of iterations")
      #### Error calculation ####
      error = (norm((Yn-Z %*% B),type='E') /norm(Yn,type='E')) + (norm((Z - Fn %*% U),type='E')/norm(Fn%*%U,type='E'))
      error1 = norm(Yn-Z%*%B,type='E')/norm(Yn,type='E')
      error2 = norm(Z-Fn%*%U,type='E')/norm(Fn%*%U,type='E')
      #### Delta change calculation ####
      dz = Z %*% (B %*% t(B) + lambda * diag(dim(B)[1])) - (Yn %*% t(B) + lambda * Fn %*% U)
      db = (t(Z) %*% Z + lambda2 * diag(dim(Z)[2])) %*% B - t(Z) %*% Yn
      du = (t(Fn) %*% Fn *lambda + lambda3 * diag(dim(U)[1])) %*% U - lambda * (t(Fn) %*% Z)
      #### Residual calculation and covariate adjustment ####
      R = (Yn - Z %*% B) * t(replicate(dim(Y)[1],attr(Yn,'scaled:scale'))) + t(replicate(dim(Y)[1],attr(Yn,'scaled:center')))
      #### Final Outputs ####
      return(list(R=R,Z=Z,B=B,U=U,o=o,error=error,error1=error1,error2=error2,dz=dz,db=db,du=du))
}


ppiNet <- function(molecularIDs = NULL,file = NULL,speciesID = 9606,Score = 100, evidence = c("neighborhood","neighborhood_transferred",
            "fusion","cooccurence","homology","coexpression","coexpression_transferred",
            "experiments","experiments_transferred","database","database_transferred","textmining",
            "textmining_transferred","combined_score")){
  
  # Detecting the input type
  if(is.null(file) && molecularIDs > 0){
    # Creating the vector to store the unique identifiers
    for_gen <- unlist(sapply(molecularIDs,function(i){
      if(grepl("-",i) > 0){
        return(unlist(strsplit(i,"-")))
      }else{
        return(i)
      }
    }))
    
    # Remove duplicated identifiers
    for_gen <- unique(sort(for_gen))
    # Transform the vector into data frame 
    new_genes <- as.data.frame(for_gen,stringsAsFactors = FALSE)
    # The column name must be gene
    names(new_genes) <- "gene"
    # Loading the STRING database
    database <- STRINGdb$new(version="10",species=speciesID,score_threshold=Score,input_directory="")
    # Obtaining the STRING ID to each identifier
    mapped <- database$map(new_genes,"gene",removeUnmappedRows = TRUE)
    # Removing the STRING ID duplicated
    mapped <- mapped[!duplicated(mapped$STRING_id),]
    # Obtaining the interactions among STRING IDs according to different types of evidence
    interactions <- database$get_interactions(mapped$STRING_id)
    # Extract the relations deleting the evidence information columns
    graph_relations <- data.frame(interactions$from,interactions$to,stringsAsFactors = FALSE)
    # Read each evidence given
    for(i in evidence){
      # Take each column with evidence information
      for(j in names(interactions)){
        # If the evidence in the "interactions" variable corresponds with one of the request evidence
        if(i == j){
          # Appends the column with the requested evidence in the data frame with interactions among STRING IDs
          graph_relations <- cbind(graph_relations,interactions[j])
        }
      }
    }
    # This data frame will be fill up with interactions with any evidence value greater than zero
    graph_ppi <- graph_relations[rowSums(graph_relations[,seq(3,ncol(graph_relations))]) > 0,]
    # This loop replace the STRING IDs with the original identifiers in the first column 
    for(n in seq_len(nrow(graph_ppi))){
      graph_ppi$interactions.from[n] <- mapped[
        graph_ppi$interactions.from[n] == mapped$STRING_id,][1]
    }
    # This loop replace the STRING IDs with the original identifiers in the second column
    for(n in seq_len(nrow(graph_ppi))){
      graph_ppi$interactions.to[n] <- mapped[
        graph_ppi$interactions.to[n] == mapped$STRING_id,][1]
    }
    
    # Is necessary transform the data frame into matrix
    edge_list <- as.matrix(as.vector(graph_ppi[,1],mode = "character"))
    edge_list <- cbind(edge_list,as.vector(graph_ppi[,2],mode = "character"))
    # Create a network based on the interactions in the matrix
    final_graph <- graph.edgelist(edge_list,directed = FALSE)
  }else if(is.null(molecularIDs) && file.exists(file)){
    # Read and create the igraph object
    final_graph <- read.graph(file = file,format = "ncol")
  }else{
    stop("A valid input was not found: Please, be sure to use an edge list or a vector of IDs as input")
  }
  # Return the PPI network
  return(final_graph)
}

