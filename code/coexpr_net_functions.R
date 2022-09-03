################################################################################
# Helper functions redes_coexpr_main
################################################################################



loadProteins <- function(){
  
  #########
  # VACV #
  #########
  prots_VACV<- read.table(file = "data_VACV_edgeR_2reps.txt", 
                          header = TRUE,sep="\t",dec=".")
  
  # Clean protein names and remove duplicates
  prot_clean_VACV <- cleanProtNames(prots_VACV)
  
  # keep only average expression time series data from infected cells, 
  # FDR and LogFC
  df_prot_filtered_VACV <- prot_clean_VACV[,c(1,2,3,11:17,20,23,26,29,32,35,18,21,24,27,30,33)]
  n_points_VACV <- 7
  
  #########
  # EBV   #
  #########
  prots_EBV <- read.table(file = "data_EBV_edgeRCtrl.txt", 
                          header = TRUE,sep="\t",dec=".")
  
  # Clean protein names and remove duplicates
  prot_clean_EBV <- cleanProtNames(prots_EBV)
  
  # keep only average expression time series data from infected cells, 
  # FDR and LogFC
  df_prot_filtered_EBV <- prot_clean_EBV[,c(1,2,3,9:13,16,19,22,25,14,17,20,23)]
  n_points_EBV <- 5
  
  #########
  # HSV1  #
  #########
  prots_HSV1 <- read.table(file = "data_HSV1_edgeR.txt", 
                           header = TRUE,sep="\t",dec=".")
  
  # Clean protein names and remove duplicates
  prot_clean_HSV1 <- cleanProtNames(prots_HSV1)
  
  # keep only average expression time series data from infected cells
  df_prot_filtered_HSV1 <- prot_clean_HSV1[,c(1,2,3,11:17,20,23,26,29,32,35,18,21,24,27,30,33)]
  n_points_HSV1<- 7
  
  
  save(df_prot_filtered_VACV, n_points_VACV, df_prot_filtered_EBV, 
       n_points_EBV, df_prot_filtered_HSV1, 
       n_points_HSV1, file = "allProteins.RData")
  
  
}


# Clean protein names and remove duplicates
# prots: data frame with a column named UniProtID and a column named GeneID
cleanProtNames <- function(protdf){
  # Clean protein names
  geneIDs <- protdf$GeneID
  cleannames <- gsub(";.*","",geneIDs)
  geneIDsc <- gsub("^$|^$", NA, cleannames)
  protdf$GeneID <- as.factor(geneIDsc)
  protIDs <- protdf$UniProtID
  cleannames <- gsub(";.*","",protIDs)
  protIDsc <- gsub("^$|^$", NA, cleannames)
  protdf$UniProtID <- as.factor(protIDsc)
  # For proteins with several IDs, keep just the first one
  protdf <- protdf[!duplicated(protdf[,1]),]
  return(protdf)
}
 


# ++++++++++++++++++++++++++++
# flattenCorrMatrix formats correlation matrix into a table with four columns: 
# Var1, Var2, rho, pval
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  
  # correlation values are passed in cormat (symmetric matrix)
  # p-values are passed in pmat. Lower triangle: Unadjusted p-values. 
  # Upper triangle: adjusted p-values
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


# Pruning adjacency matrix based on interaction values 
# transforms as a list of non redundant correlations
flattenCorrMatrix2 <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}


# Proportionality
#2*cov(logx,logy)/(var(logx)+var(logy))
concordanceCoorCoef<- function(prot1,prot2){
  x <- log10(prot1)
  y <- log10(prot2)
  cc <- 2*cov(x,y)/(var(x)+var(y))
  
}

# net is an igraph object
networkAnalysis<-function(net){

  # Centrality
  print("---------------Hubs-----------------")
  print("Degree centrality")
  print(head(sort(degree(g),decreasing = TRUE),n=10))
  print("Betweenness centrality")
  print(head(sort(betweenness(net, directed = FALSE, normalized = TRUE),
                  decreasing = TRUE),n=10))
  print("Eigenvector centrality")
  print(head(sort(eigen_centrality(net)$vector,decreasing = TRUE),n=10))
  print("Closeness centrality")
  print(head(sort(closeness(net,normalized = TRUE),decreasing = TRUE),n=10))
  
  # Distances
  cat("Edge density: ",edge_density(net),"\n")
  cat("Diameter: ",  diameter(net),"\n")
  cat("Average shortest path: ", mean_distance(g, directed = FALSE),"\n")
  
  # Other
  cat("Assortativity Degree: ", assortativity_degree(g),"\n")
  cat("Assortativity Origin: ", 
      assortativity_nominal(g, as.factor(V(g)$Origin),directed = FALSE),"\n")
  cat("Smallworldness: ", smallworldness(g)[1],"\n")
  cat("Global Clustering coefficient: ", transitivity(g,"globalundirected"),"\n")
  cat("Average per-node clustering coefficient: ", 
      transitivity(g,"localaverageundirected"),"\n")
  
}


# net is an igraph object
subnetworkAnalysis<-function(g,host_linked_virus,virus_nodes,host_nodes){
  deg <- degree(g, mode="all")
  cat('Average degree network',mean(deg),"\n")
  cat('Average degree virus proteins',mean(deg[V(g)$Origin=='Virus']),"\n")
  cat('Average degree host proteins',mean(deg[V(g)$Origin=='Host']),"\n")
  cat('Average degree host proteins linked to virus',
      mean(degree(g,host_linked_virus)),"\n")
  
  cat("Average betweenness network",
      mean(betweenness(g, directed = FALSE, normalized = TRUE)),"\n")
  cat("Average betweenness virus proteins",
      mean(betweenness(g,virus_nodes, directed = FALSE, normalized = TRUE)),"\n")
  cat("Average betweenness host proteins",
      mean(betweenness(g,host_nodes, directed = FALSE, normalized = TRUE)),"\n")
  cat("Average betweenness host proteins linked to virus",
      mean(betweenness(g,host_linked_virus, directed = FALSE, normalized = TRUE)),"\n")
  
  
  cat("Average closeness network",mean(closeness(g)),"\n")
  cat("Average closeness virus proteins",mean(closeness(g,virus_nodes)),"\n")
  cat("Average closeness host proteins",mean(closeness(g,host_nodes)),"\n")
  cat("Average closeness host proteins linked to virus",
      mean(closeness(g,host_linked_virus)),"\n")
  
  # Eigenvector centrality scores correspond to the values of the 
  # first eigenvector of the graph adjacency matrix
  e <- eigen_centrality(g)$vector
  cat("Average eigen_centrality network",mean(e),"\n")
  cat("Average eigen_centrality virus proteins",
      mean(e[names(e) %in% virus_nodes]),"\n")
  cat("Average eigen_centrality host proteins",
      mean(e[names(e) %in% host_nodes]),"\n")
  cat("Average eigen_centrality host proteins linked to virus",
      mean(e[names(e) %in% host_linked_virus]),"\n")
  
}



# g is an igraph network object
networkVisualization <- function(g){
  # Generate colors based on cell of origin:
  colrs <- c("#365C8DFF","#FDE725FF")
  colrs <- adjustcolor(colrs,alpha.f=0.8)
  V(g)$color <- colrs[as.integer(as.factor(V(g)$Origin))]
  shapes <- c("circle", "square")
  V(g)$shape  <- shapes[as.integer(as.factor(V(g)$Origin))]
  # Compute node degrees (#links) and use that to set node size:
  deg <- igraph::degree(g)
  V(g)$size <- 3*log(deg)
  # change edge color
  E(g)$edge.color <- "gray80"
  
  plot(g,layout=layout_nicely(g, dim=2),vertex.label=NA)
  
}


# clusters is an igraph communities object
# g is an igraph network object
clustersVisualization <- function(clusters,g){

  deg <- igraph::degree(g)
  V(g)$size <- 3*log(deg)
  shapes <- c("circle", "square")
  V(g)$shape  <- shapes[as.integer(as.factor(V(g)$Origin))]
  
  V(g)$community <- final_clusters$membership
  colrs <- viridis(length(final_clusters))
  colrs <- adjustcolor( colrs, alpha=.8) 
  
  plot(clusters,g,mark.groups=NULL,vertex.label=NA,col=colrs[V(g)$community],
       edge.color = c("grey", "sienna1")[igraph::crossing(final_clusters, g) + 1])
  
}


# net is an igraph object
communities <- function(net){
  ceb <- cluster_edge_betweenness(net)
  plot(ceb, net,vertex.label=NA)
  cat("Modularity: ", modularity(ceb),"\n")
  cat("Number of communities: ", length(ceb),"\n")
  # High modularity for a partitioning reflects dense connections within 
  # communities and sparse connections across communities.
  mem <- membership(ceb)
  for(i in 1:length(ceb)) {
    community <- names(mem[mem[]==i])
    if (length(community)>5){
      cat("Members of community ",i,": ",community,"\n")
    }
  }
}

# sf is an igraph object
pwLawFit <- function(sf) {
  deg <- degree(sf)
  deg.dist <- degree_distribution(sf, cumulative=F, mode="all")
  # frequency distribution (remove first element, which is relative frequency 
  # of zero degree)
  deg.distin <- length(deg)*deg.dist[-1] 
  xdist=1:max(deg)
  ydist=deg.distin
  # robust linear regression to a line
  powreg <- robustbase::lmrob(log(ydist[ydist > 0]) ~ log(xdist[ydist > 0])) 
  gof <- summary(powreg)$adj.r.squared 
  return(gof)
}



# g is an igraph object
communityDetection <- function(g){
  fast <- fastgreedy.community(g)
  edge <- edge.betweenness.community(g)
  par(mfrow=c(1,2))
  plot(fast, g, main="Fast Greedy")
  plot(edge, g, main="Edge Betweeness")
  
  x = fastgreedy.community(g)
  i <- membership(x)
  
  g<- set_vertex_attr(g, 'color', value= c('yellow', 'cyan', 'red', 'pink')[i])
  
  graphjs(g)
}

# corr_matrix is a matrix with three columns: Prot1,Prot2 and Correlation
pwLawFitDraw <- function(corr_matrix,show_data=FALSE){
  sf <- graph_from_data_frame(corr_matrix,directed=FALSE)
  deg <- degree(sf)
  deg.dist <- degree_distribution(sf, cumulative=F, mode="all")
  # frequency distribution (remove first element, which is relative frequency 
  #of zero degree)
  deg.distin <- length(deg)*deg.dist[-1] 
  deg.distin <- deg.dist[-1]
  xdist=1:max(deg)
  ydist=deg.distin
  # robust linear regression to a line
  powreg <- robustbase::lmrob(log(ydist[ydist > 0]) ~ log(xdist[ydist > 0])) 
  interc <- powreg$coefficients[1]
  slope <- powreg$coefficients[2]
  gof <- summary(powreg)$adj.r.squared  
  xfit=seq(2, max(deg), length.out = 20)
  yfit <- exp(interc)*(xfit^slope)
  #lines( log(xfit), log(yfit), type='l', col="black", lwd=3)
  if (show_data=="TRUE"){
    cat('Exponent of power law',slope,"\n")
    cat('Goodness-of-fit(adj. R2)',gof,"\n")
  }
  data <- data.frame(logk_i=log(xdist),logPk_i=log(ydist))
  ggplot(data, aes(x=logk_i, y=logPk_i)) +
    geom_point( color="#482576") +
    geom_smooth(method=lmrob , color="#1f968b", fill="#3dbc74", se=TRUE) +
    theme_ipsum()

}



cutOffScan <- function(cutoff,corsig){
  
  names_props <- c("cutoff", "nodes","links","gof", 
                   "average_local_transitivity", "modularity")
  props <- data.frame(matrix(ncol = length(names_props), nrow = 0))
  colnames(props) <- names_props
  
  for(c in cutoff){
    
    corf <- filter(corsig, abs(correlation) > c)
    
    # Get maximal component
    gm <- graph_from_data_frame(corf)
    components <- igraph::clusters(gm, mode="weak")
    biggest_cluster_id <- which.max(components$csize)
    vert_ids <- V(gm)[components$membership == biggest_cluster_id]
    gc <- igraph::induced_subgraph(gm, vert_ids)
    
    # Calculate properties
    gof <- pwLawFit(gc)
    links <- gsize(gc)
    nodes <- length(V(gc))
    trans <- transitivity(gc, type = "average")
    
    gcu <- as.undirected(gc,mode="collapse")
      
    clusters = cluster_louvain(gcu)
    modularity <- modularity(clusters)
    
    props_iter <- as.data.frame(t(c(c,nodes,links,gof,trans,modularity)))
    colnames(props_iter) <- names_props
    props <- rbind(props,props_iter)
    
  }

  return(props)
  
}


calcCorrelation <- function(corr_method, dataExprMat, adjpval=0.1, lambda=0){
  
  if (corr_method=="pearson"){ 
    
    corprot <- corr.test(dataExprMat,adjust = 'BH', ci=F)
    
    cortab <- flattenCorrMatrix(corprot$r, corprot$p)
    colnames(cortab) <- c('Prot1','Prot2','correlation','pval.adj')
    
    
  } else if (corr_method=="proportionality"){
    protmat <- t(dataExprMat)
    corrprop <- matrix(nrow = dim(protmat)[1], ncol = dim(protmat)[1])
    for(i in 1:dim(protmat)[1]){
      corrprop[i,] <- apply(protmat,1,concordanceCoorCoef,prot2=protmat[i,])
    }
    rownames(corrprop) <- rownames(protmat)
    colnames(corrprop) <- rownames(protmat)
    corrprop[lower.tri(corrprop)] <- corrprop[upper.tri(corrprop)]
    cortab <- flattenCorrMatrix2(corrprop)
    colnames(cortab) <- c('Prot1','Prot2','correlation')
    
  } else if (corr_method=="GGM"){
    
    # penalty parameter in Lasso regression, larger lambda results in faster 
    # calculation and a sparser graph, 
    # if don't set then use default sqrt(2*log(p/sqrt(n))/n)
    # a n*p matrix with n samples and p variables.
    if (lambda==0){
      p <- dim(dataExprMat)[2] # number of genes/proteins
      n <- dim(dataExprMat)[1] # Number of experiments
      lambda = sqrt(2*log(p/sqrt(n))/n)
      
    }
    outlist1 <- FastGGM(dataExprMat,lambda)
    partialCor <- outlist1$partialCor
    rownames(partialCor) <- colnames(dataExprMat)
    colnames(partialCor) <- colnames(dataExprMat)
    cortab<-flattenCorrMatrix(partialCor, outlist1$p_partialCor)
    colnames(cortab) <- c('Prot1','Prot2','correlation', 'pval')
    cortab$pval.adj <- p.adjust(cortab$pval,"BH")

    
  } else if(corr_method=="pcor"){
    partialCor <- pcor(dataExprMat)$estimate
    rownames(partialCor) <- colnames(dataExprMat)
    colnames(partialCor) <- colnames(dataExprMat)
    cortab$pval.adj <- p.adjust(cortab$pval.adj,"BH")
    cortab <- flattenCorrMatrix2(partialCor)
    colnames(cortab) <- c('Prot1','Prot2','correlation')
    #corsig <- cortab
    
  }else{
    stop("Correlation method not supported. 
         Supported methods are: pearson, proportionality, GGM")
  }
  
  
  return(cortab)
  
}


fiedler_vector <- function(g){
  M <- laplacian_matrix(g, sparse = TRUE)
  f <- function(x,extra = NULL){
    as.vector(M%*%x)
  }
  fvec <- arpack(f,sym = TRUE,options=list(n = vcount(g),nev = 2,ncv = 8, 
                                           which = "SM",maxiter = 2000))
  return(fvec$vectors[,2])
}


# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


getActivationTime <- function(protsig,n_points){
  t <- rep(1,dim(protsig)[1])
  initp_logfc <- ncol(protsig)-n_points+2
  for (i in 1:dim(protsig)[1]){
    logvc_v <- protsig[i,initp_logfc:ncol(protsig)]
    t[i] <- Position(function(x) x > 1 || x < -1, logvc_v )
  }
  return(t)
}

minDistToVirus <- function(g,nprops,virus_nodes){
  s.paths <- shortest.paths(g,algorithm = "dijkstra")
  sp <- flattenCorrMatrix2(s.paths)
  names(sp) <- c("Prot1", "Prot2", "Shortest_path")
  min_dist_virus <- rep(-1,dim(nprops)[1])
  for (i in 1:dim(nprops)[1]){ 
    ex <- rbind(sp[sp$Prot1==nprops[i,1],],sp[sp$Prot2==nprops[i,1],])
    m1 <- min(ex[ex$Prot1 %in% virus_nodes,3])
    m2 <- min(ex[ex$Prot2 %in% virus_nodes,3])
    min_dist_virus[i] <- min(m1,m2)
  }
  return(min_dist_virus)
}

avDistToVirus <- function(g,nprops,virus_nodes){
  s.paths <- shortest.paths(g,algorithm = "dijkstra")
  sp <- flattenCorrMatrix2(s.paths)
  names(sp) <- c("Prot1", "Prot2", "Shortest_path")
  min_dist_virus <- rep(-1,dim(nprops)[1])
  for (i in 1:dim(nprops)[1]){ 
    ex <- rbind(sp[sp$Prot1==nprops[i,1],],sp[sp$Prot2==nprops[i,1],])
    m1 <- ex[ex$Prot1 %in% virus_nodes,3]
    m2 <- ex[ex$Prot2 %in% virus_nodes,3]
    if (length(m1)==0 && length(m2)!=0){
      min_dist_virus[i] <- mean(m2)
    } 
    if(length(m1)!=0 && length(m2)==0){ 
      min_dist_virus[i] <- mean(m1)
    } 
    if (length(m1)!=0 && length(m2)!=0){
      min_dist_virus[i] <- mean(c(m1,m2))
    }
  }
  return(min_dist_virus)
}


# Create variant colors from viridis palette with less saturation
vir_lite = function(cols, ds=0.4, dv=0.7) {
  cols = rgb2hsv(col2rgb(cols))
  cols["v", ] = cols["v", ] + dv*(1 - cols["v", ])
  cols["s", ] = ds*cols["s", ]
  apply(cols, 2, function(x) hsv(x[1], x[2], x[3]))
}


plotTimeSeries <- function(n_points,protsig,origin='Virus',nprot=20){
  max_point <- 4+n_points-1
  virusprot <- filter(protsig, Origin==origin) 
  virusprot <- virusprot[,c(1,4:max_point)]
  rownames(virusprot) <- virusprot$UniProtID
  virusprot
  virusp.ts <- pivot_longer(virusprot, c(2:ncol(virusprot)), 
                            names_to = "time", values_to = "Rel_abund")
  meastimes <- gsub('\\D','',virusp.ts$time)
  virusp.ts$time <- as.numeric(meastimes)
  # virus host proteins
  vhprot <- protsig[,c(1,4:max_point)]
  rownames(vhprot) <- vhprot$UniProtID
  vhprot.ts <- pivot_longer(vhprot, c(2:ncol(virusprot)),
                            names_to = "time", values_to = "Rel_abund")
  meastimes <- gsub('\\D','',vhprot.ts$time)
  vhprot.ts$time <- as.numeric(meastimes)
  
  # example code to plot specific time series between proteins
  nrowsp <- nprot*n_points
  virusp.tsp <- as.data.frame(virusp.ts[1:nrowsp,]) # plots first nprot proteins
  ggplot(virusp.tsp, aes(x = time, y = Rel_abund)) + 
    geom_line(aes(color = UniProtID), size = 1) +
    scale_y_continuous(trans='log10')+
    xlab("Time from infection(hours)") +
    ylab("Relative abundance") +
    scale_color_viridis(discrete = TRUE)+
    theme_ipsum()
}


densityPlot <- function(cortab_pearson,corsig_proportionality,cortab_gaussian){
  
  pearson <- rep("PCC",length(cortab_pearson$correlation))
  proportionality <- rep("CCC",length(corsig_proportionality$correlation))
  gaussian <- rep("GGM",length(cortab_gaussian$correlation))
  
  combined_density <- data.frame(method=pearson,
                                 correlation=cortab_pearson$correlation)
  combined_density  <- rbind(combined_density,
                             data.frame(method=proportionality,
                                        correlation=corsig_proportionality$correlation))
  combined_density  <- rbind(combined_density,
                             data.frame(method=gaussian,
                                        correlation=cortab_gaussian$correlation))
  
  
  # Compare the 3 methods
  ggplot(data=combined_density, aes(x=correlation, group=method, fill=method)) +
    geom_density(adjust=1.5, alpha=.6) +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme_ipsum() +
    theme(legend.title = element_text(size=20))+
    theme(legend.text = element_text(size=15))
}


pruneAdjMatrix <- function(corsig,cutoff_filter){
  
  # Scan number of links, nodes and modularity in the biggest connected 
  #components for different correlation thresholds
  corr_cutoff <- seq(0.95, 0.995, by=0.001)
  props <- cutOffScan(corr_cutoff,corsig)
  
  if(cutoff_filter == "modularity"){
    
    # Select cutoff when modularity > 0.6
    
    corr_threshold <- props[min(which(props$modularity > 0.6)),1]
    cat("Selected cut-off for the correlation coeficient: ",corr_threshold,"\n")
    corsig  <- filter(corsig, abs(correlation) > corr_threshold)
    
  }
  
  if(cutoff_filter == "nlinks"){ 
    # For each protein node, get only the 3 edges with more correlation
    totpairs <- c(as.character(corsig[,1]),as.character(corsig[,2]))
    allProts <- unique(totpairs)
    corsig3 <- data.frame(Prot1=character(),Prot2=character(),correlation=double())  
    for(i in 1:length(allProts)){
      p <- corsig %>% filter(grepl(i, Prot1))
      p$abs.corr <- abs(p$correlation)
      p_ordered <- p[order(p$abs.corr),]
      t <- tail(p_ordered,n=nlinks) # Get only the n edges with more correlation
      t$abs.corr <-NULL
      corsig3 <- rbind(corsig3,t)
    }
    # We need to remove duplicates
    corsig <- unique(corsig3)
  }
  
  if(cutoff_filter == "maxconnections"){
    
    # Select cutoff when number of links <= 5100
    
    corr_threshold <- props[min(which(props$links <= 5100)),1]
    cat("Selected cut-off for the correlation coeficient: ",corr_threshold,"\n")
    corsig  <- filter(corsig, abs(correlation) > corr_threshold)
  }
  
  
  plot_nodes <- ggplot(props, aes(x=cutoff, y=nodes)) +
    geom_line(color="#440154FF", size=0.75)+
    geom_point(color="#440154FF")+
    geom_vline(xintercept=corr_threshold, color="orange", size=0.75, linetype=2)+
    theme_ipsum()+
    annotate(geom="text", x=corr_threshold, y=100, 
             label=corr_threshold, fontface = "bold") 
  plot_links <- ggplot(props, aes(x=cutoff, y=links)) +
    geom_line(color="#21908CFF", size=0.75)+
    geom_point(color="#21908CFF")+
    geom_vline(xintercept=corr_threshold, color="orange", size=0.75, linetype=2)+
    theme_ipsum()
  plot_modularity <- ggplot(props, aes(x=cutoff, y=modularity)) +
    geom_line(color="#FDE725FF", size=0.75)+
    geom_point(color="#FDE725FF")+
    geom_vline(xintercept=corr_threshold, color="orange", size=0.75, linetype=2)+
    theme_ipsum()
  final_picture <- ggarrange(plot_nodes, plot_links, plot_modularity,
            labels = c("a", "b", "c"),
            ncol = 1, nrow = 3)
  print(final_picture)
  return(corsig)
}


calculateCommunities <- function(g, iterations=100){
  modularity_score <- 0
  final_clusters = cluster_louvain(g)
  for (i in 1:iterations) {
    clusters = cluster_louvain(g)
    if (modularity(clusters) > modularity_score){
      final_clusters <- clusters
      modularity_score <- modularity(final_clusters)
      cat("Modularity score of the best community detection out of ",iterations,
          " iterations: " , modularity_score,"\n")
    }
  }
  print(sizes(final_clusters))
  return(final_clusters)
}

getVirusFrac <- function(final_clusters, virus_nodes){
  number_of_communities <- length(sizes(final_clusters))
  mem <- membership(final_clusters)
  virus_percent <- rep(-1, number_of_communities)
  for (i in 1:number_of_communities) {
    
    virus_percent[i] <- length(names(mem[mem==i])[names(mem[mem==i]) %in% virus_nodes])/sizes(final_clusters)[i]
    cat(i," community: ", length(names(mem[mem==i])[names(mem[mem==i]) %in% virus_nodes]),"/",sizes(final_clusters)[i] ,"\n")
    
    
  }
  names(virus_percent) <- names(sizes(final_clusters))
  return(virus_percent)
}



getCommunityDefectors <- function(defectors, virus_nodes, final_clusters, g){
  
  # Clusters without virus proteins
  gnv <- g - virus_nodes
  final_clusters_nv <- calculateCommunities(gnv)
  
  mem_nv <- membership(final_clusters_nv)
  mem_defectors <- numeric(length(defectors))
  names(mem_defectors) <- defectors
  
  for (d in defectors){
    community_d <- mem_nv[d]
    members_d <- names(mem_nv[mem_nv[]==community_d])
    overlapping <- 0
    # Iterate over communities calculated with the viruses
    for (i in 1:length(final_clusters)){ 
      com <- names(mem[mem[]==i])
      overlapping_com_i <- length(com[com %in% members_d])/length(members_d)
      # We assign the community with more overlapping members
      if ( overlapping_com_i > overlapping ){
        mem_defectors[d] <- i
        overlapping <- overlapping_com_i
      }
    }
  }
  
  final_mem <- c(mem[!names(mem) %in% defectors],mem_defectors)
  
  # Remove virus nodes for functional annotation
  return(final_mem[!names(final_mem) %in% virus_nodes])

}


# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
# For the ORA test (A) a subset of genes is used as input, and then a one-sided 
# version of Fisher's exact test will be performed to find enriched GO categories
ora <- function(mem_host,enrichment_threshold,virus){
  community <- names(mem_host[mem_host[]==1])
  info_community <- filter(protsig,protsig$UniProtID %in% community)
  genes_interest_symbol <- as.character(info_community$GeneID)
  genes_interest_entrez <- bitr(genes_interest_symbol, fromType="SYMBOL", 
                                toType="ENTREZID", OrgDb="org.Hs.eg.db", 
                                drop = TRUE)
  gene_clusters <- list(genes_interest_entrez$ENTREZID)
  
  i_clusters <- 1:length(final_clusters)
  name_communities <- sprintf("comm %d",i_clusters)
  
  for(i in 2:length(final_clusters)) {
    community <- names(mem_host[mem_host[]==i])
    cat(community, "\n")
    info_community <- filter(protsig,protsig$UniProtID %in% community)
    genes_interest_symbol <- as.character(info_community$GeneID)
    genes_interest_entrez <- bitr(genes_interest_symbol, fromType="SYMBOL", 
                                  toType="ENTREZID", OrgDb="org.Hs.eg.db", 
                                  drop = TRUE)
    gene_clusters <- append(gene_clusters,list(genes_interest_entrez$ENTREZID))
  }
  
  names(gene_clusters) <- name_communities
  
  
  # Remove small communities (<15 members)
  gene_clusters_filtered <- Filter(function(x){length(x)>=15},gene_clusters)
  
  
  ego <- compareCluster(gene_clusters_filtered, fun="enrichGO", 
                        OrgDb="org.Hs.eg.db", ont="BP", pAdjustMethod = "BH",
                        pvalueCutoff=enrichment_threshold,readable = TRUE)
  
  # remove redundant terms
  ego_simplified <- clusterProfiler::simplify(
    ego,
    cutoff = 0.5,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
  )
  
  
  # Further simplification of uninteresting terms 
  remove_description_EBV <- c("regulation of endopeptidase activity",
                              "plasma membrane invagination",
                              "coagulation",
                              "negative regulation of stress-activated MAPK cascade",
                              "regulation of striated muscle tissue development", 
                              "epiboly", "epiboly involved in wound healing",
                              "positive regulation of collagen biosynthetic process", 
                              "regulation of establishment of cell polarity",
                              "negative regulation of fibrinolysis")
  remove_description_VACV <- c("regulation of endopeptidase activity", 
                               "regulation of transcription from RNA polymerase II",
                               "female pregnancy", "regulation of protein binding",
                               "regulation of transcription from RNA polymerase II promoter in response to hypoxia",
                               "regulation of transcription from RNA polymerase II promoter in response to stress",
                               "cellular response to transforming growth factor beta stimulus",
                               "peptidyl-tyrosine phosphorylation", 
                               "morphogenesis of embryonic epithelium",
                               "regulation of reproductive process",
                               "sex determination", "response to vitamin")
  remove_description_HSV1 <- c("negative regulation of bone mineralization", 	
                               "negative regulation of endopeptidase activity",
                               "in utero embryonic development", 
                               "mitotic sister chromatid cohesion",
                               "meiotic cell cycle process", "aging", 
                               "response to radiation", 
                               "establishment of chromosome localization", 
                               "neural tube closure", "regulation of neurogenesis",
                               "positive regulation of mitochondrial translation",
                               "mitotic sister chromatid segregation", 
                               "mesenchyme development", 
                               "regulation of stem cell differentiation")
  
  remove_description = switch(  
    virus,  
    "VACV"= remove_description_VACV,
    "EBV"= remove_description_EBV,
    "HSV1" = remove_description_HSV1,
  ) 
  
  
  ego_non_redundant <- ego_simplified
  for(descr in remove_description) {
    go_id_to_remove <- ego_simplified@compareClusterResult$ID[ego_simplified@compareClusterResult$Description==descr]
    ego_non_redundant <- dropGO(ego_non_redundant, term=go_id_to_remove)
  }
  
  return(ego_non_redundant)
  
}



createNodesRelations <- function(corsig,protsig){
  # Add UniprotID and Origin as properties of the nodes
  totpairs <- c(as.character(corsig[,1]),as.character(corsig[,2]))
  d <- as.data.frame(unique(totpairs))
  colnames(d) <- c("UniProtID")
  n <- inner_join(protsig,d)
  nprops <- data.frame(n$UniProtID,n$Origin)
  colnames(nprops) <- c ("UniProtID","Origin")
  prot1 <- as.data.frame(as.character(corsig$Prot1))
  colnames(prot1) <- "UniProtID"
  prot2 <- as.data.frame(as.character(corsig$Prot2))
  colnames(prot2) <- "UniProtID"
  nodes1 <- left_join(prot1,protsig)
  nodes2 <- left_join(prot2,protsig)
  protsig_with_time_and_regulation <- data.frame(UniProtID=protsig$UniProtID,
                                                 activationTime=protsig$activationTime,
                                                 regulation=protsig$regulation)
  nprops <- merge(nprops,protsig_with_time_and_regulation)
  
  # Create relationships matrix
  relations <- data.frame("Prot1" = corsig$Prot1, "Prot2" = corsig$Prot2)
  # Add positive correlation as property of the  edges
  relations$positive_corr <- corsig$corr>=0
  
  # Construct Interaction types as a property of the edges 
  intertype <- paste(as.character(nodes1$Origin),as.character(nodes2$Origin),sep='')
  intertype <- gsub('VirusVirus', 'VV', intertype)
  intertype <- gsub('VirusHost', 'VH', intertype)
  intertype <- gsub('HostVirus', 'VH', intertype)
  intertype <- gsub('HostHost', 'HH', intertype)
  relations$InterType <- intertype
  
  setClass("netcomponents", representation(node_properties = "data.frame", 
                                           interactions = "data.frame"))
  nc <- new("netcomponents", node_properties = nprops, interactions = relations)
  
  return(nc)
}