library(ggpubr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(data.table, quietly = TRUE)
library(robustbase, quietly = TRUE)
library(psych, quietly = TRUE)
library(topGO, quietly = TRUE)
library(FastGGM, quietly = TRUE)
library(igraph, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(plyr, quietly = TRUE)
library(viridis, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(hrbrthemes, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(forcats, quietly = TRUE)
library(qgraph, quietly = TRUE)


source("C:/Users/sierram3/Documents/GitHub/vh_coexpression_networks/code/coexpr_net_functions.R")
setwd("C:/Users/sierram3/Documents/GitHub/vh_coexpression_networks/data/")  
options(warn=-1)
set.seed(1)

# ------------------------------------------------------------------------------
# coexpr_net_main
# ------------------------------------------------------------------------------
#  
# This code creates and analyzes Virus-Host co-expression networks from 
# proteomics data obtained from Quantitative Temporal Viromics experiments.
# Network creation can be done by means of different correlation methods.
# Analysis includes basic network properties and community detection with 
# functional annotation. 
#
# Inputs:
# * virus: one of "VACV", "EBV", "HSV1"
# * corr_method: one of "pearson", "proportionality", "gaussian"
# * FDR_threshold: to extract the significant differerentially expressed proteins 
# * LogFC_threshold: to extract the significant differerentially expressed proteins
# * enrichment_threshold: to filter functional annotation results
# * adj_pval_corr_threshold: to extract the significant correlated proteins 
#                           (only for gaussian and pearson methods)
# * lambda: Penalty parameter in Gaussian Graphical Model, larger lambda results 
#                 in faster calculation and a sparser graph
# * cutoff_filter: one of "modularity", "nlinks", "maxconnections"
# * nlinks: Only if cutoff_filter="nlinks". Number of links to retain per node
# 
# Outputs:
# * Inline basic network properties and pictures
# * Functional annotation of communities in a csv file
# * Node properties and relations in csv format
# * Network object exported in .R file
#
# ------------------------------------------------------------------------------

# Set to TRUE to re-load the raw data from all viruses
reload_proteins = FALSE
if(reload_proteins){
  loadProteins()
}


virus <- "HSV1" # one of "VACV", "EBV", "HSV1"
corr_method <- "proportionality" # one of "pearson", "proportionality", "gaussian"


# Establish thresholds to use:
# To determine the significant differentially expressed proteins between 
# infected cells and control
FDR_threshold <- 0.05
LogFC_threshold <- 1

# Correlations
adj_pval_corr_threshold <- 0.01 # For Pearson or GGM methods
lambda <- 2 # For GGM
cutoff_filter <- "modularity" # one of "modularity", "nlinks", "maxconnections"
nlinks <- 3 # Only used if cutoff:filter=nlinks

# Community-based over representation analysis
enrichment_threshold <- 0.05



# Check input parameter is correct
if (!(virus == "VACV" | virus == "EBV" | virus == "HSV1" )){
  stop(paste0("No data from virus ", virus))
}

#-----------------------------------------------------------------------------
#
# Load data from supported viruses
#
#-----------------------------------------------------------------------------

load("allProteins.RData")

protdf = switch(  
  virus,  
  "VACV"= df_prot_filtered_VACV,
  "EBV"= df_prot_filtered_EBV,
  "HSV1" = df_prot_filtered_HSV1,
)  
n_points = switch(  
  virus,  
  "VACV"= n_points_VACV,
  "EBV"= n_points_EBV,
  "HSV1"= n_points_HSV1,
)  

#-----------------------------------------------------------------------------
#
# Select only differentially expressed proteins, according to input thresholds:
# 
# * FDR_threshold 
# * LogFC_threshold
#
#-----------------------------------------------------------------------------

# Filtering:
# Select only proteins with significant different expression between control 
# and infected in at least one time point (assessed by FDR columns)
protFDR <- protdf[apply(protdf[,grepl("FDR",names(protdf))], 1, 
                        function(x) any(x < FDR_threshold)), ]
# Further filtering: 
# from those significantly expressed genes, keep only those with at least 
# two-fold expression differences relative to control (assessed by LogFC columns)
protsig <- protFDR[apply(protFDR[,grepl("logFC",names(protFDR))], 1, 
                         function(x) any(abs(x) > LogFC_threshold)), ]

# Add activation time and regulation information
protsig$activationTime <- getActivationTime(protsig,n_points)
protsig$regulation <- apply(protsig[,grepl("logFC",names(protsig))],1,
                            function(x) if(any(x <0)){return("down")}else{return("up")})



# ----------------------------------------------------------------------------
#
# Data exploratory analysis  
#
#-----------------------------------------------------------------------------

# Visualization of time series 
plotTimeSeries(n_points,protsig,origin='Virus',nprot=20)


#-----------------------------------------------------------------------------
#
# Calculate correlations with the different methods
#
#-----------------------------------------------------------------------------

max_point <- 4+n_points-1
m <- data.matrix(protsig)
protmat <- cbind(m[,4:max_point])
rownames(protmat) <- protsig$UniProtID
colnames(protmat) <- names(protsig)[4:max_point]
dataExprMat <- t(protmat)

# Check correlations (as scatterplots)
data <- as.data.frame(dataExprMat[,c(1:5)])
plot(data, pch=19 , cex=2 , col="#440154")

cortab_pearson <- calcCorrelation("pearson",dataExprMat)
corsig_proportionality <- calcCorrelation("proportionality",dataExprMat)
# cortab_gaussian <- calcCorrelation("GGM",dataExprMat) 
# This takes a long time, use saved dataframes instead
load("corr_ggm.RData")
cortab_gaussian = switch(  
  virus,  
  "VACV"= cortab_gaussian_VACV,
  "EBV"= cortab_gaussian_EBV,
  "HSV1"= cortab_gaussian_HSV1,
)

densityPlot(cortab_pearson,corsig_proportionality,cortab_gaussian)


# filter significant correlations
corsig_pearson <- filter(cortab_pearson, pval.adj < adj_pval_corr_threshold)
corsig_gaussian <- filter(cortab_gaussian_HSV1, 
                          pval.adj < adj_pval_corr_threshold)


corsig = switch(  
  corr_method,  
  "gaussian" = corsig_gaussian,
  "pearson" = corsig_pearson,
  "proportionality" = corsig_proportionality,
)




#-----------------------------------------------------------------------------
# Adjacency Matrix pruning
# Choose correlation cut-off  based on:
#
# * modularity
# * number of links per node
# * maximum number of links
# For the biggest connected component component
#-----------------------------------------------------------------------------

corsig <- pruneAdjMatrix(corsig,cutoff_filter)

#-----------------------------------------------------------------------------
#
# Construct protein co-expression network 
#
#-----------------------------------------------------------------------------

# Calculate max connected component.
gm <- graph_from_data_frame(corsig)
components <- igraph::clusters(gm, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(gm)[components$membership == biggest_cluster_id]
gc <- igraph::induced_subgraph(gm, vert_ids)
corsigmax <- filter(corsig, 
                    corsig$Prot2 %in% vert_ids$name |corsig$Prot1 %in% vert_ids$name)
corsig <- corsigmax 

# Generate node properties and relations
net_components <- createNodesRelations(corsig,protsig)
relations <- net_components@interactions
nprops <- net_components@node_properties

# Create the network
g <- graph_from_data_frame(relations, nprops,directed = FALSE)


#-----------------------------------------------------------------------------
#
# Network properties
#
#----------------------------------------------------------------------------- 

# Scale-free properties 
pwLawFit(g)

# number of edges (interactions)
nedges <- nrow(corsig)
cat("Number of edges: ", nedges, "\n")

# number of nodes
totpairs <- c(as.character(corsig[,1]),as.character(corsig[,2]))
nnodes <- length(unique(totpairs))

# Types of interactions
vhinter <- filter(relations, InterType=='VV')
nvhint <- nrow(vhinter)
vhinter <- filter(relations, InterType=='HH')
nvhint <- nrow(vhinter)
vhinter <- filter(relations, InterType=='VH')
nvhint <- nrow(vhinter)

cat("Number of nodes: ", nnodes, "\n")
cat('Number of virus-virus interactions',nvhint,"\n")
cat('Number of host-host interactions',nvhint,"\n")
cat('Number of virus-host interactions',nvhint,"\n")

# General network properties
networkAnalysis(g)

# Network properties differentiated by Origin (Virus, Host, VH connectors)
df <- data.frame(vhinter$Prot1, vhinter$Prot2)
v <- as.vector(unique(unlist(df)))
virus_nodes <- V(g)$name[V(g)$Origin=='Virus']
host_nodes <- V(g)$name[V(g)$Origin=='Host']
host_linked_virus <- v[!v %in% virus_nodes]
virus_linked_host <- v[v %in% virus_nodes]

# Add connector information to node properties
nprops$host_virus_connector <- nprops$UniProtID %in% host_linked_virus
nprops$virus_host_connector <- nprops$UniProtID %in% virus_linked_host

subnetworkAnalysis(g,host_linked_virus,virus_nodes,host_nodes)


#-----------------------------------------------------------------------------
#
# Communities and functional annotation 
#
#-----------------------------------------------------------------------------  


# Clusters with virus proteins
# Pick up the best partition out of 1000 iterations
final_clusters <- calculateCommunities(g)
mem <- membership(final_clusters)
number_of_communities <- length(sizes(final_clusters))

# Add community info to nprops
m <- data.frame(community = as.numeric(mem))
m$UniProtID <- names(mem)
nprops <- merge(nprops,m)

# Get proportion of virus per community
virus_percent <- getVirusFrac(final_clusters,virus_nodes)


# Functional annotation ------------------------------------------------------

# Generate network without viruses to see in which community defectors will be
# Defectors are those host nodes connected to the Viral Network

mem_host <- getCommunityDefectors(host_linked_virus, virus_nodes, 
                                  final_clusters, g)

# Community-Based Over Representation Analysis
ego <- ora(mem_host,enrichment_threshold,virus)


#-----------------------------------------------------------------------------
#
# Save results
#
#-----------------------------------------------------------------------------  

filename <- paste("enrichment_results_communities_",virus,"_",
                  corr_method,".csv", sep="")
write.table(ego@compareClusterResult, file = filename, sep = ";")

filename <- paste("nprops_",virus,"_",corr_method,".csv", sep="")
write.csv(relations,filename, row.names = FALSE)
filename <- paste("relations_",virus,"_",corr_method,".csv", sep="")
write.csv(nprops,filename, row.names = FALSE)

file_name <- paste (virus,"_results_",corr_method,".Rdata", sep = "")
save(protdf, n_points, corsig, nprops, relations, g, final_clusters, ego,
     file = file_name)
  