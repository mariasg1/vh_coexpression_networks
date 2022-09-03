library(viridis, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(igraph, quietly = TRUE)
library(eulerr, quietly = TRUE)
library(robustbase, quietly = TRUE)
library(hrbrthemes, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(ggpubr, quietly = TRUE)
library(forcats, quietly = TRUE)


source("C:/Users/sierram3/Documents/GitHub/vh_coexpression_networks/code/coexpr_net_functions.R")
setwd("C:/Users/sierram3/Documents/GitHub/vh_coexpression_networks/data/") 
options(warn=-1)
set.seed(1)



load("allProteins.RData")
load("allNets.RData")

virus <- "HSV1" # one of "VACV", "EBV", "HSV1"
corr_method <- "proportionality" # one of "pearson", "proportionality", "gaussian"
file_name <- paste (virus,"_results_",corr_method,".Rdata", sep = "")
load(file_name)

################################################################################
# Euler Diagram of network overlapping
################################################################################
g_proporc = switch(  
  virus,  
  "VACV"= g_proporc_VACCV,
  "EBV"= g_proporc_EBV,
  "HSV1"= g_proporc_HSV1,
)
g_pearson = switch(  
  virus,  
  "VACV"= g_pearson_VACCV,
  "EBV"= g_pearson_EBV,
  "HSV1"= g_pearson_HSV1,
)
g_ggm = switch(  
  virus,  
  "VACV"= g_ggm_VACCV,
  "EBV"= g_ggm_EBV,
  "HSV1"= g_ggm_HSV1,
)

edgespr <- as_edgelist(g_proporc)
edgespe <- as_edgelist(g_pearson)
edgesggm <- as_edgelist(g_ggm)

edgespr_concat <- paste(edgespr[,1],"-",edgespr[,2],sep="")
edgespe_concat <- paste(edgespe[,1],"-",edgespe[,2],sep="")
edgesggm_concat <- paste(edgesggm [,1],"-",edgesggm [,2],sep="")

area_all <- length(intersect(intersect(edgespr_concat,edgespe_concat),
                             edgesggm_concat))
area_pr_pe <- length(intersect(edgespr_concat,edgespe_concat))-area_all
area_pr_ggm <- length(intersect(edgespr_concat,edgesggm_concat))-area_all
area_pe_ggm <- length(intersect(edgespe_concat,edgesggm_concat))-area_all
area_pr <- length(edgespr_concat)-area_all-area_pr_ggm- area_pr_pe
area_pe <- length(edgespe_concat)-area_all-area_pr_pe-area_pe_ggm
area_ggm <- length(edgesggm_concat)-area_all-area_pr_ggm-area_pe_ggm

x <- c(Proportionality=area_pr,Pearson=area_pe,GGM=area_ggm,
       "Proportionality&Pearson"=area_pr_pe, 
       "Proportionality&GGM"=area_pr_ggm,
       "Pearson&GGM"=area_pe_ggm,
       "Proportionality&Pearson&GGM"=area_all
       
)

x <- c(CCC=area_pr,PCC=area_pe,GGM=area_ggm,
       "CCC&PCC"=area_pr_pe, 
       "CCC&GGM"=area_pr_ggm,
       "PCC&GGM"=area_pe_ggm,
       "CCC&PCC&GGM"=area_all
       
)

viridis_yellow_tp <- rgb(253, 231, 37, alpha = 0.2 * 255, maxColorValue = 255)
viridis_green_tp <- rgb(33, 145, 140, alpha = 0.22 * 255, maxColorValue = 255)
viridis_purple_tp <- rgb(68, 1, 84, alpha = 0.2 * 255, maxColorValue = 255)

plot(euler(x,shape = "ellipse"), quantities = list(type = "counts", font = 40),
     fills=c("white","white","white") , 
     edges = list(col = c("#440154ff",'#fde725ff','#21908dff')),
     labels=list(cex = 3,col = c("#440154ff",'#fde725ff','#21908dff')))


# Power-law fit
host_relations <- relations[relations$InterType=="HH",1:2]
virus_relations <- relations[relations$InterType=="VV",1:2]
corsig_virus <- merge(corsig,virus_relations)
corsig_host <- merge(corsig,host_relations)
pwLawFitDraw(corsig)
pwLawFitDraw(corsig_virus)
pwLawFitDraw(corsig_host)


# Network picture
networkVisualization(g)


################################################################################
# Dynamic analysis based on activation time
################################################################################

immediate_early_genes_EBV <- c("BZLF1", "BRLF1")
immediate_early_genes_VACV <- c("C11", "C6", "M2", "N2", "K1", "K3", "K5",
                                "K7", "E3", "E4", "E5", "O1", "I2", "I3", 
                                "G5", "L2", "D9", "D10", "A8", "A33", "A35",
                                "A37", "A44", "A47", "A48", "A52", "B2", 
                                "B3", "B12", "B13", "B15")
immediate_early_genes_HSV1 <- c("UL6","UL48","RS1", "RL2", "US1", "UL54", "US12")
immediate_early_genes = switch(  
  virus,  
  "VACV"= immediate_early_genes_VACV,
  "EBV"= immediate_early_genes_EBV,
  "HSV1"= immediate_early_genes_HSV1,
)
immediate_early <- protdf[protdf$GeneID %in% immediate_early_genes,1]

nprops$av_dist_to_virus_immediate_early <- avDistToVirus(g,nprops,immediate_early)
nprops_onlyhost <- filter(nprops, nprops$Origin=="Host")
nprops_onlyvirus <- filter(nprops, nprops$Origin=="Virus")
nprops_onlyhost %>%
  ggplot( aes(x=as.factor(activationTime), y=av_dist_to_virus_immediate_early, 
              fill=as.factor(activationTime))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab("Average shortest path to immediate-early viral proteins")+
  xlab("Activation time")

nprops$closeness <- closeness(g, normalized= TRUE)
nprops_onlyhost <- filter(nprops, nprops$Origin=="Host")
nprops %>%
  ggplot( aes(x=as.factor(activationTime), y=closeness, 
              fill=as.factor(activationTime))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab("Closeness")+
  xlab("Activation time")

colrs <- viridis(6)
colrs <- adjustcolor(colrs,alpha.f=0.8)
V(g)$color <- colrs[as.integer(as.factor(V(g)$activationTime))]
shapes <- c("circle", "square")
V(g)$shape  <- shapes[as.integer(as.factor(V(g)$Origin))]
# Compute node degrees (#links) and use that to set node size:
deg <- igraph::degree(g)
V(g)$size <- 3*log(deg)
# change edge color
E(g)$edge.color <- "gray80"
plot(g,layout=layout_nicely(g, dim=2),vertex.label=NA)
legend("topleft",bty = "n", cex=2,
       legend= names(protdf[5:(5+n_points-2)]),
       fill=colrs, border=NA)


################################################################################
# Compare Host and Viral Networks
################################################################################


# Clusters without virus proteins
virus_nodes <- V(g)[V(g)$Origin=="Virus"]
ghh <- g - virus_nodes
# Clusters without virus proteins
host_nodes <- V(g)[V(g)$Origin=="Host"]
gvv <- g - host_nodes


# DEGREE HISTOGRAM -------------------------------------------------------------
df1 <- data.frame(degree=degree(g))
df1$network <- rep("Combined_Network",dim(df1)[1])
df1$prot <- rownames(df1)
df2 <- data.frame(degree=degree(gvv))
df2$network <- rep("Virus_Network",dim(df2)[1])
df2$prot <- rownames(df2)
df3 <- data.frame(degree=degree(ghh))
df3$network <- rep("Host_Network",dim(df3)[1])
df3$prot <- rownames(df3)
df <- rbind(df1,df2,df3)
df$fraction <- df$degree/length(V(g))
df %>%
  mutate(text = fct_reorder(network, degree)) %>%
  ggplot( aes(x=degree, color=network, fill=network)) +
  geom_histogram(aes(y=..count../sum(..count..)),alpha=0.6, binwidth = 1) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 10)
  ) +
  xlab("Degree") +
  ylab("Number of nodes") +
  facet_wrap(~text)


# EIGENCENTRALITY VS DEGREE ----------------------------------------------------
e <- eigen_centrality(g)$vector
evv <- eigen_centrality(gvv)$vector
ehh <- eigen_centrality(ghh)$vector


# Viral Network
virus_linked_host <- nprops[nprops$virus_host_connector,1]
df1 <- data.frame(eigencentrality=evv)
df1$node_type <- rep("Virus",dim(df1)[1])
df1$prot <- rownames(df1)
df2 <- data.frame(eigencentrality=evv[names(evv) %in% virus_linked_host])
df2$node_type <- rep("Virus Connector",dim(df2)[1])
df2$prot <- rownames(df2)
df1_filtered <- filter(df1, !df1$prot %in% df2$prot)
dfv <- rbind(df1_filtered,df2)
dgv <- data.frame(degree=degree(gvv))
dgv$prot <- rownames(dgv)
dfvv <- merge(dfv,dgv)

# Calculate p-values
# Are degree and eigencentrality of connector nodes in the Viral network
# above or below the average?
# degree test
t.test(dfvv[dfvv$node_type=="Virus Connector",4], dfvv[dfvv$node_type=="Virus",4]) 
# eigencentrality test
t.test(dfvv[dfvv$node_type=="Virus Connector",2], dfvv[dfvv$node_type=="Virus",2]) 

ggplot(dfvv, aes(x=dfvv$degree, y=dfvv$eigencentrality, color=node_type)) + 
  geom_point(size=2)+
  scale_color_manual(values=c("#fde725","#31688e") )+
  xlab("Degree") +
  ylab("Eigencentrality") +
  theme_ipsum()+
  geom_vline(xintercept=mean(degree(gvv)),color=('sienna1'),
             linetype="dashed")+
  geom_hline(yintercept=mean(df1$eigencentrality), color=('sienna1'),
             linetype="dashed")

# Host Network
host_linked_virus <- nprops[nprops$host_virus_connector,1]
df3 <- data.frame(eigencentrality=ehh)
df3$node_type <- rep("Host",dim(df3)[1])
df3$prot <- rownames(df3)
df4 <- data.frame(eigencentrality=ehh[names(ehh) %in% host_linked_virus])
df4$node_type <- rep("Host_Virus_Connector",dim(df4)[1])
df4$prot <- rownames(df4)
df3_filtered <- filter(df3, !df3$prot %in% df4$prot)
dfh <- rbind(df3_filtered,df4)
dgh <- data.frame(degree=degree(ghh))
dgh$prot <- rownames(dgh)
dfhh <- merge(dfh,dgh)

# Calculate p-values
# Are degree and eigencentrality of connector nodes in the Host network
# above or below the average?
# degree test
t.test(dfhh[dfhh$node_type=="Host_Virus_Connector",4], dfhh[dfhh$node_type=="Host",4])
# eigencentrality test
t.test(dfhh[dfhh$node_type=="Host_Virus_Connector",2], dfhh[dfhh$node_type=="Host",2])

ggplot(dfhh, aes(x=dfhh$degree, y=dfhh$eigencentrality, color=node_type)) + 
  geom_point(size=2 )+
  scale_color_manual(values=c("#35b779","#440154") )+
  theme_ipsum()+
  xlab("Degree") +
  ylab("Eigencentrality") +
  geom_vline(xintercept=mean(degree(ghh)),color=('sienna1'),
             linetype="dashed")+
  geom_hline(yintercept=mean(df3$eigencentrality), color=('sienna1'),
             linetype="dashed")


# Impact on Viral Network eigencentrality after interconnecting to Host Network

df1 <- data.frame(eigencentrality_VV_network=evv)
df1$node_type <- rep("Virus",dim(df1)[1])
df1$prot <- rownames(df1)
df2 <- data.frame(eigencentrality_VV_network=evv[names(evv) %in% virus_linked_host])
df2$node_type <- rep("Virus Connector",dim(df2)[1])
df2$prot <- rownames(df2)
df1_filtered <- filter(df1, !df1$prot %in% df2$prot)
dfv <- rbind(df1_filtered,df2)
dgv <- data.frame(eigencentrality_HV_network=e)
dgv$prot <- rownames(dgv)
dfvv <- merge(dfv,dgv)

dfvv <- dfvv[dfvv$prot != "6FT_ORF232_New",] # remove outlier

ggplot(dfvv, aes(x=dfvv$eigencentrality_VV_network, y=dfvv$eigencentrality_HV_network,
                 color=node_type)) + 
  geom_point(size=2 )+
  scale_color_manual(values=c("#fde725","#31688e") )+
  theme_ipsum()+
  xlab("Eigencentrality VV network") +
  ylab("Eingencentrality HV network")

# Impact on Host Network eigencentrality after interconnecting to Viral Network

df1 <- data.frame(eigencentrality_HH_network=ehh)
df1$node_type <- rep("Host",dim(df1)[1])
df1$prot <- rownames(df1)
df2 <- data.frame(eigencentrality_HH_network=ehh[names(ehh) %in% host_linked_virus])
df2$node_type <- rep("Host Connector",dim(df2)[1])
df2$prot <- rownames(df2)
df1_filtered <- filter(df1, !df1$prot %in% df2$prot)
dfh <- rbind(df1_filtered,df2)
dgh <- data.frame(eigencentrality_HV_network=e)
dgh$prot <- rownames(dgh)
dfhh <- merge(dfh,dgh)


dfhh <- dfhh[!(dfhh$prot %in% "P48552"),] # remove outlier

# Scatterplot eigencentrality host
ggplot(dfhh, aes(x=dfhh$eigencentrality_HH_network, y=dfhh$eigencentrality_HV_network, 
                 color=node_type)) + 
  geom_point(size=2 )+
  scale_color_manual(values=c("#35b779","#440154") )+
  xlab("Eigencentrality HH network") +
  ylab("Eigencentrality HV network") +
  theme_ipsum()


################################################################################
# Community partitions
################################################################################

clustersVisualization(final_clusters,g)

# Simplified Community overview ------------------------------------------------
mem <- membership(final_clusters)
# Contract vertices
V(g)$weight <- 1
E(g)$weight <- 1
gcon <- contract.vertices(g, mem, 
                          vertex.attr.comb = list(weight = "sum", 
                                                  name = function(x)x[1], "ignore"))
E(gcon)$weight <- 1
# Simplify edges
gcon <- igraph::simplify(gcon, edge.attr.comb=list(weight="sum"))
#igraph::simplify(gcon, remove.multiple = FALSE, remove.loops = TRUE)
V(gcon)$degree <- unname(degree(gcon))
V(gcon)$name <- names(sizes(final_clusters))
V(gcon)$nnodes <- sizes(final_clusters)

ncommunities <- length(sizes(final_clusters))
colrs <- viridis(ncommunities)
colrs <- adjustcolor(colrs,alpha.f=0.8)
V(gcon)$color <- colrs[as.integer(as.factor(V(gcon)$name))]
# Compute node degrees (#links) and use that to set node size:

V(gcon)$size <- 5*log(V(gcon)$nnodes)
# change edge color
E(gcon)$edge.color <- "gray80"

commorigin <- rep(-1,ncommunities)
for (i in 1:ncommunities){ 
  commorigin[i] <- names(which.max(table(factor(nprops[nprops$community == i,
                                                       which( colnames(nprops)=="Origin" )]))))
}
V(gcon)$Origin <- commorigin
shapes <- c("circle", "square")
V(gcon)$shape  <- shapes[as.integer(as.factor(V(gcon)$Origin))]

plot(gcon,layout=layout_nicely(gcon, dim=2),edge.width = log(E(gcon)$weight),
     vertex.label.color=c("black"))

# Heatmap communities vs activation time ---------------------------------------

z <- c()
r <- c()
for (i in 1:ncommunities){ 
  actime <- table(factor(nprops[nprops$community == i
                                ,which( colnames(nprops)=="activationTime" )],
                         levels = levels(factor(nprops$activationTime))))
  z <- c(z,actime/sum(actime))
  regulation <- table(factor(nprops[nprops$community == i,
                                    which( colnames(nprops)=="regulation" )],
                             levels = levels(factor(nprops$regulation))))
  r <- c(r,regulation)
}
  
y <- rep(names(sizes(final_clusters)),each=length(levels(factor(nprops$activationTime))))
  
x <- rep(levels(factor(nprops$activationTime)),length(sizes(final_clusters)))
k <- rep(levels(factor(nprops$regulation)),length(sizes(final_clusters)))
  
data <- data.frame(activationTime=x,community=y,frac_nodes=z)
if (virus=="HSV1"){
  data$community <- factor(data$community, ordered=TRUE, 
                           levels = c("2", "1", "5", "3", "6", "4", "7"))
}
if (virus=="VACV"){
  data$community <- factor(data$community, ordered=TRUE, 
                           levels = c("3", "1", "2", "4"))
}
data <- data[data$community!=7,]
ggplot(data, aes(x=activationTime, y=community, fill= frac_nodes)) + 
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_gradientn(colors=vir_lite(viridis(3), ds=.7, dv=0))+
  theme_ipsum()+
  geom_text(aes(label = round(data$frac_nodes,2)), color = "white") 


# Functional Enrichment visualization ------------------------------------------
dotplot( ego , font.size = 8, showCategory=10) + scale_color_viridis()

