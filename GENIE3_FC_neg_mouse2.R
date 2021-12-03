library(tidyr)
library(dplyr)
library(igraph)
source('GENIE3_R/GENIE3.R')

#Read in data
expression_matrix <- read.csv('GSE152632_GEO_mberkchen_TRAP2_counts.csv')

#Make genes rownames
rownames(expression_matrix) <- expression_matrix[,1]
expression_matrix[,1] <- NULL

#Genes of interest from Dr. Chen
genes_memory_related = c("Arc","Atf3","Bdnf","Arpp21","Brd9","Camkv","Inf2","Inhba","Jph4","Mapk4","Penk","Plk2","Scg2","Sidt1","Slc29a4","Sorcs3","Sorcs1","Tiam2","Tiam3","Kcnq3","Il20rb","Glt8d2","Ctso","Calb1","Esd","Faim2","Itm2c","Matk","Ncald","Ptprz1","Rab3a","Vcp","Rfwd2","Rgs2","Ubc","Cdk14","Faah","Fam126b","Hars","Nav1","R3hdm1","Scn2b","Uchl1","Ttll7","Tub","Hsph1","Psmb1","Psmb3","Psmb4","Rnf4","Rbfox3","Atf4","Kcnk1","Glcci1","Psmd12","Hpca","Psmd4","Ptprn","Acap2","Ssb","Ncbp1","Epha7","Psmc6","Psmd3","Ctsd","Nars","Tpt1","Lmo7","Ttyh1","Tubala","Psmb2","Eif4g2","Zfand5","Fzd3","Rtn4rl1","Eif3c","Ypel5","Tpi1","Cltb")
memory_related_genes_expression <- expression_matrix[rownames(expression_matrix) %in% genes_memory_related, ]

###################################################################################################################
#FC_neg mouse2

FC_neg_cols_mouse2_GOI <- grep("^FC_neg.*mouse2$", colnames(memory_related_genes_expression))
FC_neg_mouse2_GOI <- select(memory_related_genes_expression,all_of(FC_neg_cols_mouse2_GOI))
FC_neg_mouse2_GOI_weight.matrix <- GENIE3(FC_neg_mouse2_GOI)
FC_neg_mouse2_GOI_linked.list <- get.link.list(FC_neg_mouse2_GOI_weight.matrix, report.max = 200)
write.csv(FC_neg_mouse2_GOI_linked.list, file = "Remote_memory_related_genes_results/FC_neg_mouse2_GENIE3_linkedlist_MemoryRelatedOnly_Top200.csv", row.names=F)

#Visualize in cytoscape
g_prep <- graph.data.frame(FC_neg_mouse2_GOI_linked.list, directed = F)
g <- get.adjacency(g_prep)
g_post <- graph.adjacency(g, mode = 'undirected', weighted = T)
g.cyto <-igraph.to.graphNEL(g_post)

#Cytoscape must be open before running this step
network_g <- createNetworkFromGraph("FC_neg_mouse2", graph=g.cyto)