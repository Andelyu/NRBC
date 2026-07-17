

################################################################################################################################################################
#---------------------------------------------------------part1: prepare data---------------------------------------------------------------#
################################################################################################################################################################

 
library(Seurat)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dbplyr)
library(VennDiagram)
library(CellChat)
library(nichenetr)
library(tibble)
library(tidyr )

library(RColorBrewer )
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]

setwd('~/2021_NRBC_chlyu/')
filt_NBRC_altas_seu=readRDS('20251125_filt_NBRC_altas_seu.rds')


# modified the function netAnalysis_signalingRole_heatmap 
computeCentralityLocal <- function(net) {
  centr <- vector("list")
  G <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  centr$outdeg_unweighted <- rowSums(net > 0)
  centr$indeg_unweighted <- colSums(net > 0)
  centr$outdeg <- igraph::strength(G, mode="out")
  centr$indeg <- igraph::strength(G, mode="in")
  centr$hub <- igraph::hub_score(G)$vector
  centr$authority <- igraph::authority_score(G)$vector # A node has high authority when it is linked by many other nodes that are linking many other nodes.
  centr$eigen <- igraph::eigen_centrality(G)$vector # A measure of influence in the network that takes into account second-order connections
  centr$page_rank <- igraph::page_rank(G)$vector
  igraph::E(G)$weight <- 1/igraph::E(G)$weight
  centr$betweenness <- igraph::betweenness(G)
  #centr$flowbet <- try(sna::flowbet(net)) # a measure of its role as a gatekeeper for the flow of communication between any two cells; the total maximum flow (aggregated across all pairs of third parties) mediated by v.
  #centr$info <- try(sna::infocent(net)) # actors with higher information centrality are predicted to have greater control over the flow of information within a network; highly information-central individuals tend to have a large number of short paths to many others within the social structure.
  centr$flowbet <- tryCatch({
    sna::flowbet(net)
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  centr$info <- tryCatch({
    sna::infocent(net, diag = T, rescale = T, cmode = "lower")
    # sna::infocent(net, diag = T, rescale = T, cmode = "weak")
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  return(centr)
}
netAnalysis_signalingRole_heatmap2 <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                               color.use = NULL, color.heatmap = "BuGn",title = NULL, width = 10, height = 8,
                                               font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat=mat[,colSums(mat)>0]
  mat=mat[rowSums(mat)>0,]
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA
  
  
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = c('#FFFFFF',grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 2, name = color.heatmap)))(100))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  
  #row_annotation_df=data.frame(row.names = slot(object, slot.name)$pathways,annotation=object@DB$interaction[match(slot(object, slot.name)$pathways,object@DB$interaction$pathway_name),'annotation'])
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)
  
  
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  mat[is.na(mat)]=0
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.cols,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = legend.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}


netAnalysis_signalingRole_heatmap3 <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                               color.use = NULL, color.heatmap = "BuGn",title = NULL, width = 10, height = 8,
                                               font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat=mat[,colSums(mat)>0]
  # mat=mat[rowSums(mat)>0,]
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA
  
  
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = c('#FFFFFF',grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 2, name = color.heatmap)))(100))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  
  #row_annotation_df=data.frame(row.names = slot(object, slot.name)$pathways,annotation=object@DB$interaction[match(slot(object, slot.name)$pathways,object@DB$interaction$pathway_name),'annotation'])
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)
  
  
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  mat[is.na(mat)]=0
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.cols,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = legend.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}


#-------------------------------- read data----------------------------------#

YS_subcelltype_ccd=readRDS('NRBC_YS_altas/YS_subcelltype_cellchat.rds')
FL_subcelltype_ccd=readRDS('NRBC_FL_altas/FL_subcelltyp_cellchat.rds')
FBM_subcelltype_ccd=readRDS('NRBC_BM_altas/FBM_subcelltype_cellchat.rds')
ABM_subcelltype_ccd=readRDS('NRBC_BM_altas/ABM_subcelltype_cellchat.rds')


#------------------------------------统计 incoming & outgoing signals source:target, Other celltype: NRBC celltye--------------------------------------------#
# 准备，构建所有stage list 
object.list <- list(YS = YS_subcelltype_ccd,FL=FL_subcelltype_ccd,FBM=FBM_subcelltype_ccd,
                    ABM=ABM_subcelltype_ccd)


stage_type_list=list()
levels(YS_subcelltype_ccd@idents)
levels(FL_subcelltype_ccd@idents)
levels(FBM_subcelltype_ccd@idents)
levels(ABM_subcelltype_ccd@idents)
#考虑了FL/FBM中YS——NRBC 接受到的信息
stage_type_list[['YS']][['sources.use']]=c(1:2,7:length(levels(YS_subcelltype_ccd@idents))) # no Ery
stage_type_list[['YS']][['targets.use']]=3:6    # Ery celltype
stage_type_list[['FL']][['sources.use']]=c(1:8,16:length(levels(FL_subcelltype_ccd@idents)))
stage_type_list[['FL']][['targets.use']]=9:15
stage_type_list[['FBM']][['sources.use']]=c(1:6,14:length(levels(FBM_subcelltype_ccd@idents)))
stage_type_list[['FBM']][['targets.use']]=7:13
stage_type_list[['ABM']][['sources.use']]=c(1:3,9:length(levels(ABM_subcelltype_ccd@idents)))
stage_type_list[['ABM']][['targets.use']]=4:8

rm(YS_subcelltype_ccd,FL_subcelltype_ccd,FBM_subcelltype_ccd,ABM_subcelltype_ccd);gc()


NRBC_LR_df=data.frame()
Other2Ery_df=data.frame()
Ery2Other_df=data.frame()
Ery2Ery_df=data.frame()


for ( stage in names(object.list)){
  tmp_Other2Ery_df=netVisual_bubble(object.list[[stage]],sources.use = stage_type_list[[stage]][['sources.use']],targets.use =stage_type_list[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Other2Ery_df$stage=stage;rownames(tmp_Other2Ery_df)=NULL
  tmp_Other2Ery_df=tmp_Other2Ery_df[!is.na(tmp_Other2Ery_df$prob),]
  
  Other2Ery_df=rbind(Other2Ery_df,tmp_Other2Ery_df)
  
  tmp_Ery2Other_df=netVisual_bubble(object.list[[stage]],sources.use = stage_type_list[[stage]][['targets.use']],targets.use =stage_type_list[[stage]][['sources.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Other_df$stage=stage;rownames(tmp_Ery2Other_df)=NULL
  Ery2Other_df=rbind(Ery2Other_df,tmp_Ery2Other_df)
  Ery2Other_df=Ery2Other_df[!is.na(Ery2Other_df$prob),]
  
  
  tmp_Ery2Ery_df=netVisual_bubble(object.list[[stage]],sources.use = stage_type_list[[stage]][['targets.use']],targets.use =stage_type_list[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Ery_df$stage=stage;rownames(tmp_Ery2Ery_df)=NULL
  Ery2Ery_df=rbind(Ery2Ery_df,tmp_Ery2Ery_df)
  Ery2Ery_df=Ery2Ery_df[!is.na(Ery2Ery_df$prob),]
  
}

Other2Ery_df$target_type='Other2Ery'
Ery2Other_df$target_type='Ery2Other'
Ery2Ery_df$target_type='Ery2Ery'

NRBC_altas_LR_df=rbind(rbind(Other2Ery_df,Ery2Other_df),Ery2Ery_df)
NRBC_altas_LR_df=NRBC_altas_LR_df[!is.na(NRBC_altas_LR_df$prob),]

NRBC_altas_LR_df$stage=factor(NRBC_altas_LR_df$stage,levels = c('YS','FL','FBM','ABM'))
table(NRBC_altas_LR_df$stage)
Other2Ery_df$stage=factor(Other2Ery_df$stage,levels = c('YS','FL','FBM','ABM'))
NRBC_altas_LR_df=NRBC_altas_LR_df[!is.na(NRBC_altas_LR_df$pathway_name),]
NRBC_altas_LR_df=NRBC_altas_LR_df[!is.na(NRBC_altas_LR_df$interaction_name),]

#---------------------在更细的细胞亚类水平subcelltype-----------------------------------#

#-------------------------------- read data----------------------------------#
YS_allsubcelltype_ccd=readRDS('NRBC_YS_altas/YS_all_subcelltype_cellchat.rds')
FL_allsubcelltype_ccd=readRDS('NRBC_FL_altas/FL_all_subcelltyp_cellchat.rds')
FBM_allsubcelltype_ccd=readRDS('NRBC_BM_altas/FBM_allsubcelltype_cellchat.rds')
ABM_allsubcelltype_ccd=readRDS('NRBC_BM_altas/ABM_allsubcelltype_cellchat.rds')

#------------------------------------重新统计 source:target, Other celltype: NRBC celltye--------------------------------------------#
# 准备，构建所有stage list 
object.list2 <- list(YS = YS_allsubcelltype_ccd,FL=FL_allsubcelltype_ccd,FBM=FBM_allsubcelltype_ccd,
                    ABM=ABM_allsubcelltype_ccd)



stage_type_list2=list()
levels(YS_allsubcelltype_ccd@idents)
levels(FL_allsubcelltype_ccd@idents)
levels(FBM_allsubcelltype_ccd@idents)
levels(ABM_allsubcelltype_ccd@idents)

# 未考虑YS NRBC 在FL、FBM中接收到的细胞互作
stage_type_list2[['YS']][['sources.use']]=c(5:length(levels(YS_allsubcelltype_ccd@idents))) # no Ery
stage_type_list2[['YS']][['targets.use']]=1:4  # Ery celltype
stage_type_list2[['FL']][['sources.use']]=c(6:length(levels(FL_allsubcelltype_ccd@idents)))
stage_type_list2[['FL']][['targets.use']]=1:5
stage_type_list2[['FBM']][['sources.use']]=c(6:length(levels(FBM_allsubcelltype_ccd@idents)))
stage_type_list2[['FBM']][['targets.use']]=1:5
stage_type_list2[['ABM']][['sources.use']]=c(6:length(levels(ABM_allsubcelltype_ccd@idents)))
stage_type_list2[['ABM']][['targets.use']]=1:5

rm(YS_allsubcelltype_ccd,FL_allsubcelltype_ccd,FBM_allsubcelltype_ccd,ABM_allsubcelltype_ccd);gc()


NRBC_LR_df2=data.frame()
Other2Ery_df2=data.frame()
Ery2Other_df2=data.frame()
Ery2Ery_df2=data.frame()


for ( stage in names(object.list2)){
  tmp_Other2Ery_df=netVisual_bubble(object.list2[[stage]],sources.use = stage_type_list2[[stage]][['sources.use']],targets.use =stage_type_list2[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Other2Ery_df$stage=stage;rownames(tmp_Other2Ery_df)=NULL
  tmp_Other2Ery_df=tmp_Other2Ery_df[!is.na(tmp_Other2Ery_df$prob),]
  
  Other2Ery_df2=rbind(Other2Ery_df2,tmp_Other2Ery_df)
  
  tmp_Ery2Other_df=netVisual_bubble(object.list2[[stage]],sources.use = stage_type_list2[[stage]][['targets.use']],targets.use =stage_type_list2[[stage]][['sources.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Other_df$stage=stage;rownames(tmp_Ery2Other_df)=NULL
  Ery2Other_df2=rbind(Ery2Other_df2,tmp_Ery2Other_df)
  Ery2Other_df2=Ery2Other_df2[!is.na(Ery2Other_df2$prob),]
  
  
  tmp_Ery2Ery_df=netVisual_bubble(object.list2[[stage]],sources.use = stage_type_list2[[stage]][['targets.use']],targets.use =stage_type_list2[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Ery_df$stage=stage;rownames(tmp_Ery2Ery_df)=NULL
  Ery2Ery_df2=rbind(Ery2Ery_df2,tmp_Ery2Ery_df)
  Ery2Ery_df2=Ery2Ery_df2[!is.na(Ery2Ery_df2$prob),]
  
}

Other2Ery_df2$target_type='Other2Ery'
Ery2Other_df2$target_type='Ery2Other'
Ery2Ery_df2$target_type='Ery2Ery'

NRBC_altas_LR_df2=rbind(rbind(Other2Ery_df2,Ery2Other_df2),Ery2Ery_df2)
NRBC_altas_LR_df2=NRBC_altas_LR_df2[!is.na(NRBC_altas_LR_df2$prob),]

NRBC_altas_LR_df2$stage=factor(NRBC_altas_LR_df2$stage,levels = c('YS','FL','FBM','ABM'))
table(NRBC_altas_LR_df2$stage)
Other2Ery_df2$stage=factor(Other2Ery_df2$stage,levels = c('YS','FL','FBM','ABM'))


setwd('NRBC_altas_CC')
dir.create('res_pic')
dir.create('res_data')
rownames(NRBC_altas_LR_df)=NULL
write.table(Other2Ery_df,file = 'res_data/Other2Ery_df_new.csv',quote = F,sep = "\t")
write.table(NRBC_altas_LR_df,file = 'res_data/NRBC_altas_LR_df_new.csv',quote = F,sep = "\t")

write.table(Other2Ery_df2,file = 'res_data/Other2Ery_df2.csv',quote = F,sep = "\t")
write.table(NRBC_altas_LR_df2,file = 'res_data/NRBC_altas_LR_df2.csv',quote = F,sep = "\t")


# 排除FL/FBM YS-NRBC 接收到的细胞互作 

NRBC_altas_LR_df=NRBC_altas_LR_df[-grep(pattern ='YS' ,NRBC_altas_LR_df$target),]

all_incoming_LRs=unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Other2Ery'])
all_outgoing_LRs=unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other'])
all_Ery2Ery_LRs=unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Ery'])

all_incoming_LRs2=unique(NRBC_altas_LR_df2$interaction_name[NRBC_altas_LR_df2$target_type=='Other2Ery'])
all_outgoing_LRs2=unique(NRBC_altas_LR_df2$interaction_name[NRBC_altas_LR_df2$target_type=='Ery2Other'])
all_Ery2Ery_LRs2=unique(NRBC_altas_LR_df2$interaction_name[NRBC_altas_LR_df2$target_type=='Ery2Ery'])

length(all_incoming_LRs) # 145
table(all_incoming_LRs2 %in% all_incoming_LRs) # F/T:6/144 
length(all_outgoing_LRs)# 164
table(all_outgoing_LRs2 %in% all_outgoing_LRs) # F/T:14/164
length(all_Ery2Ery_LRs) # 24 
table(all_Ery2Ery_LRs2 %in% all_Ery2Ery_LRs) # F/T:0/24 ? why, ABM altas 略有差异导致，更新ABM ALTAS细胞
rm(object.list2);gc()

# 粗颗粒细胞分类在更细颗粒细胞分类中的与NRBC细胞互作绝大部分都存在，少量未有发现的，基于表达水平,可以不考虑
all_incoming_LRs[!all_incoming_LRs %in%  all_incoming_LRs2]#  "WNT5A_FZD3"# 保存YS NRBC在其他niche中接收到的信号，则多 "RBP4_STRA6"，但是STRA6表达较低

NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% c("WNT5A_FZD3", "RBP4_STRA6"),]# MESOTHELIUM target BFUE/CFUE, HEPATOCYTE target YS_Bas/Poly
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features = c('FZD3','STRA6'),stack=T)
NRBC_altas_LR_df[NRBC_altas_LR_df$receptor=='FZD3',]# FZD3 主要接受ligand WNT5B 的作用

rm(stage_type_list,object.list);gc()


# 更细的类可以得到哪类细胞是主要ligand的提供者抛弃粗分类，直接采用细分类进行
NRBC_altas_LR_df=NRBC_altas_LR_df2
NRBC_altas_LR_df$main_type='no'
NRBC_altas_LR_df$main_type[NRBC_altas_LR_df$target_type=='Other2Ery' & NRBC_altas_LR_df$interaction_name %in% all_incoming_LRs]='yes'
NRBC_altas_LR_df$main_type[NRBC_altas_LR_df$target_type=='Ery2Other' & NRBC_altas_LR_df$interaction_name %in% all_outgoing_LRs]='yes'
NRBC_altas_LR_df$main_type[NRBC_altas_LR_df$target_type=='Ery2Ery' & NRBC_altas_LR_df$interaction_name %in% all_Ery2Ery_LRs]='yes'

NRBC_altas_LR_df$celltype=NRBC_altas_LR_df$target
NRBC_altas_LR_df$celltype[NRBC_altas_LR_df$target_type=='Ery2Other']=NRBC_altas_LR_df$source[NRBC_altas_LR_df$target_type=='Ery2Other']
NRBC_altas_LR_df$celltype=factor(NRBC_altas_LR_df$celltype,levels = NRBC_subcelltype)

NRBC_celltype_mexp_df=data.frame(AverageExpression(filt_NBRC_altas_seu,  layer = 'data',group.by = 'source_celltype')$RNA)
# 去掉不同niche中 NRBC LR 表达均值低于0.1的基因 对应的LR
for (stage in c('YS','FL','FBM','ABM')) {
  temp_LRs=NRBC_altas_LR_df[NRBC_altas_LR_df$stage==stage,]
  LR_genes=unique(c(CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name  %in% unique(temp_LRs[temp_LRs$target_type=='Other2Ery','interaction_name']),'receptor.symbol'],
                    CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name  %in% unique(temp_LRs[temp_LRs$target_type=='Ery2Other','interaction_name']),'ligand.symbol']))
  
  LR_genes=unique(unlist(strsplit(LR_genes,split=', ')))
  LR_gene_exp_test=rowSums(NRBC_celltype_mexp_df[LR_genes,grep(stage,colnames(NRBC_celltype_mexp_df))] >0.1)
  print(stage)
  print(names(LR_gene_exp_test)[LR_gene_exp_test <1])
}
 
p=DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',scale = F,features = c("IGF1R", "EPHB4", "PDGFC", "DAGLB","TNFRSF17" , "TNFRSF13B", "SDC1","PECAM1","BMP6","VEGFA","IGF1","MDK","GAS6"))+RotatedAxis()
# 去掉ABM 中"TNFRSF17" , "TNFRSF13B", "SDC1",BMP6,IGF1,MDK,GAS6 LR, 由于ABM NRBC中去掉了IGK1 表达的NRBC，导致这些基因低表达；

exp_pct_df=data.frame(p$data)
exp_pct_df=dcast(exp_pct_df[,c('id','features.plot','pct.exp')],features.plot~id)
rownames(exp_pct_df)=exp_pct_df$features.plot;exp_pct_df=exp_pct_df[,-1]

# 去掉表达均值低于0.1 且表达细胞比例低于0.1的基因
NRBC_altas_LR_df1=NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='ABM' & NRBC_altas_LR_df$target_type=='Ery2Ery',]
temp_NRBC_altas_inLR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='ABM' & NRBC_altas_LR_df$target_type=='Other2Ery' ,]
temp_NRBC_altas_outLR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='ABM' & NRBC_altas_LR_df$target_type=='Ery2Other' ,]

temp_NRBC_altas_inLR_df=temp_NRBC_altas_inLR_df[ -grep("TNFRSF17|TNFRSF13B|SDC1|BMP6|IGF1|MDK|GAS6",temp_NRBC_altas_inLR_df$receptor),]
temp_NRBC_altas_outLR_df=temp_NRBC_altas_outLR_df[ -grep("TNFRSF17|TNFRSF13B|SDC1|BMP6|IGF1|MDK|GAS6",temp_NRBC_altas_outLR_df$ligand),]
NRBC_altas_LR_df1=NRBC_altas_LR_df1[ -grep("TNFRSF17|TNFRSF13B|SDC1|BMP6|IGF1|MDK|GAS6",NRBC_altas_LR_df1$interaction_name),]
NRBC_altas_LR_df1=rbind(NRBC_altas_LR_df1,rbind(temp_NRBC_altas_inLR_df,temp_NRBC_altas_outLR_df))

NRBC_altas_LR_df=rbind(NRBC_altas_LR_df[NRBC_altas_LR_df$stage!='ABM',],NRBC_altas_LR_df1)

NRBC_altas_LR_df$celltype=NRBC_altas_LR_df$target
NRBC_altas_LR_df$celltype[NRBC_altas_LR_df$target_type=='Ery2Other']=NRBC_altas_LR_df$source[NRBC_altas_LR_df$target_type=='Ery2Other']
NRBC_altas_LR_df$celltype=factor(NRBC_altas_LR_df$celltype,levels = c('BFUE/CFUE','ProE','Bas','Poly','Orth'))

Other2Ery_df=NRBC_altas_LR_df[NRBC_altas_LR_df$target_type=='Other2Ery',]

write.table(NRBC_altas_LR_df,file = 'res_data/filt_NRBC_altas_LR_df.csv',quote = F,sep="\t")
write.table(Other2Ery_df,file = 'res_data/filt_Other2Ery_df_new.csv',quote = F,sep="\t")




##############################################################################################################################################################################
#------------------------------------part2: analysis   interactioin  across niches & primitive vs definitiive----------------------#
##############################################################################################################################################################################

NRBC_altas_LR_df=read.csv('res_data/filt_NRBC_altas_LR_df.csv',sep="\t")
Other2Ery_df=read.csv('res_data/filt_Other2Ery_df_new.csv',sep="\t")


# ----------------------incoming LR--------------------------#
the_incoming_LR_list=list('YS'=sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='YS'  & NRBC_altas_LR_df$target_type =='Other2Ery','interaction_name'])),
                            'FL'= sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='FL'  & NRBC_altas_LR_df$target_type =='Other2Ery','interaction_name'])),
                            'FBM'=sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='FBM'& NRBC_altas_LR_df$target_type =='Other2Ery','interaction_name'])),
                            'ABM'=sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='ABM'& NRBC_altas_LR_df$target_type =='Other2Ery','interaction_name']))
)
  
library(VennDiagram)
p=venn.diagram(x =the_incoming_LR_list,filename = NULL,fill=cols[1:4], scaled = T ,main ='incoming LR')
dev.off()
grid.draw(p)
  
veen_res=get.venn.partitions(the_incoming_LR_list)
incoming_LR_type_df=data.frame(LR=unique(unlist(the_incoming_LR_list)))
incoming_LR_type_df$type='un'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[1]]  ]='01_shared' # 4,3,2,1
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[2]]  ]='02_definitive'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[3]]  ]='04_no_FL'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[4]]  ]='10_BM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[5]]  ]='05_no_FBM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[6]]  ]='09_FL_ABM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[7]]  ]='15_YS_ABM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[8 ]] ]='14_ABM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[9 ]] ]='03_fetal_all'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[10]] ]='08_fetal_FLFBM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[11]] ]='07_YS_FBM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[12]] ]='13_FBM'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[13]] ]='06_YS_FL'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[14]] ]='12_FL'
incoming_LR_type_df$type[incoming_LR_type_df$LR %in% veen_res$..values..[[15]] ]='11_YS'
incoming_LR_type_df_list=split(incoming_LR_type_df$LR,incoming_LR_type_df$type)

  
  # outgoing LR
the_outgoing_LR_list=list('YS'= sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='YS'  & NRBC_altas_LR_df$target_type =='Ery2Other','interaction_name'])),
                            'FL'= sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='FL'  & NRBC_altas_LR_df$target_type =='Ery2Other','interaction_name'])),
                            'FBM'=sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='FBM'& NRBC_altas_LR_df$target_type =='Ery2Other','interaction_name'])),
                            'ABM'=sort(unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='ABM'& NRBC_altas_LR_df$target_type =='Ery2Other','interaction_name']))
)
  
library(VennDiagram)
p=venn.diagram(x =the_outgoing_LR_list,filename = NULL,fill=cols[1:4], scaled = T ,main = 'outgoning LR')
dev.off()
grid.draw(p)

veen_res=get.venn.partitions(the_outgoing_LR_list)
outgoing_LR_type_df=data.frame(LR=unique(unlist(the_outgoing_LR_list)))
outgoing_LR_type_df$type='un'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[1]]  ]='01_shared' # 4,3,2,1
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[2]]  ]='02_definitive'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[3]]  ]='04_no_FL'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[4]]  ]='10_BM'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[5]]  ]='05_no_FBM'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[6]]  ]='09_FL_ABM'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[8 ]] ]='14_ABM'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[9 ]] ]='03_fetal_all'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[10]] ]='08_fetal_FLFBM'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[11]] ]='07_YS_FBM'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[12]] ]='13_FBM'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[13]] ]='06_YS_FL'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[14]] ]='12_FL'
outgoing_LR_type_df$type[outgoing_LR_type_df$LR %in% veen_res$..values..[[15]] ]='11_YS'
outgoing_LR_type_df_list=split(outgoing_LR_type_df$LR,outgoing_LR_type_df$type)
  

saveRDS(incoming_LR_type_df_list,file = 'res_data/incoming_LR_type_df_list.rds')
saveRDS(outgoing_LR_type_df_list,file = 'res_data/outgoing_LR_type_df_list.rds')



################################################################################################################################################################
#################-----------------------check the expression level in conserved LR -------------------------------------------------------#
################################################################################################################################################################


setwd('NRBC_altas_CC')

NRBC_subcelltype=c("BFUE/CFUE","ProE","Bas","Poly","Orth" )

if(F){
    NRBC_altas_LR_df=read.csv('res_data/filt_NRBC_altas_LR_df.csv',sep="\t")
    NRBC_altas_LR_df$stage=factor(NRBC_altas_LR_df$stage,levels = c('YS','FL','FBM','ABM'))
    NRBC_altas_LR_df$celltype=factor(NRBC_altas_LR_df$celltype,levels = NRBC_subcelltype)
    Other2Ery_df=NRBC_altas_LR_df[NRBC_altas_LR_df$target_type=='Other2Ery',]
    incoming_LR_type_df_list=readRDS('res_data/incoming_LR_type_df_list.rds')
    outgoing_LR_type_df_list=readRDS('res_data/outgoing_LR_type_df_list.rds')
}

filt_NBRC_altas_seu=readRDS('../20251125_filt_NBRC_altas_seu.rds')
Idents(filt_NBRC_altas_seu)='source_celltype'
filt_NBRC_altas_seu$pd_celltype='definitive'
filt_NBRC_altas_seu$pd_celltype[grep('YS',filt_NBRC_altas_seu$tissue_stage)]='primitive'

cho_cells=c(sample(rownames(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='YS',]),size = 15000,replace = F),
            sample(rownames(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='FL',]),size = 15000,replace = F),
            sample(rownames(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='FBM',]),size = 15000,replace = F),
            sample(rownames(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='ABM',]),size = 15000,replace = F))
filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,cells=cho_cells)
rm(cho_cells)


# whole level
pd_whole_level_markers_degs_df=readRDS('../Protein_NRBC_marker/res_data/pd_whole_level_markers.rds')
fetal_adult_NRBC_whole_marker=readRDS('../Protein_NRBC_marker/res_data/fetal_adult_NRBC_whole_marker.rds')
# subcelltype level
sub_pd_all_Ery_tissue_markers=read.csv('../Protein_NRBC_marker/DE_marker/primitive_definitive_all_Ery_RNA_markers.csv')
sub_pd_all_Ery_tissue_markers=sub_pd_all_Ery_tissue_markers[,-1]
sub_fetal_adult_all_Ery_tissue_markers=read.csv('../Protein_NRBC_marker/DE_marker/fetal_adult_all_Ery_RNA_markers.csv')
sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[,-1]


all_NRBC_receptor_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']),'receptor.symbol'] # 取receptor gene
all_NRBC_ligand_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in%unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']),'ligand.symbol']# 取ligand gene
all_NRBC_receptor_genes=sort(unique(unlist(strsplit(all_NRBC_receptor_genes,split=','))))
all_NRBC_ligand_genes=sort(unique(unlist(strsplit(all_NRBC_ligand_genes,split=','))))
length(all_NRBC_receptor_genes);length(all_NRBC_ligand_genes)# 49， 80

all_NRBC_LR_genes=unique(c(all_NRBC_receptor_genes,all_NRBC_ligand_genes));length(all_NRBC_LR_genes) # 120


# conserved LR genes: expression level in sub_celltype >0.1 & in detected in all niche
LR_celltype_mexp_df2=DotPlot(filt_NBRC_altas_seu,features =unique(gsub(pattern = ' ',replacement = '',unique(c(all_NRBC_ligand_genes,all_NRBC_receptor_genes)))))
LR_celltype_mexp_df2=LR_celltype_mexp_df2$data

conserved_LR_genes=c()
for(gene in rownames(LR_celltype_mexp_df2)){
  temp=LR_celltype_mexp_df2[LR_celltype_mexp_df2$features.plot==gene,]
  temp=temp[temp$avg.exp>=0.1 &  temp$pct.exp>=10,'id' ]
  if(length(temp) <1){next}
  temp=unique(as.character(t(data.frame(strsplit(as.character(temp),split='_')))[,1]))
  if(length(temp) >3){
    conserved_LR_genes=c(conserved_LR_genes,gene)
  }
}
length(conserved_LR_genes)# 46 
saveRDS(conserved_LR_genes,file = 'conserved_LR_genes.rds')
conserved_LR_genes;length(conserved_LR_genes)# 48,46 CD74 表达细胞比例在临介附近


conserved_LR_genes


shared_receptor_genes=conserved_LR_genes[conserved_LR_genes %in% unlist(strsplit(CellChatDB.human$interaction$receptor.symbol,', '))]
shared_receptor_genes
# ## CD74 与其他receptor CD44，CXCR4/2 ，非保守基因一块才能，接收外界信息,注意要删除,
CellChatDB.human$interaction$receptor[grep('CD74',CellChatDB.human$interaction$interaction_name)]
# ITGA4_ITGB1 要成组合，"ITGA4, ITGB1",  "IFNGR2" "IFNGR1" 要组合"IFNGR1, IFNGR2"
CellChatDB.human$interaction$receptor.symbol[grep('ITGA4_ITGB1', CellChatDB.human$interaction$interaction_name)]
shared_NRBC_cLR_inLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$receptor.symbol %in% c(shared_receptor_genes, "ITGA4, ITGB1","IFNGR1, IFNGR2") ])
shared_NRBC_cLR_inLRs=shared_NRBC_cLR_inLRs[shared_NRBC_cLR_inLRs %in% NRBC_altas_LR_df$interaction_name]
shared_all_NRBC_cLR_LRs_df=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% shared_NRBC_cLR_inLRs & NRBC_altas_LR_df$target_type!='Ery2Other',]

shared_ligand_genes=conserved_LR_genes[conserved_LR_genes %in% unlist(strsplit(CellChatDB.human$interaction$ligand.symbol,', '))]
shared_ligand_genes # ITGA4_ITGB1 要成组合，"ITGA4, ITGB1"

shared_NRBC_cLR_outLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$ligand.symbol %in% c(shared_ligand_genes,"ITGA4, ITGB1") ])
shared_NRBC_cLR_outLRs=shared_NRBC_cLR_outLRs[shared_NRBC_cLR_outLRs %in% NRBC_altas_LR_df$interaction_name]
shared_all_NRBC_cLR_LRs_df=rbind(shared_all_NRBC_cLR_LRs_df,NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% shared_NRBC_cLR_outLRs & NRBC_altas_LR_df$target_type!='Other2Ery',])

left_genes=c(shared_receptor_genes,shared_ligand_genes)[!c(shared_receptor_genes,shared_ligand_genes) %in% c(CellChatDB.human$interaction$ligand.symbol,CellChatDB.human$interaction$receptor.symbol) ]
unique(left_genes)
# "ITGA4" "ITGB1"  "IFNGR2" "IFNGR1" 
CellChatDB.human$interaction$receptor.symbol[grep('ITGA4', CellChatDB.human$interaction$interaction_name)]
CellChatDB.human$interaction$receptor.symbol[grep('IFNGR2', CellChatDB.human$interaction$interaction_name)]


temp_df=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% shared_NRBC_cLR_inLRs & NRBC_altas_LR_df$target_type!='Ery2Other' ,]
temp_df=temp_df[order(temp_df$receptor,temp_df$pathway_name),]
temp_df$interaction_name=factor(temp_df$interaction_name,levels = unique(temp_df$interaction_name))
temp_df$pathway_name=factor(temp_df$pathway_name,levels = unique(temp_df$pathway_name))

#temp_df=temp_df[-grep('DHCR7|DHCR24|AKR1C3|TBXAS1',temp_df$interaction_name),]
p11=ggplot(temp_df,aes(x=celltype,y=interaction_name,size=prob,color=pathway_name,shape=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'incoming interaction based on conserved receptors')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])

fre_df=data.frame(table(unique(temp_df[,c('stage','celltype','interaction_name','annotation')])[,c('stage','celltype','annotation')]))
fre_df$celltype=factor(fre_df$celltype,levels = NRBC_subcelltype)
p12=ggplot(fre_df,aes(x=celltype,y=Freq,fill=annotation))+geom_bar(stat = 'identity')+facet_grid(~stage )+theme_classic()+ ylab(label = 'count')+ggtitle('incoming interaction')+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+NoLegend()+ylim(c(0,60))

temp_df=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in%shared_NRBC_cLR_outLRs  & NRBC_altas_LR_df$target_type=='Ery2Other' ,]
temp_df=temp_df[order(temp_df$receptor,temp_df$pathway_name),]
temp_df$interaction_name=factor(temp_df$interaction_name,levels = unique(temp_df$interaction_name))
temp_df$pathway_name=factor(temp_df$pathway_name,levels = unique(temp_df$pathway_name))
#temp_df=temp_df[-grep('DHCR7|DHCR24|AKR1C3|TBXAS1',temp_df$interaction_name),]
p22=ggplot(temp_df,aes(x=celltype,y=interaction_name,size=prob,color=pathway_name,shape=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'outgoing interaction based on conserved ligands')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])

fre_df=data.frame(table(unique(temp_df[,c('stage','celltype','interaction_name','annotation')])[,c('stage','celltype','annotation')]))
fre_df$celltype=factor(fre_df$celltype,levels = NRBC_subcelltype)
p21=ggplot(fre_df,aes(x=celltype,y=Freq,fill=annotation))+geom_bar(stat = 'identity')+facet_grid(~stage )+theme_classic()+ ylab(label = 'count')+ggtitle(label = 'outgoing interaction')+
  scale_fill_manual(values = cols[c(1:2,4)])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ylim(c(0,60))

p=p11+p22;p
ggsave(p,file='res_pic/main_figure4/pd_conserved_incoming_outing_LR_based_conserved_LR.pdf',width =20 ,height =12 ,dpi = 300)

p=p12+p21;p
ggsave(p,file='res_pic/main_figure4/pd_frequency_conserved_incoming_outing_LR_based_conserved_LR.pdf',width =15 ,height =6 ,dpi = 300)


temp_deg_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$avg_log2FC>0,]
temp_deg_df=temp_deg_df[temp_deg_df$gene %in% shared_receptor_genes, ]
temp_deg_df[shared_receptor_genes[!shared_receptor_genes %in% temp_deg_df$gene] ,1:7]=0
temp_deg_df[shared_receptor_genes[!shared_receptor_genes %in% temp_deg_df$gene] ,'gene']=shared_receptor_genes[!shared_receptor_genes %in% temp_deg_df$gene]  
temp_deg_df=temp_deg_df[order(temp_deg_df$avg_log2FC,decreasing = T),]
temp_deg_df$gene=factor(temp_deg_df$gene,levels = temp_deg_df$gene)
p1=ggplot(temp_deg_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+NoLegend()+
  ggtitle('conserved receptors')
p1
p3=VlnPlot(filt_NBRC_altas_seu,features = temp_deg_df$gene,group.by = 'final_celltype',stack = T,split.by = 'pd_celltype')
p=p1+p3+plot_layout(ncol =1,heights = c(0.6,1.2))  ;p
ggsave(p,file='res_pic/main_figure4/pd_conserved_incoming_receptors_expression.pdf',width =8 ,height =3 ,dpi = 300)


sub_pd_all_Ery_tissue_markers$celltype=factor(sub_pd_all_Ery_tissue_markers$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))
sub_temp_deg_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC>0 &  sub_pd_all_Ery_tissue_markers$gene %in% shared_receptor_genes,]
sub_temp_deg_df$gene=factor(sub_temp_deg_df$gene,levels = temp_deg_df$gene)
p2=ggplot(sub_temp_deg_df,aes(x=gene,y=celltype,col=cluster,size=avg_log2FC,alpha=avg_log2FC))+geom_point(stat = 'identity')+theme_classic()+
  ggtitle('conserved receptors')+xlab(label = '') + theme(axis.text.x = element_text(angle = 45,hjust = 1))#axis.text.x = element_blank()

ggsave(p2,file='res_pic/main_figure4/pd_conserved_incoming_receptors_DE_substage.pdf',width =8 ,height =3 ,dpi = 300)



temp_deg_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$avg_log2FC>0,]
temp_deg_df=temp_deg_df[temp_deg_df$gene %in% shared_ligand_genes, ]
temp_deg_df[shared_ligand_genes[!shared_ligand_genes %in% temp_deg_df$gene] ,1:7]=0
temp_deg_df[shared_ligand_genes[!shared_ligand_genes %in% temp_deg_df$gene] ,'gene']=shared_ligand_genes[!shared_ligand_genes %in% temp_deg_df$gene]  
temp_deg_df=temp_deg_df[order(temp_deg_df$avg_log2FC,decreasing = T),]
temp_deg_df$gene=factor(temp_deg_df$gene,levels = temp_deg_df$gene)
p1=ggplot(temp_deg_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+NoLegend()+ggtitle('conserved ligands')

p3=VlnPlot(filt_NBRC_altas_seu,features = temp_deg_df$gene,group.by = 'final_celltype',stack = T,split.by = 'pd_celltype')
p=p1+p3+plot_layout(ncol =1,heights = c(0.6,1.2))  ;p

ggsave(p,file='res_pic/main_figure4/pd_conserved_outgoing_ligands_expression.pdf',width =8 ,height =6 ,dpi = 300)


sub_pd_all_Ery_tissue_markers$celltype=factor(sub_pd_all_Ery_tissue_markers$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))
sub_temp_deg_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC>0 &  sub_pd_all_Ery_tissue_markers$gene %in% levels(temp_deg_df$gene),]
sub_temp_deg_df$gene=factor(sub_temp_deg_df$gene,levels = levels(temp_deg_df$gene))
p2=ggplot(sub_temp_deg_df,aes(x=gene,y=celltype,col=cluster,size=avg_log2FC,alpha=avg_log2FC))+geom_point(stat = 'identity')+theme_classic()+
  ggtitle('conserved ligands')+xlab(label = '')+theme(axis.text.x = element_text(angle = 45,hjust = 1))#axis.text.x = element_blank()
p2
ggsave(p2,file='res_pic/main_figure4/pd_conserved_incoming_ligands_DE_substage.pdf',width =8 ,height =3 ,dpi = 300)



################################################################################################################################################################
#################-----------------------check the expression level of LR in DEGs -------------------------------------------------------#
###############################################################################################################################################################


#  defintiive conserved LR gene: expression level in sub_celltype >0.1 & pct>10 in detected in definitive niche
LR_celltype_mexp_df=DotPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),features =unique(gsub(pattern = ' ',replacement = '',unique(c(all_NRBC_ligand_genes,all_NRBC_receptor_genes)))))
LR_celltype_mexp_df=LR_celltype_mexp_df$data
LR_celltype_mexp_df$avg.exp=round(LR_celltype_mexp_df$avg.exp,2)
saveRDS(LR_celltype_mexp_df,file = 'res_data/defintive_LR_celltype_mexp_df.rds')

definitive_conserved_LR_genes=c()
for(gene in rownames(LR_celltype_mexp_df)){
  temp=LR_celltype_mexp_df[LR_celltype_mexp_df$features.plot==gene,]
  temp=temp[temp$avg.exp>=0.1 &  temp$pct.exp>=10,'id' ]
  if(length(temp) <1){next}
  temp=unique(as.character(t(data.frame(strsplit(as.character(temp),split='_')))[,1]))
  if(length(temp) >2){
    definitive_conserved_LR_genes=c(definitive_conserved_LR_genes,gene)
  }
}
length(definitive_conserved_LR_genes)# 73
saveRDS(definitive_conserved_LR_genes,file = 'definitive_conserved_LR_genes.rds')


table(conserved_LR_genes %in% definitive_conserved_LR_genes)#T: 46
length(definitive_conserved_LR_genes[!definitive_conserved_LR_genes %in% conserved_LR_genes])# other:27 


# YS specific LR 
incoming_LR_type_df_list[['11_YS']] # 都在hDEG-LR中
# ---------------------重头分析挑选 primitive vs definitive-----------------------------------#

cadidated_pd_markers=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$avg_log2FC >1 & pd_whole_level_markers_degs_df$pct.1>0.1,]
dim(cadidated_pd_markers)# 3463

all_NRBC_receptor_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']),'receptor.symbol'] # 取receptor gene
all_NRBC_ligand_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in%unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']),'ligand.symbol']# 取ligand gene
all_NRBC_receptor_genes=sort(unique(unlist(strsplit(all_NRBC_receptor_genes,split=','))))
all_NRBC_ligand_genes=sort(unique(unlist(strsplit(all_NRBC_ligand_genes,split=','))))
length(all_NRBC_receptor_genes);length(all_NRBC_ligand_genes)# 49， 80

# Ery2Ery lignad 全部在ligand基因中
Ery2Ery_NRBC_LR_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Ery']),'ligand.symbol']# 取ligand gene
Ery2Ery_NRBC_LR_genes=sort(unique(unlist(strsplit(Ery2Ery_NRBC_LR_genes,split=','))))
table(Ery2Ery_NRBC_LR_genes %in% all_NRBC_ligand_genes )# 15

# 得到所有的高变化差异LR 基因
all_NRBC_LR_genes=unique(c(all_NRBC_receptor_genes,all_NRBC_ligand_genes));length(all_NRBC_LR_genes) # 只有NRBC LR 120个基因
table(cadidated_pd_markers$gene %in% all_NRBC_LR_genes) # T:43

cadidated_fa_markers=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$avg_log2FC >2 & fetal_adult_NRBC_whole_marker$pct.2<0.1 & fetal_adult_NRBC_whole_marker$pct.1>0.05,]
table(cadidated_fa_markers$gene %in% all_NRBC_LR_genes) # T:34
candidated_lr_genes=all_NRBC_LR_genes[ all_NRBC_LR_genes %in% c(cadidated_fa_markers$gene,cadidated_pd_markers$gene)]
length(candidated_lr_genes)# 68 mDEGs in  NRBC-LR
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =candidated_lr_genes,stack = T,split.by = 'tissue_stage',cols = cols)

#--------先分析 primitive vs definitive -------------
# 得到高特异LR 基因；
pd_candidated_lr_genes=cadidated_pd_markers[cadidated_pd_markers$gene %in% all_NRBC_LR_genes  & cadidated_pd_markers$pct.2 <0.2,]
pd_candidated_lr_genes[pd_candidated_lr_genes$gene %in% conserved_LR_genes,] # SELPLG TGFB1  NAMPT  HLA-B  ITGA4 
pd_candidated_lr_genes=pd_candidated_lr_genes[!pd_candidated_lr_genes$gene %in% conserved_LR_genes ,]
pd_candidated_lr_genes$gene[pd_candidated_lr_genes$cluster=='definitive'][!pd_candidated_lr_genes$gene[pd_candidated_lr_genes$cluster=='definitive'] %in%definitive_conserved_LR_genes ]

pd_candidated_lr_genes=pd_candidated_lr_genes[order(pd_candidated_lr_genes$cluster,pd_candidated_lr_genes$avg_log2FC,decreasing = T),]
pd_candidated_lr_genes$gene=factor(pd_candidated_lr_genes$gene,levels = pd_candidated_lr_genes$gene)
pd_candidated_hDEG_lr_genes=pd_candidated_lr_genes

p1=ggplot(pd_candidated_hDEG_lr_genes,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()
p1

sub_temp_deg_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC>0 &  sub_pd_all_Ery_tissue_markers$gene %in% pd_candidated_lr_genes$gene,]
sub_temp_deg_df$gene=factor(sub_temp_deg_df$gene,levels = pd_candidated_hDEG_lr_genes$gene)
sub_temp_deg_df$type='ligand'
sub_temp_deg_df$type[sub_temp_deg_df$gene %in% all_NRBC_receptor_genes]='receptor'
sub_temp_deg_df$type[sub_temp_deg_df$gene %in% intersect(all_NRBC_receptor_genes,all_NRBC_ligand_genes)]='Both'
sub_temp_deg_df$gene=factor(sub_temp_deg_df$gene,levels =pd_candidated_hDEG_lr_genes$gene )
p2=ggplot(sub_temp_deg_df,aes(x=gene,y=celltype,col=cluster,size=avg_log2FC,alpha=avg_log2FC,shape=type))+geom_point(stat = 'identity')+theme_classic()+
  ggtitle('hDEG-LR genes between primitive and defintive NRBC')+xlab(label = '')+theme(axis.text.x = element_text(angle = 45,hjust = 1))#axis.text.x = element_blank()
p2
ggsave(p2,file='res_pic/main_figure4/pd_hDEG_LR_expression_DE_substage.pdf',width =8 ,height =3 ,dpi = 300)


p3=VlnPlot(filt_NBRC_altas_seu,features = levels(sub_temp_deg_df$gene),group.by = 'final_celltype',stack = T,split.by = 'pd_celltype')
p=p1+p3+plot_layout(ncol =1,heights = c(0.6,1.2)) ;p
ggsave(p,file='res_pic/main_figure4/pd_hDEG_LR_expression.pdf',width =8 ,height =6 ,dpi = 300)



pd_candidated_lr_genes$gene_type='receptor'
pd_candidated_lr_genes$gene_type[pd_candidated_lr_genes$gene %in% all_NRBC_ligand_genes ]='ligand'
pd_candidated_lr_genes$gene[!pd_candidated_lr_genes$gene %in% c(CellChatDB.human$interaction$ligand.symbol ,CellChatDB.human$interaction$receptor.symbol )]# NO ,所有基因都在数据库中

pd_NRBC_mDEGs_inLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$receptor.symbol %in% pd_candidated_lr_genes$gene[pd_candidated_lr_genes$gene_type=='receptor']])
pd_NRBC_mDEGs_inLRs=pd_NRBC_mDEGs_inLRs[pd_NRBC_mDEGs_inLRs %in% NRBC_altas_LR_df$interaction_name]
pd_all_NRBC_mDEGs_LRs_df=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% pd_NRBC_mDEGs_inLRs & NRBC_altas_LR_df$target_type!='Ery2Other',]

pd_all_NRBC_mDEGs_LRs_df$target_type=factor(pd_all_NRBC_mDEGs_LRs_df$target_type,levels = c('Other2Ery','Ery2Ery','Ery2Other'))
pd_all_NRBC_mDEGs_LRs_df$cluster='definitive'
pd_all_NRBC_mDEGs_LRs_df$cluster[pd_all_NRBC_mDEGs_LRs_df$receptor %in% pd_candidated_lr_genes$gene[pd_candidated_lr_genes$cluster=='primitive'] | pd_all_NRBC_mDEGs_LRs_df$ligand %in% pd_candidated_lr_genes$gene[pd_candidated_lr_genes$cluster=='primitive']  ]='primitive'

pd_all_NRBC_mDEGs_LRs_df=pd_all_NRBC_mDEGs_LRs_df[order(pd_all_NRBC_mDEGs_LRs_df$cluster,pd_all_NRBC_mDEGs_LRs_df$target_type,pd_all_NRBC_mDEGs_LRs_df$annotation,pd_all_NRBC_mDEGs_LRs_df$interaction_name),]
pd_all_NRBC_mDEGs_LRs_df$interaction_name=factor(pd_all_NRBC_mDEGs_LRs_df$interaction_name,levels =unique( pd_all_NRBC_mDEGs_LRs_df$interaction_name))
ggplot(pd_all_NRBC_mDEGs_LRs_df,aes(x=celltype ,y=interaction_name,size=prob ,color=annotation,shape=target_type))+geom_point(alpha=0.7)+facet_grid(~stage )+theme_classic()+
  scale_color_manual(values = cols[-2])+scale_shape_manual(values = c(0:2,5))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+facet_grid(~stage )+ggtitle('primitive vs definitive LR based on hDEGs' )

order_pd_LRs=c()
for (gene in levels(pd_candidated_lr_genes$gene)) {
  order_pd_LRs=c(order_pd_LRs,unique(as.character(pd_all_NRBC_mDEGs_LRs_df$interaction_name[grep(gene,pd_all_NRBC_mDEGs_LRs_df$interaction_name)])))
    
};order_pd_LRs=unique(order_pd_LRs)

pd_all_NRBC_mDEGs_LRs_df$interaction_name=factor(pd_all_NRBC_mDEGs_LRs_df$interaction_name,levels =order_pd_LRs[length(order_pd_LRs):1])

grep('Cholesterol|Desmosterol',order_pd_LRs)
order_pd_LRs=c(order_pd_LRs[1:5],order_pd_LRs[c(12:14,16)],order_pd_LRs[ !1:length(order_pd_LRs) %in% c(1:5,12:14,16)])
grep('HLA',order_pd_LRs)
order_pd_LRs=c(order_pd_LRs[c(1:35,37)],order_pd_LRs[ !1:length(order_pd_LRs) %in% c(1:35,37)])

pd_all_NRBC_mDEGs_LRs_df$interaction_name=factor(pd_all_NRBC_mDEGs_LRs_df$interaction_name,levels =order_pd_LRs[length(order_pd_LRs):1])
pd_all_NRBC_mDEGs_LRs_df=pd_all_NRBC_mDEGs_LRs_df[!(pd_all_NRBC_mDEGs_LRs_df$interaction_name %in% c("Cholesterol-Cholesterol-LIPA_RORA" ,"Cholesterol-Cholesterol-DHCR7_RORA")  & pd_all_NRBC_mDEGs_LRs_df$stage=='ABM'),]
p22=ggplot(pd_all_NRBC_mDEGs_LRs_df,aes(x=celltype ,y=interaction_name,size=prob ,color=annotation,shape=target_type))+geom_point(alpha=0.7)+facet_grid(~stage )+
       theme_classic()+scale_color_manual(values = cols[-2])+scale_shape_manual(values = c(0:2,5))+
     theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle('primitive vs definitive LR based on hDEGs' )
p22
ggsave(p22,file='res_pic/main_figure4/pd_LR_based_hDEG.pdf',width =10 ,height =12 ,dpi = 300)

# definitive NRBC specific vs primitiive, other definitive conserved LR gene, other primitive specific genes 

new_definitive_conserved_LR_genes=definitive_conserved_LR_genes[!definitive_conserved_LR_genes %in% c(conserved_LR_genes,as.character(pd_candidated_lr_genes$gene)) ]
length(new_definitive_conserved_LR_genes) # 14 
saveRDS(new_definitive_conserved_LR_genes,file = 'res_data/new_definitive_conserved_LR_genes.rds')

# ealry_stage  expression level
sub_pd_specific_definitive_genes=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$gene %in% new_definitive_conserved_LR_genes & sub_pd_all_Ery_tissue_markers$avg_log2FC>0,]
sub_pd_specific_definitive_genes=sub_pd_specific_definitive_genes[order(sub_pd_specific_definitive_genes$avg_log2FC,decreasing = T),]
sub_pd_specific_definitive_genes=sub_pd_specific_definitive_genes[sub_pd_specific_definitive_genes$celltype=='early_Ery',]
sub_pd_specific_definitive_genes$gene=factor(sub_pd_specific_definitive_genes$gene,levels = unique(sub_pd_specific_definitive_genes$gene))
p=ggplot(sub_pd_specific_definitive_genes,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()
ggsave(p,file='res_pic/main_figure4/definitive_vs_primtive_specific_conserved_LR_early_DE.pdf',width =6 ,height =3 ,dpi = 300)


# whole  expression level
pd_specific_definitive_genes=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$gene %in% new_definitive_conserved_LR_genes & pd_whole_level_markers_degs_df$avg_log2FC>0,]
pd_specific_definitive_genes[new_definitive_conserved_LR_genes[!new_definitive_conserved_LR_genes %in% pd_specific_definitive_genes$gene],c('gene')]=new_definitive_conserved_LR_genes[!new_definitive_conserved_LR_genes %in% pd_specific_definitive_genes$gene]
pd_specific_definitive_genes[new_definitive_conserved_LR_genes[!new_definitive_conserved_LR_genes %in% pd_specific_definitive_genes$gene],c('avg_log2FC')]=NA
pd_specific_definitive_genes=pd_specific_definitive_genes[order(pd_specific_definitive_genes$avg_log2FC,decreasing = T),]
pd_specific_definitive_genes$gene=factor(pd_specific_definitive_genes$gene,levels = pd_specific_definitive_genes$gene)
p11=ggplot(pd_specific_definitive_genes,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()
p11


filt_NBRC_altas_seu$pd_celltype1=paste(filt_NBRC_altas_seu$pd_celltype,filt_NBRC_altas_seu$final_celltype,sep='_')
filt_NBRC_altas_seu$pd_celltype1=factor(filt_NBRC_altas_seu$pd_celltype1,levels = c("primitive_ProE","primitive_Bas","primitive_Poly" ,"primitive_Orth", "definitive_BFUE/CFUE","definitive_ProE" , "definitive_Bas","definitive_Poly" ,"definitive_Orth" ))
p22=DotPlot(filt_NBRC_altas_seu,features = pd_specific_definitive_genes$gene,group.by = 'pd_celltype1',scale = F)+scale_color_gradient(low = 'gray',high = 'firebrick3')+RotatedAxis()
p=p11+p22+plot_layout(ncol = 1,heights = c(0.6,1.2));p
#IL7、 NOTCH2、PPARA，F2R
ggsave(p,filename='res_pic/main_figure4/new_definitive_conserved_LR_genes_pd_expression.pdf',width = 8,height = 6)


# 获得new_definitive_conserved_LR 信息
new_conserved_fa_NRBC_mDEGs_inLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$receptor.symbol %in% new_definitive_conserved_LR_genes])
new_conserved_fa_NRBC_mDEGs_inLRs=new_conserved_fa_NRBC_mDEGs_inLRs[new_conserved_fa_NRBC_mDEGs_inLRs %in% NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']]
new_conserved_fa_all_NRBC_mDEGs_LRs_df=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% new_conserved_fa_NRBC_mDEGs_inLRs & NRBC_altas_LR_df$target_type!='Ery2Other',]

new_conserved_fa_NRBC_mDEGs_outLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$ligand.symbol %in% new_definitive_conserved_LR_genes])
new_conserved_fa_NRBC_mDEGs_outLRs=new_conserved_fa_NRBC_mDEGs_outLRs[new_conserved_fa_NRBC_mDEGs_outLRs %in% NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']]
new_conserved_fa_all_NRBC_mDEGs_LRs_df=rbind(new_conserved_fa_all_NRBC_mDEGs_LRs_df,NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% new_conserved_fa_NRBC_mDEGs_outLRs & NRBC_altas_LR_df$target_type=='Ery2Other',])

new_conserved_fa_all_NRBC_mDEGs_LRs_df=new_conserved_fa_all_NRBC_mDEGs_LRs_df[order(new_conserved_fa_all_NRBC_mDEGs_LRs_df$annotation,new_conserved_fa_all_NRBC_mDEGs_LRs_df$annotation),]
new_conserved_fa_all_NRBC_mDEGs_LRs_df$interaction_name=factor(new_conserved_fa_all_NRBC_mDEGs_LRs_df$interaction_name,levels =unique( new_conserved_fa_all_NRBC_mDEGs_LRs_df$interaction_name))
new_conserved_fa_all_NRBC_mDEGs_LRs_df$target_type=factor(new_conserved_fa_all_NRBC_mDEGs_LRs_df$target_type,levels = c('Other2Ery','Ery2Ery','Ery2Other'))
saveRDS(new_conserved_fa_all_NRBC_mDEGs_LRs_df,file = 'new_conserved_fa_all_NRBC_mDEGs_LRs_df.rds')

ggplot(new_conserved_fa_all_NRBC_mDEGs_LRs_df,aes(x=celltype ,y=interaction_name,size=prob ,color=annotation,shape=target_type))+geom_point(alpha=0.7)+facet_grid(~stage )+
  theme_classic()+scale_color_manual(values = cols[-2])+scale_shape_manual(values = c(0:2))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+facet_grid(~stage )+ggtitle('fetal vs adult LR based on hDEGs' )

defintive_specific_conserved_LR=unique(new_conserved_fa_all_NRBC_mDEGs_LRs_df$interaction_name)[unique(new_conserved_fa_all_NRBC_mDEGs_LRs_df$interaction_name) %in% c(incoming_LR_type_df_list[['02_definitive']],outgoing_LR_type_df_list[['02_definitive']])]
defintive_specific_conserved_LR=as.character(defintive_specific_conserved_LR);definitive_shared_LRs

saveRDS(defintive_specific_conserved_LR,file = 'pd_defintive_specific_conserved_LR.rds')


temp_df=new_conserved_fa_all_NRBC_mDEGs_LRs_df[new_conserved_fa_all_NRBC_mDEGs_LRs_df$interaction_name %in%defintive_specific_conserved_LR, ]
temp_df$cluster='primitive'
pd_all_NRBC_mDEGs_LRs_df1=rbind(pd_all_NRBC_mDEGs_LRs_df,temp_df)
pd_all_NRBC_mDEGs_LRs_df1$interaction_name=factor(pd_all_NRBC_mDEGs_LRs_df1$interaction_name,levels = c(as.character(unique(temp_df$interaction_name)),levels(pd_all_NRBC_mDEGs_LRs_df$interaction_name)))
p22=ggplot(pd_all_NRBC_mDEGs_LRs_df1,aes(x=celltype ,y=interaction_name,size=prob ,color=annotation,shape=target_type))+geom_point(alpha=0.7)+facet_grid(~stage )+
  theme_classic()+scale_color_manual(values = cols[-2])+scale_shape_manual(values = c(0:2,5))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+facet_grid(~stage )+ggtitle('primitive vs definitive specific LR ' )
p22
ggsave(p22,filename='res_pic/main_figure4/pd_specific_LR_dotplot.pdf',height = 12,width = 10)



# the expression of all focoused ligand or receptor genes 
focused_conserved_receptor_genes=c('TMEM219','CD55','NR1H2','EPOR','ITGA4','ITGB1','CD47','CD36','CD99')
focused_conserved_ligand_genes=c('PSAP','HLA-A','HLA-B')
pd_hDEGs_LR_genes=levels(pd_candidated_hDEG_lr_genes$gene)
pd_specific_definitive_LR_genes=levels(pd_specific_definitive_genes$gene)

p=VlnPlot(filt_NBRC_altas_seu,features = c(focused_conserved_receptor_genes, focused_conserved_ligand_genes,pd_hDEGs_LR_genes,pd_specific_definitive_LR_genes),
          group.by = 'final_celltype',stack = T,split.by = 'pd_celltype',flip = T)
p
ggsave(p,filename='res_pic/main_figure4/pd_focoused_LR_dotplot.pdf',height = 10,width = 6)


saveRDS(pd_candidated_lr_genes,file = 'res_data/pd_candidated_lr_genes.rds')
saveRDS(pd_all_NRBC_mDEGs_LRs_df1,file = 'res_data/pd_all_NRBC_mDEGs_LRs_df.rds')


# ----------------------------------------查看重点关注通路中的信号流-----------------------------------------------------#

# -------------primitive in YS -------------#
head(pd_all_NRBC_mDEGs_LRs_df[grep('COL',pd_all_NRBC_mDEGs_LRs_df$interaction_name),])

netAnalysis_signalingRole_network(object.list2[['YS']], signaling = 'GDF', width = 24, height = 2.5, font.size = 10)/
netAnalysis_signalingRole_network(object.list2[['YS']], signaling = 'Cholesterol', width = 24, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(object.list2[['YS']], signaling = 'Desmosterol', width = 24, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(object.list2[['YS']], signaling = 'EPO', width = 32, height = 2.5, font.size = 10)# ENDODERM
netAnalysis_signalingRole_network(object.list2[['FL']], signaling = 'EPO', width = 32, height = 2.5, font.size = 10) # HEPATOCYTE_II

netAnalysis_signalingRole_network(object.list2[['FL']], signaling = 'RA', width = 36, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(object.list2[['FBM']], signaling = 'RA', width = 40, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(object.list2[['ABM']], signaling = 'RA', width = 32, height = 2.5, font.size = 10)

# 展示source to target celltype 
netAnalysis_signalingRole_network(object.list2[['FL']], signaling = 'MHC-I', width = 36, height = 2.5, font.size = 10) # ,CD8 T & NK 
netAnalysis_signalingRole_network(object.list2[['FBM']], signaling = 'MHC-I', width = 40, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(object.list2[['ABM']], signaling = 'MHC-I', width = 32, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(object.list2[['FL']], signaling = 'MHC-II', width = 36, height = 2.5, font.size = 10) # target ：DC & Mac & Mono & T & Treg，
netAnalysis_signalingRole_network(object.list2[['FBM']], signaling = 'MHC-II', width = 40, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(object.list2[['ABM']], signaling = 'MHC-II', width = 32, height = 2.5, font.size = 10)

MHC_df=pd_all_NRBC_mDEGs_LRs_df[grep('HLA-',pd_all_NRBC_mDEGs_LRs_df$interaction_name),]
MHC_df$target=as.character(MHC_df$target)
MHC_df$source=factor(MHC_df$source,levels = NRBC_subcelltype)
MHC_df$target[MHC_df$target=='PDC']='pDC'
MHC_df$target[MHC_df$target=='CD8+ T-Cell']='CD8+T'
MHC_df$target[MHC_df$target=='CD4+ T-Cell']='CD4+T'
MHC_df$target[MHC_df$target=='TREG']='Treg'
MHC_df$target[grep('MONOCYTE',MHC_df$target)]='MONOCYTE'
MHC_df$target[grep('MACROPHAGE',MHC_df$target)]='MACROPHAGE'
MHC_df$target[grep('^DC|AS_DC',MHC_df$target)]='DC'
MHC_df=MHC_df[-grep('^CYCLING',MHC_df$target),]
MHC_df$target=factor(MHC_df$target,levels =c("CD4+T","GZMK cytotoxic CD4 T", "Memory CD4 T" ,"CD8+T","Naive CD8 T","TYPE_1_INNATE_T"  , "GZMB CD8 T","GZMK CD8 T", "Treg",
                                             "MONOCYTE","MACROPHAGE","NK","DC","pDC","OSTEOCLAST","THY1+ MSC","Fibro-MSC","Adipo-MSC","VSMC","ENDOTHELIUM_V"  ) )
ggplot(MHC_df,aes(x=source,y=target,color=ligand,size=prob))+geom_point(alpha=0.6)+theme_classic()+facet_grid(~stage)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))

# test the pathway_name
pd_all_NRBC_mDEGs_LRs_df$pathway_name=factor(pd_all_NRBC_mDEGs_LRs_df$pathway_name,levels =unique( pd_all_NRBC_mDEGs_LRs_df$pathway_name))
ggplot(pd_all_NRBC_mDEGs_LRs_df,aes(x=celltype ,y=pathway_name,size=prob ,color=annotation,shape=target_type))+geom_point(alpha=0.7)+facet_grid(~stage )+theme_classic()+scale_color_manual(values = cols[-2])+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+facet_grid(~stage )


# 可以考虑排除那种都高表达的基因
ggplot(NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% pd_all_NRBC_mDEGs_LRs & NRBC_altas_LR_df$target_type=='Other2Ery',],aes(x=celltype ,y=interaction_name,size=prob ,color=pathway_name))+geom_point()+facet_grid(~stage )+theme_classic()+scale_color_manual(values = cols)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+facet_grid(~stage )+ggtitle('incoming LRs')
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =all_NRBC_receptor_genes[all_NRBC_receptor_genes %in% candidated_lr_genes],stack = T,split.by = 'tissue_stage',cols = cols,flip = T)

mDEGs_NRBC_altas_outLR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% all_NRBC_mDEGs_LRs & NRBC_altas_LR_df$target_type=='Ery2Other',]
mDEGs_NRBC_altas_outLR_df=mDEGs_NRBC_altas_outLR_df[mDEGs_NRBC_altas_outLR_df$receptor!='CD8B2',] #  B 细胞分类导致的特异性存在，可以剔除

rm(list=ls());gc()


############################################################################################################################################################
#-----------------------------part3: using the nichnet to analysis the key LR to target genes-------------------------#
############################################################################################################################################################
library(nichenetr)
library(network)
library(ggnetwork)

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")

sig_network <- readRDS('signaling_network_human_21122021.rds')
gr_network <- readRDS('gr_network_human_21122021.rds')
ligand_tf_matrix <- readRDS('ligand_tf_matrix_nsga2r_final.rds')

signaling_pic_func=function(ligands_oi,targets_oi,weighted_networks=weighted_networks,ligand_tf_matrix=ligand_tf_matrix,top_n_regulators=4){
  
  active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,targets_all = targets_oi,
                                                        weighted_networks = weighted_networks,ligand_tf_matrix = ligand_tf_matrix,
                                                        top_n_regulators = top_n_regulators,minmax_scaling = TRUE) 
  
  graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                     ligands_all = ligands_oi, targets_all = targets_oi,
                                                     sig_color = "indianred", gr_color = "steelblue")
  
  
  #构建 network 对象
  # 使用 edges_df 中的 from 和 to 列构建有向图
  edges_df <- graph_min_max$edges_df
  nodes_df <- graph_min_max$nodes_df
  net_obj <- network(as.matrix(edges_df[, c("from", "to")]), directed = TRUE)
  
  # 将节点属性合并到 network 对象中
  # 这一步至关重要，它让 ggnetwork 能读取节点的颜色、类型等信息
  # 确保节点名称匹配
  net_vertices <- network.vertex.names(net_obj)
  
  # 将属性赋值给 network 对象的顶点
  # 使用 match 确保顺序一致
  net_obj %v% "label" <- nodes_df$label[match(net_vertices, nodes_df$id)]
  net_obj %v% "type" <- nodes_df$type[match(net_vertices, nodes_df$id)] # 节点类型（配体/靶点/中间分子）
  net_obj %v% "fillcolor" <- nodes_df$fillcolor[match(net_vertices, nodes_df$id)] # 节点颜色
  
  # 为了映射边的粗细，我们需要将边的权重信息也加入到 network 对象中
  # 注意：network 对象构建时边的顺序可能与 edges_df 不完全一致，需要重新匹配
  # 获取 network 对象中的边列表
  net_edges <- as.data.frame(net_obj, edge.names = TRUE)
  # 匹配 from 和 to 来合并权重
  net_edges$weight <- edges_df$penwidth[match(paste(net_edges$.tail, net_edges$.head, sep = "-"), 
                                              paste(edges_df$from, edges_df$to, sep = "-"))]
  # 将权重赋值给 network 对象的边属性
  net_obj %e% "weight" <- net_edges$weight
  
  # 开始绘图
  p <- ggplot(net_obj, aes(x = x, y = y, xend = xend, yend = yend))+
    # 绘制边 使用 weight 映射线条粗细
    geom_edges(aes(size = weight),alpha = 0.4,  color = "grey70", arrow = arrow(length = unit(1, "mm"), type = "closed")) +
    # 绘制节点 使用 fillcolor 映射填充颜色
    geom_nodes(aes(color = fillcolor), size = 5) +
    # 添加节点标签
    geom_nodetext(aes(label = label),   fontface = "bold", color = "black", size = 3, vjust = 1.5) + # vjust 调整标签位置
    # 主题调整, 移除背景网格
    theme_blank() + labs(title = "NicheNet Signaling Pathway") +
    # 手动设置颜色（保持 NicheNet 原有的配色）
    scale_color_identity() + scale_size_continuous(range = c(0.5, 2)) # 调整线条粗细范围
  
  print(unique(net_obj$label))
  return(p)
}


LR_df=CellChatDB.human$interaction
LR_df$ligand.symbol[LR_df$ligand.symbol==""]=LR_df$ligand[LR_df$ligand.symbol==""]

# --------compare the LR on primitive and definitve NRBC -------------#
primtive_NBRC_ligands=c('DHCR24','LIPA','DHCR7') # 重点分析RORA下游通路导致的改变，但是遗憾的是，以上ligand不在ligand_target_matrix中
primtive_NBRC_ligands=unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='YS' & NRBC_altas_LR_df$target_type!='Ery2Other','ligand'])
primtive_NBRC_ligands= unique(unlist(strsplit(as.character(t(data.frame(strsplit(primtive_NBRC_ligands,split = '_')))[,1]),split = '-')))
primtive_NBRC_ligands=primtive_NBRC_ligands[primtive_NBRC_ligands %in% rownames(filt_NBRC_altas_seu)]  



nichenet_predict_func=function(potential_ligands,target_genes,ligand_target_matrix,background_expressed_genes,legend_title='primitive vs definitve score'){
  predicted_pd_ligand_activites=predict_ligand_activities(geneset =target_genes,ligand_target_matrix = ligand_target_matrix,
                                                          potential_ligands = potential_ligands, background_expressed_genes =background_expressed_genes )
  
  predicted_ligand_activites_mtx=predicted_pd_ligand_activites[order(predicted_pd_ligand_activites$aupr_corrected),]  %>% 
    column_to_rownames('test_ligand') %>% dplyr::select(aupr_corrected) %>% as.matrix(ncol = 1)                                                        
  p=make_heatmap_ggplot(matrix =predicted_ligand_activites_mtx,y_name ='Ligand activaty' ,x_name = "Prioritized ligands",legend_title = "AUPR", color = "darkorange")+
    theme(axis.text.x.top = element_blank())  
  
  # target the DEGs
  ligand_target_gene_link_df <-predicted_pd_ligand_activites$test_ligand%>% lapply(get_weighted_ligand_target_links,geneset =target_genes,ligand_target_matrix = ligand_target_matrix, n = 100) %>%bind_rows() %>% drop_na()
  ligand_target_gene_link_vis=prepare_ligand_target_visualization(ligand_target_df =ligand_target_gene_link_df,ligand_target_matrix = ligand_target_matrix,cutoff = 0.25)
  p2=make_heatmap_ggplot(matrix =t(ligand_target_gene_link_vis[,rownames(predicted_ligand_activites_mtx)[rownames(predicted_ligand_activites_mtx) %in% colnames(ligand_target_gene_link_vis)]]),y_name = 'ligand',x_name = 'target',legend_title =legend_title )
  
  return(list(p1=p,p2=p2,link_vis=ligand_target_gene_link_vis))
  
}
filt_NBRC_altas_seu$pd_celltype2=paste(as.character(filt_NBRC_altas_seu$pd_celltype),as.character(filt_NBRC_altas_seu$final_celltype),sep = "_")
all_mexp_df=as.matrix(AverageExpression(filt_NBRC_altas_seu,group.by = 'pd_celltype2',features =rownames(filt_NBRC_altas_seu) )$RNA)

#primtive_NBRC_ligands=primtive_NBRC_ligands[-grep('CXCL|CCL',primtive_NBRC_ligands)]   

# cho the target geneset-#
# primitive vs definitive
top_pd_degs=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$avg_log2FC >1 & pd_whole_level_markers_degs_df$pct.1 >0.05, ] %>% group_by(cluster) %>% top_n(wt =avg_log2FC,n = 500 ) #& pd_whole_level_markers$pct.2 <0.3
table(top_pd_degs[,c('cluster')])
#definitive  primitive 
#346                500
sub_pd_all_Ery_tissue_markers=read.csv('../Protein_NRBC_marker/DE_marker/primitive_definitive_all_Ery_RNA_markers.csv')
sub_pd_all_Ery_tissue_maconserved_LR_genesrkers=sub_pd_all_Ery_tissue_markers[,-1]
top_pd_sub_degs=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >1 & sub_pd_all_Ery_tissue_markers$pct.1>0.1  ,]  %>% group_by(cluster,celltype) %>% top_n(wt =avg_log2FC,n = 500 )
table(top_pd_sub_degs$cluster)
table(top_pd_sub_degs[,c('cluster','celltype')])

primitive_positive_target_genes=unique(c(top_pd_degs$gene[top_pd_degs$cluster=='primitive'],top_pd_sub_degs$gene[top_pd_sub_degs$cluster=='primitive']))
definitive_positive_target_genes=unique(c(top_pd_degs$gene[top_pd_degs$cluster=='definitive'],top_pd_sub_degs$gene[top_pd_sub_degs$cluster=='definitive']))
definitive_positive_target_genes=definitive_positive_target_genes[-grep('^HIS',definitive_positive_target_genes)]


#known_key_regulator_erythropoiesis_genes=c('GATA1','GATA2','KLF1','TAL1','LMO2','FOG1','BCL11A','EPOR','HIF1A','HIF2A','STAT5A','ALAS2','SLC11A2','FOXO3')
#known_key_regulator_erythropoiesis_genes[known_key_regulator_erythropoiesis_genes %in%primitive_positive_target_genes ] #  "KLF1"  "ALAS2"
#known_key_regulator_erythropoiesis_genes[known_key_regulator_erythropoiesis_genes %in%definitive_positive_target_genes ] #  "GATA2"  "BCL11A" "STAT5A"

potential_ligands=primtive_NBRC_ligands;potential_ligands=potential_ligands[potential_ligands %in% colnames(ligand_target_matrix)]
primtive_target_LR_genes=conserved_LR_genes[conserved_LR_genes %in% pd_whole_level_markers_degs_df$gene[pd_whole_level_markers_degs_df$cluster=='primitive' & pd_whole_level_markers_degs_df$avg_log2FC>1 ]]
primtive_target_LR_genes=c(primtive_target_LR_genes,as.character(pd_candidated_lr_genes$gene[pd_candidated_lr_genes$cluster=='primitive']))

# 排除低表达基因
target_genes=unique(c(primitive_positive_target_genes,primtive_target_LR_genes)) # conserved_LR_genes 
target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('primi',colnames(all_mexp_df))])>0.5]

background_expressed_genes=rownames(filt_NBRC_altas_seu)[!rownames(filt_NBRC_altas_seu)%in%  target_genes]

primitive_ligand_target_res_list=nichenet_predict_func(potential_ligands = potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
p=primitive_ligand_target_res_list[[2]]
ggsave(p,filename='res_pic/main_figure4/pd_primitive_ligand_target_gene_heatmap_score.pdf',width = 6,height = 4)
saveRDS(primitive_ligand_target_res_list,file = 'res_data/pd_primitive_ligand_target_res_list.rds')
primitive_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(filt_primitive_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(primitive_LR_targetgene_enrichGO_res,file = 'res_data/primitive_LR_targetgene_enrichGO_res.rds')
cnetplot(primitive_LR_targetgene_enrichGO_res,20)+ggtitle('the  ligands of primitive NRBC to target genes')


filt_primitive_ligand_target_res_list=nichenet_predict_func(potential_ligands = potential_ligands[-grep('CXCL|CCL',potential_ligands)],target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
p=filt_primitive_ligand_target_res_list[[2]]
ggsave(p,filename='res_pic/main_figure4/pd_primitive_filt_ligand_target_gene_heatmap_score.pdf',width = 6,height = 4)

saveRDS(filt_primitive_ligand_target_res_list,file = 'res_data/pd_primitive_ligand_target_res_list.rds')


temp_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$gene %in% rownames(filt_primitive_ligand_target_res_list[[3]]) & pd_whole_level_markers_degs_df$avg_log2FC>1,]
#temp_df$gene=factor(temp_df$gene,levels = rownames(fetal_ligand_target_res_list[[3]]))
temp_df=temp_df[temp_df$pct.2<0.1 & temp_df$pct.1>0.1,]
temp_df=temp_df[order(temp_df$avg_log2FC,decreasing = T),]
temp_df$gene=factor(temp_df$gene,levels = temp_df$gene)
p1=ggplot(temp_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()+ggtitle('the hDEGs of the target genes by primitive NRBC ligand')

p2=DotPlot(filt_NBRC_altas_seu,features =as.character(temp_df$gene),scale = F,group.by = 'pd_celltype1')+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
p=p1+p2+plot_layout(heights = c(0.4,1.2),ncol = 1);p
ggsave(p,filename='res_pic/main_figure4/pd_primitive_hDEG_targetgene_expression.pdf',width = 4,height = 6)

filt_primitive_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(filt_primitive_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(filt_primitive_LR_targetgene_enrichGO_res,file = 'res_data/primitive_LR_targetgene_enrichGO_res.rds')
cnetplot(filt_primitive_LR_targetgene_enrichGO_res,20)+ggtitle('the  ligands of primitive NRBC to target genes')
dotplot(filt_primitive_LR_targetgene_enrichGO_res)

p=draw_gonetcwork_pic_func(res =filt_primitive_LR_targetgene_enrichGO_res,showCategory =20,xlimits = c(-0.2,2.2) );p
ggsave(p,filename='res_pic/main_figure4/pd_primitive_target_gene_enrichGO_res.pdf',width =6,height = 6)

primitive_LR_targetgene_enrichGO_MFres=enrichGO(gene =rownames(filt_primitive_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'MF' )
dotplot(primitive_LR_targetgene_enrichGO_MFres)


#
#ligands_oi <- c('EPO') # this can be a list of multiple ligands  if required
#targets_oi <- names(filt_primitive_ligand_target_res_list[[3]][,'EPO'][filt_primitive_ligand_target_res_list[[3]][,'EPO'] >0.05])#"IFITM2
#EPO_signaling=signaling_pic_func(ligands_oi =ligands_oi,targets_oi = targets_oi,weighted_networks =weighted_networks ,ligand_tf_matrix =ligand_tf_matrix ,top_n_regulators =4  )
# primitive_EPO_signaling_network.pdf,5 X5


ligands_oi <- c('IGFBP3') # this can be a list of multiple ligands  if required
targets_oi <- c("CDKN1A")#"IFITM2
IGFBP3_signaling=signaling_pic_func(ligands_oi =ligands_oi,targets_oi = targets_oi,weighted_networks =weighted_networks ,ligand_tf_matrix =ligand_tf_matrix ,top_n_regulators =4  )
#  IGFBP3_signaling_network.pdf,4 X4 


key_network_genes=list();targets_gene_list=list()
key_network_genes[['IGFBP3']]= unique(IGFBP3_signaling$data$label)[!unique(IGFBP3_signaling$data$label) %in% c(targets_oi,ligands_oi)]
targets_gene_list[['IGFBP3']]=targets_oi

#---------------------------------------defintive ligand to target genes------------------------------#
definitive_NBRC_ligands=unique(c(incoming_LR_type_df_list[['01_shared']],incoming_LR_type_df_list[['02_definitive']]))
definitive_NBRC_ligands
definitive_NBRC_ligands[grep('ITGA4',definitive_NBRC_ligands)]
definitive_NBRC_ligands=unique(LR_df[LR_df$interaction_name %in% definitive_NBRC_ligands,'ligand.symbol' ])
definitive_NBRC_ligands=definitive_NBRC_ligands[definitive_NBRC_ligands %in% rownames(filt_NBRC_altas_seu)]  


potential_ligands=c(definitive_NBRC_ligands,'EPO');potential_ligands=potential_ligands[potential_ligands %in% colnames(ligand_target_matrix)]
definitive_target_LR_genes=conserved_LR_genes[conserved_LR_genes %in% pd_whole_level_markers_degs_df$gene[pd_whole_level_markers_degs_df$cluster=='definitive' & pd_whole_level_markers_degs_df$avg_log2FC>1 ]]
definitive_target_LR_genes=c(definitive_target_LR_genes,as.character(pd_candidated_lr_genes$gene[pd_candidated_lr_genes$cluster=='definitive']))
target_genes=unique(c(definitive_positive_target_genes,definitive_target_LR_genes))
target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('defi',colnames(all_mexp_df))])>0.5]
background_expressed_genes=rownames(filt_NBRC_altas_seu)[!rownames(filt_NBRC_altas_seu)%in%  target_genes]

definitive_ligand_target_res_list=nichenet_predict_func(potential_ligands =potential_ligands[-grep('CXCL|CCL',potential_ligands)],target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes,legend_title = 'definitve ligand to target genes' )
p=definitive_ligand_target_res_list[[2]];p
ggsave(p,filename='res_pic/main_figure4/pd_definitive_ligand_target_gene_heatmap_score.pdf',width = 6,height = 8)


#
an_df=data.frame(definitive_ligand_target_res_list[[1]]$data[definitive_ligand_target_res_list[[1]]$data$y %in% c('EPO','MDK','APP'),])
an_df=data.frame(row.names = an_df$y,aupr_score=an_df$score)

p=pheatmap(t(definitive_ligand_target_res_list[[3]][,c('EPO','APP','MDK')]),annotation_row = an_df,border_color = 'white',cluster_rows = F,cluster_cols =F ,color =colorRampPalette(colors =  c("white" ,"#EF6548"))(100))
ggsave(as.ggplot(p),file='res_pic/main_figure4/definitive_key_ligand_target_genes_heatmap.pdf',width = 12,height = 3,dpi = 300)

saveRDS(definitive_ligand_target_res_list,file = 'res_data/pd_definitive_ligand_target_res_list.rds')

DotPlot(filt_NBRC_altas_seu,features = rownames(definitive_ligand_target_res_list[[3]]),scale = F)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')

definitive_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(definitive_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(definitive_LR_targetgene_enrichGO_res,file = 'res_data/definitive_LR_targetgene_enrichGO_res.rds')
cnetplot(definitive_LR_targetgene_enrichGO_res,20)+ggtitle('the  ligands of definitive NRBC to target genes')
dotplot(definitive_LR_targetgene_enrichGO_res)
p=draw_gonetcwork_pic_func(res =definitive_LR_targetgene_enrichGO_res,showCategory =20,xlimits = c(-0.2,2) );p
ggsave(p,filename='res_pic/main_figure4/pd_definitive_target_gene_enrichGO_res.pdf',width =6,height = 6)


temp_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$gene %in% rownames(definitive_ligand_target_res_list[[3]]) & pd_whole_level_markers_degs_df$avg_log2FC>1,]
#temp_df$gene=factor(temp_df$gene,levels = rownames(fetal_ligand_target_res_list[[3]]))
temp_df=temp_df[temp_df$pct.2<0.1 & temp_df$pct.1>0.08,]
temp_df=temp_df[order(temp_df$avg_log2FC,decreasing = T),]
temp_df$gene=factor(temp_df$gene,levels = temp_df$gene)
p1=ggplot(temp_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()+ggtitle('the hDEGs of the target genes by definitive NRBC ligand')

p2=DotPlot(filt_NBRC_altas_seu,features =as.character(temp_df$gene),scale = F)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
p=p1+p2+plot_layout(heights = c(0.4,1.2),ncol = 1);p
ggsave(p,filename='res_pic/main_figure4/pd_definitive_hDEG_targetgene_expression.pdf',width = 8,height = 12)



ligands_oi <- "MDK" # this can be a list of multiple ligands if required
targets_oi <- c( "CCNA2",'CDK6','CDK4')#"IFITM2
MDK_signaling=signaling_pic_func(ligands_oi =ligands_oi,targets_oi = targets_oi,weighted_networks =weighted_networks ,ligand_tf_matrix =ligand_tf_matrix ,top_n_regulators =4  )
#MDK_signaling_network.pdf,5x5

key_network_genes=list('MDK'= unique(MDK_signaling$data$label)[!unique(MDK_signaling$data$label) %in% c(targets_oi,ligands_oi)])
targets_gene_list=list('MDK'=targets_oi)


ligands_oi <- c('APP') # this can be a list of multiple ligands if required
targets_oi <- names(definitive_ligand_target_res_list[[3]][,'APP'][definitive_ligand_target_res_list[[3]][,'APP'] >0.04])#"IFITM2
APP_signaling=signaling_pic_func(ligands_oi =ligands_oi,targets_oi = targets_oi,weighted_networks =weighted_networks ,ligand_tf_matrix =ligand_tf_matrix ,top_n_regulators =4  )
#APP_signaling_network.pdf,5x5
key_network_genes[['APP']]= unique(APP_signaling$data$label)[!unique(APP_signaling$data$label) %in% c(targets_oi,ligands_oi)]
targets_gene_list[['APP']]=targets_oi


ligands_oi <- c('EPO') # this can be a list of multiple ligands  if required
targets_oi <- c( names(filt_primitive_ligand_target_res_list[[3]][,'EPO'][filt_primitive_ligand_target_res_list[[3]][,'EPO'] >0.05]),
                 names(definitive_ligand_target_res_list[[3]][,'EPO'][definitive_ligand_target_res_list[[3]][,'EPO'] >0.05]))
EPO_signaling=signaling_pic_func(ligands_oi =ligands_oi,targets_oi = targets_oi,weighted_networks =weighted_networks ,ligand_tf_matrix =ligand_tf_matrix ,top_n_regulators =4  )
# EPO_signaling_network.pdf,5 X6
key_network_genes[['EPO']]= unique(EPO_signaling$data$label)[!unique(EPO_signaling$data$label) %in% c(targets_oi,ligands_oi)]
targets_gene_list[['EPO']]=targets_oi

saveRDS(key_network_genes,file = 'res_data/key_network_genes.rds')
saveRDS(targets_gene_list,file = 'res_data/targets_gene_list.rds')


p=DotPlot(filt_NBRC_altas_seu,group.by = 'pd_celltype1',cols = c('gray','firebrick3'),scale = F,
        features =unique(c(key_network_genes[['EPO']],key_network_genes[['IGFBP3']],
                           key_network_genes[['MDK']],key_network_genes[['APP']])
                    ) )+RotatedAxis()
p
ggsave(p,filename='res_pic/main_figure4/key_signaling_gene_expression_dotplot.pdf',width =12 ,height = 5)

p=DotPlot(filt_NBRC_altas_seu,group.by = 'pd_celltype1',cols = c('gray','firebrick3'),scale = F,
          features = c('KRT19','HMOX1', 'NDRG1','SOX4','MYB','LMNA','CMBL','ISG20','ITGA2B'))+RotatedAxis()
p
ggsave(p,filename='res_pic/main_figure4/pd_specific_target_gene_expression_dotplot.pdf',width =8 ,height = 5)


################################################################################################################################
#--------------------------part4: expression of key ligand in altas---------------------------------------#
################################################################################################################################
cho_feature=c('IGFBP3',"EPO",'APP','MDK')
FL_altas_seu=readRDS('../NRBC_FL_altas/tmp_FL_altas_seu.rds')


p1=DotPlot(FL_altas_seu,features = cho_feature,group.by = 'subcelltype',cols = c('gray','firebrick3'),scale = F)+RotatedAxis()
p1
FL_cho_cells=as.character(unique(p1$data[p1$data$avg.exp >0.2 & p1$data$pct.exp>10 ,'id']))
FL_cho_cells
FL_cho_cells=FL_cho_cells[!FL_cho_cells %in% c( "MACROPHAGE_ERY", "MYELOCYTE", "NK/T CELLS" , "MONOCYTE" )]
p1=VlnPlot(subset(FL_altas_seu,subcelltype %in% FL_cho_cells),group.by = 'subcelltype',cols = cols,features =cho_feature,stack = T)+NoLegend()+ggtitle('FL ALTAS')
p1
ggsave(p1,filename='res_pic/main_figure4/FL_celltype_key_ligand_expression_vlnplot.pdf',height =6,width = 4)

FL_altas_seu$age=factor(FL_altas_seu$age,levels = c("CS14_4PCW" ,"CS15_5PCW", "CS17_6PCW", "CS18","CS22","CS23",  "8PCW","8.1PCW","9.1PCW","9.7PCW",  "11PCW","11.4PCW" ,"12PCW","13.9PCW","14.4PCW","15PCW","16.3PCW","16PCW","17PCW"  ))
p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c('SMOOTH MUSCLE')),
          group.by = 'age',features =c('IGFBP3'),cols = col,pt.size = 0)+NoLegend()+ggtitle('FL  SMOOTH MUSCLE')
p11=VlnPlot(subset(FL_altas_seu,subcelltype %in% c('HEPATOCYTE')),
          group.by = 'age',features =c('EPO'),cols = col,pt.size = 0,split.by = 'anno_lvl_2_final_clean')+NoLegend()+ggtitle('FL  HEPATOCYTE')
p=p+p11

ggsave(p,filename='res_pic/main_figure4/FL_IGFBP3_EPO_key_ligand_expression_vlnplot.pdf',height = 6,width = 5)
rm(p);gc()


FL_altas_seu$id=FL_altas_seu$donor
FL_altas_seu$id[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]]=paste(FL_altas_seu$donor[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]],
                FL_altas_seu$age[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]],sep='_')

FL_cho_gene_aggregated_exp=AggregateExpression(FL_altas_seu,features = cho_feature,group.by = 'id')$RNA
FL_cho_gene_aggregated_exp=log2(FL_cho_gene_aggregated_exp+1)
colnames(FL_cho_gene_aggregated_exp)=gsub(pattern = 'wk','WPC',colnames(FL_cho_gene_aggregated_exp))

order_sampleid=c( "FL-4WPC" ,"FL-5WPC" ,"FL-6WPC", "F61-CS18","F35-CS22","F32-CS22","F34-CS23", "FL-8WPC"  ,"F16-8.1PCW","F17-9.1PCW","F22-9.7PCW","F33-9.7PCW","FL-11WPC" ,
                  "F23-11.4PCW","F30-14.4PCW","F38-12PCW" , "F45-13.9PCW" ,"F30-14.4PCW", "FL-15WPC"  ,"F41-16PCW","F21-16.3PCW","F29-17PCW")

colnames(FL_cho_gene_aggregated_exp)[!colnames(FL_cho_gene_aggregated_exp) %in% order_sampleid]
FL_cho_gene_aggregated_exp=FL_cho_gene_aggregated_exp[,order_sampleid]
p=pheatmap(FL_cho_gene_aggregated_exp,cluster_rows = F,cluster_cols = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure4/key_ligand_expression_in_FL_sample_heatmap.pdf',width =6,height = 3)

FL_cho_gene_aggregated_exp_df=t(FL_cho_gene_aggregated_exp)
FL_cho_gene_aggregated_exp_df=data.frame(FL_cho_gene_aggregated_exp_df/FL_cho_gene_aggregated_exp_df[,'MDK'])
FL_cho_gene_aggregated_exp_df=FL_cho_gene_aggregated_exp_df[,cho_feature]
FL_cho_gene_aggregated_exp_df[FL_cho_gene_aggregated_exp_df==0]=NA
boxplot(FL_cho_gene_aggregated_exp_df,rm.na=T) # FL_key_ligand_ratio_refMDK_boxplot.pdf, 5 x5 
round(colMedians(as.matrix(FL_cho_gene_aggregated_exp_df),na.rm = T),digits = 1)
#IGFBP3    EPO    APP    MDK 
#0.9    0.5    0.9    1.0 

rm(FL_altas_seu);gc()

BM_altas_seu=readRDS('../NRBC_BM_altas/BM_altas_seu_v2.rds')#
VlnPlot(BM_altas_seu,group.by = 'age',features = cho_feature,stack = T)

# MDK 几乎在所有细胞中都有表达
p0=DotPlot(subset(BM_altas_seu,stage=='EBM'),features = cho_feature,group.by = 'new_celltype',cols = c('gray','firebrick3'),scale = F)+RotatedAxis()
p0
p0=VlnPlot(subset(BM_altas_seu,stage=='EBM'),features = cho_feature[cho_feature!='EPO'],group.by = 'new_celltype',stack = T,cols = cols)+NoLegend()+ggtitle('EBM ALTAS')
ggsave(p0,filename='res_pic/main_figure4/EBM_celltype_key_ligand_expression_vlnplot.pdf',height =6,width = 3)


levels(BM_altas_seu$new_celltype)

p2=DotPlot(subset(BM_altas_seu,stage=='FBM'),features = cho_feature,group.by = 'new_celltype',cols = c('gray','firebrick3'),scale = F)+RotatedAxis()
p2
FBM_cho_cells=as.character(unique(p2$data[p2$data$avg.exp >0.2 & p2$data$pct.exp>10 ,'id']))
FBM_cho_cells
FBM_cho_cells=FBM_cho_cells[!FBM_cho_cells %in%c("EO/BASO/MAST", "NEUTROPHIL","MONOCYTE","Mac_Ery" ,"NK/T CELLS","MACROPHAGE" )]
p2=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% FBM_cho_cells ),cols=cols,features = cho_feature[cho_feature!='EPO'],stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('FBM ALTAS')
p2
ggsave(p2,filename='res_pic/main_figure4/FBM_celltype_key_ligand_expression_vlnplot.pdf',height =6,width = 3)

p3=DotPlot(subset(BM_altas_seu,stage=='ABM'),features = cho_feature,group.by = 'new_celltype',cols = c('gray','firebrick3'),scale = F)+RotatedAxis()
p3
ABM_cho_cells=as.character(unique(p3$data[p3$data$avg.exp >0.2 & p3$data$pct.exp>10 ,'id']))
ABM_cho_cells=ABM_cho_cells[!ABM_cho_cells %in% c("EO/BASO/MAST","MACROPHAGE", "Treg", "CD4 T" )]
p3=VlnPlot(subset(subset(BM_altas_seu,stage=='ABM'),new_celltype %in% ABM_cho_cells),cols = cols,features = cho_feature[cho_feature!='EPO'],stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('ABM ALTAS')
p3
ggsave(p3,filename='res_pic/main_figure4/ABM_key_ligand_expression_celltype_vlnplot.pdf',width = 3,height = 6)



BM_altas_seu$donor[is.na(BM_altas_seu$donor)]=BM_altas_seu$sample[is.na(BM_altas_seu$donor)]
BM_altas_seu$id=paste(BM_altas_seu$donor,BM_altas_seu$age,sep="_")
table(BM_altas_seu$id)
BM_cho_gene_aggregated_exp=AggregateExpression(BM_altas_seu,group.by = 'id',features =cho_feature[cho_feature!='EPO'] )$RNA
BM_cho_gene_aggregated_exp=log2(BM_cho_gene_aggregated_exp+1)
BM_cho_gene_aggregated_exp=BM_cho_gene_aggregated_exp[,]
p=pheatmap(BM_cho_gene_aggregated_exp[,grep('H|F|CS',colnames(BM_cho_gene_aggregated_exp))],cluster_cols = F,cluster_rows = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure4/key_ligand_expression_in_BM_sample_heatmap.pdf',width =6,height = 3)

EBM_cho_gene_aggregated_exp_df=t(BM_cho_gene_aggregated_exp[,grep('-CS',colnames(BM_cho_gene_aggregated_exp))])
EBM_cho_gene_aggregated_exp_df=data.frame(EBM_cho_gene_aggregated_exp_df/EBM_cho_gene_aggregated_exp_df[,'MDK'])
EBM_cho_gene_aggregated_exp_df=EBM_cho_gene_aggregated_exp_df[,cho_feature[cho_feature!='EPO']]
EBM_cho_gene_aggregated_exp_df[EBM_cho_gene_aggregated_exp_df==0]=NA
boxplot(EBM_cho_gene_aggregated_exp_df,rm.na=T) ## EBM_key_ligand_ratio_refMDK_boxplot.pdf, 4 x4 
 #IGFBP3    APP    MDK 
#0.6    1.0    1.0 

FBM_cho_gene_aggregated_exp_df=t(BM_cho_gene_aggregated_exp[,grep('PCW',colnames(BM_cho_gene_aggregated_exp))])
FBM_cho_gene_aggregated_exp_df=data.frame(FBM_cho_gene_aggregated_exp_df/FBM_cho_gene_aggregated_exp_df[,'MDK'])
FBM_cho_gene_aggregated_exp_df=FBM_cho_gene_aggregated_exp_df[,cho_feature[cho_feature!='EPO']]
FBM_cho_gene_aggregated_exp_df[FBM_cho_gene_aggregated_exp_df==0]=NA
boxplot(FBM_cho_gene_aggregated_exp_df[grep('CW',rownames(FBM_cho_gene_aggregated_exp_df)),],rm.na=T) ## FBM_key_ligand_ratio_refMDK_boxplot.pdf, 4 x4 
round(colMedians(as.matrix(FBM_cho_gene_aggregated_exp_df[grep('CW',rownames(FBM_cho_gene_aggregated_exp_df)),]),na.rm = T),digits = 1)
#IGFBP3    APP    MDK 
#0.6    1.0    1.0 

ABM_cho_gene_aggregated_exp_df=t(BM_cho_gene_aggregated_exp[,grep('y',colnames(BM_cho_gene_aggregated_exp))])
ABM_cho_gene_aggregated_exp_df=data.frame(ABM_cho_gene_aggregated_exp_df/ABM_cho_gene_aggregated_exp_df[,'MDK'])
ABM_cho_gene_aggregated_exp_df=ABM_cho_gene_aggregated_exp_df[,cho_feature[cho_feature!='EPO']]
ABM_cho_gene_aggregated_exp_df[ABM_cho_gene_aggregated_exp_df==0]=NA
boxplot(ABM_cho_gene_aggregated_exp_df,rm.na=T) ## ABM_key_ligand_ratio_refMDK_boxplot.pdf, 4 x4 
round(colMedians(as.matrix(ABM_cho_gene_aggregated_exp_df),na.rm = T),digits = 1)
#IGFBP3    APP    MDK 
#0.9    1.1    1.0 

YS_altas_seu= readRDS('../NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' )
YS_altas_seu=NormalizeData(YS_altas_seu)
levels(YS_altas_seu$subcelltype)

p4=DotPlot(YS_altas_seu,features = cho_feature,group.by = 'subcelltype',cols = c('gray','firebrick3'),scale = F)+RotatedAxis()
p4
YS_cho_celltype=as.character(unique(p4$data[p4$data$avg.exp >0.2 & p4$data$pct.exp>10 ,'id']))
YS_cho_celltype
YS_cho_celltype=YS_cho_celltype[!YS_cho_celltype %in% c("ERYTHROID","MONOCYTE","MACROPHAGE")]
   
ggsave(p4,filename='res_pic/main_figure4/key_ligand_expression_celltype_YS_vlnplot1.pdf',width = 4,height = 6)

YS_cho_celltype=YS_cho_celltype[YS_cho_celltype %in% c("ENDODERM","FIBROBLAST","ENDOTHELIUM", "MESOTHELIUM", "SMOOTH_MUSCLE" )]
p4=VlnPlot(subset(YS_altas_seu,subcelltype %in% YS_cho_celltype),group.by = 'subcelltype',features =c('IGFBP3','EPO'),stack = T)+NoLegend()+ggtitle('YS ALTAS')
p4
ggsave(p4,filename='res_pic/main_figure4/key_ligand_expression_celltype_YS_vlnplot.pdf',width = 3,height = 6)

p41=VlnPlot(subset(YS_altas_seu,subcelltype %in% c("FIBROBLAST","ENDOTHELIUM", "MESOTHELIUM", "SMOOTH_MUSCLE" )),pt.size = 0,cols = cols,group.by = 'stage',features =c('IGFBP3'))+NoLegend()
p42=VlnPlot(subset(YS_altas_seu,subcelltype %in% c("ENDODERM" )),cols = cols,pt.size = 0.1,group.by = 'stage',features =c('EPO'))+NoLegend()
p=p41/p42;p
ggsave(p,filename='res_pic/main_figure4/key_ligand_expression_stage_YS_vlnplot.pdf',width = 4,height = 5)

table(YS_altas_seu@meta.data[YS_altas_seu$subcelltype %in% c("FIBROBLAST","ENDOTHELIUM", "MESOTHELIUM", "SMOOTH_MUSCLE" ),c('stage')])
table(YS_altas_seu@meta.data[YS_altas_seu$subcelltype=='ENDODERM',c('stage')])



YS_cho_gene_aggregated_exp=AggregateExpression(YS_altas_seu,features = cho_feature,group.by = 'id')$RNA
YS_cho_gene_aggregated_exp=log2(YS_cho_gene_aggregated_exp+1)
p=pheatmap(YS_cho_gene_aggregated_exp,cluster_rows = F,cluster_cols = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure4/key_ligand_expression_in_YS_sample_heatmap.pdf',height = 3,width = 6)

YS_cho_gene_aggregated_exp_df=t(YS_cho_gene_aggregated_exp)
YS_cho_gene_aggregated_exp_df=data.frame(YS_cho_gene_aggregated_exp_df/YS_cho_gene_aggregated_exp_df[,'MDK'])
YS_cho_gene_aggregated_exp_df=YS_cho_gene_aggregated_exp_df[,cho_feature]
YS_cho_gene_aggregated_exp_df[YS_cho_gene_aggregated_exp_df==0]=NA
boxplot(YS_cho_gene_aggregated_exp_df,rm.na=T) # YS_key_ligand_ratio_refMDK_boxplot.pdf, 5 x5 
round(colMedians(as.matrix(YS_cho_gene_aggregated_exp_df),na.rm = T),digits = 1)

#APP    EPO IGFBP3    MDK 
#1.0    0.5    0.9    1.0 



