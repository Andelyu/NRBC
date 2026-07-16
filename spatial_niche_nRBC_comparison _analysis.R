#---------------------------------------------------tissue RNA marker ----------------------------------------------------#

################################################################################################################################################################
#---------------------------------------------------prepare pkg and function ----------------------------------------------------#
################################################################################################################################################################
library(ggrepel)
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(SingleCellExperiment)
library(RColorBrewer )
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tibble)
library(pheatmap)
library(clusterProfiler)
library(ggplotify)
library(ggplot2)
library(reshape2)
library(corrplot)
library(APL)
library(patchwork)
library(yulab.utils)

#install.packages("patchwork")
library(RColorBrewer )
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
getwd()
setwd('/home/gibh/2021_NRBC_chlyu')

dir.create('Protein_NRBC_marker')

find_mDEGs_func=function(seu,group,sfile,mgenes=new_all_merged_plasma_protein){
  early_Ery_tissue_markers=FindAllMarkers(subset(seu,Ery_stage=='early_Ery'))
  mid_Ery_tissue_markers=  FindAllMarkers(subset(seu,Ery_stage=='mid_Ery'))
  late_Ery_tissue_markers= FindAllMarkers(subset(seu,Ery_stage=='late_Ery'))
  
  early_Ery_tissue_markers$celltype='early_Ery'
  mid_Ery_tissue_markers$celltype='mid_Ery'
  late_Ery_tissue_markers$celltype='late_Ery'
  all_Ery_tissue_markers=rbind(early_Ery_tissue_markers,rbind(mid_Ery_tissue_markers,late_Ery_tissue_markers))
  rownames(all_Ery_tissue_markers)=NULL
  all_Ery_tissue_markers$group=group
  all_Ery_tissue_markers$membrane='no'
  all_Ery_tissue_markers$membrane[all_Ery_tissue_markers$gene %in% mgenes]='yes'
  
  write.csv(all_Ery_tissue_markers,file = sfile,quote = F)
  rm(early_Ery_tissue_markers,mid_Ery_tissue_markers,late_Ery_tissue_markers);gc()
  
  pos_all_Ery_tissue_markers=all_Ery_tissue_markers[all_Ery_tissue_markers$avg_log2FC >0  & all_Ery_tissue_markers$p_val_adj <0.01,]
  if(dim(pos_all_Ery_tissue_markers)[2] <1 ){ stop('no DEGs are found!!!')}
  
  count_tmp_df=cbind( data.frame(table(pos_all_Ery_tissue_markers[,c('cluster','celltype')])[1,]), data.frame(table(pos_all_Ery_tissue_markers[,c('cluster','celltype')])[2,]))
  colnames(count_tmp_df)= strsplit(group,'_')[[1]]
  type1=strsplit(group,split = '_')[[1]][1]
  type2=strsplit(group,split = '_')[[1]][2]
  
  count_tmp_df['all',type1]=length(unique(pos_all_Ery_tissue_markers$gene[pos_all_Ery_tissue_markers$cluster==type1]))
  count_tmp_df['all',type2]=length(unique(pos_all_Ery_tissue_markers$gene[pos_all_Ery_tissue_markers$cluster==type2]))
  count_tmp_df['add_all',]=length(unique(pos_all_Ery_tissue_markers$gene))
  
  mall_Ery_tissue_markers=pos_all_Ery_tissue_markers[pos_all_Ery_tissue_markers$membrane=='yes',]
  count_tmp_df2=cbind( data.frame(table(mall_Ery_tissue_markers[,c('cluster','celltype')])[1,]), data.frame(table(mall_Ery_tissue_markers[,c('cluster','celltype')])[2,]))
  colnames(count_tmp_df2)= paste('m',strsplit(group,'_')[[1]],sep = "")
  count_tmp_df2['all',paste('m',type1,sep='')]=length(unique(mall_Ery_tissue_markers$gene[mall_Ery_tissue_markers$cluster==type1]))
  count_tmp_df2['all',paste('m',type2,sep='')]=length(unique(mall_Ery_tissue_markers$gene[mall_Ery_tissue_markers$cluster==type2]))
  count_tmp_df2['add_all',]=length(unique(mall_Ery_tissue_markers$gene))
  
  count_tmp_df=cbind(count_tmp_df,count_tmp_df2)
  
  return(list(all_Ery_tissue_markers,count_tmp_df))
}
get_cmmarkers_func=function(Ery_tissue_markers,pct_diff=0.2,pct.1=0.1,pct.2=0.3,avg_log2FC=1,exp_value=1,celltype_mexp_df=celltype_mexp_df,max_score=20,CC_gene_df=CC_gene_df,CD_maker_genes = all_CD_maker_genes){
  
  # get hDEGs
  Ery_tissue_markers=Ery_tissue_markers[Ery_tissue_markers$avg_log2FC >0  & Ery_tissue_markers$p_val_adj <0.01,] 
  Ery_tissue_markers$pct_diff=Ery_tissue_markers$pct.1-Ery_tissue_markers$pct.2
  Ery_tissue_markers=Ery_tissue_markers[Ery_tissue_markers$pct_diff >pct_diff & Ery_tissue_markers$pct.1 > pct.1 & Ery_tissue_markers$pct.2< pct.2 & Ery_tissue_markers$avg_log2FC >avg_log2FC,]
  print(dim(Ery_tissue_markers))
  
  if(dim(Ery_tissue_markers)[2] <1){ stop('no hDEGs are found!!')}
  
  group=unique(Ery_tissue_markers$group)
  type1=strsplit(group,split = '_')[[1]][1]
  type2=strsplit(group,split = '_')[[1]][2]
  count_tmp_df=cbind( data.frame(table(Ery_tissue_markers[,c('cluster','celltype')])[1,]), data.frame(table(Ery_tissue_markers[,c('cluster','celltype')])[2,]))
  colnames(count_tmp_df)= c(type1,type2)
  count_tmp_df['all',type1]=length(unique(Ery_tissue_markers$gene[Ery_tissue_markers$cluster==type1]))
  count_tmp_df['all',type2]=length(unique(Ery_tissue_markers$gene[Ery_tissue_markers$cluster==type2]))
  count_tmp_df['add_all',]=length(unique(Ery_tissue_markers$gene))
  
  mEry_tissue_markers=Ery_tissue_markers[Ery_tissue_markers$membrane=='yes',]
  #mEry_tissue_markers=Ery_tissue_markers
  count_tmp_df2=cbind( data.frame(table(mEry_tissue_markers[,c('cluster','celltype')])[1,]), data.frame(table(mEry_tissue_markers[,c('cluster','celltype')])[2,]))
  colnames(count_tmp_df2)= paste('m',strsplit(group,'_')[[1]],sep = "")
  count_tmp_df2['all',paste('m',type1,sep='')]=length(unique(mEry_tissue_markers$gene[mEry_tissue_markers$cluster==type1]))
  count_tmp_df2['all',paste('m',type2,sep='')]=length(unique(mEry_tissue_markers$gene[mEry_tissue_markers$cluster==type2]))
  count_tmp_df2['add_all',]=length(unique(mEry_tissue_markers$gene))
  
  count_tmp_df=cbind(count_tmp_df,count_tmp_df2)
  
  # get cmmarkers
  
  # filtout the control group celltype_mexp > 1
  
  filt_type1_mmakers=unique(mEry_tissue_markers$gene[mEry_tissue_markers$cluster==type1])[rowSums(celltype_mexp_df[unique(mEry_tissue_markers$gene[mEry_tissue_markers$cluster==type1]),grep(type2,colnames(celltype_mexp_df))] >exp_value) >0]
  filt_type2_mmakers=unique(mEry_tissue_markers$gene[mEry_tissue_markers$cluster==type2])[rowSums(celltype_mexp_df[unique(mEry_tissue_markers$gene[mEry_tissue_markers$cluster==type2]),grep(type1,colnames(celltype_mexp_df))] >exp_value) >0]
  if(length(c(filt_type1_mmakers,filt_type2_mmakers)) >0){
    mEry_tissue_markers=mEry_tissue_markers[!mEry_tissue_markers$gene %in% c(filt_type1_mmakers,filt_type2_mmakers),]
    
  }
  print(dim(mEry_tissue_markers))
  
  # caculate the marker_score
  temp1=mEry_tissue_markers[mEry_tissue_markers$cluster==type1,]
  if(dim(temp1)[1] >1 ){
    temp1$pct1_maxmexp=rowMax(as.matrix(celltype_mexp_df[temp1$gene,grep(type1,colnames(celltype_mexp_df))]))
    temp1$pct2_maxmexp=rowMax(as.matrix(celltype_mexp_df[temp1$gene,grep(type2,colnames(celltype_mexp_df))]))
    
  }else if(dim(temp1)[1] ==1){
    temp1$pct1_maxmexp=max(as.matrix(celltype_mexp_df[temp1$gene,grep(type1,colnames(celltype_mexp_df))]))
    temp1$pct2_maxmexp=max(as.matrix(celltype_mexp_df[temp1$gene,grep(type2,colnames(celltype_mexp_df))]))
    
  }
  
  temp2=mEry_tissue_markers[mEry_tissue_markers$cluster==type2,]
  if(dim(temp2)[1] >1 ){
    temp2$pct1_maxmexp=rowMax(as.matrix(celltype_mexp_df[temp2$gene,grep(type2,colnames(celltype_mexp_df))]))
    temp2$pct2_maxmexp=rowMax(as.matrix(celltype_mexp_df[temp2$gene,grep(type1,colnames(celltype_mexp_df))]))
  }else if(dim(temp1)[1] ==1){
    temp2$pct1_maxmexp=max(as.matrix(celltype_mexp_df[temp2$gene,grep(type2,colnames(celltype_mexp_df))]))
    temp2$pct2_maxmexp=max(as.matrix(celltype_mexp_df[temp2$gene,grep(type1,colnames(celltype_mexp_df))]))
    
  }
  
  mEry_tissue_markers=rbind(temp1,temp2)
  if(table(is.na( mEry_tissue_markers$pct2_maxmexp==0)) >0){mEry_tissue_markers$pct2_maxmexp[ mEry_tissue_markers$pct2_maxmexp==0]=0.001}
  mEry_tissue_markers$marker_score1=log2(mEry_tissue_markers$pct_diff*mEry_tissue_markers$avg_log2FC*mEry_tissue_markers$pct1_maxmexp/(mEry_tissue_markers$pct2_maxmexp^2))
  mEry_tissue_markers$marker_score1[mEry_tissue_markers$marker_score1 >max_score]=max_score
  
  mEry_tissue_markers=mEry_tissue_markers[order(mEry_tissue_markers$cluster,mEry_tissue_markers$marker_score1,decreasing = T),]
  mEry_tissue_markers=mEry_tissue_markers[!duplicated(mEry_tissue_markers$gene),]
  mEry_tissue_markers$marker_type= CC_gene_df[match(x = mEry_tissue_markers$gene,table = CC_gene_df$CC_genes),'type']
  mEry_tissue_markers$marker_type[is.na( mEry_tissue_markers$marker_type)]='no'
  mEry_tissue_markers$marker_type[mEry_tissue_markers$gene %in% CD_maker_genes]='CD_maker'
  mEry_tissue_markers$marker_type[mEry_tissue_markers$gene %in% intersect(CC_gene_df$CC_genes[CC_gene_df$CC_type %in% c('receptor','receptor_ligand')],CD_maker_genes)]='CC_receptor_CD_marker'
  mEry_tissue_markers$marker_type[mEry_tissue_markers$gene %in% intersect(CC_gene_df$CC_genes[CC_gene_df$CC_type %in% c('ligand')],CD_maker_genes)]='CC_ligand_CD_marker'
  
  
  top10_ordered_mEry_tissue_markers=mEry_tissue_markers %>% group_by(cluster) %>% top_n(n = 10,wt = marker_score1)
  other_df=mEry_tissue_markers[mEry_tissue_markers$marker_type!='no',]
  top10_ordered_mEry_tissue_markers=rbind(top10_ordered_mEry_tissue_markers,other_df[!other_df$gene %in% top10_ordered_mEry_tissue_markers$gene,])
  p=ggplot( mEry_tissue_markers,aes(x=avg_log2FC,y=marker_score1,color=marker_type))+geom_point()+facet_grid(~cluster)+theme_bw()+
    geom_text_repel(data = top10_ordered_mEry_tissue_markers,mapping = aes(x=avg_log2FC-0.1,y=marker_score1+0.1,label =gene),segment.color = "black",box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"),max.overlaps =10)#+ylim(c(-5,25))
  p
  return(list(mEry_tissue_markers,count_tmp_df,p))
}
subcelltype_gseGO_func=function(RNA_markers,keyType='SYMBOL'){
  subcelltype_gseGO_glist=list()
  subcelltype_gseGO_list=list()
  tmp_df=data.frame()
  for (type in unique(RNA_markers$celltype)) {
    print(type)
    subcelltype_gseGO_glist[[type]]=RNA_markers[RNA_markers$celltype==type,'avg_log2FC']
    names(subcelltype_gseGO_glist[[type]])=RNA_markers[RNA_markers$celltype==type,'gene']
    subcelltype_gseGO_glist[[type]]=sort(subcelltype_gseGO_glist[[type]],decreasing = T)
    subcelltype_gseGO_list[[type]] =gseGO(geneList = subcelltype_gseGO_glist[[type]],ont ='all' ,OrgDb = org.Hs.eg.db,keyType = keyType)
    
    if(dim(subcelltype_gseGO_list[[type]]@result)[1] >0){
      subcelltype_gseGO_list[[type]]@result$celltype=type
      tmp_df=rbind(tmp_df,subcelltype_gseGO_list[[type]]@result)}
  }
  return(list(subcelltype_gseGO_list,tmp_df))
}
subcelltype_enrichGO_func=function(RNA_markers){
  subcelltype_enrichGO_glist=list()
  subcelltype_enrichGO_list=list()
  RNA_markers=RNA_markers[RNA_markers$avg_log2FC>0,]
  for (type in unique(RNA_markers$celltype)) {
    cluster=as.character(unique(RNA_markers$cluster))
    subcelltype_enrichGO_glist[[type]][[ cluster[1] ]]=unique(RNA_markers[RNA_markers$cluster==cluster[1] & RNA_markers$celltype==type ,'gene'])
    subcelltype_enrichGO_glist[[type]][[ cluster[2] ]]=unique(RNA_markers[RNA_markers$cluster==cluster[2] & RNA_markers$celltype==type ,'gene'])
    subcelltype_enrichGO_list[[type]][[ cluster[1] ]] =enrichGO(gene =subcelltype_enrichGO_glist[[type]][[ cluster[1] ]] ,OrgDb =org.Hs.eg.db ,keyType = 'SYMBOL',ont = 'all')
    subcelltype_enrichGO_list[[type]][[ cluster[2] ]] =enrichGO(gene =subcelltype_enrichGO_glist[[type]][[ cluster[2] ]] ,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'all')
  }
  
  return(subcelltype_enrichGO_list)
}



# -----------------prepare membrane gene info ---------------------

if(F){
  all_merged_plasma_protein=readRDS('Protein_NRBC_marker/all_summary_membrane.RDS')
  length(all_merged_plasma_protein) # 3222
  
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # GO:0005886 
  pm_genes <- getBM(
    attributes = c('entrezgene_id', 'external_gene_name', 'uniprotswissprot', 
                   "go_id",                 # GO编号
                   "name_1006",             # GO术语名称（如"plasma membrane"）
                   "namespace_1003",        # GO命名空间（如"cellular_component"）
                   "go_linkage_type" ),  # **证据代码**（关键！如IDA、IEA）
    filters = 'go_parent_term', # 使用go_parent_term可以捕获所有下级term
    values = 'GO:0005886', # 范围宽泛。指整个质膜结构，包含内在膜蛋白、外周膜蛋白（仅附着在膜表面）、脂类以及任何与膜物理相关的蛋白
    mart = ensembl
  )
  
  #，GO:0005886（质膜）、GO:0097524（膜微区）和GO:0001533（眼睑脂质滴） 这三个术语之间的关系可以概括为：GO:0005886包含GO:0097524，而GO:0001533是一个在功能和组织上完全独立、具有组织特异性的细胞器
  # 查看去重后的基因数量
  length(unique(pm_genes$external_gene_name))
  pm_genes=pm_genes[pm_genes$external_gene_name!='',]
  head(pm_genes)
  dim(pm_genes)
  length(unique(pm_genes$external_gene_name))#5684
  table(all_merged_plasma_protein %in% pm_genes$external_gene_name) # F/T: 558/2671 
  length(unique(pm_genes$uniprotswissprot))#5605
  
  new_all_merged_plasma_protein=unique(c(all_merged_plasma_protein,unique(pm_genes$external_gene_name)))
  length(new_all_merged_plasma_protein)
  saveRDS(new_all_merged_plasma_protein,file = 'Protein_NRBC_marker/new_all_merged_plasma_protein.rds')
  
  
  
  # GO:0005886
 cell_surface <- getBM(
    attributes = c('entrezgene_id', 'external_gene_name', 'uniprotswissprot', 
                   "go_id",                 # GO编号
                   "name_1006",             # GO术语名称（如"plasma membrane"）
                   "namespace_1003",        # GO命名空间（如"cellular_component"）
                   "go_linkage_type" ),  # **证据代码**（关键！如IDA、IEA）
    filters = 'go', # 使用go_parent_term可以捕获所有下级term
    values = c('GO:0009986'), # 质膜
    mart = ensembl
  )
  
 table(unique(cell_surface$external_gene_name) %in% all_merged_plasma_protein ) # T/F: 469/215
 table(unique(cell_surface$external_gene_name) %in% new_all_merged_plasma_protein ) # T/F: 469/215 # F/T:84/600
 
 saveRDS(cell_surface,file = 'Protein_NRBC_marker/cell_surface.rds')
 
}else{
  new_all_merged_plasma_protein=readRDS('Protein_NRBC_marker/new_all_merged_plasma_protein.rds')
}


#--------------------------prepare TF & plasam membrane genes & CD makers data-----------------------#

TF_df_TFtargetdb=read.csv('ref_data/TF/TF-Target-information.txt',sep="\t")
TF_df_trrustdb=read.csv('ref_data/TF/trrust_rawdata.human.tsv',sep="\t")


library(CellChat)
CellChatDB <- CellChatDB.human
recepotors=unlist(sapply(unique(CellChatDB$interaction$receptor.symbol),function(x){strsplit(x,split = ', ')}))
recepotors=unique(as.character(recepotors));length(recepotors) #717

ligand_genes=unlist(sapply(unique(CellChatDB$interaction$ligand.symbol),function(x){strsplit(x,split = ', ')}))
ligand_genes=unique(as.character(ligand_genes));length(ligand_genes) #777

table(recepotors %in% ligand_genes) # shared 129 

#table(recepotors %in% all_merged_plasma_protein ) # F/T:59  /658
table(recepotors %in% new_all_merged_plasma_protein ) # F/T:23  /694

# 部分ligand 是脂溶性的，可以够自由穿过细胞膜，无需传统意义上的跨膜受体来识别并结合。如RetinoicAcid-RA

CC_gene_df=data.frame(CC_genes=unique(c(recepotors,ligand_genes)))
CC_gene_df$CC_type='receptor';CC_gene_df$CC_type[CC_gene_df$CC_genes %in% ligand_genes]='ligand';CC_gene_df$CC_type[CC_gene_df$CC_genes %in% intersect(ligand_genes,recepotors)]='receptor_ligand'
CC_gene_df$transmembrane='no' ; CC_gene_df$transmembrane[CC_gene_df$CC_genes %in% c(unlist(strsplit(unique(CellChatDB$interaction$receptor.symbol[CellChatDB$interaction$receptor.transmembrane==TRUE]),split = ', ')),
                                                                                    unlist(strsplit(unique(CellChatDB$interaction$ligand.symbol[CellChatDB$interaction$ligand.transmembrane==TRUE]),split = ', ')) )]='transm'
CC_gene_df$type=paste(CC_gene_df$CC_type,CC_gene_df$transmembrane,sep="_")

out_membrane_receptors=recepotors[!recepotors %in% new_all_merged_plasma_protein]
out_membrane_receptors_inf=data.frame()
for(  left_symbol in out_membrane_receptors){
  temp_df=unique(CellChatDB$interaction[ grep(left_symbol,CellChatDB$interaction$receptor.symbol),c('receptor.symbol','receptor.location')])
  temp_df$receptor=left_symbol
  out_membrane_receptors_inf=rbind(out_membrane_receptors_inf,temp_df)
}
#out_membrane_receptors_inf
rm(CellChatDB);gc()

CD_antibody_genes=mapIds(x = org.Hs.eg.db,keys = rownames(filt_NBRC_altas_seu),column = 'ALIAS',keytype = 'SYMBOL',multiVals ='list' )
CD_antibody_genes=sapply(CD_antibody_genes, function(x){x[grep('^CD[1-9]',x)]})
CD_antibody_genes=CD_antibody_genes[lengths(CD_antibody_genes)>0];length(CD_antibody_genes)# 431
CD_antibody_genes=CD_antibody_genes[-grep('AS1',names(CD_antibody_genes))]
length(CD_antibody_genes)# 425
length(as.character(CD_antibody_genes))
CD_antibody_genes[['SCARB1']] # 非CD marker 
all_CD_maker_genes=unique(c(names(CD_antibody_genes),as.character(unlist(CD_antibody_genes))) )



#################################################################################################################################################################
#--------------------------------------------------prepare data ----------------------------------------------------#
#################################################################################################################################################################

#--------------------------prepare NRBC data-----------------------#
if(F){
  YS_altas_Ery_seu=readRDS('NRBC_YS_altas/YS_altas_Ery_seu.rds')
  FL_altas_Ery_seu=readRDS('NRBC_FL_altas/tmp_FL_altas_Ery_seu.rds')
  BM_NRBC_altas_seu=readRDS('NRBC_BM_altas/res_data/BM_NRBC_altas_seu.rds')
  
  (DimPlot(YS_altas_Ery_seu,group.by = 'final_celltype',reduction = 'ref.umap',cols = cols)+ggtitle('YS'))+
    (DimPlot(FL_altas_Ery_seu,group.by = 'subcelltype',reduction = 'umap2',cols = cols)+ggtitle('FL'))
  DimPlot(BM_NRBC_altas_seu,group.by = 'new_celltype',reduction = 'ref.umap',cols = cols,split.by = 'stage')
  

  # early Ery: BFUE/CFUE, ProE
  # mid Ery: Bas
  # late Ery: Poly, Orth
  
  
  #统计membrane protein 对variable genes的贡献
  YS_altas_Ery_seu=FindVariableFeatures(YS_altas_Ery_seu,nfeatures = length(rownames(YS_altas_Ery_seu)))
  var_feature_meta_df=YS_altas_Ery_seu[["RNA"]]@meta.data[,61:68]
  var_feature_meta_df=var_feature_meta_df[order(var_feature_meta_df$var.features.rank,decreasing = F),]
  var_feature_meta_df$membrane='no'
  var_feature_meta_df$membrane[var_feature_meta_df$var.features %in% new_all_merged_plasma_protein]='yes'
  YS_top2000_var_mgenes=var_feature_meta_df[1:2000,];YS_top2000_var_mgenes=YS_top2000_var_mgenes[YS_top2000_var_mgenes$membrane=='yes','var.features']
  
  YS_mgene_invgene_ratio_df=data.frame(c(prop.table(table(var_feature_meta_df$membrane[1:500]))[2],prop.table(table(var_feature_meta_df$membrane[1:1000]))[2],
                                         prop.table(table(var_feature_meta_df$membrane[1:1500]))[2],prop.table(table(var_feature_meta_df$membrane[1:2000]))[2],
                                         prop.table(table(rownames(YS_altas_Ery_seu) %in%  new_all_merged_plasma_protein))[2]))
  
  rownames(YS_mgene_invgene_ratio_df)=c('mgene_ratio_top500vgene','mgene_ratio_top1000vgene','mgene_ratio_top1500vgene','mgene_ratio_top2000vgene','mgene_ratio_allgene')
  colnames(YS_mgene_invgene_ratio_df)='ratio';YS_mgene_invgene_ratio_df$ratio=round(YS_mgene_invgene_ratio_df$ratio,3)
  YS_mgene_invgene_ratio_df$top=c(500,1000,1500,2000,2100)
  YS_mgene_invgene_ratio_df$membrane='yes'
  
  p1=ggplot(var_feature_meta_df[1:2000,],aes(x=vf_vst_counts_rank,y=vf_vst_counts_variance.standardized,color=membrane))+geom_point()+ geom_vline(xintercept=c(500,1000,1500,2000),lty=6,col="black",lwd=0.5)+
    geom_text_repel(data = YS_mgene_invgene_ratio_df,mapping =aes(x=top,y=4,label=as.character(ratio)))+theme_classic()+ggtitle('YS')
  
  
  
  head(FL_altas_Ery_seu[["RNA"]]@meta.data)
  var_feature_meta_df=FL_altas_Ery_seu[["RNA"]]@meta.data
  var_feature_meta_df=var_feature_meta_df[order(var_feature_meta_df$vf_vst_counts.FL_rank,decreasing = F),]
  var_feature_meta_df$membrane='no'
  var_feature_meta_df$membrane[var_feature_meta_df$var.features %in% new_all_merged_plasma_protein]='yes'
  FL_top2000_var_mgenes=var_feature_meta_df[1:2000,];FL_top2000_var_mgenes=FL_top2000_var_mgenes[FL_top2000_var_mgenes$membrane=='yes','var.features']
  
  FL_mgene_invgene_ratio_df=data.frame(c(prop.table(table(var_feature_meta_df$membrane[1:500]))[2],prop.table(table(var_feature_meta_df$membrane[1:1000]))[2],
                                         prop.table(table(var_feature_meta_df$membrane[1:1500]))[2],prop.table(table(var_feature_meta_df$membrane[1:2000]))[2],
                                         prop.table(table(rownames(FL_altas_Ery_seu) %in%  new_all_merged_plasma_protein))[2]))
  
  rownames(FL_mgene_invgene_ratio_df)=c('mgene_ratio_top500vgene','mgene_ratio_top1000vgene','mgene_ratio_top1500vgene','mgene_ratio_top2000vgene','mgene_ratio_allgene')
  colnames(FL_mgene_invgene_ratio_df)='ratio';FL_mgene_invgene_ratio_df$ratio=round(FL_mgene_invgene_ratio_df$ratio,3)
  FL_mgene_invgene_ratio_df$top=c(500,1000,1500,2000,2100)
  FL_mgene_invgene_ratio_df$membrane='yes'
  
  p2=ggplot(var_feature_meta_df[1:2000,],aes(x=vf_vst_counts.FL_rank,y=vf_vst_counts.FL_variance.standardized,color=membrane))+geom_point()+ geom_vline(xintercept=c(500,1000,1500,2000),lty=6,col="black",lwd=0.5)+
    geom_text_repel(data = FL_mgene_invgene_ratio_df,mapping =aes(x=top,y=60,label=as.character(ratio)))+theme_classic()+ggtitle('FL')
  
  # 统计一下FL 中的双阳性细胞
  HBB_HBE1_Ery_seu=subset(FL_altas_Ery_seu,HBE1 >1 & HBB >1);HBB_HBE1_Ery_seu# 1924, FL 存在HBE1—>HBB转变的过程, 以后可以从这个数据中研究转变过程
  subset(FL_altas_Ery_seu,HBE1 >0 & HBB <1) # 2691
  table(HBB_HBE1_Ery_seu$id)# 主要来自CS18_F61样本，该样本也是YS——Ery的主要来源
  # 由高转低，在小鼠中发现FL 存在YS-EMP来源的Ery，但是人里面是缺乏的，我们之前的FL数据显示。这是跟小鼠的不同点，YS-EMP主要产生了mon-derived Max。
  table( subset(FL_altas_Ery_seu,HBE1 >0 & HBB <1)$id)/table(FL_altas_Ery_seu$id)

  head(BM_NRBC_altas_seu)
  FBM_NRBC_altas_seu=subset(BM_NRBC_altas_seu,stage=='FBM')
  FBM_NRBC_altas_seu=FindVariableFeatures(FBM_NRBC_altas_seu,nfeatures =length(rownames(FBM_NRBC_altas_seu)))
  head(FBM_NRBC_altas_seu[["RNA"]]@meta.data)
  var_feature_meta_df=FBM_NRBC_altas_seu[["RNA"]]@meta.data
  var_feature_meta_df=var_feature_meta_df[order(var_feature_meta_df$var.features.rank,decreasing = F),]
  var_feature_meta_df=var_feature_meta_df[var_feature_meta_df$vf_vst_counts_mean>0,]
  var_feature_meta_df$membrane='no'
  var_feature_meta_df$membrane[var_feature_meta_df$var.features %in% new_all_merged_plasma_protein]='yes'
  FBM_top2000_var_mgenes=var_feature_meta_df[1:2000,];FBM_top2000_var_mgenes=FBM_top2000_var_mgenes[FBM_top2000_var_mgenes$membrane=='yes','var.features']
  
  FBM_mgene_invgene_ratio_df=data.frame(c(prop.table(table(var_feature_meta_df$membrane[1:500]))[2],prop.table(table(var_feature_meta_df$membrane[1:1000]))[2],
                                          prop.table(table(var_feature_meta_df$membrane[1:1500]))[2],prop.table(table(var_feature_meta_df$membrane[1:2000]))[2],
                                          prop.table(table(rownames(FBM_NRBC_altas_seu) %in%  new_all_merged_plasma_protein))[2]))
  
  rownames(FBM_mgene_invgene_ratio_df)=c('mgene_ratio_top500vgene','mgene_ratio_top1000vgene','mgene_ratio_top1500vgene','mgene_ratio_top2000vgene','mgene_ratio_allgene')
  colnames(FBM_mgene_invgene_ratio_df)='ratio';FBM_mgene_invgene_ratio_df$ratio=round(FBM_mgene_invgene_ratio_df$ratio,3)
  FBM_mgene_invgene_ratio_df$top=c(500,1000,1500,2000,2100)
  FBM_mgene_invgene_ratio_df$membrane='yes'
  
  p3=ggplot(var_feature_meta_df[1:2000,],aes(x=var.features.rank,y=vf_vst_counts_variance.standardized,color=membrane))+geom_point()+ geom_vline(xintercept=c(500,1000,1500,2000),lty=6,col="black",lwd=0.5)+
    geom_text_repel(data = FBM_mgene_invgene_ratio_df,mapping =aes(x=top,y=30,label=as.character(ratio)))+theme_classic()+ggtitle('FBM')
  
  
  
  ABM_NRBC_altas_seu=subset(BM_NRBC_altas_seu,stage=='ABM')
  ABM_NRBC_altas_seu[['RNA']]=JoinLayers(ABM_NRBC_altas_seu[['RNA']])
  ABM_NRBC_altas_seu=FindVariableFeatures(ABM_NRBC_altas_seu,nfeatures =length(rownames(ABM_NRBC_altas_seu)))
  
  head(ABM_NRBC_altas_seu[["RNA"]]@meta.data)
  var_feature_meta_df=ABM_NRBC_altas_seu[["RNA"]]@meta.data
  var_feature_meta_df=var_feature_meta_df[order(var_feature_meta_df$var.features.rank,decreasing = F),]
  var_feature_meta_df=var_feature_meta_df[var_feature_meta_df$vf_vst_counts_mean>0,]
  var_feature_meta_df$membrane='no'
  var_feature_meta_df$membrane[var_feature_meta_df$var.features %in% new_all_merged_plasma_protein]='yes'
  ABM_top2000_var_mgenes=var_feature_meta_df[1:2000,];ABM_top2000_var_mgenes=ABM_top2000_var_mgenes[ABM_top2000_var_mgenes$membrane=='yes','var.features']
  
  ABM_mgene_invgene_ratio_df=data.frame(c(prop.table(table(var_feature_meta_df$membrane[1:500]))[2],prop.table(table(var_feature_meta_df$membrane[1:1000]))[2],
                                          prop.table(table(var_feature_meta_df$membrane[1:1500]))[2],prop.table(table(var_feature_meta_df$membrane[1:2000]))[2],
                                          prop.table(table(rownames(ABM_NRBC_altas_seu) %in%  new_all_merged_plasma_protein))[2]))
  
  rownames(ABM_mgene_invgene_ratio_df)=c('mgene_ratio_top500vgene','mgene_ratio_top1000vgene','mgene_ratio_top1500vgene','mgene_ratio_top2000vgene','mgene_ratio_allgene')
  colnames(ABM_mgene_invgene_ratio_df)='ratio';ABM_mgene_invgene_ratio_df$ratio=round(ABM_mgene_invgene_ratio_df$ratio,3)
  ABM_mgene_invgene_ratio_df$top=c(500,1000,1500,2000,2100)
  ABM_mgene_invgene_ratio_df$membrane='yes'
  
  p4=ggplot(var_feature_meta_df[1:2000,],aes(x=var.features.rank,y=vf_vst_counts_variance.standardized,color=membrane))+geom_point()+ geom_vline(xintercept=c(500,1000,1500,2000),lty=6,col="black",lwd=0.5)+
    geom_text_repel(data = ABM_mgene_invgene_ratio_df,mapping =aes(x=top,y=30,label=as.character(ratio)))+theme_classic()+ggtitle('ABM')
  
  
  p=p1+p2+p3+p4;p
  ggsave(p,filename = 'Protein_NRBC_marker/mgene_invariable_gene_analysis.pdf',width = 12,height = 10)
  
  rm(ABM_NRBC_altas_seu,FBM_NRBC_altas_seu)
  
  library(VennDiagram)
  dev.off()
  p=venn.diagram(x =list('YS'=YS_top2000_var_mgenes,'FL'=FL_top2000_var_mgenes,'FBM'=FBM_top2000_var_mgenes,'ABM'=ABM_top2000_var_mgenes),
                 filename = NULL,fill=cols[1:4], scaled = T,main='membrane genes in top2000 variable genes',force.unique = T,main.cex = 2,sub.cex = 2,total.population = T)
  grid.draw(p)
  ggsave(p,file='Protein_NRBC_marker/NRBC_altas_mgene_in_top2000vgene_venn.pdf',width = 6,height = 6)
  inter <- get.venn.partitions(list('YS'=YS_top2000_var_mgenes,'FL'=FL_top2000_var_mgenes,'FBM'=FBM_top2000_var_mgenes,'ABM'=ABM_top2000_var_mgenes))
  saveRDS(list('YS'=YS_top2000_var_mgenes,'FL'=FL_top2000_var_mgenes,'FBM'=FBM_top2000_var_mgenes,'ABM'=ABM_top2000_var_mgenes),file = 'Protein_NRBC_marker/tissue_NRBC_mvar_genes_list.rds')
  
  
  
  NBRC_altas_seu=merge(YS_altas_Ery_seu,c(FL_altas_Ery_seu,BM_NRBC_altas_seu));
  NBRC_altas_seu[['prediction.score.celltype']]=NULL
  NBRC_altas_seu[['RNA']]=JoinLayers(NBRC_altas_seu[['RNA']])
  NBRC_altas_seu[['umap']]=merge(YS_altas_Ery_seu[['ref.umap']],c(FL_altas_Ery_seu[['umap2']],BM_NRBC_altas_seu[['ref.umap']]))
  
  NBRC_altas_seu$tissue_stage=NBRC_altas_seu$stage
  NBRC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('tissue_stage')]='YS'
  NBRC_altas_seu@meta.data[rownames(FL_altas_Ery_seu@meta.data),c('tissue_stage')]='FL'
  
  NBRC_altas_seu$final_celltype=as.character(NBRC_altas_seu$new_celltype)
  NBRC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('final_celltype')]=as.character(YS_altas_Ery_seu$final_celltype)
  NBRC_altas_seu@meta.data[rownames(FL_altas_Ery_seu@meta.data),c('final_celltype')]=as.character(FL_altas_Ery_seu$subcelltype)
  table(NBRC_altas_seu$final_celltype)
  
  # 不考虑 FL/FBM的YS_NRBC
  filt_cells=rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$tissue_stage %in% c('FL','FBM') & NBRC_altas_seu$final_celltype %in% c("YS_Bas/Poly","YS_Orth"),])
  cho_cells=rownames(NBRC_altas_seu@meta.data)[!rownames(NBRC_altas_seu@meta.data) %in%  filt_cells ];length(cho_cells)
  filt_NBRC_altas_seu=subset(NBRC_altas_seu,cells=cho_cells)
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu, tissue_stage!='EBM')
  filt_NBRC_altas_seu$tissue_stage=factor(filt_NBRC_altas_seu$tissue_stage,levels = c('YS','FL','FBM','ABM'))
  filt_NBRC_altas_seu$final_celltype=factor(filt_NBRC_altas_seu$final_celltype,levels = c("BFUE/CFUE","ProE","Bas","Poly" ,"Orth"))
  DimPlot(filt_NBRC_altas_seu,group.by = 'final_celltype',split.by = 'tissue_stage',cols = cols,raster=FALSE)/
    FeaturePlot(filt_NBRC_altas_seu,features = 'HBE1',split.by = 'tissue_stage')
  rm(NBRC_altas_seu);gc()
  
  
  filt_NBRC_altas_seu$Ery_stage='late_Ery'
  filt_NBRC_altas_seu$Ery_stage[filt_NBRC_altas_seu$final_celltype %in% c('BFUE/CFUE','ProE')]='early_Ery'
  filt_NBRC_altas_seu$Ery_stage[filt_NBRC_altas_seu$final_celltype %in% c('Bas')]='mid_Ery'
  filt_NBRC_altas_seu$Ery_stage=factor(filt_NBRC_altas_seu$Ery_stage,levels = c('early_Ery','mid_Ery','late_Ery'))
  filt_NBRC_altas_seu$tissue_stage=factor(filt_NBRC_altas_seu$tissue_stage,levels = c('YS','FL','FBM','ABM'))
  
  table(filt_NBRC_altas_seu@meta.data[,c('Ery_stage','final_celltype')])
  
  rm(YS_altas_Ery_seu,FL_altas_Ery_seu,FBM_NRBC_altas_seu,BM_NRBC_altas_seu);gc()
  
  Idents(filt_NBRC_altas_seu)='tissue_stage'
  
  
  # FL/FBM 存在HBE1->HBB的转换
  temp_seu= subset(filt_NBRC_altas_seu, tissue_stage %in% c('FL','FBM'))
  filt_cells2=WhichCells(temp_seu,expression =  HBE1 >2,slot = 'data')
  length(filt_cells2) # HBE1 >1 : 2245, HBE1 >2: 841 
  DimPlot(subset(filt_NBRC_altas_seu,cells=filt_cells2),group.by = 'Ery_stage',split.by = 'tissue_stage',cols = cols,raster=FALSE)+
    FeaturePlot(subset(filt_NBRC_altas_seu,cells=filt_cells2),features = 'HBE1',split.by = 'tissue_stage')
  cho_cells2=rownames(filt_NBRC_altas_seu@meta.data)[!rownames(filt_NBRC_altas_seu@meta.data) %in%  filt_cells2 ];length(cho_cells2)
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,cells=cho_cells2)
  rm(temp_seu,cho_cells,cho_cells2,filt_cells,filt_cells2);gc()
  
  sort(unique(filt_NBRC_altas_seu$final_celltype))
  filt_NBRC_altas_seu$age[is.na(filt_NBRC_altas_seu$age)]=filt_NBRC_altas_seu$stage[is.na(filt_NBRC_altas_seu$age)]
  
  colnames(filt_NBRC_altas_seu@meta.data)=gsub(pattern = 'resource',replacement ='source' ,colnames(filt_NBRC_altas_seu@meta.data))
  filt_NBRC_altas_seu$source[is.na(filt_NBRC_altas_seu$source)]=filt_NBRC_altas_seu$orig.dataset[is.na(filt_NBRC_altas_seu$source)]
  filt_NBRC_altas_seu$age[filt_NBRC_altas_seu$tissue_stage=='YS']=filt_NBRC_altas_seu$stage[filt_NBRC_altas_seu$tissue_stage=='YS']
  
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,source !='GSE253355')
  
  filt_NBRC_altas_seu_meta=filt_NBRC_altas_seu@meta.data
  filt_NBRC_altas_seu_meta=filt_NBRC_altas_seu_meta[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'component', 'age',  'sex', 'sort.ids', 'fetal.ids','orig.dataset',
                                                       'sequencing.type','id','tissue_stage', 'final_celltype', 'Ery_stage','source', 'source_celltype')]
  
  
  # 计算不同组织微环境稳定造血时期
  YS_CS15_CS18_sample_df=filt_NBRC_altas_seu_meta[filt_NBRC_altas_seu_meta$tissue_stage=='YS' & filt_NBRC_altas_seu_meta$age %in% c('CS15', "CS17" ,"CS18" ),]
  YS_celltype_mratio_df=data.frame( colMeans(prop.table(table(YS_CS15_CS18_sample_df[,c('id','final_celltype')]),margin = 1)))
  colnames(YS_celltype_mratio_df)='YS'
  
  FL_CS18post_sample_df=filt_NBRC_altas_seu_meta[filt_NBRC_altas_seu_meta$tissue_stage=='FL' & !filt_NBRC_altas_seu_meta$age %in% c('CS14_4PCW', "CS15_5PCW" ,"CS17_6PCW" ),]
  FL_celltype_mratio_df=data.frame( colMeans(prop.table(table(FL_CS18post_sample_df[,c('id','final_celltype')]),margin = 1)))
  colnames(FL_celltype_mratio_df)='FL'
  
  FBM_sample_df=filt_NBRC_altas_seu_meta[filt_NBRC_altas_seu_meta$tissue_stage=='FBM',]
  FBM_celltype_mratio_df=data.frame( colMeans(prop.table(table(FBM_sample_df[,c('id','final_celltype')]),margin = 1)))
  colnames(FBM_celltype_mratio_df)='FBM'
  
  ABM_sample_df=filt_NBRC_altas_seu_meta[filt_NBRC_altas_seu_meta$tissue_stage=='ABM',]
  table(ABM_sample_df$source)# CD34- CD235+: GSE150774:CD34+,GSE133181,GSE135194,GSE169426,  Mononuclear cells:GSE165645, Mononuclear:CD34+=4:1:GSE181989
  ABM_celltype_mratio_CD34pos_df=data.frame( round(colMeans(table(ABM_sample_df[!ABM_sample_df$source %in% c('GSE150774','GSE181989','GSE165645'),c('id','final_celltype')])),digits = 0) )
  colnames(ABM_celltype_mratio_CD34pos_df)='mcount'
  ABM_celltype_mratio_GYPAposCD34neg_df=data.frame( round(colMeans(table(ABM_sample_df[ABM_sample_df$source=='GSE150774',c('id','final_celltype')])),digits = 0))
  colnames(ABM_celltype_mratio_GYPAposCD34neg_df)='mcount'
  ABM_celltype_mratio_GYPAposCD34neg_df['BFUE/CFUE','mcount']=0
  
  ABM_NRBC_ratio_df=prop.table(ABM_celltype_mratio_CD34pos_df+ ABM_celltype_mratio_GYPAposCD34neg_df[rownames(ABM_celltype_mratio_CD34pos_df),])
  colnames(ABM_NRBC_ratio_df)='ABM'
  
  YS_celltype_mratio_df['BFUE/CFUE','YS']=0
  NRBC_ratio_df= cbind(YS_celltype_mratio_df,cbind(cbind(FL_celltype_mratio_df,FBM_celltype_mratio_df),ABM_NRBC_ratio_df)[rownames(YS_celltype_mratio_df),])
  NRBC_ratio_df$celltytpe=rownames(NRBC_ratio_df)
  NRBC_ratio_df= melt(NRBC_ratio_df)
  
  NRBC_ratio_df$celltytpe=factor(NRBC_ratio_df$celltytpe,levels = c("BFUE/CFUE" ,"ProE","Bas","Poly","Orth"  ))
  p=ggplot(NRBC_ratio_df ,aes(x=variable,y=value,fill=celltytpe))+geom_bar(position = "fill",stat = 'identity')+
    scale_fill_manual(values =  cols)+theme_classic()+theme(axis.text.x = element_text(hjust =1 ,angle = 45,face = 'bold'))
  p
  ggsave(p,file='res_pic/main_figure1/tissue_celltype_ratio_barplot.pdf',width =6 ,height = 6)
  
  
  
  # 之前由于symbo不一致导致的注释别名和基因symbol重复，临时纠正情况
  if(F){
    temp_data=GetAssayData(filt_NBRC_altas_seu,assay = 'RNA',layer = 'counts')
    cho_temp_data=temp_data[c('KIAA1524','KIAA0922', 'C1orf27', 'C17orf62','CIP2A','TMEM131L','ODR4','CYBC1'), ]
    cho_temp_data['CIP2A',]=cho_temp_data['CIP2A',]+cho_temp_data['KIAA1524',]
    cho_temp_data['TMEM131L',]=cho_temp_data['TMEM131L',]+cho_temp_data['KIAA0922',]
    cho_temp_data['ODR4',]=cho_temp_data['ODR4',]+cho_temp_data['C1orf27',]
    cho_temp_data['CYBC1',]=cho_temp_data['CYBC1',]+cho_temp_data['C17orf62',]
    cho_temp_data=cho_temp_data[c('CIP2A','TMEM131L','ODR4','CYBC1'),]
    
    temp_data=temp_data[!rownames(temp_data) %in% c('KIAA1524','KIAA0922', 'C1orf27', 'C17orf62','CIP2A','TMEM131L','ODR4','CYBC1') ,]
    temp_data=rbind(temp_data,cho_temp_data)
    
    filt_NBRC_altas_seu<-CreateSeuratObject(counts =temp_data,meta.data = filt_NBRC_altas_seu@meta.data )
    filt_NBRC_altas_seu=NormalizeData(filt_NBRC_altas_seu) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
    rm(temp_data);gc()
    
    # ATP5 家族基因注释不一致,做差异分析时候得到，校正后发现是因为基因注释信息导致的，故校正
    ATP5_genes=rownames(filt_NBRC_altas_seu)[grep('^ATP5',rownames(filt_NBRC_altas_seu))]
    ATP5_genes=mapIds(x = org.Hs.eg.db,keys = ATP5_genes,keytype = 'ALIAS',column = 'SYMBOL' )
    cor_ATP5_genes=ATP5_genes[names(ATP5_genes)!=ATP5_genes]
    table(as.character(cor_ATP5_genes)  %in% rownames(filt_NBRC_altas_seu))
    
    cor_ATP5_genes_count=GetAssayData(object = filt_NBRC_altas_seu,assay = 'RNA',layer = 'count')[names(cor_ATP5_genes),]
    rownames(cor_ATP5_genes_count)=as.character(cor_ATP5_genes)
    
    cor_ATP5_genes_count2=GetAssayData(object = filt_NBRC_altas_seu,assay = 'RNA',layer = 'count')[as.character(cor_ATP5_genes),]
    cor_ATP5_genes_count2=cor_ATP5_genes_count2+cor_ATP5_genes_count
    
    out_ATP5_genes_count2=GetAssayData(object = filt_NBRC_altas_seu,assay = 'RNA',layer = 'count')
    out_ATP5_genes_count2=out_ATP5_genes_count2[rownames(out_ATP5_genes_count2)[!rownames(out_ATP5_genes_count2) %in% c(as.character(cor_ATP5_genes),names(cor_ATP5_genes))],]
    count=rbind(cor_ATP5_genes_count2,out_ATP5_genes_count2);rm(out_ATP5_genes_count2,cor_ATP5_genes_count2,cor_ATP5_genes_count);gc()
    table(rownames(count) %in% names(cor_ATP5_genes))
    
    new_filt_NBRC_altas_seu=CreateSeuratObject(counts = GetAssayData(filt_NBRC_altas_seu,assay = 'RNA',layer = 'counts'),meta.data = filt_NBRC_altas_seu@meta.data,min.cells = 10)
    new_filt_NBRC_altas_seu[['umap']]=filt_NBRC_altas_seu[['umap']]
    rm(filt_NBRC_altas_seu);gc()
    
    filt_NBRC_altas_seu=new_filt_NBRC_altas_seu;rm(new_filt_NBRC_altas_seu);gc()
    
    rm(count);gc()
    saveRDS(filt_NBRC_altas_seu,file = 'temp_filt_NBRC_altas_seu.rds')
    
    org_hsa_gene_symbols=mapIds(x = org.Hs.eg.db,keys = rownames(filt_NBRC_altas_seu),keytype = 'ALIAS',column = 'SYMBOL')
    table(as.character(org_hsa_gene_symbols)==names(org_hsa_gene_symbols))
    cor_org_hsa_gene_symbols=org_hsa_gene_symbols[as.character(org_hsa_gene_symbols)!=names(org_hsa_gene_symbols)]
    cor_org_hsa_gene_symbols=cor_org_hsa_gene_symbols[!is.na(cor_org_hsa_gene_symbols)]
    length(cor_org_hsa_gene_symbols)
    table(duplicated(as.character(cor_org_hsa_gene_symbols))) # 75 TRUE, org.Hs.eg.db注释也不完全准确， CGB5/CGB7/CGB8:是不同的基因，但是symbol信息都是CGB3，同一个蛋白不同的亚基组成 ,出现重复的symbol，排除
    cor_org_hsa_gene_symbols=cor_org_hsa_gene_symbols[!cor_org_hsa_gene_symbols %in% as.character(cor_org_hsa_gene_symbols[duplicated(cor_org_hsa_gene_symbols)])]
    length(cor_org_hsa_gene_symbols)
    
    table(as.character(cor_org_hsa_gene_symbols) %in%  rownames(filt_NBRC_altas_seu))# 1114 TRUE
    cor_org_hsa_gene_symbols=cor_org_hsa_gene_symbols[cor_org_hsa_gene_symbols %in%  rownames(filt_NBRC_altas_seu)]
    du_tmp_assay1=GetAssayData(filt_NBRC_altas_seu,layer='counts')[names(cor_org_hsa_gene_symbols),]
    rownames(du_tmp_assay1)=as.character(cor_org_hsa_gene_symbols)
    du_tmp_assay2=GetAssayData(filt_NBRC_altas_seu,layer='counts')[as.character(cor_org_hsa_gene_symbols),]
    
    du_tmp_assay=du_tmp_assay1+du_tmp_assay2;rm(du_tmp_assay2,du_tmp_assay1);gc()
    
    left_tmp_assay=GetAssayData(filt_NBRC_altas_seu,layer='counts')
    left_tmp_assay=left_tmp_assay[!rownames(left_tmp_assay) %in%  c(names(cor_org_hsa_gene_symbols),as.character(cor_org_hsa_gene_symbols)),]
    left_tmp_assay=rbind(left_tmp_assay,du_tmp_assay)
    new_filt_NBRC_altas_seu=CreateSeuratObject(counts = left_tmp_assay,meta.data = filt_NBRC_altas_seu@meta.data,min.cells = 10)
    new_filt_NBRC_altas_seu[['umap']]=filt_NBRC_altas_seu[['umap']]
    Idents(new_filt_NBRC_altas_seu)='tissue_stage'
    rm(left_tmp_assay,du_tmp_assay);gc()
    filt_NBRC_altas_seu=new_filt_NBRC_altas_seu;rm(new_filt_NBRC_altas_seu)
    
    counts=GetAssayData(filt_NBRC_altas_seu,layer = 'counts')
    counts['MT-ATP8',]=counts['MT-ATP8',]+counts['ATP8',];counts['MT-ATP6',]=counts['MT-ATP6',]+counts['ATP6',]
    counts=counts[!rownames(counts) %in%  c('ATP6','ATP8'),]
    SetAssayData(filt_NBRC_altas_seu)
    
  }
    
  VlnPlot(filt_NBRC_altas_seu,features =as.character(cor_HIST_genes),stack = T )+NoLegend()
  VlnPlot(filt_NBRC_altas_seu,features =names(cor_HIST_genes)[!names(cor_HIST_genes) %in% c('HIST1H2BA','HIST1H4G')],stack = T )+NoLegend()
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  filt_NBRC_altas_seu <- CellCycleScoring(filt_NBRC_altas_seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  filt_NBRC_altas_seu$source_celltype=paste(filt_NBRC_altas_seu$tissue_stage,filt_NBRC_altas_seu$final_celltype,sep = '_')
  filt_NBRC_altas_seu$source_celltype=factor(filt_NBRC_altas_seu$source_celltype,levels = c( "YS_ProE","YS_Bas","YS_Poly","YS_Orth", "FL_BFUE/CFUE" ,"FL_ProE", "FL_Bas","FL_Poly","FL_Orth",
                                                                                             "FBM_BFUE/CFUE" ,"FBM_ProE","FBM_Bas","FBM_Poly","FBM_Orth", "ABM_BFUE/CFUE","ABM_ProE","ABM_Bas","ABM_Poly","ABM_Orth"))
  filt_NBRC_altas_seu$final_celltype=factor(filt_NBRC_altas_seu$final_celltype,levels = c('BFUE/CFUE','ProE','Bas','Poly','Orth'))
  filt_NBRC_altas_seu$sourceid=paste(filt_NBRC_altas_seu$tissue_stage,filt_NBRC_altas_seu$id,sep = '_')
  
  # 少了 GSE253355 Ery数据,其中Ery 包含达IGKC表达细胞：2970,GSE253355 样本的年龄都是50-80岁的老人为主，占比很高,剔除这部分样本数据
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,subset=IGKC <=1 ) #因为要比较的是健康青壮年，与老人的差异显著，暂时不做此分析研究比较，专门做一个论文课题分析
  IGKC_nRBC_seu==subset(filt_NBRC_altas_seu,subset=IGKC <=1 )
  saveRDS(filt_NBRC_altas_seu,file = '20251125_filt_NBRC_altas_seu.rds')
  
  #saveRDS(filt_NBRC_altas_seu,file = 'Protein_NRBC_marker/merged_NRBC_data/filt_NBRC_altas_seu.rds') 
}else{
  #filt_NBRC_altas_seu=readRDS('Protein_NRBC_marker/merged_NRBC_data/filt_NBRC_altas_seu.rds')# 
  filt_NBRC_altas_seu=readRDS('20251125_filt_NBRC_altas_seu.rds')
}


rownames(filt_NBRC_altas_seu)[grep('^HB',rownames(filt_NBRC_altas_seu))]
HBGenes=c("HBA1","HBA2","HBB", "HBD" ,"HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
VlnPlot(filt_NBRC_altas_seu,features =HBGenes,group.by = 'source_celltype',stack = T,cols = cols)+NoLegend()




#################################################################################################################################################################
#--------------------------------------------------figure 1: hematopoietic niche shape the differenation of nRBC---------------------------------------------------#
#################################################################################################################################################################
# nRBC from different niche, are quite different and 异质性的

filt_NBRC_altas_seu[['prediction.score.celltype']]=NULL
p=pheatmap(cor(data.frame(AverageExpression(filt_NBRC_altas_seu,group.by = 'tissue_stage',layer = 'data' )$RNA)),main = 'based on all genes')# ,6 X 6,Protein_NRBC_marker/res_pic/main_figure1/tissue_basedonallgene_cor_heatmap.pdf
ggsave(as.ggplot(p),filename = 'res_pic/main_figure1/tissue_basedonallgene_cor_heatmap.pdf',width = 6,height = 6,dpi = 300)

p=pheatmap(cor(data.frame(AverageExpression(filt_NBRC_altas_seu,group.by = 'tissue_stage',layer = 'data',features =VariableFeatures(filt_NBRC_altas_seu) )$RNA)),main = 'based on variable genes')# ,6 X 6,Protein_NRBC_marker/res_pic/main_figure1/tissue_basedonvariable_cor_heatmap.pdf
ggsave(as.ggplot(p),filename = 'res_pic/main_figure1/tissue_basedonvariable_cor_heatmap.pdf',width = 6,height = 6,dpi = 300)


YS_celltype_mexp_df=data.frame(AverageExpression(subset(filt_NBRC_altas_seu,tissue_stage=='YS'),  layer = 'data',group.by = 'final_celltype')$RNA)
FL_celltype_mexp_df=data.frame(AverageExpression(subset(filt_NBRC_altas_seu,tissue_stage=='FL'),  layer = 'data',group.by = 'final_celltype')$RNA)
FBM_celltype_mexp_df=data.frame(AverageExpression(subset(filt_NBRC_altas_seu,tissue_stage=='FBM'),layer = 'data',group.by = 'final_celltype')$RNA)
ABM_celltype_mexp_df=data.frame(AverageExpression(subset(filt_NBRC_altas_seu,tissue_stage=='ABM'),layer = 'data',group.by = 'final_celltype')$RNA)

colnames(YS_celltype_mexp_df)=paste0('YS_',colnames(YS_celltype_mexp_df))
colnames(FL_celltype_mexp_df)=paste0('FL_',colnames(FL_celltype_mexp_df))
colnames(FBM_celltype_mexp_df)=paste0('FBM_',colnames(FBM_celltype_mexp_df))
colnames(ABM_celltype_mexp_df)=paste0('ABM_',colnames(ABM_celltype_mexp_df))

celltype_mexp_df=cbind(cbind(cbind(YS_celltype_mexp_df,FL_celltype_mexp_df),FBM_celltype_mexp_df),ABM_celltype_mexp_df)
rm(YS_celltype_mexp_df,FL_celltype_mexp_df,FBM_celltype_mexp_df,ABM_celltype_mexp_df);gc()
colnames(celltype_mexp_df)=c(paste0('fetal_',colnames(celltype_mexp_df)[1:14]),paste0('adult_',colnames(celltype_mexp_df)[15:19]))

cor_celltype_mexp_df=cor(celltype_mexp_df)
p=pheatmap(cor_celltype_mexp_df,cluster_rows = T,cluster_cols = T,color =colorRampPalette(colors = c('white','#C31E1F'))(100))
ggsave(as.ggplot(p),filename = 'res_pic/main_figure1/tissue_finalcelltype_basedonallgene_cor_heatmap.pdf',width = 10,height = 10,dpi = 300)

p=pheatmap(cor(celltype_mexp_df[VariableFeatures(filt_NBRC_altas_seu),]),cluster_rows = T,main = 'based on top 3000 var-genes',cluster_cols = T,color =colorRampPalette(colors = c('white','#C31E1F'))(100))
ggsave(as.ggplot(p),filename = 'res_pic/main_figure1/tissue_finalcelltype_basedonvarablegene_cor_heatmap.pdf',width = 10,height = 10,dpi = 300)

filt_NBRC_altas_seu$source_celltype=paste(filt_NBRC_altas_seu$tissue_stage,filt_NBRC_altas_seu$final_celltype,sep = ":")

# NRBC reference marker
p=VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =c('KIT','GATA1','KLF1','TFRC','GYPA','CCNB1','NCL'),stack = T,cols = col,split.by = 'tissue_stage')+NoLegend();p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure1/classical_nRBC_marker_expression_in_niches_subceltlype_vlnplot.pdf',width = 8,height = 8)



#################################################################################################################################################################
##----------------------------------------------------Figure2 primitive vs definitive nRBC-------------------------# 
#################################################################################################################################################################

filt_NBRC_altas_seu$NRBC_type='definitive'
filt_NBRC_altas_seu$NRBC_type[filt_NBRC_altas_seu$tissue_stage=='YS']='primitive'
Idents(filt_NBRC_altas_seu)='NRBC_type'
pd_tissue_NRBC_DE_res=FindAllMarkers(filt_NBRC_altas_seu,only.pos = T)
saveRDS(pd_tissue_NRBC_DE_res,file = 'Protein_NRBC_marker/res_data/main_figure2/primary_pd_tissue_NRBC_DE_res.rds')

pd_enrichgo_res=compareCluster(geneClusters =list('primitive'=unique(pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$cluster=='primitive','gene']),
                                                'definitive'=unique(pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$cluster=='definitive','gene'])),
                                               keyType = 'SYMBOL',ont = "BP",fun = 'enrichGO',  OrgDb='org.Hs.eg.db' )
                            
saveRDS(pd_enrichgo_res,file = 'Protein_NRBC_marker/res_data/main_figure2/whole_tissue_pd_NRBC_marker_enrichgo_res.rds')

p=dotplot(pd_enrichgo_res,showCategory=10);p


# ggplot美化, 不同niche来源nRBC具有不同的特征状态
top10_enrichgo_res_df=pd_enrichgo_res@compareClusterResult %>% group_by(Cluster) %>% do(head(.,10))
top10_enrichgo_res_df$ratio=as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[1,])/as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[2,])
top10_enrichgo_res_df$Description=factor(top10_enrichgo_res_df$Description,levels =unique(top10_enrichgo_res_df$Description) )
p=ggplot(top10_enrichgo_res_df,aes(x=Cluster,y=Description,color=-log10(p.adjust),size=Count))+geom_point()+theme_bw()+scale_color_gradient(low = '#4387B5',high = 'firebrick3')+theme(text = element_text(face = 'bold'))
p 

# -----------采用emapplot 全局比较 展示------------#
library(ggtree);library(enrichplot)
pd_enrichgo_res2=pairwise_termsim(pd_enrichgo_res)
emapplot(pd_enrichgo_res2, node_label = "group")
p=emapplot(pd_enrichgo_res2, node_label = "category");p
ggsave(p,filename = 'Protein_NRBC_marker/res_pic/main_figure2/pd_whole_enrichGO_emapplot.pdf',width = 6,height = 8)

p=emapplot(pd_enrichgo_res2, node_label = "group")
table(p$data$color2) # group information :positive cell-cell proliferation adhesion ,cellular macroautophagy autophagy component ,endosomal vesicle Golgi modificatio, splicing DNA replication via 

Idents(filt_NBRC_altas_seu)='source_celltype'

# primitive: autophagy related pathways, 囊泡运输系统，RNA转录剪切，获取其核心驱动基因:高差异基因
#（1）RNA转录剪切
#c('RNA splicing, via transesterification reactions','RNA splicing','RNA splicing, via transesterification reactions with bulged adenosine as nucleophile','mRNA splicing, via spliceosome')
cho_enrich_pathways=pd_enrichgo_res2@compareClusterResult$Description[grep('RNA splicing',pd_enrichgo_res2@compareClusterResult$Description)]
splicing_genes=unique(unlist(strsplit(pd_enrichgo_res@compareClusterResult[pd_enrichgo_res@compareClusterResult$Description %in% cho_enrich_pathways,'geneID'],split = '/')))
top_splicing_genes=pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$gene %in% splicing_genes & pd_tissue_NRBC_DE_res$pct.1 >0.1, ] %>% top_n(wt =avg_log2FC,n = 30 )
top_splicing_genes=top_splicing_genes[order(top_splicing_genes$pct.1,decreasing = T),]
top_splicing_genes=top_splicing_genes[1:10,]
VlnPlot(filt_NBRC_altas_seu,features =top_splicing_genes$gene,group.by = 'source_celltype',stack = T ) +NoLegend()

#核心剪接机器的“特种部队”：确保精准与效率：LSM10, SCNM1, FAM50A,TSSC4 
# RNA质量控制与“清洁工”：维护转录组完整性：ZCCHC8, RBM38, ZMAT5
#3. 转录-剪接耦合：协调生产流水线：TAF10， 它是通用转录因子 TFIID 的亚基，负责启动转录。
#生物学意义：它的出现证明了“转录与剪接是耦合的”。在 Primitive nBRC 中，为了追求速度，转录机器（TAF10）和剪接机器（SF3B1等）必须紧密配合。
# TAF10 的高表达确保了转录起始的高效，同时可能通过招募剪接因子，实现“边转录、边剪接”。
#4. 代谢支持：为剪接提供燃料
# 基因：SLC38A2（末期）, TRPT1
# 含有U12型内含子的基因虽少（人类仅约700-800个），但其功能高度集中，主要编码与DNA复制、RNA转录与剪接、蛋白质翻译、细胞周期调控及MAPK信号转导等基础生命活动密切相关的蛋白质
# 次要剪接体：SCNM1,ZMAT5:U12型内含子

order_top_splicing_genes=c('LSM10', 'SCNM1', 'FAM50A','TSSC4','ZCCHC8', 'RBM38', 'ZMAT5','TAF10','SLC38A2', 'TRPT1')
VlnPlot(filt_NBRC_altas_seu,features =order_top_splicing_genes,group.by = 'final_celltype',stack = T,split.by = 'NRBC_type',cols = cols ) 


#基础剪接	SF3B1, U2AF1, SRSF2	剪接体核心组分，决定剪接效率
#调控因子	RBM25, RAC1, MYC	通过蛋白-蛋白/蛋白-RNA互作调控剪接
#机制：RAC1通过肌动蛋白动态重组促进剪接因子（如SRSF2）在核内的快速转运，形成剪接-细胞骨架协同系统，显著提升剪接效率（Nature Cell Biology, 2024）。
#SRSF2：调控血红蛋白基因（HBB, HBA）的外显子包含，使mRNA成熟速度快2.5倍
#RBM25	高 (log2FC=2.1)	低 (log2FC=-1.0)	剪接-迁移协同因子：连接RNA剪接与细胞运动
# MYC转录因子，驱动骨髓中造血干细胞增殖
#特殊机制	AKAP17A, CAAP1, RTCB	SOS剪接，应对转座子威胁,暂不考虑
# 除了MYC，其他基因在YS NRBC 中明显表达更高
known_RNAsplicing_key_genes=c('SF3B1', 'U2AF1', 'SRSF2','RBM25', 'RAC1') # 除了MYC和U2AF1， U2AF1可能是测序导致
known_RNAsplicing_key_genes[known_RNAsplicing_key_genes %in%splicing_genes ] #"SF3B1" "SRSF2" "RBM25"
VlnPlot(filt_NBRC_altas_seu,features =c(known_RNAsplicing_key_genes,order_top_splicing_genes),group.by = 'source_celltype',stack = T ,split.by = 'tissue_stage') +NoLegend()
VlnPlot(filt_NBRC_altas_seu,features =c(order_top_splicing_genes,known_RNAsplicing_key_genes),group.by = 'NRBC_type',stack = T,split.by = 'final_celltype',cols = cols ) 

# 1.1 主要剪接体核心基因
major_spliceosome_genes <- c(
  # U1 snRNP
  "SNRNP70", "SNRPA", "SNRPC",
  # U2 snRNP
  "SF3A1", "SF3A2", "SF3A3", "SF3B1", "SF3B2", "SF3B3", "SF3B4", "SF3B5",
  # U4/U6.U5 tri-snRNP
  "PRPF3", "PRPF4", "PRPF6", "PRPF8", "PRPF31", "SNRNP200", "EFTUD2", "SART1",
  # 催化复合物
  "CDC40", "CWC15", "CWC22", "CWC27", "BCAS2",
  # EJC复合物
  "RBM8A", "MAGOH", "EIF4A3",
  # SR蛋白家族
  "SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6", "SRSF7", "SRSF9", "SRSF10",
  # hnRNP家族
  "HNRNPA1", "HNRNPK", "HNRNPL", "HNRNPU"
)

# 1.2 次要剪接体核心基因
minor_spliceosome_genes <- c(
  # U11/U12 snRNP特有蛋白
  "ZCRB1", "SNRNP35", "SNRNP48", "SNRNP25", "RNPC3", "ZMAT5", "ZRSR1", "ZRSR2",
  # U4atac/U6atac特有蛋白
  "CENATAC", "TXNL4B",  # DIM2 = TXNL4B
  # 次要剪接体特异结构组分
  "SCNM1", "CRIPT", "RBM48", "ARMC7", "PPIL2"
)

# 1.3 确保基因存在于Seurat对象中（避免报错）
major_genes_avail <- major_spliceosome_genes[major_spliceosome_genes %in% rownames(filt_NBRC_altas_seu)]
minor_genes_avail <- minor_spliceosome_genes[minor_spliceosome_genes %in% rownames(filt_NBRC_altas_seu)]

an_df=data.frame(row.names = c(major_genes_avail,minor_genes_avail),type=rep(c('major_spliceosome_genes','minor_spliceosome_genes'),c(length(major_genes_avail),length(minor_genes_avail))))
temp_df=AverageExpression(filt_NBRC_altas_seu,features =c(major_genes_avail,minor_genes_avail))$RNA
pheatmap(temp_df[rownames(an_df),],cluster_cols = F,cluster_rows = F,scale = 'row',annotation_row =an_df ,color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))


#autophagy_pathways=c( "macroautophagy","regulation of autophagy","regulation of macroautophagy","vacuole organization")#,"proteasome-mediated ubiquitin-dependent protein catabolic process"
cho_enrich_pathways=pd_enrichgo_res2@compareClusterResult$Description[grep('autophagy',pd_enrichgo_res2@compareClusterResult$Description)]
autophagy_genes=unique(unlist(strsplit(pd_enrichgo_res@compareClusterResult[pd_enrichgo_res@compareClusterResult$Description %in% cho_enrich_pathways,'geneID'],split = '/')))
top_autophagy_genes=pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$gene %in% autophagy_genes & pd_tissue_NRBC_DE_res$pct.1 >0.1, ] %>% top_n(wt =avg_log2FC,n = 30 )
top_autophagy_genes=top_autophagy_genes[order(top_autophagy_genes$pct.1,decreasing = T),]
top_autophagy_genes=top_autophagy_genes[1:10,]
VlnPlot(filt_NBRC_altas_seu,features =top_autophagy_genes$gene,group.by = 'source_celltype',stack = T ) +NoLegend()

# macroautophage 自噬相关基因 
known_autophage_genes=c('MAP1LC3B','BNIP3L','SQSTM1','OPTN') # "ULK1"     "PIK3C3" 几乎不表达，'FUNDC1','LAMP1'低表达，无明显差别
top_autophagy_genes$gene[top_autophagy_genes$gene %in% known_autophage_genes ]# SQSTM1,OPTN
p=VlnPlot(filt_NBRC_altas_seu ,features= known_autophage_genes,stack = T,group.by = 'source_celltype')+NoLegend();p


# 囊泡运输
cho_enrich_pathways=pd_enrichgo_res2@compareClusterResult$Description[grep('vesicle|vecuolar|endosomal|cytosolic',pd_enrichgo_res2@compareClusterResult$Description)]
cho_enrich_pathways=cho_enrich_pathways[-grep(pattern ='synaptic|cytosolic transport' ,cho_enrich_pathways)] # 移除不相关通路
vesicle_transport_genes=unique(unlist(strsplit(pd_enrichgo_res@compareClusterResult[pd_enrichgo_res@compareClusterResult$Description %in% cho_enrich_pathways,'geneID'],split = '/')))
vesicle_transport_genes=vesicle_transport_genes[vesicle_transport_genes %in% pd_tissue_NRBC_DE_res$gene[pd_tissue_NRBC_DE_res$avg_log2FC>1  ]];length(vesicle_transport_genes)
vesicle_transport_genes=vesicle_transport_genes[!vesicle_transport_genes %in% top_autophagy_genes$gene ]
top_vesicle_transport_genes=pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$gene %in% vesicle_transport_genes & pd_tissue_NRBC_DE_res$pct.1 >0.1, ] %>% top_n(wt =avg_log2FC,n = 30 )
top_vesicle_transport_genes=top_vesicle_transport_genes[order(top_vesicle_transport_genes$pct.1,decreasing = T),]
top_vesicle_transport_genes=top_vesicle_transport_genes[1:10,]
# 去除其中自噬中获得的基因
#top_vesicle_transport_genes=top_vesicle_transport_genes[top_vesicle_transport_genes$gene %in% top_autophagy_genes$gene, ]
VlnPlot(filt_NBRC_altas_seu ,features= top_vesicle_transport_genes$gene,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()

# 货物识别与囊泡出芽	ARF家族 (如 ARF1, ARF6)	启动囊泡形成，招募被膜蛋白。ARF6 直接调控细胞膜重塑和迁移。
#适配蛋白	AP2复合体 (AP2M1) 等	连接被膜蛋白与特定货物受体。	可能调控粘附分子（如ITGB1）的内吞与循环，影响迁移。
# RAB家族 GTP酶	囊泡的“导航系统”和“身份标识”。不同RAB定位不同细胞区室：
# RAB5: 早期内体，RAB5A
# RAB7: 晚期内体→溶酶体，RAB7A
# RAB11: 回收内体→质膜（与迁移前沿关系密切）几乎不表达
# RAB27: 分泌颗粒/溶酶体胞吐,RAB27A 几乎不表达
#SNARE蛋白家族	囊泡与靶膜的“分子拉链”： v-SNARE (囊泡膜): 如 VAMP2, VAMP7, 差异不明显
#VAMP3和STX4是SNARE蛋白家族的两个关键成员，它们常常作为一对“搭档”，在囊泡运输的最后一步——囊泡与靶膜融合中发挥核心作用
known_vesicle_transport=c('ARF1', 'ARF6','AP2M1','RAB5A','RAB7A','VAMP3','STX4');vesicle_transport_genes[vesicle_transport_genes %in% top_vesicle_transport_genes$gene ]#  No
VlnPlot(filt_NBRC_altas_seu ,features= vecle_transport,stack = T,group.by = 'source_celltype')+NoLegend()


#------------------------------ 排除nRBC亚类细胞比例影响 ，统一细胞亚类比例-------------------------------------#

if(T){
  # fetal NRBC: YS、FL、FBM，三个阶段各抽取1:1:1的细胞，构成early、mid、late NRBC， 与整体不抽样分析，结果差异很小
  
  YS_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='YS',downsample =3000)
  YS_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='YS',downsample =3000)
  YS_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='YS',downsample =3000)
  YS_filt_NBRC_altas_seu=merge(YS_early_Ery_filt_NBRC_altas_seu,c(YS_mid_Ery_filt_NBRC_altas_seu,YS_late_Ery_filt_NBRC_altas_seu))
  rm(YS_early_Ery_filt_NBRC_altas_seu,YS_mid_Ery_filt_NBRC_altas_seu,YS_late_Ery_filt_NBRC_altas_seu)
  
  FL_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FL',downsample =1000)
  FL_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FL',downsample =1000)
  FL_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FL',downsample =1000)
  FL_filt_NBRC_altas_seu=merge(FL_early_Ery_filt_NBRC_altas_seu,c(FL_mid_Ery_filt_NBRC_altas_seu,FL_late_Ery_filt_NBRC_altas_seu))
  rm(FL_early_Ery_filt_NBRC_altas_seu,FL_mid_Ery_filt_NBRC_altas_seu,FL_late_Ery_filt_NBRC_altas_seu)
  
  FBM_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FBM',downsample =1000)
  FBM_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FBM',downsample =1000)
  FBM_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FBM',downsample =1000)
  FBM_filt_NBRC_altas_seu=merge(FBM_early_Ery_filt_NBRC_altas_seu,c(FBM_mid_Ery_filt_NBRC_altas_seu,FBM_late_Ery_filt_NBRC_altas_seu))
  rm(FBM_early_Ery_filt_NBRC_altas_seu,FBM_mid_Ery_filt_NBRC_altas_seu,FBM_late_Ery_filt_NBRC_altas_seu)
  
  ABM_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='ABM',downsample =1000)
  ABM_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='ABM',downsample =1000)
  ABM_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='ABM',downsample =1000)
  ABM_filt_NBRC_altas_seu=merge(ABM_early_Ery_filt_NBRC_altas_seu,c(ABM_mid_Ery_filt_NBRC_altas_seu,ABM_late_Ery_filt_NBRC_altas_seu))
  rm(ABM_early_Ery_filt_NBRC_altas_seu,ABM_mid_Ery_filt_NBRC_altas_seu,ABM_late_Ery_filt_NBRC_altas_seu)
  
  pd_subset_filt_NBRC_altas_seu=merge(ABM_filt_NBRC_altas_seu,c(FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu,YS_filt_NBRC_altas_seu))
  rm(ABM_filt_NBRC_altas_seu,FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu,YS_filt_NBRC_altas_seu)
  gc()
  
  pd_subset_filt_NBRC_altas_seu <- JoinLayers(pd_subset_filt_NBRC_altas_seu)
  pd_subset_filt_NBRC_altas_seu$type_stage='definitive'
  pd_subset_filt_NBRC_altas_seu$type_stage[pd_subset_filt_NBRC_altas_seu$tissue_stage=='YS']='primitive'
  Idents(pd_subset_filt_NBRC_altas_seu)='type_stage'
  saveRDS(pd_subset_filt_NBRC_altas_seu,file = 'Protein_NRBC_marker/res_data/temp_pd_subset_filt_NBRC_altas_seu.rds')
  
}else{
  pd_subset_filt_NBRC_altas_seu=readRDS('Protein_NRBC_marker/res_data/temp_pd_subset_filt_NBRC_altas_seu.rds')
  pd_subset_filt_NBRC_altas_seu$source_celltype=factor(pd_subset_filt_NBRC_altas_seu$source_celltype,levels = levels(filt_NBRC_altas_seu$source_celltype))
}

Idents(pd_subset_filt_NBRC_altas_seu)='type_stage'
pd_whole_level_markers=FindAllMarkers(pd_subset_filt_NBRC_altas_seu)
saveRDS(pd_whole_level_markers,file = 'Protein_NRBC_marker/res_data/pd_whole_level_markers.rds')

top_pd_whole_level_markers=pd_whole_level_markers[pd_whole_level_markers$avg_log2FC >1 & pd_whole_level_markers$pct.2<0.2,] %>% group_by(cluster) %>% do(head(.,10))
p=DotPlot(filt_NBRC_altas_seu,features = c(top_pd_whole_level_markers$gene[11:20],top_pd_whole_level_markers$gene[1:10]),cols = c('gray','firebrick3'),scale = F,group.by = 'source_celltype')+RotatedAxis() ;p
ggsave(p,filename = 'Protein_NRBC_marker/res_pic/main_figure2/pd_top10_marker_dotplot.pdf',width = 8,height = 10)

genelist=pd_whole_level_markers[pd_whole_level_markers$cluster=='primitive','avg_log2FC'];names(genelist)=pd_whole_level_markers[pd_whole_level_markers$cluster=='primitive','gene']
genelist=sort(genelist,decreasing = T)
pd_whole_level_marker_gseGO_res=data.frame(gseGO(geneList =genelist,ont = 'BP',OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',pvalueCutoff = 0.05 ))
pd_whole_level_marker_gseGO_res$res='up';pd_whole_level_marker_gseGO_res$res[pd_whole_level_marker_gseGO_res$NES<0]='down'
table(pd_whole_level_marker_gseGO_res$res)#down203,up:43 
saveRDS(pd_whole_level_marker_gseGO_res,file = 'Protein_NRBC_marker/res_data/main_figure1/pd_whole_level_marker_gseGO_res.rds')


top_pd_whole_level_marker_gseGO_res=pd_whole_level_marker_gseGO_res %>% group_by( res)  %>% do(head(.,20))
top_pd_whole_level_marker_gseGO_res=top_pd_whole_level_marker_gseGO_res[order(top_pd_whole_level_marker_gseGO_res$res,top_pd_whole_level_marker_gseGO_res$NES,decreasing = T),]
top_pd_whole_level_marker_gseGO_res$Description=factor(top_pd_whole_level_marker_gseGO_res$Description,levels = unique(top_pd_whole_level_marker_gseGO_res$Description))
p1=ggplot(top_pd_whole_level_marker_gseGO_res,aes(x =res ,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'gray',high = 'firebrick3')+
  theme(text = element_text(face = 'bold'))+ggtitle("primitive vs definitive  in whole level")+scale_y_discrete(labels=function(x)str_wrap(x,width = 50))
p1

ggsave(p1,filename = 'Protein_NRBC_marker/res_pic/main_figure2/pd_whole_deg_gsego_dotplot.pdf',width = 6,height = 8)

#--------与原始比较，存在差异，但是P；巨噬，RNA剪接都有，definitive：基本都是免疫与细胞粘附，都存在--------------#
if(F){
  table(pd_whole_level_markers[pd_whole_level_markers$avg_log2FC >0,'cluster'])# definitive / primitive:1427/8985 
  primitive_Ery_pos_markers=pd_whole_level_markers[pd_whole_level_markers$avg_log2FC >0 & pd_whole_level_markers$cluster=='primitive','gene']
  definitive_Ery_pos_markers=pd_whole_level_markers[pd_whole_level_markers$avg_log2FC >0 & pd_whole_level_markers$cluster=='definitive','gene']
  pd_enrichgo_res2_new=compareCluster(geneClusters =list('primitive'=unique(primitive_Ery_pos_markers),
                                                         'definitive'=unique(definitive_Ery_pos_markers)),
                                      keyType = 'SYMBOL',ont = "BP",fun = 'enrichGO',  OrgDb='org.Hs.eg.db' )
  
  dotplot(pd_enrichgo_res2_new,showCategory=10) # "proteasome-mediated ubiquitin-dependent protein catabolic process" ,"macroautophagy"  ,也包含子在其中，只是没有那么明显
  
  pd_enrichgo_res2_new=pairwise_termsim(pd_enrichgo_res2_new)
  emapplot(pd_enrichgo_res2_new, node_label = "category",showCategory = 30) # 结果相似
  
}


group='primitive_definitive'
sfile='Protein_NRBC_marker/DE_marker/primitive_definitive_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = pd_subset_filt_NBRC_altas_seu,group = group,sfile = sfile)
sub_pd_all_Ery_tissue_markers=res[[1]]
sub_pd_count_df=res[[2]]
rm(res);gc()

length(unique(sub_pd_all_Ery_tissue_markers$gene) )
table(sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >0,c('celltype','cluster')])
length(unique(sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >0 & sub_pd_all_Ery_tissue_markers$cluster=='primitive' ,'gene']))
length(unique(sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >0 & sub_pd_all_Ery_tissue_markers$cluster=='definitive' ,'gene']))

table(unique(sub_pd_all_Ery_tissue_markers$gene) %in% unique(pd_whole_level_markers$gene)) 
#FALSE  TRUE 
#2577 10225
table(unique(pd_whole_level_markers$gene) %in% unique(sub_pd_all_Ery_tissue_markers$gene)) 
#FALSE  TRUE 
#187 10225
# 查看在亚群中鉴定到的新marker中top marker 表达情况，效果比较差
new_sub_pd_all_Ery_tissue_markers=sub_pd_all_Ery_tissue_markers[!sub_pd_all_Ery_tissue_markers$gene %in% unique(pd_whole_level_markers$gene) & sub_pd_all_Ery_tissue_markers$avg_log2FC >0 ,]
new_sub_pd_all_Ery_tissue_markers=new_sub_pd_all_Ery_tissue_markers[new_sub_pd_all_Ery_tissue_markers$avg_log2FC >1 & new_sub_pd_all_Ery_tissue_markers$pct.1 >0.1 & new_sub_pd_all_Ery_tissue_markers$pct.2 <0.3,]
length(unique(new_sub_pd_all_Ery_tissue_markers$gene))# 435 
top5_new_sub_pd_all_Ery_tissue_markers=new_sub_pd_all_Ery_tissue_markers %>% group_by(celltype) %>%top_n(wt = avg_log2FC,n = 5) 
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =top5_new_sub_pd_all_Ery_tissue_markers$gene,stack = T)+NoLegend()

#查看在整体群中鉴定到的marker，而在亚群中未鉴定到，中top marker 表达情况，效果比较好,ABM nBRC pct低，而 表达大水平接近
out_pd_whole_level_markers=pd_whole_level_markers[!pd_whole_level_markers$gene %in% unique(sub_pd_all_Ery_tissue_markers$gene),]
table(out_pd_whole_level_markers$avg_log2FC >1 &  out_pd_whole_level_markers$pct.2 <0.3 & out_pd_whole_level_markers$pct.1 >0.1)
#FALSE  TRUE 
#179     8 ,几乎不存在特异性，除了KRCC1， 候选筛选还是应该考虑sub_pd_all_Ery_tissue_markers
out_pd_whole_level_markers=out_pd_whole_level_markers[out_pd_whole_level_markers$avg_log2FC >1 &  out_pd_whole_level_markers$pct.2 <0.3 & out_pd_whole_level_markers$pct.1 >0.1,]
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =out_pd_whole_level_markers$gene,stack = T,cols = cols)
DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =out_pd_whole_level_markers$gene,scale = T)+RotatedAxis()

sub_pd_all_Ery_tissue_markers_gseGO_res=subcelltype_gseGO_func(RNA_markers = sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$cluster=='primitive',])
sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$res='up';sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$res[sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$NES <0]='down'
saveRDS(sub_pd_all_Ery_tissue_markers_gseGO_res,file = 'Protein_NRBC_marker/res_data/main_figure2/sub_pd_all_Ery_tissue_markers_gseGO_res.rds')

top_sub_pd_all_Ery_tissue_markers_gseGO_res=sub_pd_all_Ery_tissue_markers_gseGO_res[[2]][sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$ONTOLOGY=='BP',] %>% group_by(celltype,res) %>% do(head(.,10))
top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype=factor(top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))
top_sub_pd_all_Ery_tissue_markers_gseGO_res=top_sub_pd_all_Ery_tissue_markers_gseGO_res[order(top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype,top_sub_pd_all_Ery_tissue_markers_gseGO_res$res,top_sub_pd_all_Ery_tissue_markers_gseGO_res$NES),]
top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description=factor(top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description,levels = unique(top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description))
p2=ggplot(top_sub_pd_all_Ery_tissue_markers_gseGO_res,aes(x =celltype ,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+
  scale_color_gradient2(low = 'blue',mid = 'gray',high = 'firebrick3')+theme(text = element_text(face = 'bold'))+ggtitle(' primitive vs definitive in substage')
  #scale_y_discrete(label=function(x)str_wrap(string = x,width = 50))
  
p2
ggsave(p2,filename = 'Protein_NRBC_marker/res_pic/main_figure2/pd_erysubstage_gsetop10_dotplot.pdf',width = 8,height = 8)


# 利用emapplot 查看top enriched pathway 之间的关系
test=pairwise_termsim(sub_pd_all_Ery_tissue_markers_gseGO_res[[1]]$early_Ery)
emapplot(test,color='NES',showCategory=30)+ggtitle('pd_early_Ery_top30')
# 细胞迁移通路与血管的管发育通路共享部分驱动共享基因
p_core_genes=top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype=='early_Ery' & top_sub_pd_all_Ery_tissue_markers_gseGO_res$NES >0,]$core_enrichment
names(p_core_genes)=as.character(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype=='early_Ery' & top_sub_pd_all_Ery_tissue_markers_gseGO_res$NES >0,]$Description)
p_core_genes_list=sapply(p_core_genes,function(x){strsplit(x,split = '/')})
dev.off()
p=venn.diagram(p_core_genes_list[c('positive regulation of cell motility','tissue morphogenesis','blood vessel development')],
               filename = NULL,fill=cols[1:3], scaled = T,main='test',force.unique = T,main.cex = 2,sub.cex = 2,total.population = T)
grid.draw(p)

test=pairwise_termsim(sub_pd_all_Ery_tissue_markers_gseGO_res[[1]]$mid_Ery)
emapplot(test,color='NES',showCategory=30)+ggtitle('pd_mid_Ery_top30')
test=pairwise_termsim(sub_pd_all_Ery_tissue_markers_gseGO_res[[1]]$late_Ery)
emapplot(test,color='NES')+ggtitle('pd_late_Ery_top30')


# 比较enrichGO 得到的结果
pd_all_Ery_tissue_markers_enrichGO_list=subcelltype_enrichGO_func(RNA_markers = sub_pd_all_Ery_tissue_markers[ sub_pd_all_Ery_tissue_markers$avg_log2FC>0,] )#no term enriched under specific pvalueCutoff.

pd_all_Ery_tissue_markers_enrichGO_df=data.frame()
for (celltype in c('early_Ery','mid_Ery','late_Ery')) {
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['primitive']]@result$celltype=celltype
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['definitive']]@result$celltype=celltype
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['primitive']]@result$type='primitive'
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['definitive']]@result$type='definitive'
  temp_df=rbind(pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['primitive']]@result,pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['definitive']]@result)
  pd_all_Ery_tissue_markers_enrichGO_df=rbind(pd_all_Ery_tissue_markers_enrichGO_df,temp_df)
}

pd_all_Ery_tissue_markers_enrichGO_df$ratio=as.numeric(t(data.frame(strsplit(pd_all_Ery_tissue_markers_enrichGO_df$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit( pd_all_Ery_tissue_markers_enrichGO_df$GeneRatio,split = '/')))[,2])
pd_all_Ery_tissue_markers_enrichGO_df$celltype=factor(pd_all_Ery_tissue_markers_enrichGO_df$celltype,levels =c('early_Ery','mid_Ery','late_Ery') )
saveRDS(pd_all_Ery_tissue_markers_enrichGO_df,file = 'Protein_NRBC_marker/res_data/main_figure2/pd_all_Ery_tissue_markers_enrichGO_df.rds')

top_pd_all_Ery_tissue_markers_enrichGO_df=pd_all_Ery_tissue_markers_enrichGO_df %>%group_by(type,celltype) %>% do(head(.,10))
top_pd_all_Ery_tissue_markers_enrichGO_df=top_pd_all_Ery_tissue_markers_enrichGO_df[order(top_pd_all_Ery_tissue_markers_enrichGO_df$type,top_pd_all_Ery_tissue_markers_enrichGO_df$celltype),]
top_pd_all_Ery_tissue_markers_enrichGO_df$Description =factor(top_pd_all_Ery_tissue_markers_enrichGO_df$Description,levels = unique(top_pd_all_Ery_tissue_markers_enrichGO_df$Description))
p=ggplot(top_pd_all_Ery_tissue_markers_enrichGO_df,aes(x=type,y=Description,color=-log10(p.adjust),size=ratio))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid='gray',high = 'firebrick3')+facet_grid(~celltype)+
  theme( axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'enrichGO')+scale_y_discrete(label=function(x)str_wrap(x,width = 50))
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure2/pd_erystage_enrichtop10_dotplot.pdf',width = 8,height = 8)

# 发现YS NRBC 显著高表达MTq1家族基因，抗氧化酶基因SOD3和FAM213A（PRXL2A），应激标记基因GDF15 ,ABAC1：囊泡运输的关键调控基因， 可能参与调控 自噬体形成 或成熟过程，BEX1基因:神经中研究比较多，
#基因富集分析也显著富集 自噬、囊泡运输和金属离子应答通路通路，且细胞复制多个通路在末期下调。

 
# ---------------细胞迁移------------------#还有细胞骨架组织等通路可以补,actin
cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep('motility|migration',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
cell_move_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
cell_move_genes=cell_move_genes[cell_move_genes %in% pd_tissue_NRBC_DE_res$gene[pd_tissue_NRBC_DE_res$avg_log2FC>1  ]];length(cell_move_genes)
top_cell_move_genes=pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$gene %in% cell_move_genes & pd_tissue_NRBC_DE_res$pct.1 >0.1, ] %>% top_n(wt =avg_log2FC,n = 30 )
top_cell_move_genes=top_cell_move_genes[order(top_cell_move_genes$pct.1,decreasing = T),]
top_cell_move_genes=top_cell_move_genes[1:10,]
VlnPlot(filt_NBRC_altas_seu ,features= top_cell_move_genes$gene,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()

known_cell_move_genes=c('VCL','PXN','ITGB1','MYH9','MYH10','MYL6','MYL12A','MYL12B', 'RAC1', 'CDC42', 'RHOA','DIAPH1','DIAPH3', 'ACTR2', 'ACTR3',
                  'ARPC2','WASF2','PTEN')
known_cell_move_genes[known_cell_move_genes %in% top_cell_move_genes$gene] # no
known_cell_move_genes[known_cell_move_genes %in% cell_move_genes] # no 

temp_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$gene %in% cell_move_genes & pd_whole_level_markers_degs_df$avg_log2FC >0, ]
temp_df$gene=factor(temp_df$gene,cell_move_genes)
p1=ggplot(temp_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)

p2=VlnPlot(object = filt_NBRC_altas_seu,features =cell_move_genes ,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols,stack = T,raster=FALSE)
p1+p2+plot_layout(ncol = 1,heights = c(0.6,1.2))


#----------------blood_development_genes------------#
# primitive 显著富集到多个血管发育相关通路，可能原因（血管母细胞分子印记）（2）参与调控血管发生

cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep('blood|vasculature|tube',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
blood_development_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
blood_development_genes=blood_development_genes[blood_development_genes %in% pd_tissue_NRBC_DE_res$gene[pd_tissue_NRBC_DE_res$avg_log2FC>1  ]];length(blood_development_genes)
top_blood_development_genes=pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$gene %in% blood_development_genes & pd_tissue_NRBC_DE_res$pct.1 >0.1, ] %>% top_n(wt =avg_log2FC,n = 30 )
top_blood_development_genes=top_blood_development_genes[order(top_blood_development_genes$pct.1,decreasing = T),]
top_blood_development_genes=top_blood_development_genes[1:10,]
VlnPlot(filt_NBRC_altas_seu ,features= top_blood_development_genes$gene,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()

table(blood_development_genes %in%cell_move_genes ) # 34/36，有一半基因与细胞迁移有关，可以关注

# HTATIP2,RORA,HDAC5,BTG1,GADD45A,ID1,SERPINF1,RORA,HMOX1,ZFAND5,DDIT3,EPN1,MKKS
# BTG1：细胞周期阻滞与分化促进，GADD45A	应激反应与基因组维稳	作为 “压力传感器” ，响应DNA损伤/氧化应激，通过阻滞细胞周期、启动DNA修复，为YS-nRBC在快速增殖中维持基因组完整性。
#HTATIP2（又名TIP30/CC3）是一种具有氧化还原酶活性和转录共激活功能的多功能蛋白
#PDCD10：程序性细胞死亡10，参与调节细胞凋亡、促进细胞增殖、调控细胞迁移，并维持高尔基体正常结构

# TAL1和LYL1是著名的同时调控早期造血和血管发育的关键因子，它们在血岛的内皮-造血转化中发挥核心作用
#TAL1	核心生血转录因子。在内皮细胞中表达，驱动其向造血命运转化，是形成血岛造血中心的绝对关键因子。
#LYL1	与TAL1功能部分冗余的基本螺旋-环-螺旋转录因子，共同调控早期造血干/祖细胞的产生。
#LMO4	TAL1转录复合物的关键衔接子/辅因子。通过组装复合物精确调控造血特异性基因的表达，是生血过程的标志。
#CITED2	转录辅激活因子，对心脏、神经嵴和血管发育至关重要。缺失导致严重血管缺陷，是早期血管母细胞功能的标志。
known_endo_derived_gene_marker=c('TAL1','LYL1','LMO4','CITED2')
VlnPlot(filt_NBRC_altas_seu,features =known_endo_derived_gene_marker,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols,stack = T )


#--------------------------DNA 复制-------------------#
#BTG2 (PC3/TIS21), BTG3 (ANA), BTG4 (PC3B)	与BTG1高度同源，功能冗余。作为转录辅调节因子，通过抑制细胞周期蛋白（如Cyclin D1）转录或与CCR4-NOT复合物互作，调控mRNA降解，从而阻滞细胞周期（常在G0/G1期）并促进分化
#p21 (CDKN1A), p27 (CDKN1B), p57 (CDKN1C)	直接结合并抑制多种Cyclin-CDK复合物，导致G1期阻滞，是多种信号通路（如p53）诱导分化、实现细胞
cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep('mitotic|spindle|chromatid|chromosome',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
DNA_replication_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
length(DNA_replication_genes)
DNA_replication_genes=DNA_replication_genes[DNA_replication_genes %in% sub_pd_all_Ery_tissue_markers$gene[sub_pd_all_Ery_tissue_markers$avg_log2FC>1 & sub_pd_all_Ery_tissue_markers$celltype=='late_Ery'  ]];length(DNA_replication_genes)
top_DNA_replication_genes=pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$gene %in% DNA_replication_genes & pd_tissue_NRBC_DE_res$pct.1 >0.1, ] %>% top_n(wt =avg_log2FC,n = 30 )
top_DNA_replication_genes=top_DNA_replication_genes[order(top_DNA_replication_genes$pct.1,decreasing = T),]
top_DNA_replication_genes=top_DNA_replication_genes[1:10,]
VlnPlot(filt_NBRC_altas_seu ,features= top_DNA_replication_genes$gene,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()


known_DNA_replication=c('CDT1','CDC6','PCNA','FEN1','CDK1','NDC80','MCM2','MCM5') # 细胞周期 
known_prodifferention=c('GATA1','KLF1','BTG1','BTG2','BTG3','CDKN1A','CDKN1B')
VlnPlot(object = filt_NBRC_altas_seu,features =c(known_DNA_replication,'known_prodifferention') ,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols,stack = T,raster=FALSE)

# 血红蛋白血红素合成关键基因
#GATA1是红系分化的主调控因子，它直接开启了整个红系分化程序。 在figure1
#KLF1是GATA1的关键下游靶基因，被称为“红系分化的效应器”。  figure1
#它直接激活您提到的ALAS2（血红素合成限速酶）、SLC25A37（线粒体铁转运蛋白） 等一系列功能基因
#因此，KLF1的更高表达，直接解释了为何YS-NRBC的铁利用和血红素合成基因会“更高更早地表达”。

# TFRC 转铁蛋白，   figure1
# NFE2L2：编码NRF2，是抗氧化反应的关键调节因子。在铁过载或氧化压力下被激活，可上调 FTH1、FTL等基因，以应对铁诱导的氧化损伤。在NRBC都存在
# 还有核心靶基因：NQO1、GCLM，NQO1主要防御外源性醌类物质，GCLM驱动的谷胱甘肽系统是对抗广泛氧化损伤的“通用货币”。

#SLC40A1：编码铁转运蛋白，是唯一已知的细胞铁输出通道。
#ALAS2：编码δ-氨基酮戊酸合酶2，是红细胞特异性、血红素合成途径的第一个限速酶。其表达受铁供应水平调控。
#SLC25A37 / SLC25A38：分别编码线粒体铁转运蛋白Mitoferrin-1和-2，负责将细胞质中的铁转运至线粒体基质，是血红素合成的关键步骤。
# FECH：编码亚铁螯合酶，催化铁插入原卟啉IX形成血红素的最后一步。该基因突变可导致红细胞原卟啉病。
#HMOX1 (血红素加氧酶1) 是细胞应对铁和氧化压力的核心枢纽，其作用可概括为“解毒者”、“铁释放者”和“压力响应中心”
#HMOX1是一种诱导酶，其主要功能是催化有毒的“游离血红素”分解，生成胆绿素（随后转化为胆红素）、一氧化碳（CO）和亚铁离子（Fe²⁺）。这一过程是铁代谢循环的关键环节：

# 1.1核心合成酶（按通路顺序）
heme_synthesis_genes <- c(
  "ALAS2",           # 红系特异性ALA合酶
  # 胞质反应酶
  "ALAD",            # ALA脱水酶
  "HMBS",            # PBG脱氨酶
  "UROS",            # 尿卟啉原III合酶
  "UROD",            # 尿卟啉原脱羧酶
  # 线粒体反应酶
  "CPOX",            # 粪卟啉原氧化酶
  "PPOX",            # 原卟啉原氧化酶
  "FECH"             # 亚铁螯合酶
)


# 1.2 血红素解毒模块 (Heme Degradation/Detoxification)
heme_detox_genes <- c(
  "HMOX1",    # 血红素加氧酶1（核心解毒酶）
  "BLVRB",    # 胆绿素还原酶B
  "BLVRA"     # 胆绿素还原酶A
)

# 1.3 血红素外排/转运模块 (Heme Export/Transport)
heme_export_genes <- c(
  "FLVCR1",   # 血红素外排泵（细胞膜），表达低，无明显差异，可以考虑排除
  "SLC48A1",  # HRG1，溶酶体血红素转运蛋白
  "ABCB6",    # ABC转运蛋白，参与原卟啉IX转运，表达低，无明显差异，可以考虑排除
  "ABCB10"    # 线粒体ABC转运蛋白，参与血红素合成调控
)

# 1.4 铁代谢与储存模块 (Iron Metabolism & Storage)
iron_genes <- c(
  "FTH1",     # 铁蛋白重链
  "FTL",      # 铁蛋白轻链
  "SLC25A37", # mitoferrin-1，线粒体铁转运
  "SLC25A38"  # 甘氨酸转运，为ALAS提供底物
)

known_hema_function_genes <- c(heme_synthesis_genes, heme_detox_genes,  heme_export_genes,  iron_genes)
VlnPlot(filt_NBRC_altas_seu,features =known_hema_function_genes,group.by = 'source_celltype',stack = T,split.by = 'tissue_stage',cols = cols )+NoLegend()

known_hema_function_genes=c('NFE2L2','FTH1','FTL','SLC40A1','ALAS2','SLC25A37','SLC25A38','FECH','HMOX1')
VlnPlot(filt_NBRC_altas_seu,features =hema_function_genes,group.by = 'source_celltype',stack = T,split.by = 'tissue_stage',cols = cols )+NoLegend()

#-------------------------iron 压力-------------------#
cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep(' ion',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
iron_stress_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
length(iron_stress_genes)
VlnPlot(filt_NBRC_altas_seu ,features=iron_stress_genes,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()


iron_stress_genes=iron_stress_genes[iron_stress_genes %in% sub_pd_all_Ery_tissue_markers$gene[sub_pd_all_Ery_tissue_markers$avg_log2FC>1 & sub_pd_all_Ery_tissue_markers$celltype=='late_Ery'  ]];length(iron_stress_genes)
top_iron_stress_genes=pd_tissue_NRBC_DE_res[pd_tissue_NRBC_DE_res$gene %in% iron_stress_genes & pd_tissue_NRBC_DE_res$pct.1 >0.1, ] %>% top_n(wt =avg_log2FC,n = 30 )
top_iron_stress_genes=top_iron_stress_genes[order(top_iron_stress_genes$pct.1,decreasing = T),]
top_iron_stress_genes=top_iron_stress_genes[1:10,] # "MT1E","MT1F","MT1G","MT1H"已经在top marker中展现
top_iron_stress_genes=top_iron_stress_genes[!top_iron_stress_genes$gene %in% c("MT1E","MT1F","MT1G","MT1H"),]
top_iron_stress_genes=top_iron_stress_genes[order(top_iron_stress_genes$gene),]
VlnPlot(filt_NBRC_altas_seu ,features=top_iron_stress_genes$gene,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()



# ---------------HSPC-derived NRBC unique signature -----------------#
erythropoesis_gene=c('BCL11A','GATA2')
# ICAM4 与巨噬细胞整合素结合，CD44 与ECM 结合
cell_adension_genes=c('CD44','ICAM4')
# 核纤层蛋白LMNA、细胞骨架马达MYH9/MYH10、激酶AURKA/AURKB及成核因子DIAPH1，收缩调节激酶ROCK1
enucleation_genes=c('LMNA','AURKA','AURKB','MYH9', 'MYH10', 'ROCK1','DIAPH1')
#TAPBP	TAP结合蛋白（亦称 Tapasin）	抗原肽装载的“质量控制员”与“桥梁”	专一服务于 MHC-I类分子。
#它在内质网中物理性桥接抗原肽转运体 和新合成的MHC-I重链-β2M复合物，确保将合适长度和亲和力的抗原肽高效、精准地装载到MHC-I分子的肽结合槽中，是MHC-I成熟的关键步骤
#B2M	β2-微球蛋白	MHC-I类分子的结构基石与稳定伴侣	是 所有MHC-I类分子（如HLA-A, B, C）的不变轻链。
#它与MHC-I重链非共价结合，对其正确折叠、稳定结构和转运至细胞表面绝对必需。没有B2M，MHC-I重链无法稳定存在。它与MHC-II类分子无关。
VlnPlot(filt_NBRC_altas_seu ,features= erythropoesis_gene,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()

cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep('MHC|antigen',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
MHC_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
sort(MHC_genes)
VlnPlot(filt_NBRC_altas_seu ,features= MHC_genes,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()
MHC_genes=c("B2M","HLA-A","HLA-B","HLA-C","HLA-E","CD74","HLA-DMA","HLA-DPA1","HLA-DPB1","HLA-DRA" ,"HLA-DRB1", "HLA-DRB5","TAPBP") # order, and cho
VlnPlot(filt_NBRC_altas_seu ,features= MHC_genes,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()

peptide_load=c('TAPBP','TAP1','TAP2','PDIA3','ERAP1')
Proteasome_genes=rownames(filt_NBRC_altas_seu)[grep('PSMB',rownames(filt_NBRC_altas_seu))][-10]

VlnPlot(filt_NBRC_altas_seu ,features =c(peptide_load,regulator_genes),stack = T,group.by = 'source_celltype',cols = cols,split.by = 'tissue_stage')


cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep('monocyte|leukocyte|T cell',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
immune_remodle_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
sort(immune_remodle_genes)
immune_remodle_genes=immune_remodle_genes[!immune_remodle_genes %in% MHC_genes ]
VlnPlot(filt_NBRC_altas_seu ,features =immune_remodle_genes,stack = T,group.by = 'source_celltype',cols = cols)+NoLegend()
#转录调控,信号转导，细胞因子以及受体，细胞粘附
cho_immune_remodle_genes=c( "SOX4",'FOXP1',"ZBTB16", "SMARCA2", "PTPRC","INPP5D", "RHOH","TESPA1","NCKAP1L", 'IL1B',"LGALS3","LGALS9","CSF2RB","IL2RG","ITGA4","CD44" ,'SPN') # PNP:嘌呤代谢酶，ATAD5：	DNA损伤应答/周期调控
VlnPlot(filt_NBRC_altas_seu ,features =cho_immune_remodle_genes,stack = T,group.by = 'source_celltype',cols = cols)+NoLegend()


#SOX4	SRY-box转录因子4	免疫细胞发育的“命运调控器”。在B细胞和T细胞的早期发育、谱系决定与存活中不可或缺。它调控多种靶基因，影响淋巴细胞的增殖、分化和功能。
#TESPA1:T细胞受体信号转导的“专用适配器”。在胸腺T细胞阳性选择阶段，特异性介导TCR信号向下游（如NFAT通路）的高效传递，对功能性T细胞库的形成至关重要。
#CD74	恒定链（Ii链）	MHC-II类分子抗原呈递的“专职伴侣与导航员”。作为伴侣蛋白，防止MHC-II在内质网中过早结合抗原；作为靶向信号，引导MHC-II至内吞系统装载外源抗原肽。
#DOCK8	胞质因子DOCK8	淋巴细胞迁移、存活与功能的“支架”。一种鸟嘌呤核苷酸交换因子（GEF），通过激活CDC42等小G蛋白，调控淋巴细胞（尤其CD8⁺ T和B细胞）的骨架重组、免疫突触形成、存活和效应功能。
#INPP5D	SHIP1磷酸酶	髓系与B细胞的负向调控“刹车”。通过降解PIP3，负向调控PI3K-Akt-mTOR这一关键的促存活、增殖和代谢通路，防止免疫细胞过度活化，在免疫耐受和稳态中起关键作用。
#LEF1	淋巴样增强结合因子1	T/B淋巴细胞发育的“命运决定者”。是Wnt/β-catenin信号通路下游的关键转录因子，对早期T细胞发育、CD4⁺ T细胞谱系决定以及B细胞生成不可或缺。
#RHOH	Rho GTP酶H	T细胞发育与信号转导的“特化调控器”。一种缺乏GTP酶活性的特殊Rho蛋白，主要在T细胞中表达。在胸腺细胞发育的β选择阶段和成熟T细胞受体（TCR）信号传导中起关键调控作用。
# SOX4,TESPA1, CD74
# 显著性检验
immune_pathway_related_genes=c('SOX4','TESPA1','CD74','DOCK8','INPP5D','LEF1','RHOH') #
VlnPlot(filt_NBRC_altas_seu ,features=c(MHC_genes,immune_pathway_related_genes) ,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols,pt.size = 0,raster=FALSE,alpha = 0.6,stack = T)+NoLegend()


# 去除HMOX1 ，放在血红素中
gene_list=list('coregene_vasculature_development'=top_blood_development_genes$gene[-10],'known_keygene_vasculature'=known_endo_derived_gene_marker,
               'coregene_RNA_splicing'=top_splicing_genes$gene,'known_keygene_RNA_splicing'=known_RNAsplicing_key_genes,
               'known_hema_metabolism_genes'=known_hema_function_genes, 'coregene_iron_tress'=top_iron_stress_genes$gene, 
               'coregene_autophagy_genes'=top_autophagy_genes$gene,'known_keygene_autophage'=known_autophage_genes,
               'coregene_vesicle_transport'=top_vesicle_transport_genes$gene,'known_keygene_vesicle_transport'=known_vesicle_transport, 
               'coregene_DNA_reppication'=top_DNA_replication_genes$gene,'known_DNA_replication'=known_DNA_replication,'known_prodifferention'=known_prodifferention,
               'coregene_cell_migration'=top_cell_move_genes$gene,'known_cell_migration'=known_cell_move_genes,
               'coregene_immune_remodeling'=immune_pathway_related_genes,'coregene_MHC'=MHC_genes
               
)
saveRDS(gene_list,file = 'Protein_NRBC_marker/res_data/main_figure2/keypathways_keycore_gene_list.rds')
# 彼此存在重复,特别是kn 与core之间
all_genes=unique(unlist(gene_list))
core_genes=unique(as.character(unlist(gene_list[names(gene_list)[grep('core',names(gene_list))]])))
kn_keygenes=unique(as.character(unlist(gene_list[names(gene_list)[grep('kn',names(gene_list))]])))

sapply(names(gene_list)[grep('core',names(gene_list))], function(x){gene_list[[x]]=gene_list[[x]][ gene_list[[x]] %in% kn_keygenes ] })

Idents(filt_NBRC_altas_seu)='source_celltype'
mexp_df=as.matrix(AverageExpression(filt_NBRC_altas_seu,features =as.character(unlist(gene_list)),group.by = 'source_celltype' )$RNA)
mexp_df=mexp_df[all_genes,]
an_df=data.frame(row.names = all_genes,type=rep('top_coregenes',length(all_genes)))
an_df[kn_keygenes[!kn_keygenes %in% core_genes],'type']='known_keygenes'

p=pheatmap(t(mexp_df),annotation_col =an_df,cluster_rows = F,annotation_colors = list('knwon_keygene'=cols[1],'top_coregenes'=cols[2]), # cols[1]:
         cluster_cols = F,scale = 'column',border_color = 'white',color = colorRampPalette(colors = c('gray','white','firebrick3'))(100))
ggsave(as.ggplot(p),filename='Protein_NRBC_marker/res_pic/main_figure2/keypathay_topcoregene_exp_heatmap.pdf',width = 20,height = 5,dpi = 300)



####-------------显著富集通路活性打分计算-----------------#
library(UCell)
gene_sets <- list( macroautophagy = unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='macroautophagy'],split = '/')),
                  'regulation_of_autophagy'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='regulation of autophagy'],split = '/')),
                  'vacuole organization'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='vacuole organization'],split = '/')),
                  'proteasome-mediated_ubiquitin-dependent_protein_catabolic_process'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='proteasome-mediated ubiquitin-dependent protein catabolic process'],split = '/')),
                  'DNA replication'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='DNA replication'],split = '/')),
                  'positive_regulation_of_T_cell_activation'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='positive regulation of T cell activation'],split = '/')),
                  'regulation_of_leukocyte_proliferation'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='regulation of leukocyte proliferation'],split = '/')),
                  'positive_regulation_of_leukocyte_cell-cell adhesion'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='positive regulation of leukocyte cell-cell adhesion'],split = '/')),
                  'stress response_to_metal_ion'= unlist(strsplit(top_pd_whole_level_marker_gseGO_res$core_enrichment[top_pd_whole_level_marker_gseGO_res$Description=='stress response to metal ion'],split = '/')),
                  'detoxification_of_copper_ion'= unlist(strsplit(top_pd_whole_level_marker_gseGO_res$core_enrichment[top_pd_whole_level_marker_gseGO_res$Description=='detoxification of copper ion'],split = '/')),
                  'intracellular_zinc_ion_homeostasis'= unlist(strsplit(top_pd_whole_level_marker_gseGO_res$core_enrichment[top_pd_whole_level_marker_gseGO_res$Description=='intracellular zinc ion homeostasis'],split = '/'))
                  )
# primitive nRBC 可能存在加速功能成熟
#erythrocyte differentiation	GO:0030218
#negative regulation of erythrocyte differentiation	GO:0045649
# heme biosynthetic process	GO:0006783
# hemoglobin biosynthetic process	GO:0042541
# oxygen transport	GO:0015671
cho_gene_list=mapIds(x = org.Hs.eg.db,keys =c( 'GO:0030218','GO:0006783','GO:0042541','GO:0015671'),column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
names(cho_gene_list)=c('erythrocyte_differentiation','heme_biosynthetic_process','hemoglobin_biosynthetic_process','oxygen_transport')
gene_sets=c(gene_sets,cho_gene_list)
set.seed(123)
filt_NBRC_altas_seu <- AddModuleScore_UCell(filt_NBRC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快
saveRDS(filt_NBRC_altas_seu@meta.data,file ='20251125_filt_NBRC_altas_seu_meta.rds' )





cho_pathways=c("heme_biosynthetic_process_UCell","hemoglobin_biosynthetic_process_UCell","oxygen_transport_UCell" ,"stress_response_to_metal_ion_UCell", "regulation_of_autophagy_UCell", "vacuole_organization_UCell","DNA_replication_UCell","regulation_of_leukocyte_proliferation_UCell"    )
p=VlnPlot(filt_NBRC_altas_seu,features = cho_pathways,group.by = 'Ery_stage',raster=FALSE,stack = T,split.by = 'NRBC_type',flip = T)
p$theme <- p$theme + theme(
  axis.text.y = element_text(hjust = 0, vjust = 0.5),
  # 如果是分面的情况，还需要修改strip.text
  strip.text.y = element_text(angle = 0, hjust = 0,face ='plain' )
)

print(p)
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure2/pd_key_pathway_score_stage_vlnplot.pdf',width = 8,height = 8)


test_p_list=list()
# 执行 Wilcoxon 检验
for( cho_pathway in cho_pathways ){
  for( celltype in levels(filt_NBRC_altas_seu$Ery_stage)){
    group1_expr <- filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu@meta.data$NRBC_type == "definitive"& filt_NBRC_altas_seu@meta.data$Ery_stage == celltype, cho_pathway]
    group2_expr <- filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu@meta.data$NRBC_type == "primitive" & filt_NBRC_altas_seu@meta.data$Ery_stage == celltype, cho_pathway]
    wilcox_result <- wilcox.test(group1_expr, group2_expr)
    test_p_list[[cho_pathway]][[celltype]] <- wilcox_result$p.value
    
  }
}

# 都是0 


# HSPC-derive NRBC 最显著特的特征就是免疫调控和末期的细胞复制
pathway=c('chromosome segregation','DNA-templated DNA replication','leukocyte proliferation','positive regulation of cell-cell adhesion','leukocyte cell-cell adhesion','B cell activation','positive regulation of T cell activation')
gene_list=unique(unlist(strsplit(pd_enrichgo_res@compareClusterResult[pd_enrichgo_res@compareClusterResult$Description %in% pathway & pd_enrichgo_res@compareClusterResult$Cluster=='definitive','geneID'],split = '/',fixed = T)))
gene_list_df=pd_whole_level_markers[pd_whole_level_markers$gene %in% gene_list,]

DoHeatmap(filt_NBRC_altas_seu,features = gene_list_df[gene_list_df$avg_log2FC >1 & gene_list_df$pct.2 <0.1,'gene'])
pheatmap(celltype_mexp_df[gene_list_df[gene_list_df$avg_log2FC >1 & gene_list_df$pct.2 <0.1,'gene'],],cluster_cols = F,scale = 'row',color =colorRampPalette(colors = c('navy','white','#C31E1F'))(100))



#-------------------------------------------------------YS vs FL/FBM/ABM--------------------------------------------#
Idents(filt_NBRC_altas_seu)='tissue_stage'
group='YS_FL'
sfile='Protein_NRBC_marker/DE_marker/YS_FL_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = subset(filt_NBRC_altas_seu,tissue_stage %in% c('YS','FL')),group = group,sfile = sfile)
YS_FL_all_Ery_tissue_markers=res[[1]]
YS_FL_count_df=res[[2]]
rm(res);gc()

group='YS_FBM'
sfile='Protein_NRBC_marker/DE_marker/YS_FBM_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = subset(filt_NBRC_altas_seu,tissue_stage %in% c('YS','FBM')),group = group,sfile = sfile)
YS_FBM_all_Ery_tissue_markers=res[[1]]
YS_FBM_count_df=res[[2]]
rm(res);gc()


group='YS_ABM'
sfile='Protein_NRBC_marker/DE_marker/YS_ABM_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = subset(filt_NBRC_altas_seu,tissue_stage %in% c('YS','ABM')),group = group,sfile = sfile)
YS_ABM_all_Ery_tissue_markers=res[[1]]
YS_ABM_count_df=res[[2]]
rm(res);gc()

p=venn.diagram(x = list('FL_YS'=YS_FL_all_Ery_tissue_markers$gene,'FBM_YS'=YS_FBM_all_Ery_tissue_markers$gene,'ABM_YS'=YS_ABM_all_Ery_tissue_markers$gene),
                  filename ='Protein_NRBC_marker/res_pic/main_figure2/HSPC_derived_FLFBMABM_YS_DEGS_venn.tiff' ,fill=cols[1:3],alpha=0.6)
grid.draw(p) # 
length(unique(unlist(list('FL_YS'=YS_FL_all_Ery_tissue_markers$gene,'FBM_YS'=YS_FBM_all_Ery_tissue_markers$gene,'ABM_YS'=YS_ABM_all_Ery_tissue_markers$gene)))) # 14205

YS_FL_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =YS_FL_all_Ery_tissue_markers[YS_FL_all_Ery_tissue_markers$cluster=='FL',] )
YS_FBM_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =YS_FBM_all_Ery_tissue_markers[YS_FBM_all_Ery_tissue_markers$cluster=='FBM',] )
YS_ABM_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =YS_ABM_all_Ery_tissue_markers[YS_ABM_all_Ery_tissue_markers$cluster=='ABM',] )
YS_FL_subcelltype_gseGO_list[[2]]$group='YS_FL'
YS_FBM_subcelltype_gseGO_list[[2]]$group='YS_FBM'
YS_ABM_subcelltype_gseGO_list[[2]]$group='YS_ABM'
degs_gseGO_res_df=YS_FL_subcelltype_gseGO_list[[2]]
degs_gseGO_res_df=rbind(degs_gseGO_res_df,rbind(YS_FBM_subcelltype_gseGO_list[[2]],YS_ABM_subcelltype_gseGO_list[[2]]))
write.csv(degs_gseGO_res_df,file = 'Protein_NRBC_marker/res_data/main_figure1/tissue_nRBC_degs_gseGO_res_df1.csv')

degs_gseGO_res_df=degs_gseGO_res_df[degs_gseGO_res_df$ONTOLOGY=='BP',]
degs_gseGO_res_df$group=factor(degs_gseGO_res_df$group,levels = c('YS_FL','YS_FBM','YS_ABM'))
degs_gseGO_res_df$res='up'
degs_gseGO_res_df$res[degs_gseGO_res_df$NES <0]='dn'
degs_gseGO_res_df$celltype=factor(degs_gseGO_res_df$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))

top_degs_gseGO_bp_res_df=degs_gseGO_res_df %>% group_by(group,celltype,res)  %>%do(head(.,10)) # top_n(wt =-log10(p.adjust),n = 10 ) #
head(sort(table(top_degs_gseGO_bp_res_df$Description),decreasing = T),40)
head(top_degs_gseGO_bp_res_df[,c(1:10,13:14)])
top_degs_gseGO_bp_res_df$Description=factor(top_degs_gseGO_bp_res_df$Description,levels = unique(top_degs_gseGO_bp_res_df$Description))
# ABM NRBC对比fetal NRBC具有免疫调节作用，而fetal NRBC侧重于器官的形态发育，特别是 对血管的发育调控
#ggplot(top_degs_gseGO_bp_res_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+facet_grid(~group)+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+
#  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'gseGO of DEGs')

temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group=='YS_FL',]
temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
p11=ggplot(temp_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'FL vs YS nRBC gseGO of DEGs')+theme(text = element_text(face = 'bold'))

temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group=='YS_FBM',]
temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
p12=ggplot(temp_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'FBM vs YS  nRBC gseGO of DEGs')+theme(text = element_text(face = 'bold'))

temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group=='YS_ABM',]
temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
p13=ggplot(temp_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'ABM vs YS nRBC gseGO of DEGs')+theme(text = element_text(face = 'bold'))

p=p11+p12+p13;p
ggsave(p,width =24 ,height =10,filename='Protein_NRBC_marker/res_pic/main_figure2/DEGS_gseGOBP_substate_YS_HSPC_derived_dotplot.pdf' )

#--------------top20-------------------#
top_degs_gseGO_bp_res_df2=degs_gseGO_res_df %>% group_by(group,celltype,res)  %>%do(head(.,20)) # top_n(wt =-log10(p.adjust),n = 10 ) #
all_shared_pathway_inDefinitive=intersect(intersect(top_degs_gseGO_bp_res_df2[top_degs_gseGO_bp_res_df2$group=='YS_FL',]$Description,
                                                    top_degs_gseGO_bp_res_df2[top_degs_gseGO_bp_res_df2$group=='YS_FBM',]$Description),
                                          top_degs_gseGO_bp_res_df2[top_degs_gseGO_bp_res_df2$group=='YS_ABM',]$Description)
temp_df=degs_gseGO_res_df[degs_gseGO_res_df$Description %in% c(all_shared_pathway_inDefinitive), ]
temp_df=temp_df[order(temp_df$group,temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))

p=ggplot(temp_df,aes(x=group,y=Description ,color=NES,size=-log10(p.adjust)))+geom_point()+facet_grid(~celltype)+
  theme_classic()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+RotatedAxis()
p
ggsave(p,width =8 ,height =8,filename='Protein_NRBC_marker/res_pic/main_figure2/top20_DEGS_gseGOBP_substate_YS_HSPC_derived_shared_dotplot.pdf' )



#  富集分析发现YS nRBC 晚期显著下调细胞有丝分裂多个相关通路，分析 不同niche 中nRBC 细胞周期分布
YS_id_phase_ratio=prop.table(table(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='YS',c('sourceid','Phase')]),margin = 1)
FL_id_phase_ratio=prop.table(table(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='FL',c('sourceid','Phase')]),margin = 1)
FBM_id_phase_ratio=prop.table(table(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='FBM',c('sourceid','Phase')]),margin = 1)
ABM_id_phase_ratio=prop.table(table(filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage=='ABM',c('sourceid','Phase')]),margin = 1)

fetal_id_phase_ratio=rbind(data.frame(YS_id_phase_ratio),rbind(data.frame(FL_id_phase_ratio),data.frame(FBM_id_phase_ratio)))
fetal_id_phase_ratio$sourceid=factor(fetal_id_phase_ratio$sourceid,levels = c(rownames(YS_id_phase_ratio),rownames(FL_id_phase_ratio)[15:21],rownames(FL_id_phase_ratio)[10:14],rownames(FL_id_phase_ratio)[1:9],rownames(FBM_id_phase_ratio)))
p=ggplot(data.frame(fetal_id_phase_ratio),aes(x=sourceid,y=Freq,fill=Phase))+geom_bar(stat = 'identity')+theme_bw()+scale_fill_manual(values =  cols)+theme(axis.text.x = element_text(angle = 90,hjust = 1))
ggsave(p,width =8 ,height =6,filename='Protein_NRBC_marker/res_pic/main_figure2/NRBC_altas_CellCycle_phase_ratio_sourceid_barplot.pdf' )

p=DimPlot(filt_NBRC_altas_seu,split.by = 'tissue_stage',group.by = 'Phase',cols = cols,ncol = 1,raster = F)
ggsave(p,width =4 ,height =8,filename='Protein_NRBC_marker/res_pic/main_figure2/NRBC_altas_CellCycle_phase_umap.pdf' )



p=as.ggplot(pheatmap(t(prop.table(table(filt_NBRC_altas_seu@meta.data[,c('source_celltype','Phase')]),margin = 1)),display_numbers = T,cluster_rows = F,cluster_cols = F,main = 'ratio of Phase',
                     color=colorRampPalette(colors = brewer.pal(11, "PiYG")[c(8,6,3)])(100)))
p
ggsave(p,width =8 ,height =4,filename='Protein_NRBC_marker/res_pic/main_figure2/niche_NRBC_altas_CellCycle_phase_ratio_heatmap.pdf' )



# 查看已il家族受体表达情况
VlnPlot(filt_NBRC_altas_seu ,features =c('IL2RA','IL2RB','IL2RG','IL4R','IL10RA','IL10RB','IL6R','IL6ST','IL10RA'),stack = T,group.by = 'source_celltype',cols = cols)+NoLegend()


# --------------defninitive vs primitive nRBC 在细胞脱核方面的差异, 后期考虑-------------------#


# 1. 准备阶段

#核心浓缩复合体
#NCAPD2	凝缩蛋白I复合物亚基	早期有丝分裂中起始染色质凝集，介导染色质环挤压和轴向压缩
#NCAPG	凝缩蛋白I复合物亚基	与NCAPD2协同，参与染色质纤维的折叠与稳定
#NCAPH	凝缩蛋白I复合物亚基	调节复合物ATP酶活性，控制凝集速率和程度
#SMC2	染色体结构维持蛋白2	凝缩蛋白复合物核心ATP酶，驱动染色质环挤压
#SMC4	染色体结构维持蛋白4	形成SMC2-SMC4异二聚体，是凝缩蛋白复合物的结构骨架
#TOP2A	DNA拓扑异构酶IIα	解除DNA超螺旋和连环化，是染色质凝集的先决条件
# CASP3 Caspase介导的核膜重塑
key_genes1=c(c('NCAPD2','NCAPG','NCAPH','SMC2','SMC4','TOP2A'))
VlnPlot(filt_NBRC_altas_seu,stack = T,features =key_genes1 ,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)

#组蛋白修饰与表观调控		
#H3S10ph	组蛋白H3 Ser10磷酸化	有丝分裂标志性修饰，由Aurora B激酶催化，促进染色质凝集起始,H3C1 
#H3S28ph	组蛋白H3 Ser28磷酸化	与H3S10ph协同，在有丝分裂前期和中期增强染色质压缩
#AURKB	Aurora激酶B	催化H3S10ph和H3S28ph，是凝集启动的关键激酶
#AURKA	Aurora激酶A	调控中心体成熟和纺锤体组装，间接影响凝集空间模式
#HAT1	组蛋白乙酰转移酶1	调控组蛋白H4 Lys5/K12乙酰化，影响凝集起始和核小体间接触
#HP1α (CBX5)	异染色质蛋白1α	识别H3K9me3，维持异染色质区域在凝集中的紧密度
#组蛋白去乙酰化	HDAC2、HDAC6
#NA甲基化调控	TET2、TET3

#染色质结构调控		
#HMGN1/2	高迁移率族核小体结合蛋白	促进染色质去凝集；在终末分化中其下调允许染色质凝集
#LMNA	核纤层蛋白A/C	核纤层解体是染色质凝集的前提；LMNA与染色质相互作用
#BANF1	屏障自整合因子1	在凝集后期调控核膜破裂和染色质-核膜连接解离
key_genes2=c(c('AURKA','AURKB','HAT1','CBX5',	'HDAC2','HDAC6','TET2','TET3','HMGN1','HMGN2','LMNA','BANF1'))
VlnPlot(filt_NBRC_altas_seu,stack = T,features =key_genes2,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)

# 2. 启动执行阶段
#细胞骨架重排与收缩环形成通路
# 收缩环相关：MYH9/MYH10、DIAPH1、ROCK
#	MYH9 / MYH10	非肌肉肌球蛋白重链，构成收缩环，提供核挤出所需的机械力	敲除或抑制MYH10导致脱核严重受损
# Rho GTPase，驱动收缩环组装：RAC1,RAC2,	Rho激酶，调节肌动蛋白骨架收缩及细胞极性	ROCK抑制剂处理显著降低脱核效率
# DIAPH1 / DIAPH3	成核因子，调控肌动蛋白聚合与收缩环组装	在红系终末分化中高表达，失活导致脱核失败
## ACTB / ACTG1	肌动蛋白，构成收缩环及细胞皮层骨架	终末期红系细胞中维持高表达
key_genes3=c(c('MYH9','MYH10','DIAPH1','DIAPH3','RAC1','RAC2','ACTB','ACTG1'))
VlnPlot(filt_NBRC_altas_seu,stack = T,features = key_genes3,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)

# LINC复合物将细胞骨架的力量传递给细胞核
key_genes32=c(c('SYNE1','SYNE2','SUN1','SUN2'))



#三、钙离子信号与收缩力激活通路
#蛋白名称	基因Symbol	基因ID	主要功能
#钙调蛋白1	CALM1	801	钙感受器，激活下游激酶
#钙调蛋白2	CALM2	805	与CALM1功能冗余
#钙调蛋白3	CALM3	808	与CALM1功能冗余
#钙/钙调蛋白依赖性蛋白激酶IIα	CAMK2A	815	调控肌动蛋白重组
#钙/钙调蛋白依赖性蛋白激酶IIβ	CAMK2B	816	调控收缩信号
#钙/钙调蛋白依赖性蛋白激酶IIδ	CAMK2D	817	参与肌球蛋白激活
#钙/钙调蛋白依赖性蛋白激酶IIγ	CAMK2G	818	调控细胞骨架
#肌球蛋白轻链激酶	MYLK	4638	磷酸化肌球蛋白轻链，激活收缩
#细胞周期蛋白依赖性激酶1	CDK1	983	与Cyclin B1形成复合物，启动有丝分裂样程序
#细胞周期蛋白B1	CCNB1	891	CDK1的调节亚基
#核转运蛋白β1	KPNB1	3837	与CDK9互作，在钙信号上游发挥作用
#细胞周期蛋白依赖性激酶9	CDK9	1025	与Importin β复合物互作调控脱核
# 钙调蛋白 在primtiive NRBC 中持续高表达至Orth，而在definiitve NRBC 在Orrh表达出现明显下调，钙/钙调蛋白依赖性蛋白激酶几乎不表达，
# 细胞周期蛋白 相关基因几乎都在在primtiive 早期表达，中期下调，而在definitive NRBC 多持续表达至晚期
key_genes4=c(c('CALM1','CALM2','CALM3','CAMK2A',	'CAMK2B','CAMK2D','MYLK','CDK1','CCNB1','KPNB1','CDK9'))
VlnPlot(filt_NBRC_altas_seu,stack = T,features = key_genes4,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)

VlnPlot(filt_NBRC_altas_seu,stack = T,features = c(key_genes1,key_genes2,key_genes3,key_genes32,key_genes4),group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)

temp_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$celltype=='early_Ery' & sub_pd_all_Ery_tissue_markers$ge %in%  c(key_genes1,key_genes2,key_genes3,'VRK1') ,]
temp_df=temp_df[temp_df$avg_log2FC >0,];temp_df$gene=factor(temp_df$gene,level=c(key_genes1,key_genes2,key_genes3,'VRK1'))
p=ggplot(temp_df,aes(x=gene,y=avg_log2FC,fill= cluster))+geom_bar(stat = 'identity')+theme_classic()+RotatedAxis()
ggsave(p,width =12 ,height =4,filename='Protein_NRBC_marker/res_pic/main_figure2/pd_nuclear_keygene_DE_expression.pdf' )

p=VlnPlot(filt_NBRC_altas_seu,stack = T,features = c(key_genes1,key_genes2,key_genes3,'VRK1'),group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols);p
ggsave(p,width =12 ,height =6,filename='Protein_NRBC_marker/res_pic/main_figure2/pd_nuclear_keygene_VlnPlot_expression.pdf' )

# 脱核关键调控因子
# ZMPSTE24	参与核纤层蛋白A前体的加工成熟	缺失导致核异常及脱核障碍
# 	Rap1 GTPase激活蛋白，缺失导致小鼠胚胎致死及红系脱核缺陷，在primitive NRBC 中显著高表达
# HMGN1	
#红系终末分化信号（EPO / GATA1 / KLF1）
#↓
#TSPO2 / AURKA/B 激活 → 核浓缩及核膜重塑
#↓
#DIAPH / RAG / MYH10 → 收缩环组装与激活
#↓
#ACTIN聚合 → 机械力产生 → 核挤出
#↓
#自噬（ATG7 / BNIP3L / LC3B）清除多余细胞器
#↓
#网织红细胞释放

VlnPlot(filt_NBRC_altas_seu,stack = T,features = c('HMGN1','HMGN2','ZMPSTE24'),group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)
DotPlot(filt_NBRC_altas_seu,features = c('LMNA','MYH9','MYH10','ROCK1','ROCK2','AURKA','AURKB','DIAPH1','DIAPH3'),group.by = 'source_celltype')+RotatedAxis()

# 涉及GO：BP 富集通路
# GO:0043131	红细胞脱核 (erythrocyte enucleation）
# GO:0043353	去核红细胞分化 (enucleate erythrocyte differentiation)
# GO:0061930	红细胞脱核调节 (regulation of erythrocyte enucleation)
# "GO:0030036",          # actin cytoskeleton organization
#GO:0061931    positive regulation of erythrocyte enucleation	positively_regulates
#GO:0061930    regulation of erythrocyte enucleation	regulates
#GO:0061932    negative regulation of erythrocyte enucleation	negatively_regulates

GOBP_enucleated_pahtways=c('erythrocyte enucleation','enucleate erythrocyte differentiation','regulation of erythrocyte enucleation','actin cytoskeleton organization',
                           'positive regulation of erythrocyte enucleation	positively_regulates',
                           'negative regulation of erythrocyte enucleation	negatively_regulates','enucleation','enucleate erythrocyte maturation','enucleate erythrocyte development')
names(GOBP_enucleated_pahtways)=c('GO:0043131','GO:0043353','GO:0061930',"GO:0030036",'GO:0061931','GO:0061932','GO:0090601','GO:0043354','GO:0048822')
gene_sets=mapIds(x = org.Hs.eg.db,keys =names(GOBP_enucleated_pahtways),column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
names(gene_sets)=as.character(GOBP_enucleated_pahtways)

gene_sets$`GO:0043131`=c('NEMP1','PLEK2')
# TSPO2 最明显
VlnPlot(filt_NBRC_altas_seu,stack = T,features = gene_sets$`enucleate erythrocyte differentiation`,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)
VlnPlot(filt_NBRC_altas_seu,stack = T,features = gene_sets$`enucleate erythrocyte development`,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols)




#################################################################################################################################################################
#------------------------------------------figure 3 ：the comparision of HSPC_derived_NRBC -------------------------------------------#
#################################################################################################################################################################

filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,NRBC_type=='definitive')
Idents(filt_NBRC_altas_seu)='tissue_stage'
HSPC_derived_NRBC_DE_res=FindAllMarkers(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'))
saveRDS(HSPC_derived_NRBC_DE_res,file ='Protein_NRBC_marker/DE_marker/HSPC_derived_nRBC_wholelevel_RNA_markers.rds' )

# 基因信息中挑选具有特异性的marker

top_gene_markers_HSPC_derived_NRBC=HSPC_derived_NRBC_DE_res[-grep('^AC[0-9]|CH507-|^LIN',HSPC_derived_NRBC_DE_res$gene),] %>% filter(avg_log2FC>1 & pct.2 < 0.2 & pct.1 > 0.2) %>% group_by(cluster)  %>%top_n(wt=avg_log2FC,10)# %>%  do(head(., n = 10))
top_gene_markers_HSPC_derived_NRBC=top_gene_markers_HSPC_derived_NRBC[order(top_gene_markers_HSPC_derived_NRBC$cluster,top_gene_markers_HSPC_derived_NRBC$avg_log2FC,decreasing = T),]
p=DotPlot(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'),group.by = 'source_celltype',features =unique(top_gene_markers_HSPC_derived_NRBC$gene)[c(21:30,11:20,1:10)],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p

FL_top10_markers=unique(top_gene_markers_HSPC_derived_NRBC$gene)[21:30]

#top10 marker 中， FBM 明显高表达DNA复制依赖的组蛋白，这显然不合生物学事实
# 分析原因：10X 测序采用的依赖ployA尾巴磁珠捕获mRNA的技术，而依赖复制性的组蛋白基因，是无polyA尾巴，采用3’末端茎环结构+茎环结合蛋白替代polyA尾巴
VlnPlot(filt_NBRC_altas_seu,group.by = 'source',split.by = 'tissue_stage',cols = cols,stack = T,features = c( "HIST1H4K","HIST2H2AA4" ,"HIST2H2AA3", "HIST1H2AI"))
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols,stack = T,features = c("H2AFJ","H2AFV","H2AFX","H2AFY","H2AFY2","H2AFZ" )) # 非复制性型全部捕获，表达高
HIS_genes=mapIds(x = org.Hs.eg.db,keys = rownames(filt_NBRC_altas_seu)[grep('^HIS',rownames(filt_NBRC_altas_seu))],keytype = 'ALIAS',column = 'SYMBOL') 
HIS_genes=HIS_genes[as.character(HIS_genes)!=names(HIS_genes)] # 全部为复制依赖行组蛋白，仅仅在FBM中检测高，
#排除复制组蛋白基因的干扰
HSPC_derived_NRBC_DE_res=HSPC_derived_NRBC_DE_res[!HSPC_derived_NRBC_DE_res$gene %in% names(HIS_genes),]

table(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$avg_log2FC >0,'cluster'])
#   FL   FBM   ABM 
# 10593   969  2186
# FBM 差异基因少

# FBM NRBC 比较缺乏特异表达基因

# 查看top20
top_gene_markers_HSPC_derived_NRBC=HSPC_derived_NRBC_DE_res[-grep('^AC[0-9]|CH507-|^LIN',HSPC_derived_NRBC_DE_res$gene),] %>% filter(avg_log2FC>1 & pct.2 < 0.2 & pct.1 > 0.2) %>% group_by(cluster)  %>%top_n(wt=avg_log2FC,20)
p=DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =top_gene_markers_HSPC_derived_NRBC$gene[top_gene_markers_HSPC_derived_NRBC$cluster=='FBM'],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
# 发现FBM 缺乏特异标志基因，EIF3C 是翻译起始因子，ATP6V0C是	V-ATP酶关键亚基，都属于关键基因，BOLA2/B铁硫簇组装相关，对细胞基本能量代谢至关重要也倾向为管家基因 
# EIF3CL	否	EIF3C的假基因同源物。其功能本身不明确，且假基因的表达通常具有组织或阶段特异性，并非持续表达。
# DDTL	未知	功能未明确基因。缺乏足够的功能和表达谱数据，无法归类为管家基因
#U2AF1L5	否	剪接因子类似物。RNA剪接是基础过程，但该特定因子功能研究少，可能具有组织或阶段特异性，非通用
# GET4	蛋白质靶向与质量控制，参与所有真核细胞必需的、持续进行的蛋白质定向运输与质量控制，表达广泛且稳定。
#SMN1	RNA加工：编码“运动神经元存活蛋白”，是剪接体小核核糖核蛋白复合物组装的核心 其功能对细胞生存绝对基础，但作为已分化细胞中RNA加工机器的调节核心，其表达水平可能因细胞类型和状态而有适应性变化，非绝对恒定。
#SMN2	RNA加工（功能缺陷）：是SMN1的旁系同源基因，因关键位点变异，主要产生截短的不稳定蛋白，仅贡献少量全长功能蛋白。
#U2AF1 是一个在RNA剪接中扮演核心角色的基因，其编码的蛋白是剪接体组装过程中的一个关键通用剪接因子。
# SLX1A/SLX1B DNA结构特异性修复。DNA完整性的维持是基础，但其激活具有应激性（如复制压力），表达可能不像核糖体基因那样绝对恒定，在其他来源FL中也表达
VlnPlot(UCB_NRBC_altas,group.by = 'celltype',features = top_gene_markers_HSPC_derived_NRBC$gene[top_gene_markers_HSPC_derived_NRBC$cluster=='FBM'],stack = T,split.by = 'type')

del_gene=c('EIF3C','EIF3CL','U2AF1','ATP6V0C','DDTL','BOLA2','BOLA2B','U2AF1L5','GET4','SMN1','SMN2','SLX1A','SLX1B', "SBF2-AS1",'PHOSPHO1')#  "SBF2-AS1" anti-mRNA, RGS6 表达太低
FBM_unique_genes=top_gene_markers_HSPC_derived_NRBC[!top_gene_markers_HSPC_derived_NRBC$gene %in% del_gene & top_gene_markers_HSPC_derived_NRBC$cluster=='FBM', ]$gene

if(F){
  fl_fbm_pos_marker=FindMarkers(filt_NBRC_altas_seu,ident.1 = 'FBM',ident.2 = 'FL')
  fbm_pos_marker=fl_fbm_pos_marker[fl_fbm_pos_marker$avg_log2FC >1 & fl_fbm_pos_marker$pct.1>0.1 & fl_fbm_pos_marker$pct.2 <0.3,]
  fbm_pos_marker=fbm_pos_marker[!rownames(fbm_pos_marker) %in% c(names(HIS_genes),del_gene),]
  fbm_pos_marker=fbm_pos_marker[-grep('^CH507|AC00',rownames(fbm_pos_marker)),]
  dim(fbm_pos_marker)# 69
  DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =rownames(fbm_pos_marker),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
  # FBM NRBC缺乏特异
  
  
  FL_pos_marker=fl_fbm_pos_marker[abs(fl_fbm_pos_marker$avg_log2FC )>1 & fl_fbm_pos_marker$pct.2>0.2 & fl_fbm_pos_marker$pct.1 <0.3,]
  FL_pos_marker=FL_pos_marker[!rownames(FL_pos_marker) %in% c(names(HIS_genes),del_gene),]
  FL_pos_marker=FL_pos_marker[-grep('^CH507|AC00',rownames(FL_pos_marker)),]
  dim(FL_pos_marker) # 1563
  FL_pos_marker=FL_pos_marker[order(FL_pos_marker$avg_log2FC),]
  DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =rownames(FL_pos_marker)[1:60],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
  # FL NRBC缺乏特异高表达基因, 差异表达
}


candidated_ABM_specific_genes1  =HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='ABM' &HSPC_derived_NRBC_DE_res$avg_log2FC >2 & HSPC_derived_NRBC_DE_res$pct.1 >0.1 & HSPC_derived_NRBC_DE_res$pct.2<0.2 ,]
candidated_ABM_specific_genes1=candidated_ABM_specific_genes1[order(candidated_ABM_specific_genes1$avg_log2FC,decreasing = T),];length(candidated_ABM_specific_genes1$gene)
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =candidated_ABM_specific_genes1$gene[1:20],stack = T,split.by = 'tissue_stage')
#"AC005943.2"    "RP11-411B6.6"     "RP11-111K18.1" , "LINC00570"， ， "RP11-354E11.2"非编码RNA特异高表达，暂时不考虑这些基因,KIAA0125 : USP45,泛素化酶,DNA损伤修复，这是USP45最明确的功能之一。
ABM_cho_topmarkers=c( 'HBD',"CA1","PDZK1IP1", "ANXA1","NECAB1","ANKRD28","TSC22D3","IFIT1B",'LGALS9',"HLA-B","HLA-DRA","HLA-DRB1") # top10 中挑选基因,
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =ABM_cho_topmarkers,stack = T,split.by = 'tissue_stage')



whole_celltype_mexp=AverageExpression(filt_NBRC_altas_seu,group.by = 'tissue_stage')$RNA
fetal_unique_marker_genes=HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='ABM' & HSPC_derived_NRBC_DE_res$avg_log2FC <0,]
fetal_unique_marker_genes=fetal_unique_marker_genes[fetal_unique_marker_genes$avg_log2FC < -2 & fetal_unique_marker_genes$pct.1<0.3 & fetal_unique_marker_genes$pct.2 >0.2, ]
fetal_unique_marker_genes$score=-1*fetal_unique_marker_genes$avg_log2FC *fetal_unique_marker_genes$pct.2/(fetal_unique_marker_genes$pct.1+0.001)*whole_celltype_mexp[fetal_unique_marker_genes$gene,'FL']/(whole_celltype_mexp[fetal_unique_marker_genes$gene,'ABM']+0.001)*(fetal_unique_marker_genes$pct.2-fetal_unique_marker_genes$pct.1)
fetal_unique_marker_genes$score=log2(fetal_unique_marker_genes$score+1)
fetal_unique_marker_genes=fetal_unique_marker_genes[order(fetal_unique_marker_genes$score,-1*fetal_unique_marker_genes$avg_log2FC,decreasing = T),]
fetal_unique_marker_genes=fetal_unique_marker_genes[!fetal_unique_marker_genes$gene %in% top_gene_markers_HSPC_derived_NRBC$gene , ]
dim(fetal_unique_marker_genes)# 133
saveRDS(fetal_unique_marker_genes,file = 'Protein_NRBC_marker/res_data/main_figure3/fetal_unique_marker_genes.rds')

DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =unique(fetal_unique_marker_genes$gene)[1:40],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =unique(fetal_unique_marker_genes$gene)[1:20],stack = T)
temp_df=fetal_unique_marker_genes[1:20,]
temp_df=temp_df[order(temp_df$pct.2,decreasing = T),]
temp_df

VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =unique(temp_df$gene)[1:20],stack = T)

top_fetal_unique_marker_genes=c('HBG1','HBG2','HBZ', "TUBB6","HSPA1A","HSPA1B",'IGF2BP1','IGF2BP3','DLK1','CISH','HIF3A')# 还可以考虑CHD7,TIMP3,CISH,后面好像有获得

p=DotPlot(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'),group.by = 'source_celltype',features =c(FL_top10_markers,FBM_unique_genes,
          top_fetal_unique_marker_genes,ABM_cho_topmarkers),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure3/top_marker_HSPC_derived_nRBC.pdf',width = 16,height = 8)


# FL/FBM nRBC 在整体水平上缺乏特异maker，进一步在亚类水平比较,已经在整体水平找到少量FL 分子印记marker
if(F){
  early_Ery_markers=FindAllMarkers(subset(subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM','ABM')),Ery_stage=='early_Ery'))
  mid_Ery_markers=FindAllMarkers(subset(subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM','ABM')),Ery_stage=='mid_Ery'))
  late_Ery_markers=FindAllMarkers(subset(subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM','ABM')),Ery_stage=='late_Ery'))
  
  early_Ery_markers$membrane='no'
  early_Ery_markers$membrane[early_Ery_markers$gene %in% new_all_merged_plasma_protein ]='yes'
  mid_Ery_markers$membrane='no'
  mid_Ery_markers$membrane[mid_Ery_markers$gene %in% new_all_merged_plasma_protein ]='yes'
  late_Ery_markers$membrane='no'
  late_Ery_markers$membrane[late_Ery_markers$gene %in% new_all_merged_plasma_protein ]='yes'
  
  early_Ery_markers$celltype='early_Ery'
  mid_Ery_markers$celltype='mid_Ery'
  late_Ery_markers$celltype='late_Ery'
  
  subcelltype_Ery_markers=rbind(early_Ery_markers,rbind(mid_Ery_markers,late_Ery_markers))
  subcelltype_Ery_markers=subcelltype_Ery_markers[subcelltype_Ery_markers$avg_log2FC >0,]
  FL_subcelltype_Ery_markers=subcelltype_Ery_markers[subcelltype_Ery_markers$cluster=='FL',]
  high_FL_subcelltype_Ery_markers=FL_subcelltype_Ery_markers[FL_subcelltype_Ery_markers$avg_log2FC>1 & FL_subcelltype_Ery_markers$pct.1 >0.3 & FL_subcelltype_Ery_markers$pct.2 <0.2,]
  high_FL_subcelltype_Ery_markers=high_FL_subcelltype_Ery_markers[order(high_FL_subcelltype_Ery_markers$avg_log2FC,decreasing = T),]
  length(unique(high_FL_subcelltype_Ery_markers$gene))# 272
  
  left_high_FL_subcelltype_Ery_markers=unique(high_FL_subcelltype_Ery_markers$gene)[rowMaxs(as.matrix(celltype_mexp_df[unique(high_FL_subcelltype_Ery_markers$gene),10:19])) <1]
  high_FL_subcelltype_Ery_markers=high_FL_subcelltype_Ery_markers[high_FL_subcelltype_Ery_markers$gene %in% left_high_FL_subcelltype_Ery_markers, ]
  top_high_FL_subcelltype_Ery_markers=names(sort(rowMaxs(as.matrix(celltype_mexp_df[unique(high_FL_subcelltype_Ery_markers$gene),10:19])),decreasing = F))
  
  # 除了前三个基因，其他都有表达，FL 面对FBM、ABM缺乏特异表达marker，有特异表达蛋白有三个：AFP，SERPINA1，APOA1
  DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =top_high_FL_subcelltype_Ery_markers,cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
  top_high_FL_subcelltype_Ery_markers=top_high_FL_subcelltype_Ery_markers[1:3]
  DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =top_high_FL_subcelltype_Ery_markers,cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
  
  p=DotPlot(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'),group.by = 'source_celltype',features =c(unique(top_gene_markers_HSPC_derived_NRBC$gene),'AFP','SERPINA1','APOA1'),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
  p
}

HSPC_nRBC_enrichgo_res=compareCluster(geneClusters =list('FL'=unique(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='FL' & HSPC_derived_NRBC_DE_res$avg_log2FC >0,'gene']),
                                               'FBM'=unique(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='FBM'& HSPC_derived_NRBC_DE_res$avg_log2FC >0,'gene']),
                                               'ABM'=unique(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='ABM'& HSPC_derived_NRBC_DE_res$avg_log2FC >0,'gene']) ),keyType = 'ALIAS',ont = "BP",
                            fun = 'enrichGO',  OrgDb='org.Hs.eg.db' )
saveRDS(HSPC_nRBC_enrichgo_res,file = 'Protein_NRBC_marker/res_data/main_figure3/HSPC_nRBC_enrichgo_res.rds')

p=dotplot(HSPC_nRBC_enrichgo_res,showCategory=10);p
# ggplot美化
top10_enrichgo_res_df=HSPC_nRBC_enrichgo_res@compareClusterResult %>% group_by(Cluster) %>% do(head(.,15))
top10_enrichgo_res_df$ratio=as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[1,])/as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[2,])
top10_enrichgo_res_df$Description=factor(top10_enrichgo_res_df$Description,levels =unique(top10_enrichgo_res_df$Description) )
p=ggplot(top10_enrichgo_res_df,aes(x=Cluster,y=Description,color=-log10(p.adjust),size=ratio))+geom_point()+theme_bw()+scale_color_gradient(low = '#4387B5',high = 'firebrick3')
p 

ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure3/comapreenrichGO_HSPC_derived_nRBC.pdf',width = 8,height = 10)



#------------------------------------------------------ HSPC_derived fetal vs adult DE analysis---------------------------------------------------#
if(F){
  # fetal NRBC: YS、FL、FBM，三个阶段各抽取1:1:1的细胞，构成early、mid、late NRBC， 与整体不抽样分析，结果差异很小
  
  # YS_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='YS',downsample =1000)
  # YS_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='YS',downsample =1000)
  # YS_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='YS',downsample =1000)
  # YS_filt_NBRC_altas_seu=merge(YS_early_Ery_filt_NBRC_altas_seu,c(YS_mid_Ery_filt_NBRC_altas_seu,YS_late_Ery_filt_NBRC_altas_seu))
  # rm(YS_early_Ery_filt_NBRC_altas_seu,YS_mid_Ery_filt_NBRC_altas_seu,YS_late_Ery_filt_NBRC_altas_seu)
  
  FL_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FL',downsample =1500)
  FL_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FL',downsample =1500)
  FL_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FL',downsample =1500)
  FL_filt_NBRC_altas_seu=merge(FL_early_Ery_filt_NBRC_altas_seu,c(FL_mid_Ery_filt_NBRC_altas_seu,FL_late_Ery_filt_NBRC_altas_seu))
  rm(FL_early_Ery_filt_NBRC_altas_seu,FL_mid_Ery_filt_NBRC_altas_seu,FL_late_Ery_filt_NBRC_altas_seu)
  
  FBM_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FBM',downsample =1500)
  FBM_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FBM',downsample =1500)
  FBM_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FBM',downsample =1500)
  FBM_filt_NBRC_altas_seu=merge(FBM_early_Ery_filt_NBRC_altas_seu,c(FBM_mid_Ery_filt_NBRC_altas_seu,FBM_late_Ery_filt_NBRC_altas_seu))
  rm(FBM_early_Ery_filt_NBRC_altas_seu,FBM_mid_Ery_filt_NBRC_altas_seu,FBM_late_Ery_filt_NBRC_altas_seu)
  
  ABM_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='ABM',downsample =3000)
  ABM_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='ABM',downsample =3000)
  ABM_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='ABM',downsample =3000)
  ABM_filt_NBRC_altas_seu=merge(ABM_early_Ery_filt_NBRC_altas_seu,c(ABM_mid_Ery_filt_NBRC_altas_seu,ABM_late_Ery_filt_NBRC_altas_seu))
  rm(ABM_early_Ery_filt_NBRC_altas_seu,ABM_mid_Ery_filt_NBRC_altas_seu,ABM_late_Ery_filt_NBRC_altas_seu)
  
  #subset_filt_NBRC_altas_seu=merge(ABM_filt_NBRC_altas_seu,c(FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu,YS_filt_NBRC_altas_seu))
  #rm(ABM_filt_NBRC_altas_seu,FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu,YS_filt_NBRC_altas_seu)
  
  subset_filt_NBRC_altas_seu=merge(ABM_filt_NBRC_altas_seu,c(FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu))
  rm(ABM_filt_NBRC_altas_seu,FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu)
  
  gc()
  
  subset_filt_NBRC_altas_seu <- JoinLayers(subset_filt_NBRC_altas_seu)
  subset_filt_NBRC_altas_seu$type_stage='fetal'
  subset_filt_NBRC_altas_seu$type_stage[subset_filt_NBRC_altas_seu$tissue_stage=='ABM']='adult'
  Idents(subset_filt_NBRC_altas_seu)='type_stage'
  subset_filt_NBRC_altas_seu$source_celltype=factor(subset_filt_NBRC_altas_seu$source_celltype,levels =c( "FL_BFUE/CFUE","FL_ProE","FL_Bas","FL_Poly", "FL_Orth","FBM_BFUE/CFUE","FBM_ProE","FBM_Bas","FBM_Poly",  "FBM_Orth",
                                                                                                          "ABM_BFUE/CFUE", "ABM_ProE","ABM_Bas","ABM_Poly","ABM_Orth") )
  saveRDS(subset_filt_NBRC_altas_seu,file = 'Protein_NRBC_marker/res_data/temp_subset_filt_NBRC_altas_seu.rds')
  
}else{
  subset_filt_NBRC_altas_seu=readRDS('Protein_NRBC_marker/res_data/temp_subset_filt_NBRC_altas_seu.rds')
}


# 或者采用FL、FBM以及ABM中ABM vs FL/FBM 得到的degs list
fetal_adult_NRBC_whole_marker=FindAllMarkers(subset_filt_NBRC_altas_seu,group.by = 'type_stage')
     
fetal_marker_list=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$cluster=='fetal','avg_log2FC']
names(fetal_marker_list)=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$cluster=='fetal','gene']
fetal_marker_list=sort(fetal_marker_list,decreasing = T)
fetal_adult_gsego_res=gseGO(geneList = fetal_marker_list,OrgDb = org.Hs.eg.db,ont = 'BP',keyType = 'ALIAS')
dotplot(fetal_adult_gsego_res)
saveRDS(fetal_adult_gsego_res,file='Protein_NRBC_marker/res_data/main_figure3/wholelevel_fetal_adult_gsego_res.rds')

fetal_adult_gsego_res_df=fetal_adult_gsego_res@result;fetal_adult_gsego_res_df$res='up';fetal_adult_gsego_res_df$res[fetal_adult_gsego_res_df$NES <0]='down'
fetal_adult_gsego_res_df=fetal_adult_gsego_res_df[fetal_adult_gsego_res_df$p.adjust <0.05,] %>% group_by(res) %>% do(head(.,10))
fetal_adult_gsego_res_df=fetal_adult_gsego_res_df[order(fetal_adult_gsego_res_df$res,fetal_adult_gsego_res_df$NES),]
fetal_adult_gsego_res_df$Description=factor(fetal_adult_gsego_res_df$Description,levels = fetal_adult_gsego_res_df$Description)
p=ggplot(fetal_adult_gsego_res_df,aes(x=res,y=Description,size=-log(p.adjust),color=NES))+geom_point()+theme_bw()+scale_color_gradient(low = '#4387B5',high = 'firebrick3')+ggtitle(label = 'Fetal vs Adult')
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure3/top10_fetal_adult_gseGO.pdf',width = 6,height = 6)



Idents(subset_filt_NBRC_altas_seu)='type_stage'
group='fetal_adult'
sfile='Protein_NRBC_marker/DE_marker/fetal_adult_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = subset(subset_filt_NBRC_altas_seu,type_stage %in% c('fetal','adult')),group = group,sfile = sfile)
sub_fetal_adult_all_Ery_tissue_markers=res[[1]]
sub_fetal_addult_count_df=res[[2]]
rm(res);gc()
sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[!sub_fetal_adult_all_Ery_tissue_markers$gene %in% names(HIS_genes),]

late_positive_sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$celltype=='late_Ery' & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >0,]
late_positive_sub_fetal_adult_all_Ery_tissue_markers=late_positive_sub_fetal_adult_all_Ery_tissue_markers[order(late_positive_sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC,decreasing = T),]
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=late_positive_sub_fetal_adult_all_Ery_tissue_markers[1:100,]
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=top_late_positive_sub_fetal_adult_all_Ery_tissue_markers[order(top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC*top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$pct.1,decreasing = T),]
DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',late_positive_sub_fetal_adult_all_Ery_tissue_markers$gene[1:50])+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=top_late_positive_sub_fetal_adult_all_Ery_tissue_markers[1:8,]
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=top_late_positive_sub_fetal_adult_all_Ery_tissue_markers[!top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$gene %in% ABM_cho_topmarkers,]
DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$gene)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
# LGALS3 
ABM_cho_topmarkers=c(ABM_cho_topmarkers[1:9],'LGALS3',ABM_cho_topmarkers[10:12])
p=DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =c(FL_top10_markers,FBM_unique_genes,top10_fetal_unique_marker_genes,ABM_cho_topmarkers),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure3/top_marker_HSPC_derived_nRBC.pdf',width = 16,height = 8)


fetal_positive_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$cluster=='fetal',]
fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[fetal_positive_all_Ery_tissue_markers$avg_log2FC >0,]
fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[fetal_positive_all_Ery_tissue_markers$pct.1>0.1,]
fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[order(fetal_positive_all_Ery_tissue_markers$avg_log2FC,decreasing = T),]
top_fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[1:50,]
top_fetal_positive_all_Ery_tissue_markers=top_fetal_positive_all_Ery_tissue_markers[order(top_fetal_positive_all_Ery_tissue_markers$pct.1,decreasing = T),]
DotPlot(filt_NBRC_altas_seu,features = unique(top_fetal_positive_all_Ery_tissue_markers$gene),cols = c('gray','firebrick3'),scale=F)+RotatedAxis()
# 胚胎印记基因："IGF2BP1"，"DLK1"，"HIF3A"，"GATA5"[调控血管以及心脏发育关键转录因子] ，"MEG3", DLK1-DIO3印记簇核心的一个长链非编码RNA（lncRNA），具有明确的母源表达特征
#LIN28B	RNA结合蛋白（上游调控者）	核心“开关”，抑制let-7，维持胎儿程序
# HIF3A 广泛意义上的胎儿时期环境影响标志,
# HMGA2 , LIN28B 下游靶基因
top_fetal_unique_marker_genes=c(top_fetal_unique_marker_genes[1:9],'MEG3','GATA5','LIN28B','HMGA2',top_fetal_unique_marker_genes[10:11])

defintive_markers=list(FL_top10_markers=FL_top10_markers,FBM_unique_genes=FBM_unique_genes,top_fetal_unique_marker_genes=top_fetal_unique_marker_genes,ABM_cho_topmarkers=ABM_cho_topmarkers)
p=DotPlot(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'),group.by = 'source_celltype',features =as.character(unlist(defintive_markers)),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure3/top_marker_HSPC_derived_nRBC.pdf',width = 16,height = 8)

saveRDS(defintive_markers,file = 'Protein_NRBC_marker/res_data/main_figure3/defintive_markers.rds')


sub_fetal_adult_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$cluster=='fetal',] )
sub_fetal_adult_subcelltype_gseGO_list[[2]]$group='fetal_adult'
saveRDS(sub_fetal_adult_subcelltype_gseGO_list,file = 'Protein_NRBC_marker/res_data/main_figure2/sub_fetal_adult_subcelltype_gseGO_list.rds')

top_sub_fetal_adult_subcelltype_gseGO_df=sub_fetal_adult_subcelltype_gseGO_list[[2]]
top_sub_fetal_adult_subcelltype_gseGO_df=top_sub_fetal_adult_subcelltype_gseGO_df[top_sub_fetal_adult_subcelltype_gseGO_df$ONTOLOGY=='BP', ]
top_sub_fetal_adult_subcelltype_gseGO_df$ratio=sapply(strsplit(top_sub_fetal_adult_subcelltype_gseGO_df$core_enrichment,split = '/'), length)/top_sub_fetal_adult_subcelltype_gseGO_df$setSize

top_sub_fetal_adult_subcelltype_gseGO_df$res='up'
top_sub_fetal_adult_subcelltype_gseGO_df$res[top_sub_fetal_adult_subcelltype_gseGO_df$NES <0]='down'
top_sub_fetal_adult_subcelltype_gseGO_df=top_sub_fetal_adult_subcelltype_gseGO_df %>% group_by(celltype,res) %>% do(head(.,10))
top_sub_fetal_adult_subcelltype_gseGO_df$celltype=factor(top_sub_fetal_adult_subcelltype_gseGO_df$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))

top_sub_fetal_adult_subcelltype_gseGO_df=top_sub_fetal_adult_subcelltype_gseGO_df[order(top_sub_fetal_adult_subcelltype_gseGO_df$celltype,top_sub_fetal_adult_subcelltype_gseGO_df$NES),]
top_sub_fetal_adult_subcelltype_gseGO_df$Description=factor(top_sub_fetal_adult_subcelltype_gseGO_df$Description,levels = unique(top_sub_fetal_adult_subcelltype_gseGO_df$Description))
p5=ggplot(top_sub_fetal_adult_subcelltype_gseGO_df,aes(x=celltype,y=Description,size=ratio,color=NES))+geom_point()+scale_color_gradient2(low = 'navy',mid = 'white',high = 'firebrick3')+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,face = 'bold'))+ggtitle('fetal vs adult nRBC:gseGO of DEGs')
p5

ggsave(p5,width =6 ,height =8,filename='Protein_NRBC_marker/res_pic/main_figure3/DEGS_gseGOBP_fetal_adult_dotplot.pdf' )


sort(unique(top_sub_fetal_adult_subcelltype_gseGO_df$Description))
temp_df=top_sub_fetal_adult_subcelltype_gseGO_df
gene_sets= list('positive regulation of lymphocyte activation'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='positive regulation of lymphocyte activation'],'/')),
                'antigen processing and presentation of peptide antigen'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='antigen processing and presentation of peptide antigen'],'/')),
                'adaptive immune response'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='adaptive immune response'],'/')),
                'leukocyte mediated cytotoxicity'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='leukocyte mediated cytotoxicity'],'/')))

filt_NBRC_altas_seu <- AddModuleScore_UCell(filt_NBRC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快
saveRDS(filt_NBRC_altas_seu@meta.data,file ='20251125_filt_NBRC_altas_seu_meta2.rds' )

VlnPlot(subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM','ABM')),group.by = 'source_celltype',stack = T,
        features = paste0(names(gene_sets),'_UCell'))



#------------------------分组单独比较-----------------------------------#
gc()
group='FL_FBM'
sfile='Protein_NRBC_marker/DE_marker/FL_FBM_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM')),group = group,sfile = sfile)
FL_FBM_all_Ery_tissue_markers=res[[1]]
FL_FBM_count_df=res[[2]]
rm(res);gc()
# FL 与FBM吧nRBC 更为相似可能是导致FL nRBC 相较于BM nRBC无显著差异基因的原因 

group='FL_ABM'
sfile='Protein_NRBC_marker/DE_marker/FL_ABM_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','ABM')),group = group,sfile = sfile)
FL_ABM_all_Ery_tissue_markers=res[[1]]
FL_ABM_count_df=res[[2]]
rm(res);gc()


group='FBM_ABM'
sfile='Protein_NRBC_marker/DE_marker/ABM_FBM_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = subset(filt_NBRC_altas_seu,tissue_stage %in% c('FBM','ABM')),group = group,sfile = sfile)
FBM_ABM_all_Ery_tissue_markers=res[[1]]
FBM_ABM_count_df=res[[2]]
rm(res);gc()

# FBM中组蛋白表达太高，影响富集结果 
FL_FBM_all_Ery_tissue_markers=FL_FBM_all_Ery_tissue_markers[!FL_FBM_all_Ery_tissue_markers$gene %in% names(HIS_genes),]
FBM_ABM_all_Ery_tissue_markers=FBM_ABM_all_Ery_tissue_markers[!FBM_ABM_all_Ery_tissue_markers$gene %in% names(HIS_genes),]

table(FL_FBM_all_Ery_tissue_markers[FL_FBM_all_Ery_tissue_markers$avg_log2FC >1 & FL_FBM_all_Ery_tissue_markers$pct.1 >0.1 & FL_FBM_all_Ery_tissue_markers$pct.2 <0.3,c('cluster','celltype')])
table(FL_ABM_all_Ery_tissue_markers[FL_ABM_all_Ery_tissue_markers$avg_log2FC >1 & FL_ABM_all_Ery_tissue_markers$pct.1 >0.1 & FL_ABM_all_Ery_tissue_markers$pct.2 <0.3,c('cluster','celltype')])
table(FBM_ABM_all_Ery_tissue_markers[FBM_ABM_all_Ery_tissue_markers$avg_log2FC >1 & FBM_ABM_all_Ery_tissue_markers$pct.1 >0.1 & FBM_ABM_all_Ery_tissue_markers$pct.2 <0.3,c('cluster','celltype')])


FL_FBM_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =FL_FBM_all_Ery_tissue_markers[FL_FBM_all_Ery_tissue_markers$cluster=='FL',],keyType = 'ALIAS' )
FL_ABM_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =FL_ABM_all_Ery_tissue_markers[FL_FBM_all_Ery_tissue_markers$cluster=='FL',],keyType = 'ALIAS' )
FL_FBM_subcelltype_gseGO_list[[2]]$group='FL_FBM'
FL_ABM_subcelltype_gseGO_list[[2]]$group='FL_ABM'
degs_gseGO_res_df2=data.frame()
degs_gseGO_res_df2=rbind(degs_gseGO_res_df2,rbind(FL_FBM_subcelltype_gseGO_list[[2]],FL_ABM_subcelltype_gseGO_list[[2]]))

ABM_FBM_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =FBM_ABM_all_Ery_tissue_markers[FBM_ABM_all_Ery_tissue_markers$cluster=='FBM',] ,keyType = 'ALIAS' )
ABM_FBM_subcelltype_gseGO_list[[2]]$group='FBM_ABM'
degs_gseGO_res_df2=rbind(degs_gseGO_res_df2,ABM_FBM_subcelltype_gseGO_list[[2]])
#write.csv(degs_gseGO_res_df2,file = 'Protein_NRBC_marker/res_data/main_figure1/tissue_nRBC_degs_gseGO_res_df2.csv') # symbol
write.csv(degs_gseGO_res_df2,file = 'Protein_NRBC_marker/res_data/main_figure3/HSPC_derived_nRBC_degs_gseGO_res_df2_ALIAS.csv')

degs_gseGO_res_df2=degs_gseGO_res_df2[degs_gseGO_res_df2$ONTOLOGY=='BP',]   
degs_gseGO_res_df2$group=factor(degs_gseGO_res_df2$group,levels = c('FL_FBM','FL_ABM','FBM_ABM'))
degs_gseGO_res_df2$res='up'
degs_gseGO_res_df2$res[degs_gseGO_res_df2$NES <0]='down'

top_degs_gseGO_bp_res_df=degs_gseGO_res_df2[degs_gseGO_res_df2$p.adjust <0.01,] %>% group_by(group,celltype,res) %>%top_n(wt = -log10(p.adjust),n=10)  %>% do(head(.,10))
head(sort(table(top_degs_gseGO_bp_res_df$Description),decreasing = T),40)
head(top_degs_gseGO_bp_res_df[,c(1:10,13:14)])
top_degs_gseGO_bp_res_df$Description=factor(top_degs_gseGO_bp_res_df$Description,levels = unique(top_degs_gseGO_bp_res_df$Description))
top_degs_gseGO_bp_res_df$celltype=factor(top_degs_gseGO_bp_res_df$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))


temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group %in% c('FL_FBM'),]
temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
p2=ggplot(temp_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+facet_grid(~group)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'FL vs FBM nRBC gseGO of DEGs')+theme(text = element_text(face = 'bold'),panel.grid = element_blank())
p2

temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group %in% c('FL_ABM'),]
temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
p3=ggplot(temp_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+facet_grid(~group)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'FL vs ABM nRBC gseGO of DEGs')+theme(text = element_text(face = 'bold'),panel.grid = element_blank())
p3

temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group=='FBM_ABM',]
temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
p4=ggplot(temp_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face = 'bold'))+ggtitle(label = 'FBM vs ABM nRBC gseGO of DEGs')
p4

p=p2+p3+p4;p
ggsave(p,width =22 ,height =12,filename='../Protein_NRBC_marker/res_pic/main_figure3/DEGS_gseGOBP_substate_FL_FBM_ABM_dotplot.pdf' )

###################### 分析通路中的核心驱动基因的表达，发现血管发育核心驱动基因以及凝血调控核心驱动基因中大部分（103/130，75/95）在FL 痕量表达（pct <0.1）, 
#剩余基因基因中除了NR4A1、GPX1和SERPINEA1，其他基因在FL 与FBM NRBC 的表达模式基本一致#
# vessel development
#  查看core gene expression 发现FL和FBM中表达模式类同，无明显差异 
if(F){
  
  temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group %in% c('FL_FBM'),]
  temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
  temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
  coree_enrichment_genesets=unique(unlist(strsplit(temp_df$core_enrichment[grep('vas|vesssel|angio',temp_df$Description)],split = '/')))
  coree_enrichment_genesets_df=FL_FBM_all_Ery_tissue_markers[FL_FBM_all_Ery_tissue_markers$gene %in%  coree_enrichment_genesets,]
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[coree_enrichment_genesets_df$cluster=='FL',]
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[order(coree_enrichment_genesets_df$avg_log2FC,-log10(coree_enrichment_genesets_df$p_val_adj),decreasing = T),]
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[abs(coree_enrichment_genesets_df$avg_log2FC) >1,]
  coree_enrichment_genesets_df$lg.p_val_adj=-log10(coree_enrichment_genesets_df$p_val_adj);coree_enrichment_genesets_df$lg.p_val_adj[coree_enrichment_genesets_df$lg.p_val_adj >200]=200
  length(unique(coree_enrichment_genesets_df$gene))# 130
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[coree_enrichment_genesets_df$pct.1 >0.1,]
  length(unique(coree_enrichment_genesets_df$gene)) # 27
  
  label_gene=unique(c(unique(coree_enrichment_genesets_df$gene[coree_enrichment_genesets_df$lg.p_val_adj>100]),unique(coree_enrichment_genesets_df$gene[coree_enrichment_genesets_df$avg_log2FC>5])))
  label_gene_df=coree_enrichment_genesets_df[coree_enrichment_genesets_df$gene %in% label_gene, ]
  label_gene_df=label_gene_df[label_gene_df$avg_log2FC >5| label_gene_df$lg.p_val_adj >100,]
  label_gene_df=label_gene_df[!duplicated(label_gene_df$gene),]
  
  p=ggplot(coree_enrichment_genesets_df[!duplicated(coree_enrichment_genesets_df$gene),],aes(x=avg_log2FC,y=lg.p_val_adj,size=pct.1,color=avg_log2FC,shape=celltype))+geom_point()+theme_classic()+
    geom_label(aes(label=gene,x=avg_log2FC+0.5,y=lg.p_val_adj-1),alpha=0.6)+ggtitle('core gene of vessel development')
  p
  
  p/DotPlot(filt_NBRC_altas_seu,features =unique(coree_enrichment_genesets_df$gene),group.by = 'source_celltype' ,scale = F)+RotatedAxis()
  
  # coagulation
  coree_enrichment_genesets=unique(unlist(strsplit(temp_df$core_enrichment[grep('coagulation|wound',temp_df$Description)],split = '/')))
  coree_enrichment_genesets_df=FL_FBM_all_Ery_tissue_markers[FL_FBM_all_Ery_tissue_markers$gene %in%  coree_enrichment_genesets,]
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[coree_enrichment_genesets_df$cluster=='FL',]
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[order(coree_enrichment_genesets_df$avg_log2FC,-log10(coree_enrichment_genesets_df$p_val_adj),decreasing = T),]
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[abs(coree_enrichment_genesets_df$avg_log2FC) >1,]
  coree_enrichment_genesets_df$lg.p_val_adj=-log10(coree_enrichment_genesets_df$p_val_adj);coree_enrichment_genesets_df$lg.p_val_adj[coree_enrichment_genesets_df$lg.p_val_adj >200]=200
  length(unique(coree_enrichment_genesets_df$gene))# 95
  coree_enrichment_genesets_df=coree_enrichment_genesets_df[coree_enrichment_genesets_df$pct.1 >0.1,]
  length(unique(coree_enrichment_genesets_df$gene)) # 20
  
  label_gene=unique(c(unique(coree_enrichment_genesets_df$gene[coree_enrichment_genesets_df$lg.p_val_adj>100]),unique(coree_enrichment_genesets_df$gene[coree_enrichment_genesets_df$avg_log2FC>5])))
  label_gene_df=coree_enrichment_genesets_df[coree_enrichment_genesets_df$gene %in% label_gene, ]
  label_gene_df=label_gene_df[label_gene_df$avg_log2FC >5| label_gene_df$lg.p_val_adj >100,]
  label_gene_df=label_gene_df[!duplicated(label_gene_df$gene),]
  
  p=ggplot(coree_enrichment_genesets_df[!duplicated(coree_enrichment_genesets_df$gene),],aes(x=avg_log2FC,y=lg.p_val_adj,size=pct.1,color=avg_log2FC,shape=celltype))+geom_point()+theme_classic()+
    geom_label(aes(label=gene,x=avg_log2FC+0.5,y=lg.p_val_adj-1),alpha=0.6)+ggtitle('core gene of coagulation')
  
  
  p/DotPlot(filt_NBRC_altas_seu,features =unique(coree_enrichment_genesets_df$gene),group.by = 'source_celltype' ,scale = F)+RotatedAxis()
  
}



temp_df=top_degs_gseGO_bp_res_df[top_degs_gseGO_bp_res_df$group %in% c('FL_FBM'),]
temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
filt_NBRC_altas_seu=AddModuleScore_UCell(filt_NBRC_altas_seu,features = list(angiogenesis=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='angiogenesis'],'/')),
                                                                             'blood_vessel_morphogenesis'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='blood vessel morphogenesis'],'/')),
                                                                             'endothelial_cell_migration'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='endothelial cell migration'],'/')),
                                                                             'coagulation'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='coagulation'],'/')),
                                                                             ncores = 6))



saveRDS(filt_NBRC_altas_seu@meta.data[,c('angiogenesis_UCell','blood_vessel_morphogenesis_UCell','endothelial_cell_migration_UCell','coagulation_UCell')],file = 'NRBC_angiogenesis_pathway.rds')

p=VlnPlot(subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM','ABM')),group.by = 'source_celltype',stack = T,cols = cols,
        features = c('angiogenesis_UCell','blood_vessel_morphogenesis_UCell','endothelial_cell_migration_UCell','coagulation_UCell', paste0(names(gene_sets),'_UCell')))

ggsave(p,filename = 'Protein_NRBC_marker/res_pic/main_figure3/key_pathway_UCell_score_definitive_nRBC_vlnplot.pdf',width = 10,height = 10)


#---------------------------the ratio of cell cycle phase--------------------------------# 
phase_stats_by_sample <- filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage!='YS',c('id','Phase','final_celltype','tissue_stage')] %>%
  group_by(tissue_stage, id, final_celltype, Phase) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue_stage, id, final_celltype) %>%
  mutate(
    total = sum(count),
    proportion = count / total * 100,
    proportion_label = sprintf("%.1f%%", proportion)
  ) %>%
  ungroup()
# 选则样本中至少10个NRBC 存在
p_final <- ggplot(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='G1',], aes(x = tissue_stage, y = proportion, fill = tissue_stage)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1,alpha=0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +facet_wrap(~final_celltype, scales = "free_y", nrow = 1) +theme_classic() +
  scale_fill_manual(values = cols, name = "Tissue Stage") +labs(title = "Cell Cycle G1 Phase Distribution in subcelltype",x = "",y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),strip.text = element_text(size = 11, face = "bold")) +ylim(c(0,115))+
  stat_compare_means(comparisons = list(c("FL", "FBM"),c("FL", "ABM"), c("FBM", "ABM")),method = "wilcox.test",label = "p.signif" )

print(p_final)

# 计算样本数目
for( celltype in unique(phase_stats_by_sample$final_celltype)){
  print(celltype)
  print(table(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='G1' &phase_stats_by_sample$final_celltype==celltype ,'tissue_stage']))
}

ggsave(p_final,filename = 'Protein_NRBC_marker/res_pic/main_figure3/definitive_G1_phase_wilcox_test.pdf',height = 6,width = 18)

p_final <- ggplot(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='S',], aes(x = tissue_stage, y = proportion, fill = tissue_stage)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1,alpha=0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +facet_wrap(~final_celltype, scales = "free_y", nrow = 1) +theme_classic() +
  scale_fill_manual(values = cols, name = "Tissue Stage") +labs(title = "Cell Cycle S Phase Distribution in subcelltype",x = "",y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),strip.text = element_text(size = 11, face = "bold")) +ylim(c(0,100))+
  stat_compare_means(comparisons = list(c("FL", "FBM"),c("FL", "ABM"), c("FBM", "ABM")),method = "wilcox.test",label = "p.signif" )

print(p_final)

ggsave(p_final,filename = 'Protein_NRBC_marker/res_pic/main_figure3/definitive_S_phase_wilcox_test.pdf',height = 6,width = 18)


p_final <- ggplot(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='G2M',], aes(x = tissue_stage, y = proportion, fill = tissue_stage)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1,alpha=0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +facet_wrap(~final_celltype, scales = "free_y", nrow = 1) +theme_classic() +
  scale_fill_manual(values = cols, name = "Tissue Stage") +labs(title = "Cell Cycle G2M Phase Distribution in subcelltype",x = "",y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),strip.text = element_text(size = 11, face = "bold")) +ylim(c(0,115))+
  stat_compare_means(comparisons = list(c("FL", "FBM"),c("FL", "ABM"), c("FBM", "ABM")),method = "wilcox.test",label = "p.signif" )

print(p_final)

ggsave(p_final,filename = 'Protein_NRBC_marker/res_pic/main_figure3/definitive_G2M_phase_wilcox_test.pdf',height = 6,width = 18)

Idents(filt_NBRC_altas_seu)='source_celltype'
VlnPlot(filt_NBRC_altas_seu,features = c('CCNB1','CDK1'),pt.size = 0,split.by = 'tissue_stage',stack = T)


######################################################################################################################################################
#----------------bulk RNAseq validated the marker expression
######################################################################################################################################################
#-----------------------体外期脐带血Ery1：Ery4，Ery5：成人外周血
defintive_markers=readRDS('Protein_NRBC_marker/res_data/main_figure3/defintive_markers.rds')
load('/home/gibh/2021_NRBC_chlyu/ref_data/hema_ref_bullk_RNAseq_se.Rdata',verbose = T)
Ery_annotation_df=data.frame(hema.se@colData[grep('Ery',hema.se@colData$celltype),])
Ery_count_df=data.frame(hema.se@assays@data$logcounts[,rownames(Ery_annotation_df)])
Ery_count_df[as.character(unlist(defintive_markers))[ !as.character(unlist(defintive_markers))%in%  rownames(Ery_count_df)],colnames(Ery_count_df)]=0

symbol_an_df=data.frame(row.names = as.character(unlist(defintive_markers)),type=rep(names(sapply(defintive_markers, length)),as.numeric(sapply(defintive_markers, length))))
p=pheatmap(Ery_count_df[as.character(unlist(defintive_markers)),],cluster_rows = F,cluster_cols = F,show_colnames = F,annotation_row = symbol_an_df,
         color = colorRampPalette(colors = c('white','firebrick3'))(100),annotation_col = Ery_annotation_df[,c('label.fine','celltype')])
ggsave(as.ggplot(p),filename = 'Protein_NRBC_marker/res_pic/main_figure3/definitive_marker_expression_UCB_PBMC_bulkRNAseq_heatmap.pdf',width = 8,height = 8)

FL_ABM_nRBC_df_refer_control_df=readRDS('ref_data/bulk_RNAseq/nr_FL_ABM_nRBC_ref_ACTB.rds')
FL_ABM_nRBC_df_df=readRDS('ref_data/bulk_RNAseq/FL_ABM_nRBC_df_df_nr_exp.rds')

p=pheatmap(FL_ABM_nRBC_df_refer_control_df[c('ACTB', as.character(unlist(defintive_markers))),-grep('MM|thy_ery',colnames(FL_ABM_nRBC_df_refer_control_df))],cluster_rows = F,cluster_cols = F,annotation_row = symbol_an_df,
         color = colorRampPalette(colors = c('white','firebrick3'))(100))
ggsave(as.ggplot(p),filename = 'Protein_NRBC_marker/res_pic/main_figure3/definitive_marker_expression_FL_ABM_bulkRNAseq_heatmap.pdf',width = 8,height = 8)

MHC_genes=c("B2M","HLA-A","HLA-B","HLA-C","HLA-E","CD74","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1", "HLA-DQB1" ,"HLA-DQB2" ,"HLA-DRA" ,"HLA-DRB1","HLA-DRB4", "HLA-DRB5","HLA-DRB6","TAPBP") # order, and cho
peptide_load=c('TAPBP','TAP1','TAP2','PDIA3','ERAP1')
Proteasome_genes=rownames(Ery_count_df)[grep('PSMB',rownames(Ery_count_df))]
Proteasome_genes
Proteasome_genes=c(Proteasome_genes[1],Proteasome_genes[3:10],Proteasome_genes[2])
# CIITA（MHC-II和共刺激分子主调控因子）
co_stimulatory =c('CIITA','CD86','CD80','ICOSLG','CD40','OX40L', 'TNFSF4', 'TNFSF9', 'TNFSF14', 'TNFSF18')
co_inhibtor=c('CD274', 'PDCD1LG2', 'CD276', 'VTCN1', 'HHLA2', 'IDO1')

all_antigen_process_genes=unique(c(MHC_genes,peptide_load,Proteasome_genes,co_stimulatory,co_inhibtor))
p=pheatmap(Ery_count_df[all_antigen_process_genes[all_antigen_process_genes %in% rownames(Ery_count_df)],],
         cluster_rows = F,cluster_cols = F,show_colnames = F,annotation_col = Ery_annotation_df[,c('label.fine','celltype')],color = colorRampPalette(colors = c('white','firebrick3'))(100))
ggsave(as.ggplot(p),filename = 'Protein_NRBC_marker/res_pic/main_figure3/MHC_peptide_load_gene_expression_UCB_bulkRNAseq_heatmap.pdf',width = 8,height = 8)

p=pheatmap(FL_ABM_nRBC_df_refer_control_df[c('ACTB',all_antigen_process_genes[all_antigen_process_genes %in% rownames(FL_ABM_nRBC_df_refer_control_df)]),-grep('MM|thy_ery',colnames(FL_ABM_nRBC_df_refer_control_df))],
           cluster_rows = F,cluster_cols = F,show_colnames = T,color = colorRampPalette(colors = c('white','firebrick3'))(100))

ggsave(as.ggplot(p),filename = 'Protein_NRBC_marker/res_pic/main_figure3/MHC_peptide_load_gene_expression_FL_ABM_bulkRNAseq_heatmap.pdf',width = 8,height = 8)





# ---------------------查看之前购买marker的表达情况-------------------------#
if(F){
  buy_antiobody_genes=c('CD81','CD109','DLK1','PDZK1IP1','CXCR4','IGF1R','ALCAM','IL3RA','ITGAV','SEMA4D') #CXCR4,CD221,CD166,CD123,CD51,CD100
  DotPlot(filt_NBRC_altas_seu,features =buy_antiobody_genes,scale = F,group.by = 'source_celltype')+RotatedAxis()+scale_color_gradient2(low = 'blue',mid = 'gray',high = 'firebrick3')
  
  ABM_altas_seu=readRDS('NRBC_BM_altas/tmp_ABM_altas_seu_new.rds')
  ABM_altas_NRBC_seu=subset(ABM_altas_seu,subcelltype %in% c('Bas','BFUE/CFUE','Orth','Poly','ProE'))
  ABM_altas_NRBC_seu$subcelltype=factor(ABM_altas_NRBC_seu$subcelltype,levels = c('BFUE/CFUE','ProE','Bas','Poly','Orth'))
  subset(ABM_altas_NRBC_seu,subset=IGKC >1) #3737 
  
  table(subset(ABM_altas_NRBC_seu,subset=IGKC >1)$subcelltype)# 主要是Bas和BFUE/CFUE
  #Bas BFUE/CFUE      Orth      Poly      ProE 
  #2068      1162       152        57       298 
  
  ABM_altas_NRBC_seu$IGHC='negtive'
  ABM_altas_NRBC_seu$IGHC[rownames(ABM_altas_NRBC_seu@meta.data) %in% colnames(subset(ABM_altas_NRBC_seu,subset=IGKC >1))]='positive'
  VlnPlot(ABM_altas_NRBC_seu,features =buy_antiobody_genes,group.by = 'subcelltype',split.by = 'IGHC',stack = T)
  
  
  
  pheatmap(GSE301441_BM_NRBC_df[buy_antiobody_genes[buy_antiobody_genes %in% rownames(GSE301441_BM_NRBC_df)],],cluster_cols = F,cluster_rows = F,main = 'ABM')
  pheatmap(FL_primary_NRBC_df[buy_antiobody_genes[buy_antiobody_genes %in% rownames(FL_primary_NRBC_df) ],],cluster_cols = F,cluster_rows = F,main = 'FL')
  
  pheatmap(FL_ABM_nRBC_df_refer_control_df[c(buy_antiobody_genes[buy_antiobody_genes %in% rownames(FL_ABM_nRBC_df_refer_control_df)] ,'ACTB'),], 
           labels_row = c("CD81","CD109","DLK1","PDZK1IP1",'CXCR4','CD221','CD166','CD123','CD51','CD100','ACTB'),
           cluster_cols = F,cluster_rows = F,main = 'FL/ABM refer_ACTB_as_control')
  
  
  
}






