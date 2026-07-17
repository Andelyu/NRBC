
############################################################################################################################### 
#------------------------------------YS Ery altas data analysis, version2 ------------------------------------
############################################################################################################################### 
source('/home/gibh/2021_NRBC_chlyu/zx_lab_NRBC/scripts/scRNAseq_pipline/scRNAseq_analysis_model.R')

library(SeuratDisk)
library(anndata)
library(Seurat)
library(RColorBrewer)
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
setwd('/home/gibh/2021_NRBC_chlyu/')

#----------------------singleR reference bulk RNAseq----------------------#
if(F){
  #----------------------- refer gene symbol ------------------------#
  all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  
  all_shared_ensembl_id_info$org_symbol=as.character(mapIds(x = org.Hs.eg.db,keys = all_shared_ensembl_id_info$X,keytype = 'ENSEMBL',column = 'SYMBOL'))
  all_shared_ensembl_id_info$org_symbol[is.na(all_shared_ensembl_id_info$org_symbol)]=all_shared_ensembl_id_info$ref_symbol[is.na(all_shared_ensembl_id_info$org_symbol)]
  table(all_shared_ensembl_id_info$ref_symbol==all_shared_ensembl_id_info$org_symbol) # F:2482,T:33248 /23824
  
  all_shared_ensembl_id_info$ref_symbol1=all_shared_ensembl_id_info$fl_symbol
  all_shared_ensembl_id_info$ref_symbol1[is.na(all_shared_ensembl_id_info$ref_symbol1)]=all_shared_ensembl_id_info$org_symbol[is.na(all_shared_ensembl_id_info$ref_symbol1)]
  table(is.na(all_shared_ensembl_id_info$ref_symbol1))
  table(all_shared_ensembl_id_info$ref_symbol1==all_shared_ensembl_id_info$org_symbol) # F:11906 ,T: 23824
  rownames(all_shared_ensembl_id_info)=all_shared_ensembl_id_info$X
  
  write.csv(all_shared_ensembl_id_info,file = 'ref_data/all_shared_ensembl_id_info.csv')
  
  # ---------------------prepare refer bulk data data---------------#
  load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T)  
  table(rownames(YS_altas_Ery_seu@assays$RNA) %in% rownames(nrbc_ref_se2@assays@data$logcounts)) #FALSE/TRUE :12958 16594 
  
  logcounts=nrbc_ref_se2@assays@data$logcounts
  symbol_ToEn=as.character(mapIds(x = org.Hs.eg.db,keys = rownames(nrbc_ref_se2@assays@data$logcounts),keytype = 'ALIAS',column = 'ENSEMBL' ))
  table(symbol_ToEn %in% all_shared_ensembl_id_info$X) # F/T: 3487 18407
  symbol_ToEn=all_shared_ensembl_id_info[match(x =as.character(symbol_ToEn) ,table = all_shared_ensembl_id_info$X),'ref_symbol1']
  table(is.na(symbol_ToEn))# TRUE:3487, FALSE:18407
  symbol_ToEn[is.na(symbol_ToEn)]=rownames(logcounts)[is.na(symbol_ToEn)]
  du_num=which(symbol_ToEn %in% symbol_ToEn[duplicated(symbol_ToEn)] )
  du_logcounts=logcounts[du_num,]
  du_logcounts=aggregate(du_logcounts,list(symbol_ToEn[du_num]),sum);du_logcounts=column_to_rownames(du_logcounts,'Group.1')
  logcounts=rbind(logcounts[-du_num,],du_logcounts)
  
  pheatmap(cor(logcounts))
  nrbc_ref_se2=SummarizedExperiment(assays =list(logcounts=logcounts),colData = data.frame(nrbc_ref_se2@colData) )
  
  
  logcounts=nrbc_ref_se2@assays@data$logcounts[,-1:-4]
  pheatmap(cor(logcounts))
  nrbc_ref_se=SummarizedExperiment(assays =list(logcounts=logcounts),colData = data.frame(nrbc_ref_se@colData) )
  
  rm(logcounts);gc()
  
  logcounts=nrbc_ref_se2@assays@data$logcounts[,-grep('HSPC|CD34',colnames(nrbc_ref_se2@assays@data$logcounts))]
  colnames(logcounts)=c("BFUE","CFUE",colnames(logcounts)[3:length(colnames(logcounts))])
  colData=data.frame(celltype=as.character(t(data.frame(strsplit(colnames(logcounts),split = '_')))[,1]),
                     sourceis=c(rep('GSE128268',3),rep('GSE53983',15),rep('GSE61566',6)))
  rownames(colData)= colnames(logcounts)
  
  colnames(logcounts)=c(colnames(logcounts)[1:2],'ProE',colnames(logcounts)[4:length(colnames(logcounts))])
  annotation_df=data.frame(row.names =  colnames(logcounts),sourceid=c(rep('GSE128268',3),rep('GSE53983',15),rep('GSE61566',6)))
  p=pheatmap(annotation_row = annotation_df,cor(logcounts),color = colorRampPalette(colors = c('navy','white','firebrick3'))(100),main = 'singleR:refence bulkRNAseq',cluster_cols = F)
  
  # 剔除HSPC
  nrbc_ref_se2=SummarizedExperiment(assays =list(logcounts=logcounts),colData =colData )

  logcounts=nrbc_ref_se@assays@data$logcounts[,-grep('HSPC|CD34',colnames(nrbc_ref_se@assays@data$logcounts))]
  nrbc_ref_se=SummarizedExperiment(assays =list(logcounts=logcounts),colData =nrbc_ref_se@colData[-grep('HSPC|CD34',rownames(nrbc_ref_se@colData)),] )# 去掉HSCP细胞
  colnames(nrbc_ref_se@colData)='celltype'
  
  save(nrbc_ref_se,nrbc_ref_se2,file = '/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata')
  
  library(ggplotify)
  ggsave(as.ggplot(p),filename='ref_data/NRBC_ref_bullk_RNAseq_cor_heatmap.pdf',width = 10,height = 8)
  
}else{
  load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T)  
  all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  
}




############################################################################################################################### 
#-------------------------------------prepare  data------------------------------------
############################################################################################################################### 

# ----------------------------------- YS  data information ------------------------------#
# 多个数据来源，导致gene symbol并没有统一，需要统一注释不同来源gene symbol
# 这个样本数据也是集合了多篇文章的数据，包括： 
# Popescu_et_al_Nature_2019:Decoding human fetal liver haematopoiesis 的YS：F32,F35(CS23(8 PWC)，74 个Ery,太少，进行细胞比例的时候删掉了), F37，数据在验证集中已经使用
# Wang_et_al_Cell_Stem Cell_2021 ：Decoding Human Megakaryocyte Development， Wang_CS10, Wang_CS11， 验证集中已经使用，
# HM_YS：MAPPING HUMAN HAEMATOPOIETIC STEM CELLS FROM HAEMOGENIC ENDOTHELIUM TO BIRTH ，2022 nature, :mikola_1，只有一个样本:CS14
# WE_YS：F138:CS17，一个样本,
# Mapping the developing human immune system across organs*,Science,"YS_F37_4PCW"  "LI_F61_7PCW"  "YS_F61_7PCW"  "LI_F19_10PCW"；YS_F37_4PCW在上述数据中已经在验证集被使用
# De_novo：F61,F79,F80，三个样本，CS18，CS18，CS15,其中F61在 接下来的数据中已经被使用
# 总结：F79,F80，F138，mikola_1四个样本没有使用过
# YS-FL match sample id: F32,F35,F61 三个样本
# F开头的样本中，同一个样本，有不同区域来源的NRBC， 需要再次优化

if(F){
  #------------------------这个YS数据，来自多个文章-----------------#
  Convert(source ='NRBC_YS_altas/raw_ref_data/ys_portal_object_20221208.h5ad', dest="h5seurat" )
  YS_altas_seu=LoadH5Seurat('NRBC_YS_altas/raw_ref_data/ys_portal_object_20221208.h5seurat', meta.data = FALSE)
  meta_df=read.csv('NRBC_YS_altas/raw_ref_data/ys_portal_object_20221208_meta.csv',header = T)
  # 查定，修改F138-->F158
  rownames(meta_df)=meta_df$X;meta_df=meta_df[,-1]
  meta_df$id=paste(meta_df$stage,meta_df$fetal.ids,sep = '_')
  
  YS_altas_seu@meta.data=meta_df
  
  table(YS_altas_seu$LVL1);table(YS_altas_seu$LVL2);table(YS_altas_seu$LVL3)
  t(table(YS_altas_seu@meta.data[,c('LVL1','LVL2')]))
  
  
  # M. Haniffa, Decoding human fetal liver haematopoiesis. Nature574, 365–371 (2019).Dorin-Mirel Popescu 1，E-MTAB-7407
  # M. Haniffa, S. A. Teichmann, Mapping the developing human immune system across organs. Science376, eabo0510 (2022).Chenqu Suo 1，E-MTAB-11343
  # H. K. A. Mikkola, Mapping human haematopoietic stem cells from haemogenic endothelium to birth. Nature604, 534–540 (2022).Vincenzo Calvanese 1， GSE162950. Data from published：GSE135202(无feature 信息).
  # Decoding Human Megakaryocyte Development. Cell Stem Cell28, 535–549.e8 (2021). GSE144024
  # The role of the yolk sac in human fetal development and identification of a hepatocyte-like cell in the human yolk sac: E-MTAB-10552, Michael Mather 1Muzz Haniffa 
  # S. Webb, M. Haniffa, E. Stephenson, “Human fetal yolk sac scRNA-seq data (sample ID: F158 for Haniffa Lab; 16099 for HDBR),” BioStudies, E-MTAB-11673 (2022);
  # We :(CellRanger version 6.0.1 using human reference genome GRCh38-2020-A ,https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11673
  
  # Smart2-seq M. Haniffa, M. Mather, R. Botting, “The role of the yolk sac in human fetal development and identification of a hepatocyte-like cell in the human yolk sac (SS2),” BioStudies, E-MTAB-10888 (2023)
  
  c('GNB2L1','ATP5E','RACK1','ATP5F1E') %in% rownames(YS_altas_seu)
  object_list=SplitObject(YS_altas_seu,split.by = 'orig.dataset')
  object_list[[1]]=CreateSeuratObject(counts = GetAssayData(object_list[[1]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data =object_list[[1]]@meta.data );dim(object_list[[1]])
  object_list[[2]]=CreateSeuratObject(counts = GetAssayData(object_list[[2]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[2]]@meta.data);dim(object_list[[2]])
  object_list[[3]]=CreateSeuratObject(counts = GetAssayData(object_list[[3]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[3]]@meta.data);dim(object_list[[3]])
  object_list[[4]]=CreateSeuratObject(counts = GetAssayData(object_list[[4]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[4]]@meta.data);dim(object_list[[4]])
  object_list[[5]]=CreateSeuratObject(counts = GetAssayData(object_list[[5]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[5]]@meta.data);dim(object_list[[5]])
  
  GSE162950_feature_info=read.table('NRBC_YS_altas/raw_ref_data/GSE162950_features.tsv',sep="\t")
  YS_feature_info=read.table('NRBC_YS_altas/YS_10x_genes.tsv',sep="\t",header = T)
  
  We_feature_info=read.table('NRBC_YS_altas/raw_ref_data/gene_info.txt',sep="\t",header = F) # cellranger_genome_2020A
  We_feature_info=We_feature_info[,-1];We_feature_info$V2=gsub(pattern = ' ',replacement = '',We_feature_info$V2)
  table(We_feature_info$V2 %in% all_shared_ensembl_id_info$X)# F/T:3695/32905 
  table(We_feature_info[!We_feature_info$V2 %in% all_shared_ensembl_id_info$X,'V3'])
  #lncRNA  protein_coding
  # 3621              74 
  
  names(object_list)
  #  Wang_et_al_Cell_Stem_Cell_2021:zx_lab,Popescu_et_al_Nature_2019 & De_novo: FL , Popescu_跟FL数据是同一来源, HM_YS、 WE_YS注释信息各不同, 需要找到原文数据
  # 以 Popescu_et_al_Nature_2019 数据注释为参考，统一全文，需要重新注解
  
  # 统计不同来源数据注释基因差异
  # 发现CS10_Wang_CS10、CS11_Wang_CS11、CS14_mikola_1、CS17_F158: Wang_et_al_Cell_Stem_Cell_2021: zx_lab,
  # Popescu_et_al_Nature_2019+De_novo+WE_YS:gene feature 注释信息比较一致，F*样本的gene feature 注释信息比较一致，
  ## 5ad文件无feature信息, Wang_et_al_Cell_Stem: GSE144024
  
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$abm_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$zxlab_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$fl_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$ref_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$ref_symbol1)})
  sapply(object_list, function(x){table(rownames(x) %in% GSE162950_feature_info$V2)}) # 同Wang_et_al_Cell_Stem
  sapply(object_list, function(x){table(rownames(x) %in% YS_feature_info$X)}) # WE_YS:2174
  sapply(object_list, function(x){table(rownames(x) %in% We_feature_info$V4)}) # WE_YS:F/T:9/23814
  
  
  new_symbol_list=list()
  refId="Wang_et_al_Cell_Stem Cell_2021"
  new_symbol_list[[refId]]=all_shared_ensembl_id_info[match( x = rownames(object_list[[refId]]),table = all_shared_ensembl_id_info$zxlab_symbol),'ref_symbol1']
  out_re_symbols=rownames(object_list[[refId]])[! rownames(object_list[[refId]]) %in% all_shared_ensembl_id_info$zxlab_symbol ];out_re_symbols
  out_re_symbols=as.character(t(data.frame(strsplit(out_re_symbols,split = '.',fixed = T)))[,1]);out_re_symbols
  new_symbol_list[[refId]][! rownames(object_list[[refId]]) %in% all_shared_ensembl_id_info$zxlab_symbol ]=out_re_symbols
  
  
  refId="WE_YS"
  new_symbol_list[[refId]]=all_shared_ensembl_id_info[We_feature_info$V2[match( x = rownames(object_list[[refId]]),table = We_feature_info$V4)],'ref_symbol1']
  table(is.na(new_symbol_list[[refId]])) # F/T:22217/1606, lncRNA 占绝大多数，可能跟其测序方法有关
  new_symbol_list[[refId]][is.na(new_symbol_list[[refId]])]=rownames(object_list[[refId]])[is.na(new_symbol_list[[refId]])] 
  out_re_symbols=rownames(object_list[[refId]])[! rownames(object_list[[refId]]) %in% We_feature_info$V4 ];out_re_symbols
  out_re_symbols=as.character(t(data.frame(strsplit(out_re_symbols,split = '-1',fixed = T)))[,1]);out_re_symbols
  out_re_symbols[!out_re_symbols %in% all_shared_ensembl_id_info$ref_symbol1] # 3 不在
  trans_out_re_symbols=all_shared_ensembl_id_info[We_feature_info$V2[match( x = out_re_symbols,table = We_feature_info$V4)],'ref_symbol1']
  trans_out_re_symbols[is.na(trans_out_re_symbols)]=out_re_symbols[is.na(trans_out_re_symbols)]
  new_symbol_list[[refId]][! rownames(object_list[[refId]]) %in% We_feature_info$V4 ]=trans_out_re_symbols
  
  HM_to_refsymbol=mapIds(x = org.Hs.eg.db,keys =rownames( object_list[['HM_YS']]),keytype = 'ALIAS',column = 'ENSEMBL')
  HM_to_refsymbol=all_shared_ensembl_id_info[match(x = as.character(HM_to_refsymbol),table =all_shared_ensembl_id_info$X),'ref_symbol1']
  HM_to_refsymbol[is.na(HM_to_refsymbol)]=rownames( object_list[['HM_YS']])[is.na(HM_to_refsymbol)]
  new_symbol_list[['HM_YS']]=HM_to_refsymbol;rm(HM_to_refsymbol)
  table(new_symbol_list[['HM_YS']] %in% all_shared_ensembl_id_info$ref_symbol1) # F/T:368/16301
  
  new_object_list=list();library(tidyverse)
  for(refId in names(new_symbol_list)){
    du_symbol=new_symbol_list[[refId]][duplicated(new_symbol_list[[refId]])]
    tmp_assay=GetAssayData(object_list[[refId]],assay = 'RNA')
    rownames(tmp_assay)=new_symbol_list[[refId]]
    du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% du_symbol,]
    tmp_assay=tmp_assay[!rownames(tmp_assay) %in% du_symbol,]
    
    du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum) 
    du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
    du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
    
    tmp_assay=rbind(tmp_assay,du_tmp_assay)
    object_list[[refId]]=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =object_list[[refId]]@meta.data )
  }
  
  rm(du_tmp_assay,tmp_assay);gc()
  
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$ref_symbol1)})
  #Wang_et_al_Cell_Stem Cell_2021 Popescu_et_al_Nature_2019 De_novo HM_YS WE_YS
  #FALSE                              3                        13      21   371  1556
  #TRUE                           18397                     19716   23062 16032 22240
  
  YS_altas_seu_umap=YS_altas_seu[['umap']]
  YS_altas_seu=merge(object_list[[1]],c(object_list[[2]],object_list[[3]],object_list[[4]],object_list[[5]]))
  YS_altas_seu[['RNA']]=JoinLayers(YS_altas_seu[['RNA']])
  table(rownames(YS_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol1)# F/T:1939/24374 
  
  YS_altas_seu[['umap']]= YS_altas_seu_umap
  saveRDS(YS_altas_seu,file = 'NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds')
  rm(object_list,YS_altas_seu_umap);gc()
  
  
  DimPlot(YS_altas_seu,group.by = c('orig.dataset','stage'),cols = col,raster=FALSE)
  DimPlot(YS_altas_seu,group.by = c('LVL1','LVL2'),cols = col,raster=FALSE)
  
  
  
  
}else{
  YS_altas_seu=readRDS('NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds')
}


################################  ############################################################################################### 
#----------------------------------defined the YS_NRBC subcelltype --------------------------------------------------------#
############################################################################################################################### 
#----------- QC ------------#
YS_altas_Ery_seu=subset(YS_altas_seu,LVL1=='ERYTHROID');rm(YS_altas_seu);gc() 
DimPlot(YS_altas_Ery_seu,group.by = c('fetal.ids','stage','LVL3'),cols = col,raster=FALSE)

# 剔除 def_Ery
cells=subset(YS_altas_Ery_seu,subset =  HBB > 4,layer='counts' );cells # HBB>2:4858, HBB>4:599, 
FeaturePlot(cells,reduction = 'umap',features = c('HBB','HBE1'))
YS_altas_Ery_seu=subset(YS_altas_Ery_seu,cell=colnames(YS_altas_Ery_seu)[!colnames(YS_altas_Ery_seu) %in% colnames(cells)])



#------------------------singleR NRBC 亚类注释----------------------------#
YS_altas_Ery_seu[['RNA']]=split(YS_altas_Ery_seu[['RNA']],f = YS_altas_Ery_seu$id)
YS_altas_Ery_seu=FindVariableFeatures(YS_altas_Ery_seu) %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:15,reduction.name = 'raw_uamp') 

YS_altas_Ery_seu[['RNA']]=JoinLayers(YS_altas_Ery_seu[['RNA']])
YS_altas_Ery_seu=singleR_analysis_func(refdata =nrbc_ref_se2,test =as.matrix(YS_altas_Ery_seu@assays$RNA$data)[VariableFeatures(YS_altas_Ery_seu),],outdata = YS_altas_Ery_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype')
YS_altas_Ery_seu=YS_altas_Ery_seu[[1]]
YS_altas_Ery_seu=singleR_analysis_func(refdata =nrbc_ref_se,test =as.matrix(YS_altas_Ery_seu@assays$RNA$data)[VariableFeatures(YS_altas_Ery_seu),],outdata = YS_altas_Ery_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype2')
YS_altas_Ery_seu=YS_altas_Ery_seu[[1]]

table(is.na(YS_altas_Ery_seu$Pre_celltype)) # F/T:42471 19594
YS_altas_Ery_seu$Pre_celltype=factor( YS_altas_Ery_seu$Pre_celltype,levels = c('BFUE','CFUE','ProE','eBas','lBas','Poly','Orth'))
YS_altas_Ery_seu$Pre_celltype2=factor( YS_altas_Ery_seu$Pre_celltype2,levels = c('BFUE','CFUE','ProE','eBas','lBas','Poly','Orth'))
table(is.na(YS_altas_Ery_seu$Pre_celltype2)) # F/T:62064     1


# test batch effect, deepseek 
YS_altas_Ery_seu=RunHarmony(YS_altas_Ery_seu,group.by.vars=c('sequencing.type','orig.dataset','id'),reduction.save = "harmony")
YS_altas_Ery_seu=RunUMAP(YS_altas_Ery_seu,dims = 1:10,reduction = 'harmony',reduction.name='har_uamp')
YS_altas_Ery_seu=RunHarmony(YS_altas_Ery_seu,group.by.vars=c('sequencing.type','orig.dataset','id','lanes'),reduction.save = "harmony1")
YS_altas_Ery_seu=RunUMAP(YS_altas_Ery_seu,dims = 1:10,reduction = 'harmony1',reduction.name='har1_uamp')
YS_altas_Ery_seu=RunHarmony(YS_altas_Ery_seu,group.by.vars=c('orig.dataset','id'),reduction.save = "harmony2")
YS_altas_Ery_seu=RunUMAP(YS_altas_Ery_seu,dims = 1:15,reduction = 'harmony2',reduction.name='har2_uamp')

DimPlot(YS_altas_Ery_seu,group.by = c('orig.dataset','Pre_celltype2','LVL3'),cols = col,raster=FALSE,reduction = 'raw_uamp')/
DimPlot(YS_altas_Ery_seu,group.by = c('orig.dataset','Pre_celltype2','LVL3'),cols = col,raster=FALSE,reduction = 'har_uamp')/
DimPlot(YS_altas_Ery_seu,group.by =c('orig.dataset','Pre_celltype2','LVL3'),cols = col,raster=FALSE,reduction = 'har1_uamp')/
DimPlot(YS_altas_Ery_seu,group.by =c('orig.dataset','Pre_celltype2','LVL3'),cols = col,raster=FALSE,reduction = 'har2_uamp')

# choose harmony2:orig.dataset +id 
YS_altas_Ery_seu[['harmony']]=NULL;YS_altas_Ery_seu[['harmony1']]=NULL;YS_altas_Ery_seu[['har_uamp']]=NULL;YS_altas_Ery_seu[['har1_uamp']]=NULL

YS_altas_Ery_seu=RunUMAP(YS_altas_Ery_seu,dims = 1:6,reduction = 'harmony2',reduction.name='har3_uamp',return.model = T)
DimPlot(YS_altas_Ery_seu,group.by =c('orig.dataset','Pre_celltype2','LVL3'),cols = col,raster=FALSE,reduction = 'har3_uamp')/
  FeaturePlot(YS_altas_Ery_seu,features =c( 'HBE1','HBB'),reduction = 'har3_uamp',ncol = 3)
DimPlot(YS_altas_Ery_seu,group.by =c('orig.dataset','Pre_celltype2','LVL3'),cols = col,raster=FALSE,reduction = 'har3_uamp')

old_YS_altas_Ery_seu=readRDS('NRBC_YS_altas/new_YS_altas_Ery_seu.rds')
old_YS_altas_Ery_seu_meta=old_YS_altas_Ery_seu@meta.data;rm(old_YS_altas_Ery_seu);gc()
YS_altas_Ery_seu$celltype_v1=old_YS_altas_Ery_seu_meta[rownames(YS_altas_Ery_seu@meta.data),c('celltype')]
DimPlot(YS_altas_Ery_seu,group.by =c('RNA_snn_res.0.1','celltype_v1','Pre_celltype2'),cols = col,raster=FALSE,reduction = 'har3_uamp')


YS_altas_Ery_seu=FindNeighbors(YS_altas_Ery_seu,reduction = 'harmony2')
YS_altas_Ery_seu=FindClusters(YS_altas_Ery_seu,resolution = c(0.1,0.2,0.3,0.4))
DimPlot(YS_altas_Ery_seu,reduction = 'har3_uamp',cols = col,group.by =c('RNA_snn_res.0.1','RNA_snn_res.0.2','RNA_snn_res.0.3','RNA_snn_res.0.4'),label = T)
# cho RNA_snn_res.0.1 & RNA_snn_res.0.2
DimPlot(YS_altas_Ery_seu,group.by =c('RNA_snn_res.0.1','RNA_snn_res.0.2','Pre_celltype2'),cols = col,raster=FALSE,reduction = 'har3_uamp')
Idents(YS_altas_Ery_seu)='RNA_snn_res.0.1'
YS_altas_Ery_seu=FindSubCluster(YS_altas_Ery_seu,cluster = '0',graph.name = 'RNA_snn',resolution = 0.5)
DimPlot(YS_altas_Ery_seu,group.by =c('RNA_snn_res.0.1','RNA_snn_res.0.2','sub.cluster','Pre_celltype2'),cols = col,raster=FALSE,reduction = 'har3_uamp')
DimPlot(YS_altas_Ery_seu,group.by =c('sub.cluster','Pre_celltype2','celltype_v1'),cols = col,raster=FALSE,reduction = 'har3_uamp')
temp=AggregateExpression(subset(YS_altas_Ery_seu,RNA_snn_res.0.1=='0'),group.by = 'sub.cluster',features = VariableFeatures(YS_altas_Ery_seu))$RNA
pheatmap(cor(as.matrix(temp)))
YS_altas_Ery_seu$cluster=YS_altas_Ery_seu$sub.cluster
YS_altas_Ery_seu$cluster[YS_altas_Ery_seu$sub.cluster %in% c('0_0','0_2','0_5','0_6')]='0_0'
YS_altas_Ery_seu$cluster[YS_altas_Ery_seu$sub.cluster %in% c('0_1','0_3','0_4')]='0_1'

celltype_features=c('CD34','KIT',"TFRC","GYPA",'HBE1',"MALAT1",'NCL','CD63','MKI67')
DimPlot(YS_altas_Ery_seu,group.by =c('cluster','celltype_v1','Pre_celltype2'),cols = col,raster=FALSE,reduction = 'har3_uamp')/
(VlnPlot(YS_altas_Ery_seu,features = celltype_features,group.by='cluster',stack = T,cols = cols)+NoLegend()+VlnPlot(YS_altas_Ery_seu,features = celltype_features,group.by='celltype_v1',stack = T,cols = cols)+NoLegend())





############################################################################################################################################################################
#-----------------------------------------F61 as reference umap --------------------------------#
############################################################################################################################################################################



YS_altas_F61_NRBC_seu=subset(YS_altas_Ery_seu,orig.ident=='F61')

YS_altas_F61_NRBC_seu=RunHarmony(YS_altas_F61_NRBC_seu,group.by.vars=c('lanes'))
ElbowPlot(YS_altas_F61_NRBC_seu)
# test pca number
YS_altas_F61_NRBC_seu=RunUMAP(YS_altas_F61_NRBC_seu,reduction = 'harmony',dims = 1:10,reduction.name ='umap1' )
YS_altas_F61_NRBC_seu=RunUMAP(YS_altas_F61_NRBC_seu,reduction = 'harmony',dims = 1:15,reduction.name ='umap_pc15' )
YS_altas_F61_NRBC_seu=RunUMAP(YS_altas_F61_NRBC_seu,reduction = 'harmony',dims = 1:6,reduction.name ='umap_pc6')
YS_altas_F61_NRBC_seu=RunUMAP(YS_altas_F61_NRBC_seu,reduction = 'harmony',dims = 1:8,reduction.name ='umap_pc8')

DimPlot(YS_altas_F61_NRBC_seu,reduction = 'umap1',cols=cols,label = T)+
  DimPlot(YS_altas_F61_NRBC_seu,reduction = 'umap_pc15',cols=cols,label = T)+
  DimPlot(YS_altas_F61_NRBC_seu,reduction = 'umap_pc6',cols=cols,label = T)+
  DimPlot(YS_altas_F61_NRBC_seu,reduction = 'umap_pc8',cols=cols,label = T)

YS_altas_F61_NRBC_seu=RunUMAP(YS_altas_F61_NRBC_seu,reduction = 'harmony',dims = 1:6,reduction.name ='umap1',return.model = T,spread =0.5 ,seed.use = 42 )
DimPlot(YS_altas_F61_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('seurat_clusters','celltype'))

YS_altas_F61_NRBC_seu=FindNeighbors(YS_altas_F61_NRBC_seu,reduction = 'harmony',dims = 1:20)
YS_altas_F61_NRBC_seu=FindClusters(YS_altas_F61_NRBC_seu,resolution = c(0.2,0.1))
DimPlot(YS_altas_F61_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('RNA_snn_res.0.2','RNA_snn_res.0.1'))

celltype_marker=c('KIT',"TFRC","GYPA",'HBE1',"MALAT1",'CD63','MKI67')
DimPlot(YS_altas_F61_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('RNA_snn_res.0.2','celltype_v1'))/
VlnPlot(YS_altas_F61_NRBC_seu,features = celltype_marker,group.by = 'RNA_snn_res.0.2',stack = T)+NoLegend()

temp=AggregateExpression(YS_altas_F61_NRBC_seu,group.by = 'RNA_snn_res.0.2',features = VariableFeatures(YS_altas_F61_NRBC_seu))$RNA
pheatmap(cor(as.matrix(temp)))

cluster_celltype_list=list('5'='lBas2','6'='ProE','4'='eBas','3'='mBas','7'='mBas','1'='lBas','0'='Poly','2'='Orth')
YS_altas_F61_NRBC_seu$celltype=as.character(cluster_celltype_list[as.character(YS_altas_F61_NRBC_seu$RNA_snn_res.0.2)])
YS_altas_F61_NRBC_seu$celltype=factor(YS_altas_F61_NRBC_seu$celltype,levels = c('ProE','eBas','mBas','lBas','lBas2','Poly','Orth'))

YS_altas_F61_NRBC_seu$Pre_celltype2=factor( YS_altas_F61_NRBC_seu$Pre_celltype2,levels = c('BFUE','CFUE','ProE','eBas','lBas','Poly','Orth'))

p=DimPlot(YS_altas_F61_NRBC_seu,group.by = c('RNA_snn_res.0.2'),cols = col,reduction = 'umap1')+ggtitle('YS_F61.Seurat.cluster')+
  (DimPlot(YS_altas_F61_NRBC_seu,group.by = c('Pre_celltype2'),cols = col,reduction = 'umap1')+ggtitle('YS_F61.Predict.celltype'))
p
ggsave(p,filename='NRBC_YS_altas/res_pic/YS_altas_F61_NRBC_umap_cluster_predict_celltype_umap.pdf',width = 10,height = 5)

YS_altas_F61_NRBC_seu$final_celltype=as.character(YS_altas_F61_NRBC_seu$celltype)
YS_altas_F61_NRBC_seu$final_celltype[YS_altas_F61_NRBC_seu$final_celltype %in% c('eBas','mBas','lBas','lBas2')]='Bas'
YS_altas_F61_NRBC_seu$final_celltype=factor(YS_altas_F61_NRBC_seu$final_celltype,levels =c('ProE','Bas','Poly','Orth') )  

celltype_marker=c('KIT',"TFRC","GYPA",'HBE1',"NCL",'CD63','MKI67')
p=DimPlot(YS_altas_F61_NRBC_seu,group.by = c('final_celltype'),cols = col,reduction = 'umap1')+ggtitle('YS_F61.celltype')+
  VlnPlot(YS_altas_F61_NRBC_seu,group.by ='final_celltype',features =celltype_marker,stack = T ,cols = cols)+NoLegend()
p
ggsave(p,filename='NRBC_YS_altas/res_pic/YS_altas_F61_NRBC_umap_celltype_marker.pdf',width = 10,height = 5)

saveRDS(YS_altas_F61_NRBC_seu,file = 'NRBC_YS_altas/raw_ref_data/YS_altas_Ery_seu.rds')



Other_YS_altas_Ery_seu=subset(YS_altas_Ery_seu,orig.ident!='F61')
anchors=FindTransferAnchors(reference =YS_altas_F61_NRBC_seu,query = Other_YS_altas_Ery_seu, reference.reduction='pca')
Other_YS_altas_Ery_seu=MapQuery(anchorset =anchors,query = Other_YS_altas_Ery_seu,reference =YS_altas_F61_NRBC_seu,reference.reduction = 'pca',reduction.model = 'umap1',refdata =  list(celltype = "celltype")  )
DimPlot(Other_YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,label = T,group.by = c('predicted.celltype'))
FeaturePlot(Other_YS_altas_Ery_seu,features = c('TFRC','HBE1','GYPA'),reduction = 'ref.umap',cols = c('gray','firebrick3'))

Other_YS_altas_Ery_seu=FindNeighbors(Other_YS_altas_Ery_seu,reduction = 'ref.umap',dims = 1:2)
Other_YS_altas_Ery_seu=FindClusters(Other_YS_altas_Ery_seu,resolution = 0.1)
DimPlot(Other_YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,label = T,group.by = c('predicted.celltype','seurat_clusters'))
Other_YS_altas_Ery_seu$celltype=Other_YS_altas_Ery_seu$predicted.celltype
Other_YS_altas_Ery_seu$celltype[Other_YS_altas_Ery_seu$seurat_clusters %in% c('11','12')]='lBas'
Other_YS_altas_Ery_seu$celltype[Other_YS_altas_Ery_seu$seurat_clusters %in% c('8')]='ProE'


YS_altas_Ery_seu[['ref.pca']]=merge(CreateDimReducObject(embeddings = YS_altas_F61_NRBC_seu[['pca']]@cell.embeddings[,1:30],loadings = YS_altas_F61_NRBC_seu[['pca']]@feature.loadings[,1:30]),Other_YS_altas_Ery_seu[['ref.pca']])
YS_altas_Ery_seu[['ref.umap']]=merge(YS_altas_F61_NRBC_seu[['umap1']],Other_YS_altas_Ery_seu[['ref.umap']])
YS_altas_Ery_seu$celltype='UN'
YS_altas_Ery_seu@meta.data[rownames(Other_YS_altas_Ery_seu@meta.data),'celltype']=Other_YS_altas_Ery_seu$celltype
YS_altas_Ery_seu@meta.data[rownames(YS_altas_F61_NRBC_seu@meta.data),'celltype']=as.character(YS_altas_F61_NRBC_seu$celltype)
YS_altas_Ery_seu$celltype=factor(YS_altas_Ery_seu$celltype,levels =levels(YS_altas_F61_NRBC_seu$celltype) )  

#finally define the  subcelltype of NRBC 
YS_altas_Ery_seu$final_celltype=as.character(YS_altas_Ery_seu$celltype)
YS_altas_Ery_seu$final_celltype[YS_altas_Ery_seu$final_celltype %in% c('eBas','mBas','lBas','lBas2')]='Bas'
YS_altas_Ery_seu$final_celltype=factor(YS_altas_Ery_seu$final_celltype,levels =c('ProE','Bas','Poly','Orth') )  
DimPlot(YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,label = T,group.by = c('celltype_v1','celltype','final_celltype'))

p=DimPlot(YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,group.by = c('id'))
ggsave(p,filename='NRBC_YS_altas/res_pic/YS_altas_Ery_sample_info_refUmap.pdf',width = 6,height = 5)

p=DimPlot(YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,group.by = c('final_celltype'))+ggtitle('YS: 62,065')
p
ggsave(p,filename='NRBC_YS_altas/res_pic/YS_altas_Ery_celltype_refUmap.pdf',width = 6,height = 5)


celltype_df=data.frame(table(YS_altas_Ery_seu@meta.data[,c('id','final_celltype')]))
celltype_fre_Pro_df3=data.frame(table(YS_altas_Ery_seu@meta.data[,c('id')]))
colnames(celltype_fre_Pro_df3)=c('id','Count')

scale_factor=23663# 可以美化双轴
p=ggplot(celltype_df,aes(x = id,y=Freq,fill=final_celltype))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+scale_fill_manual(values =  col[1:11])+
   theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold')) +
  geom_label(data = celltype_fre_Pro_df3,aes(x=id,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
p
ggsave(p,filename = 'NRBC_YS_altas/res_pic/YS_subcelltype_freq_barplot.pdf',width = 8,height = 6,dpi = 300)
rm(p)

table(YS_altas_Ery_seu@meta.data[,c('id','component','final_celltype')])


temp_df=data.frame(table(YS_altas_Ery_seu@meta.data[,c('id','component')]))
temp_df=temp_df[temp_df$Freq >0,]
p=ggplot(temp_df,aes(x=id,y=Freq,color=component,group=component))+geom_point()+geom_line()+theme_classic()+
  RotatedAxis()+scale_color_manual(values = cols[-2])+theme(text = element_text(face = 'bold'))#  log2(Freq) 更能反应波动
ggsave(as.ggplot(p),filename = 'NRBC_YS_altas/res_pic/YS_source_count_poin_line.pdf',width = 6,height = 4,dpi = 300)

saveRDS(YS_altas_Ery_seu,file = 'NRBC_YS_altas/YS_altas_Ery_seu.rds')

rm(list=ls());gc()



#############################################################################################################################
#---------add  the subcelltype of nRBC to YS ALTAS----------#
#############################################################################################################################
DimPlot(YS_altas_seu,group.by = 'LVL2',cols = cols) # 169,494 cells
YS_altas_seu$subcelltype=YS_altas_seu$LVL2
sort(table(YS_altas_seu$subcelltype))# MONO_Mac_DC 细胞很少
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('EOSINOPHIL_BASOPHIL','MAST_CELL')]='EO/BASO/MAST'
DimPlot(subset(YS_altas_seu,subcelltype %in% c('MONOCYTE_MACROPHAGE','MONO MAC PRE DC2','PDC PRECURSOR','MONO MAC DC2','PROMONOCYTE','MONOCYTE_YS_1', 'MONOCYTE','MOP','DC')),group.by = 'LVL2',cols = cols)
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('MONO MAC PRE DC2','PDC PRECURSOR','MONO MAC DC2','DC')]='MOMO_MAC_DC'
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('PROMONOCYTE')]='MOP'
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('MONOCYTE_YS_1')]='MONOCYTE'
YS_altas_seu$subcelltype[YS_altas_seu$LVL3 %in% c('DEF_HSPC_1','DEF_HSPC_2')]='DEF_HSPC'
YS_altas_seu$subcelltype[YS_altas_seu$LVL3 %in% c('HSPC_1','HSPC_2')]='HSPC'
YS_altas_seu$subcelltype[YS_altas_seu$LVL3 %in% c('PRIM_HSPC_1','PRIM_HSPC_2')]='PRIM_HSPC'

YS_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),'subcelltype']=as.character(YS_altas_Ery_seu$final_celltype)
sort(table(YS_altas_seu$subcelltype)) # 599 ERYTHROID:HBB>4,而被剔除，YS 中有部分Ery细胞既表达HBE1又表达HBB， 可以知道，YS 中存在一波过渡造血，或者说血红蛋白转换
YS_altas_seu$subcelltype=gsub(pattern = 'MK',replacement = 'MEGAKARYOCYTE',YS_altas_seu$subcelltype)

celltype_levels=c('EO/BASO/MAST','MEGAKARYOCYTE','ProE','Bas','Poly','Orth', 'ERYTHROID','MACROPHAGE','PRIM_HSPC','DEF_HSPC','HSPC','CMP','LMPP','MEMP','MOP','MONOCYTE','MOMO_MAC_DC','ELP','ILC','NK','B_CELL','NEUTROPHIL_PRECURSOR','MONOCYTE_MACROPHAGE',                 
                  'MESOTHELIUM','SMOOTH_MUSCLE', 'ENDODERM','ENDOTHELIUM','FIBROBLAST' )
table(unique(YS_altas_seu$subcelltype )%in% celltype_levels)
YS_altas_seu$subcelltype=factor(YS_altas_seu$subcelltype,levels = celltype_levels)
DimPlot(subset(YS_altas_seu,subcelltype!='ERYTHROID'),group.by = 'subcelltype',cols = cols,raster = F) # 169,494-599=168,895 cells
saveRDS(YS_altas_seu,file ='NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' )




############################################################################################################################################################################
#-----------------------------------------cellchat analysis--------------------------------#




######################################################################################################################################################################
# ----------------------YS NRBC in embryo blood--------------#
######################################################################################################################################################################

# CD12-CS15 时期血液系统单细胞数据,此时，只有 YS 产生红细胞，FL还未产生 dNRBC
if(F){
  library(Matrix)
  early_organogenesis_count=readMM(file = gzfile('ref_data/ref_scRNAseq_data/GSE157329_early_organogenesis_in_human_embryos/GSE157329_raw_counts.mtx.gz'),header = T)
  dim(early_organogenesis_count)
  
  gene_info=read.table(gzfile('ref_data/ref_scRNAseq_data/GSE157329_early_organogenesis_in_human_embryos/GSE157329_gene_annotate.txt.gz'),header = T)
  head(gene_info)
  rownames(early_organogenesis_count)=gene_info$gene_short_name
  dim(gene_info)
  
  
  early_organogenesis_meta=read.table(file = gzfile('ref_data/ref_scRNAseq_data/GSE157329_early_organogenesis_in_human_embryos/GSE157329_cell_annotate.txt.gz'),sep="\t",header = T)
  head(early_organogenesis_meta)
  rownames(early_organogenesis_meta)=early_organogenesis_meta$cell_id;early_organogenesis_meta=early_organogenesis_meta[,-1]
  colnames(early_organogenesis_count)=rownames(early_organogenesis_meta)
  early_organogenesis_count=early_organogenesis_count[rowSums(early_organogenesis_count)>5,]
  dup_num=which(rownames(early_organogenesis_count) %in% rownames(early_organogenesis_count)[duplicated(rownames(early_organogenesis_count))])
  dup_early_organogenesis_count=early_organogenesis_count[dup_num,]
  
  dup_early_organogenesis_count=data.frame(LINC01238=colSums(dup_early_organogenesis_count[1:2,]),
                                           CYB561D2=colSums(dup_early_organogenesis_count[3:4,]),
                                           MATR3=colSums(dup_early_organogenesis_count[5:6,]),
                                           TMSB15B=colSums(dup_early_organogenesis_count[7:8,]),
                                           LINC01505=colSums(dup_early_organogenesis_count[9:10,]),
                                           GOLGA8M=colSums(dup_early_organogenesis_count[9:10,])
  )
  dup_early_organogenesis_count=t(dup_early_organogenesis_count)
  
  
  early_organogenesis_count=early_organogenesis_count[-dup_num,]
  early_organogenesis_count=rbind(early_organogenesis_count,dup_early_organogenesis_count)
  early_organogenesis_count=as(early_organogenesis_count,'dgTMatrix')
  
  early_organogenesis_seu=CreateSeuratObject(counts =early_organogenesis_count, meta.data =early_organogenesis_meta,project = 'early_organogenesis' )
  rm(early_organogenesis_count,early_organogenesis_meta)
  
  
  early_organogenesis_seu=NormalizeData(early_organogenesis_seu)
  early_organogenesis_seu[['RNA']]=split(early_organogenesis_seu[['RNA']],f = early_organogenesis_seu$sample)
  
  early_organogenesis_seu=subset(early_organogenesis_seu,developmental.system=='blood')
  early_organogenesis_seu=FindVariableFeatures(early_organogenesis_seu)
  early_organogenesis_seu=ScaleData(early_organogenesis_seu)
  early_organogenesis_seu=RunPCA(early_organogenesis_seu)
  ElbowPlot(early_organogenesis_seu)
  early_organogenesis_seu=RunUMAP(early_organogenesis_seu,dims = 1:15)
  DimPlot(early_organogenesis_seu,group.by =c('embryo','dissection_part','stage','annotation'),cols = col)
  
  FeaturePlot(early_organogenesis_seu,features = c('HBE1','HBZ','HBG1','HBG2'),cols = c('white','firebrick3'))
  
  early_organogenesis_seu=IntegrateLayers(early_organogenesis_seu,method = HarmonyIntegration,orig.reduction = "pca",new.reduction='harmony')
  early_organogenesis_seu=RunUMAP(early_organogenesis_seu,reduction = 'harmony',reduction.name = 'har_umap',dims = 1:15)
  #viscera: 内脏, 内脏下干:viscera-lowerTrunk
  DimPlot(early_organogenesis_seu,group.by =c('embryo','dissection_part','stage','annotation'),cols = col,reduction = 'har_umap')
  early_organogenesis_seu[['RNA']]=JoinLayers(early_organogenesis_seu[['RNA']])
  
  FeaturePlot(early_organogenesis_seu,features = c('HBE1','HBZ','HBG1','HBG2'),cols = c('gray','firebrick3'),reduction = 'har_umap')
  FeaturePlot(early_organogenesis_seu,features = c('HBE1','HBZ','HBG1','HBG2'),cols = c('gray','firebrick3'),reduction = 'har_umap',split.by = 'stage')
  
  early_organogenesis_seu=RunUMAP(early_organogenesis_seu,reduction = 'harmony',reduction.name = 'har_umap1',dims = 1:6)
  DimPlot(early_organogenesis_seu,group.by =c('embryo','dissection_part','stage','annotation'),cols = col,reduction = 'har_umap1')
  FeaturePlot(early_organogenesis_seu,features = c('KIT','TFRC','GYPA','CD63'),cols = c('gray','firebrick3'),reduction = 'har_umap1')
  
  
  VlnPlot(early_organogenesis_seu,features = 'HBD',split.by = 'annotation',cols = col)+
    FeaturePlot(early_organogenesis_seu,features = 'HBD',cols = c('gray','firebrick3'),reduction = 'har_umap')
  
  DotPlot(subset(early_organogenesis_seu,annotation=='erythroid'),group.by = 'stage',features = heamoglobin_genes)+RotatedAxis()
  saveRDS(early_organogenesis_seu,file = 'ref_data/early_organogenesis_seu.rds')
  
  Ery_early_organogenesis_seu=subset(early_organogenesis_seu,annotation=='erythroid')
  Ery_early_organogenesis_seu=RunUMAP(Ery_early_organogenesis_seu,reduction = 'harmony',dims = 1:8,reduction.name = 'umap1')
  Ery_early_organogenesis_seu=FindVariableFeatures(Ery_early_organogenesis_seu) %>%NormalizeData()
  Ery_early_organogenesis_seu=singleR_analysis_func(refdata =nrbc_ref_se,test =as.matrix(Ery_early_organogenesis_seu@assays$RNA$data)[VariableFeatures(Ery_early_organogenesis_seu),],outdata = Ery_early_organogenesis_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype')
  Ery_early_organogenesis_seu=Ery_early_organogenesis_seu[[1]]
  Ery_early_organogenesis_seu$Pre_celltype=factor(Ery_early_organogenesis_seu$Pre_celltype,levels =c('BFUE','CFUE','ProE','eBas','mBas','lBas','Poly','Orth'))
  FeaturePlot(Ery_early_organogenesis_seu,features = c('HBE1'),reduction = 'umap1')
  Ery_early_organogenesis_seu <- CellCycleScoring(Ery_early_organogenesis_seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  
  p=FeaturePlot(Ery_early_organogenesis_seu,features = 'HBE1',reduction = 'umap1')
  ggsave(p,filename = 'NRBC_YS_altas/res_pic/YS_embryo_blood_systerm_nRBC_HBE1_uumap.pdf',width = 4,height = 4,dpi = 300)
  
  p=DimPlot(Ery_early_organogenesis_seu,reduction = 'umap1',group.by = c('stage','Pre_celltype','Phase'),cols = col)
  p
  ggsave(p,filename = 'NRBC_YS_altas/res_pic/YS_embryo_blood_systerm_nRBC_stage_celltype_Phase_uumap.pdf',width = 12,height = 4,dpi = 300)
  
  celltype_df=data.frame(table(Ery_early_organogenesis_seu@meta.data[,c('stage','Pre_celltype')]))
  celltype_fre_Pro_df3=data.frame(table(Ery_early_organogenesis_seu@meta.data[,c('stage')]))
  colnames(celltype_fre_Pro_df3)=c('stage','Count')
  
  p=ggplot(celltype_df,aes(x = stage,y=Freq,fill=Pre_celltype))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+scale_fill_manual(values =  col[1:11])+
    theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold')) +
    geom_label(data = celltype_fre_Pro_df3,aes(x=stage,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
    ggtitle('embryo blood system')
  p
  ggsave(p,filename = 'NRBC_YS_altas/res_pic/YS_embryo_blood_systerm_nRBC_analysis_barplot.pdf',width = 4,height = 8,dpi = 300)
  
  saveRDS(Ery_early_organogenesis_seu,file = 'NRBC_YS_altas/raw_ref_data/Ery_early_organogenesis_seu.rds')
  
  
}else{
  early_organogenesis_seu=readRDS( 'ref_data/early_organogenesis_seu.rds')
  Ery_early_organogenesis_seu=readRDS('NRBC_YS_altas/raw_ref_data/Ery_early_organogenesis_seu.rds')
}





