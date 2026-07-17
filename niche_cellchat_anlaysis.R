
################################################################################################################################################
#--------------------------YS niche cellchat anlaysis-----------------------#
################################################################################################################################################

# Rscript1 level 1 粗颗粒niche 细胞互作水平进行分析
YS_altas_seu=readRDS(''NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' ')
future.seed=TRUE
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )
cho_cells=sample(colnames(YS_altas_seu),size=50000)
YS_altas_seu=subset(YS_altas_seu,cells=cho_cells)
YS_altas_seu=subset(YS_altas_seu, subcelltype %in% names(table(YS_altas_seu$subcelltype))[table(YS_altas_seu$subcelltype) >20])
YS_altas_seu$subcelltype=droplevels(YS_altas_seu$subcelltype, exclude = setdiff(levels(YS_altas_seu$subcelltype),unique(YS_altas_seu$subcelltype)))
YS_altas_seu=NormalizeData(YS_altas_seu)
cellchat=createCellChat(YS_altas_seu,assay = 'RNA',group.by='subcelltype')
rm(YS_altas_seu)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
celllchat= computeCommunProb(cellchat, type = "triMean",raw.use = FALSE,population.size = TRUE)
#celllchat=computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
saveRDS(cellchat,file='YS_subcelltype_cellchat.rds')
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat,file='YS_subcelltype_cellchat.rds')

# Rscript2 level2 更细颗粒niche 细胞互作进行分析
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )
tmp_seu=readRDS('raw_ref_data/dealt_YS_altas_seu_20251028.rds')
cho_cells=sample(colnames(tmp_seu),size=50000)
tmp_seu=subset(tmp_seu,cells=cho_cells)
tmp_seu$all_subcelltype=tmp_seu$LVL3
Ery_subcelltype=c('ProE','Bas','Poly','Orth')
tmp_seu$all_subcelltype[tmp_seu$subcelltype %in% Ery_subcelltype]=as.character(tmp_seu$subcelltype[tmp_seu$subcelltype %in% Ery_subcelltype])

tmp_seu=subset(tmp_seu, all_subcelltype %in% names(table(tmp_seu$all_subcelltype))[table(tmp_seu$all_subcelltype) >20])
tmp_seu$all_subcelltype=factor(tmp_seu$all_subcelltype,levels=c(Ery_subcelltype,sort(unique(tmp_seu$all_subcelltype))[!sort(unique(tmp_seu$all_subcelltype)) %in% Ery_subcelltype]))
#tmp_seu$all_subcelltype=droplevels(tmp_seu$all_subcelltype, exclude = setdiff(levels(tmp_seu$all_subcelltype),unique(tmp_seu$all_subcelltype)))
tmp_seu=NormalizeData(tmp_seu)
cellchat=createCellChat(GetAssayData(tmp_seu,layer='data'),meta=tmp_seu@meta.data,group.by='all_subcelltype')
rm(tmp_seu)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 6) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat= computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat,file='YS_all_subcelltype_cellchat.rds')




################################################################################################################################################
#--------------------------FL niche cellchat anlaysis-----------------------#
################################################################################################################################################
# Rscript1 level1
library(CellChat)
library(Seurat)
future.seed=TRUE
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )
tmp_seu=readRDS('tmp_FL_altas_seu.rds')
cho_random_cells=sample(rownames(tmp_seu@meta.data),size=50000)
tmp_seu=subset(tmp_seu,cells=cho_random_cells)
tmp_seu=subset(tmp_seu, subcelltype %in% names(table(tmp_seu$subcelltype))[table(tmp_seu$subcelltype) >19])
tmp_seu$subcelltype = droplevels(tmp_seu$subcelltype, exclude = setdiff(levels(tmp_seu$subcelltype),unique(tmp_seu$subcelltype)))
cellchat=createCellChat(tmp_seu,assay = 'RNA',,group.by='subcelltype')
rm(tmp_seu)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat=computeCommunProb(cellchat)
print('----------------sucess------------------')

cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat,file='FL_subcelltyp_cellchat.rds')


# Rscript2 level2
library(CellChat)
library(Seurat)
future.seed=TRUE
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )
tmp_seu=readRDS('tmp_FL_altas_seu.rds')
cho_random_cells=sample(rownames(tmp_seu@meta.data),size=50000)
tmp_seu=subset(tmp_seu,cells=cho_random_cells)
tmp_seu$all_subcelltype=tmp_seu$anno_lvl_2_final_clean
Ery_subcelltype=c("BFUE/CFUE","ProE","Bas","Poly","Orth" )
tmp_seu$all_subcelltype[tmp_seu$subcelltype %in% Ery_subcelltype]=as.character(tmp_seu$subcelltype[tmp_seu$subcelltype %in% Ery_subcelltype])
tmp_seu=subset(tmp_seu, all_subcelltype %in% names(table(tmp_seu$all_subcelltype))[table(tmp_seu$all_subcelltype) >19])
tmp_seu$all_subcelltype=factor(tmp_seu$all_subcelltype,levels=c(Ery_subcelltype,sort(unique(tmp_seu$all_subcelltype))[!sort(unique(tmp_seu$all_subcelltype)) %in% Ery_subcelltype ]) )
#tmp_seu$subcelltype = droplevels(tmp_seu$subcelltype, exclude = setdiff(levels(tmp_seu$subcelltype),unique(tmp_seu$subcelltype)))
tmp_seu=NormalizeData(tmp_seu)
cellchat=createCellChat(GetAssayData(tmp_seu,layer='data'),meta=tmp_seu@meta.data,group.by='all_subcelltype')
rm(tmp_seu)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat=computeCommunProb(cellchat)
print('----------------sucess------------------')

cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat,file='FL_all_subcelltyp_cellchat.rds')



################################################################################################################################################
#--------------------------FBM niche cellchat anlaysis-----------------------#
################################################################################################################################################
# Rscript1 level1
future.seed=TRUE
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )

BM_altas_seu=readRDS('BM_altas_seu_v2.rds')# GSE253355 少了大半细胞
FBM_altas_seu=subset(BM_altas_seu,stage=='FBM');rm(BM_altas_seu);gc()
cho_cells=sample(rownames(FBM_altas_seu@meta.data),size=50000)
FBM_altas_seu=subset(FBM_altas_seu,cells=cho_cells)
FBM_altas_seu$new_celltype=droplevels(FBM_altas_seu$new_celltype, exclude = setdiff(levels(FBM_altas_seu$new_celltype),unique(FBM_altas_seu$new_celltype)))
cellchat=createCellChat(FBM_altas_seu,assay = 'RNA',group.by='new_celltype')
rm(FBM_altas_seu);gc()

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 8) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat= computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat,file='FBM_subcelltype_cellchat.rds')


# Rscript2 level2
future.seed=TRUE
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )

BM_altas_seu=readRDS('BM_altas_seu_v2.rds')# GSE253355 少了大半细胞
FBM_altas_seu=subset(BM_altas_seu,stage=='FBM');rm(BM_altas_seu);gc()
cho_cells=sample(rownames(FBM_altas_seu@meta.data),size=50000)
FBM_altas_seu=subset(FBM_altas_seu,cells=cho_cells)
FBM_altas_seu$all_subcelltype=FBM_altas_seu$anno_lvl_2_final_clean
Ery_subcelltype=c('BFUE/CFUE','ProE','Bas','Poly','Orth','YS_Bas/Poly','YS_Orth')
FBM_altas_seu$all_subcelltype[FBM_altas_seu$new_celltype %in% Ery_subcelltype]=as.character(FBM_altas_seu$new_celltype[FBM_altas_seu$new_celltype %in% Ery_subcelltype])
FBM_altas_seu$all_subcelltype=factor(FBM_altas_seu$all_subcelltype,levels=c(Ery_subcelltype,sort(unique(FBM_altas_seu$all_subcelltype))[!sort(unique(FBM_altas_seu$all_subcelltype)) %in% Ery_subcelltype ] ))

cellchat=createCellChat(GetAssayData(FBM_altas_seu,layer='data'),meta=FBM_altas_seu@meta.data,group.by='all_subcelltype')
rm(FBM_altas_seu);gc()

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 8) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat= computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat,file='FBM_allsubcelltype_cellchat.rds')


################################################################################################################################################
#--------------------------ABM niche cellchat anlaysis-----------------------#
################################################################################################################################################
# Rscript1 level1

future.seed=TRUE
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )

ABM_altas_seu=readRDS('BM_altas_seu_v2.rds')
ABM_altas_seu=subset(ABM_altas_seu,stage=='ABM')
ABM_altas_seu=subset(ABM_altas_seu,resource !='GSE150774' & new_celltype!='B CELLs')
cho_cells=sample(rownames(ABM_altas_seu@meta.data),size=50000)# 12 个样本，抽取3/12=1/4
ABM_altas_seu=subset(ABM_altas_seu,cells=cho_cells)
ABM_altas_seu$new_celltype=droplevels(ABM_altas_seu$new_celltype, exclude = setdiff(levels(ABM_altas_seu$new_celltype),unique(ABM_altas_seu$new_celltype)))
cellchat=createCellChat(GetAssayData(ABM_altas_seu,layer='data'),meta=ABM_altas_seu@meta.data,group.by='new_celltype')
rm(ABM_altas_seu);gc()

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 8) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat= computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat,file='ABM_subcelltype_cellchat.rds')

#Rscript2 level 2
future.seed=TRUE
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )

ABM_altas_seu=readRDS('BM_altas_seu_v2.rds')
ABM_altas_seu=subset(ABM_altas_seu,stage=='ABM')
ABM_altas_seu=subset(ABM_altas_seu,resource !='GSE150774' & new_celltype!='B CELLs')# 未分类亚类且散，数目较少
cho_cells=sample(rownames(ABM_altas_seu@meta.data),size=50000)# 12 个样本，抽取3/12=1/4
ABM_altas_seu=subset(ABM_altas_seu,cells=cho_cells)

Ery_subcelltype=c('BFUE/CFUE','ProE','Bas','Poly','Orth')
ABM_altas_seu_meta=ABM_altas_seu@meta.data
ABM_altas_seu_meta$all_subcelltype=ABM_altas_seu_meta$ct
ABM_altas_seu_meta$all_subcelltype[is.na(ABM_altas_seu_meta$all_subcelltype)]=ABM_altas_seu_meta$celltype[is.na(ABM_altas_seu_meta$all_subcelltype)] # GSE253355 注释是celltype
ABM_altas_seu_meta$all_subcelltype[is.na(ABM_altas_seu_meta$all_subcelltype)]=as.character(ABM_altas_seu_meta$new_celltype[is.na(ABM_altas_seu_meta$all_subcelltype)])
ABM_altas_seu_meta$all_subcelltype[ABM_altas_seu_meta$new_celltype %in% Ery_subcelltype]=as.character(ABM_altas_seu_meta$new_celltype[ABM_altas_seu_meta$new_celltype %in% Ery_subcelltype])
df=data.frame(table(ABM_altas_seu_meta[,c('all_subcelltype','new_celltype')]))
df=df[df$Freq >0,]
df_list=split(as.character(df$all_subcelltype),as.character(df$new_celltype))
df[df$new_celltype %in% names(sapply(df_list, length)[sapply(df_list, length) >1]),]
uniue_name_list=list('Proliferation HSC'='Cycling HSPC','Megakaryocyte/erythroid progenitor'='MEP',
                     'Erythroblast'='MEGAKARYOCYTE','Megakaryocyte'='MEGAKARYOCYTE','Megakaryocyte progenitor'='MEGAKARYOCYTE','MEP'='MEGAKARYOCYTE','Platelet'='MEGAKARYOCYTE',
                     'Eosinophil'='Ba/Eo/Ma','Neutrophil'='NEUTROPHIL','Monocyte'='MONOCYTE','CD14 monocyte'='MONOCYTE','CD16 monocyte'='MONOCYTE','Common lymphoid progenitor'='CLP','Plasma Cell'='Plasma cell','Osteoblast'='OSTEOBLAST' )
ABM_altas_seu_meta$final_allsubcelltype=ABM_altas_seu_meta$all_subcelltype
ABM_altas_seu_meta$final_allsubcelltype[ ABM_altas_seu_meta$all_subcelltype %in% names(uniue_name_list)]= as.character(uniue_name_list[ABM_altas_seu_meta$all_subcelltype[ABM_altas_seu_meta$all_subcelltype %in% names(uniue_name_list)]])

ABM_altas_seu@meta.data[,'all_subcelltype']=ABM_altas_seu_meta[rownames(ABM_altas_seu@meta.data),'final_allsubcelltype']
ABM_altas_seu$all_subcelltype=factor(ABM_altas_seu$all_subcelltype,levels=c(Ery_subcelltype,sort(unique(ABM_altas_seu$all_subcelltype))[!sort(unique(ABM_altas_seu$all_subcelltype)) %in% Ery_subcelltype ]))

#ABM_altas_seu$new_celltype=droplevels(ABM_altas_seu$new_celltype, exclude = setdiff(levels(ABM_altas_seu$new_celltype),unique(ABM_altas_seu$new_celltype)))
cellchat=createCellChat(GetAssayData(ABM_altas_seu,layer='data'),meta=ABM_altas_seu@meta.data,group.by='all_subcelltype')
rm(ABM_altas_seu);gc()

print('---------------start CC analysis--------------')
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 8) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat= computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat,file='ABM_allsubcelltype_cellchat.rds')






