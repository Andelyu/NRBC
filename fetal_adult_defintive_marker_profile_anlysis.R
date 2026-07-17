
########################################################################################################################################################################
#-------------------------------------------- fetal NRBC marker#-------------------------------------------- #
########################################################################################################################################################################
setwd("/home/gibh/2021_NRBC_chlyu/NRBC_altas_CC")
load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T)
fetal_adult_NRBC_whole_marker=readRDS('../Protein_NRBC_marker/res_data/fetal_adult_NRBC_whole_marker.rds')
fetal_ligand_target_res_list=readRDS('res_data/fa_fetal_ligand_target_res_list.rds')
fa_candidated_lr_genes=readRDS('res_data/fa_candidated_lr_genes.rds')

adult_ligand_target_res_list=readRDS('res_data/fa_adult_ligand_target_res_list.rds')
fetal_ligand_target_res_list=readRDS('res_data/fa_fetal_ligand_target_res_list.rds')


NRBC_altas_LR_df=read.csv('res_data/filt_NRBC_altas_LR_df.csv',sep="\t")
Other2Ery_df=read.csv('res_data/filt_Other2Ery_df_new.csv',sep="\t")
all_NRBC_receptor_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']),'receptor.symbol'] # 取receptor gene
all_NRBC_ligand_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in%unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']),'ligand.symbol']# 取ligand gene
all_NRBC_receptor_genes=sort(unique(unlist(strsplit(all_NRBC_receptor_genes,split=','))))
all_NRBC_ligand_genes=sort(unique(unlist(strsplit(all_NRBC_ligand_genes,split=','))))
length(all_NRBC_receptor_genes);length(all_NRBC_ligand_genes)# 49， 80


temp_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% rownames(fetal_ligand_target_res_list[[3]]) & fetal_adult_NRBC_whole_marker$avg_log2FC>1,]
#temp_df$gene=factor(temp_df$gene,levels = rownames(fetal_ligand_target_res_list[[3]]))
temp_df=temp_df[temp_df$pct.2<0.1 & temp_df$pct.1>0.1,]
temp_df=temp_df[order(temp_df$avg_log2FC,decreasing = T),]

levels(fa_candidated_lr_genes$gene)[1:11]
# fetal NRBC 金标准：HbF：HBG1/2 抗体鉴定，HBZ
# C4A 可以考虑去掉，fetal 重点：DLK1（early stage），SLC6A9： late stage
# adult: ，HBD，重点关注CA1，ANXA1，TNFSF13B
top_fetal_marker=c('HBG1','HBG2','HBZ','TUBB6','HSPA1A','HSPA1B','IGF2BP1','IGF2BP3','DLK1','MEG3','GATA5','LIN28B','HMGA2','CISH','HIF3A')
fetal_markers=c(top_fetal_marker,c('ROBO2',"WNT5B","PTGR1","ADORA2B"),temp_df$gene);fetal_markers=unique(fetal_markers)# MEG3: lcnRNA
adult_markers=c('CA1','PDZK1IP1',"ANXA1",'NECAB1','TSC22D3','IFIT1B',"TNFSF13B", 'CLEC2B','CXCR4') # 'CA1', fetal bone marrow, 'IFIT1B','NECAB1','TSC22D3' 舍弃

p=DotPlot(filt_NBRC_altas_seu,features = c(fetal_markers,'HBD',adult_markers),scale = F,group.by = 'source_celltype')+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3');p
ggsave(p,filename='../Protein_NRBC_marker/res_pic/main_figure6/fa_marker_candidated_dotplot.pdf',height =6,width = 12,dpi =300 )




#-----------------------------------------------------1.1 primary bulk RNAseq level validated ---------------------------------------------#
# ------------------------bulk RNAseq to validated the unique markers-----------------------#
library(ggplot2)
library(reshape2)
library(edgeR)


#------------------------------ABM nRBC  progenitor bulkRNAseq-----------------#
GSE301441_BM_NRBC_df1=read.table('../ref_data/bulk_RNAseq/GSE301441_BM_primary_nRBC/progenitor_processeddata.txt')
GSE301441_BM_NRBC_df1[c('GAPDH','ACTB'),] # ACTB 作为内参


GSE301441_BM_NRBC_df=read.table('../ref_data/bulk_RNAseq/GSE301441_BM_primary_nRBC/terminal_processeddata.txt')
GSE301441_BM_NRBC_df[c('GAPDH','ACTB'),] # ACTB 作为内参


all_genes=unique(c(rownames(GSE301441_BM_NRBC_df1),rownames(GSE301441_BM_NRBC_df)))

GSE301441_BM_NRBC_df1=rbind(GSE301441_BM_NRBC_df1,data.frame(row.names = all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df1)]))
GSE301441_BM_NRBC_df1[all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df1)],1:dim(GSE301441_BM_NRBC_df1)[2]]=0
GSE301441_BM_NRBC_df=rbind(GSE301441_BM_NRBC_df,data.frame(row.names = all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df)]))
GSE301441_BM_NRBC_df[all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df)],1:dim(GSE301441_BM_NRBC_df)[2]]=0
GSE301441_BM_NRBC_df=cbind(GSE301441_BM_NRBC_df1,GSE301441_BM_NRBC_df);dim(GSE301441_BM_NRBC_df)


GSE301441_BM_NRBC_df=as.matrix(GSE301441_BM_NRBC_df)
GSE301441_BM_NRBC_df=log2(cpm(GSE301441_BM_NRBC_df)+1)

refer_control_GSE301441_BM_NRBC_df=data.frame(t(apply(GSE301441_BM_NRBC_df, FUN =  function(x){x/GSE301441_BM_NRBC_df['ACTB',]},MARGIN = 1)))
refer_control_GSE301441_BM_NRBC_df[c('GAPDH','ACTB'),] # ACTB 作为内参



# -----------GYPA+ NRBC bulRNAseq -----------#
# Erythroid Differentiation Enhances RNA Mis-Splicing in SF3B1-Mutant Myelodysplastic Syndromes with Ring Sideroblasts
ABM_Ery_df=read.csv('../ref_data/bulk_RNAseq/ABM_Ery_BulkRNAseq/BulkRNAseq_Counts.txt',sep='\t',header = T)
rownames(ABM_Ery_df)=ABM_Ery_df$GENE_ID
ABM_Ery_df=ABM_Ery_df[,-1]
ABM_Ery_df=ABM_Ery_df[rowMeans(ABM_Ery_df) >1,]

ABM_metadata_df=read.csv('../ref_data/bulk_RNAseq/ABM_Ery_BulkRNAseq/BulkRNAseq_metadata.tsv',sep='\t',header = T)
ABM_metadata_df=ABM_metadata_df[ABM_metadata_df$DISEASE_STATUS=='NBM' & ABM_metadata_df$SAMPLE_TYPE %in% c('GPA_NBM' , 'RET_NBM' ),]
ABM_metadata_df$id=paste(ABM_metadata_df$PATIENT_ID,ABM_metadata_df$SAMPLE_TYPE,sep = '_')
rownames(ABM_metadata_df)=ABM_metadata_df$FASTQ_FILE
ABM_Ery_df=ABM_Ery_df[,ABM_metadata_df$FASTQ_FILE]
colnames(ABM_Ery_df)=ABM_metadata_df[colnames(ABM_Ery_df),'id']
ABM_Ery_df=data.frame(ABM_Ery_df)


ABM_Ery_df$Symbol=as.character(mapIds(x = org.Hs.eg.db,keys =rownames(ABM_Ery_df),column = 'SYMBOL',keytype = 'ENSEMBL'))
ABM_Ery_df=ABM_Ery_df[!is.na(ABM_Ery_df$Symbol),]
dup_symbols=ABM_Ery_df$Symbol[duplicated(ABM_Ery_df$Symbol)]

temp_df=data.frame()
for(id in dup_symbols ){
  temp=colMaxs(as.matrix(ABM_Ery_df[ABM_Ery_df$Symbol %in% id,1:8])) 
  temp_df=rbind(temp_df,temp)
}

ABM_Ery_df=ABM_Ery_df[!ABM_Ery_df$Symbol %in% dup_symbols,]
rownames(ABM_Ery_df)=ABM_Ery_df$Symbol;ABM_Ery_df=ABM_Ery_df[,1:8]

colnames(temp_df)=colnames(ABM_Ery_df)
rownames(temp_df)=dup_symbols

ABM_Ery_df=rbind(ABM_Ery_df,temp_df)

ABM_Ery_df=ABM_Ery_df[, sort(colnames(ABM_Ery_df))]
ABM_Ery_df=ABM_Ery_df[,1:4]

ABM_Ery_df=as.matrix(ABM_Ery_df)
ABM_Ery_df=log2(cpm(ABM_Ery_df)+1)
ABM_Ery_df[c('GAPDH','ACTB'),] # ACTB 作为内参


refer_control_ABM_Ery_df=data.frame(t(apply(ABM_Ery_df,MARGIN = 1,function(x){x/ABM_Ery_df['ACTB',]})))
refer_control_ABM_Ery_df[c('GAPDH','ACTB'),] # ACTB 作为内参


#------------------------------FL nRBC bulkRNAseq-----------------#
FL_primary_NRBC_df=read.table('../ref_data/bulk_RNAseq/Anxiuli_lab_bulkRNAseq/FL_nRBC_invivo.txt',sep="\t",header = T)
rownames(FL_primary_NRBC_df)=FL_primary_NRBC_df$X;FL_primary_NRBC_df=FL_primary_NRBC_df[,-1]
FL_primary_NRBC_df=as.matrix(FL_primary_NRBC_df)
FL_primary_NRBC_df=log2(cpm(FL_primary_NRBC_df)+1)
FL_primary_NRBC_df[c('GAPDH','ACTB'),] # ACTB 表达稳定，可以作为内参

refer_control_FL_primary_NRBC_df=data.frame(t(apply(FL_primary_NRBC_df, FUN =  function(x){x/FL_primary_NRBC_df['ACTB',]+0.1},MARGIN = 1)))
refer_control_FL_primary_NRBC_df[c('GAPDH','ACTB'),] # ACTB 作为内参

#----------------------merge-------------------------#
all_genes=unique(c(rownames(GSE301441_BM_NRBC_df),rownames(FL_primary_NRBC_df),rownames(ABM_Ery_df)))
refer_control_GSE301441_BM_NRBC_df[all_genes[!all_genes %in% rownames(refer_control_GSE301441_BM_NRBC_df)],1:dim(refer_control_GSE301441_BM_NRBC_df)[2]]=0
refer_control_FL_primary_NRBC_df[all_genes[!all_genes %in% rownames(refer_control_FL_primary_NRBC_df)],1:dim(refer_control_FL_primary_NRBC_df)[2]]=0
refer_control_ABM_Ery_df[all_genes[!all_genes %in% rownames(refer_control_ABM_Ery_df)],1:dim(refer_control_ABM_Ery_df)[2]]=0

FL_ABM_nRBC_df_refer_control_df=cbind(cbind(refer_control_FL_primary_NRBC_df[all_genes,],refer_control_GSE301441_BM_NRBC_df[all_genes,]),refer_control_ABM_Ery_df[all_genes,])
FL_ABM_nRBC_df_refer_control_df[c('GAPDH','ACTB'),]
FL_ABM_nRBC_df_refer_control_df[is.na(FL_ABM_nRBC_df_refer_control_df)]=0

FL_primary_NRBC_df=data.frame(FL_primary_NRBC_df)
GSE301441_BM_NRBC_df=data.frame(GSE301441_BM_NRBC_df)
ABM_Ery_df=data.frame(ABM_Ery_df)
FL_primary_NRBC_df[all_genes[!all_genes %in%  rownames(FL_primary_NRBC_df)],1:dim(FL_primary_NRBC_df)[2]]=0

FL_ABM_nRBC_df_df=cbind(cbind(FL_primary_NRBC_df,GSE301441_BM_NRBC_df[all_genes,]),ABM_Ery_df[all_genes,])
FL_ABM_nRBC_df_df[is.na(FL_ABM_nRBC_df_df)]=0

saveRDS(FL_ABM_nRBC_df_refer_control_df,file = '../ref_data/bulk_RNAseq/nr_FL_ABM_nRBC_ref_ACTB.rds')
saveRDS(FL_ABM_nRBC_df_df,file = '../ref_data/bulk_RNAseq/FL_ABM_nRBC_df_df_nr_exp.rds')

colnames(FL_ABM_nRBC_df_refer_control_df)
# healthy-Ery ；CA1 表达低，且还较高表达HBG1，IGF2BP3 fetal nRBC marker，更接近fetal，可能是样本污染或则应景，排除，且MM 中ery 也表达fetal 结果趋势

FL_ABM_nRBC_df_refer_control_df=FL_ABM_nRBC_df_refer_control_df[,c(colnames(FL_ABM_nRBC_df_refer_control_df)[1:12],colnames(FL_ABM_nRBC_df_refer_control_df)[c(13:15,19:23)],
                                                                   colnames(FL_ABM_nRBC_df_refer_control_df)[29:32],colnames(FL_ABM_nRBC_df_refer_control_df)[c(16:18,24:28)] )]
FL_ABM_nRBC_df_refer_control_df=FL_ABM_nRBC_df_refer_control_df[,-grep('healthy_ery',colnames(FL_ABM_nRBC_df_refer_control_df))]
FL_ABM_nRBC_df_refer_control_df=as.matrix(FL_ABM_nRBC_df_refer_control_df)
col_an_df=data.frame(row.names =colnames(FL_ABM_nRBC_df_refer_control_df),sample=c(rep('FL',12),rep('ABM',7),rep('MABM','8'))) 
col_an_df$sample=factor(col_an_df$sample,levels = c('FL','ABM','MABM'))


p=pheatmap(annotation_col = col_an_df,FL_ABM_nRBC_df_refer_control_df[c('ACTB',fetal_markers,'HBD',adult_markers),],
           cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c("#00008066",'white','firebrick3'))(200))
ggsave(as.ggplot(p),filename='../Protein_NRBC_marker/res_pic/main_figure6/fa_fetal_adult_canididated_specific_markers_bulk_FL_ABM_NRBC_heatmap.pdf',width = 8,height = 8)



#-----------------------------------------------------1.2 protein level validated ---------------------------------------------#
# primary tissue data------------------#
BM_nRBC_MS_df=read.csv('../ref_data/Protein_NRBC/NC_2018_BM_hemo_celltype_LF_MS_data.tsv',header = T,sep="\t",skip = 1)
BM_nRBC_MS_df=BM_nRBC_MS_df[,c(colnames(BM_nRBC_MS_df)[1:4],'ERP.number.of.donors','ERP.normalized.LF.sum')]
BM_nRBC_MS_df=BM_nRBC_MS_df[!is.na(BM_nRBC_MS_df$ERP.normalized.LF.sum),]
BM_Ery_MS_protein_LF_genes=unique(unlist(strsplit(BM_nRBC_MS_df$gene.name[BM_nRBC_MS_df$ERP.normalized.LF.sum >0.1],split = ';')));length(BM_Ery_MS_protein_LF_genes)

adult_HSPC_Ery_MS_protein_df=read.table(file = '../Protein_NRBC_marker/recent_MS_Ery_protein_omics/dealt_adult_HSPC_Ery_MS_protein_df.tsv',sep="\t")
adult_HSPC_Ery_MS_protein_genes=unique(unlist(strsplit(adult_HSPC_Ery_MS_protein_df$Gene.names,split = ';')));length(adult_HSPC_Ery_MS_protein_genes)


#明确含跨膜结构域	DLK1, SLC6A9, PTGR1, ADORA2B, CLEC2B, PDZK1IP1, TNFSF13B
#膜相关但无经典TM	TIMP3 (GPI锚定), WNT5B (脂修饰)
#无跨膜结构域	HBG1, HBG2, CISH, GAL, C4A, IGF2BP1, IGF2BP3, MEG3, HBD, CA1, ANXA1, IFIT1B, NECAB1, TSC22D3
gene_an_df=data.frame(row.names = c(fetal_markers[!fetal_markers %in% c("HBG1","HBG2","HBZ" )],adult_markers))
gene_an_df$gene=rownames(gene_an_df)
gene_an_df$gene=factor(gene_an_df$gene,levels = gene_an_df$gene)

gene_an_df$transmembrane='no'
gene_an_df[ c('DLK1', 'SLC6A9', 'PTGR1', 'ADORA2B', 'CLEC2B', 'PDZK1IP1', 'TNFSF13B','CXCR4'),'transmembrane']='yes'
gene_an_df$type='marker'
gene_an_df$type[rownames(gene_an_df) %in% all_NRBC_receptor_genes]='receptor'
gene_an_df$type[rownames(gene_an_df) %in% all_NRBC_ligand_genes]='ligand'
gene_an_df$type[rownames(gene_an_df) %in% c(rownames(fetal_ligand_target_res_list[[3]]),rownames(adult_ligand_target_res_list[[3]]))]='target_genes'
gene_an_df$type[rownames(gene_an_df) %in% c('DLK1','ANXA1')]='ligand_target_genes'
gene_an_df$MS='no'
gene_an_df$MS[rownames(gene_an_df) %in% c(BM_Ery_MS_protein_LF_genes,adult_HSPC_Ery_MS_protein_genes)]='detected'
gene_an_df$HPA_IHC='no'
gene_an_df$HPA_IHC[rownames(gene_an_df) %in% c('SLC6A9','CLEC2B','IFIT1B')]='no_sure:low'# may low
gene_an_df$HPA_IHC[rownames(gene_an_df) %in% c('CA1','TNFSF13B','ANXA1','HSPA1A')]='high'
gene_an_df$HPA_IHC[rownames(gene_an_df) %in% c('PDZK1IP1','CISH','TUBB6','HMGA2')]='Medium'
#gene_an_df$HPA_IHC[rownames(gene_an_df) %in% c('TSC22D3')]='Low'
gene_an_df$HPA_IHC[rownames(gene_an_df) %in% c('TIMP3','WNT5B','ADORA2B','GATA5')]='No data'
gene_an_df$HPA_IHC[rownames(gene_an_df) %in% c('DLK1','GAL','PTGR1','IGF2BP1','IGF2BP3','NECAB1','HSPA1B','LIN28B','NECAB1','IFIT1B')]='Not detected'

p1=ggplot(gene_an_df,aes(x=gene,y='transmembrane',shape=transmembrane))+geom_point()+theme_classic()+RotatedAxis()+scale_shape_manual(values = c(1,3))+ylab('')
p2=ggplot(gene_an_df,aes(x=gene,y='MS',shape=MS))+geom_point()+theme_classic()+RotatedAxis()+scale_shape_manual(values = c(3, 1,4))+ylab('')+theme(axis.text.x = element_blank())

gene_an_df$HPA_IHC_res='No data/not sure'
gene_an_df$HPA_IHC_res[gene_an_df$HPA_IHC %in%  c('high','Medium','Low')]='Positive'
gene_an_df$HPA_IHC_res[gene_an_df$HPA_IHC %in%  c('Not detected')]='Negtive'

p3=ggplot(gene_an_df,aes(x=gene,y='HPA_IHC_res',shape=HPA_IHC_res))+geom_point()+theme_classic()+RotatedAxis()+scale_shape_manual(values = c(4, 1,3))+ylab('')+theme(axis.text.x = element_blank())

p=p3/p2/p1;p
ggsave(p,filename='../Protein_NRBC_marker/res_pic/main_figure6/fa_fetal_adult_canididated_specific_markers.pdf',width = 12,height = 6)




###################################################################################################################################################
#----------------------------validated the marker profile in erythroblast related cell line
###################################################################################################################################################

#  分析HEL(fetal NRBC, HBG1/HBG2) 、K562( HBE1) 、HUDEP-2 ( HBB)常用红细胞系 中marker 的表达情况
# /mnt/data/bio_program/2021_NRBC_chlyu/ref_data/bulk_RNAseq/erythroblast_cellline_data/HEL erythroleukemia_GSE203060_featureCounts.xlsx

HEL_GSE203060_df=read_excel('../ref_data/bulk_RNAseq/erythroblast_cellline_data/HEL erythroleukemia_GSE203060_featureCounts.xlsx',skip = 1)
table(duplicated(HEL_GSE203060_df$Geneid))
HEL_GSE203060_df=data.frame(HEL_GSE203060_df)
rownames(HEL_GSE203060_df)=HEL_GSE203060_df$Geneid
colnames(HEL_GSE203060_df)=gsub(pattern = '_Aligned.sortedByCoord.out.bam',replacement = '',colnames(HEL_GSE203060_df))
colnames(HEL_GSE203060_df)=gsub(pattern = 'PE01_DMSO.',replacement = 'HEL.',colnames(HEL_GSE203060_df))
HEL_GSE203060_df=HEL_GSE203060_df[,-1:-6]
HEL_GSE203060_df[is.na(HEL_GSE203060_df)]=0
HEL_GSE203060_df=HEL_GSE203060_df[rowSums(HEL_GSE203060_df) >1,]

HEL_GSE203060_df=log2(cpm(HEL_GSE203060_df)+1)
HEL_GSE203060_df=data.frame(HEL_GSE203060_df)
c(fetal_markers,'HBD',adult_markers)[!c(fetal_markers,'HBD',adult_markers) %in% rownames(HEL_GSE203060_df)]
HEL_GSE203060_df[c("HSPA1A", "MEG3"),1:dim(HEL_GSE203060_df)[2]]=0

K562_GSE311284_df=read.csv('../ref_data/bulk_RNAseq/erythroblast_cellline_data/GSE311284_HRG1_K562_counts.csv')
rownames(K562_GSE311284_df)=K562_GSE311284_df$X
K562_GSE311284_df=K562_GSE311284_df[,-1]
K562_GSE311284_df=log2(cpm(K562_GSE311284_df)+1)
#colnames(K562_GSE311284_df)=c(rep(c('WT','HRG1_KO','WT_dealt','KO_dealt'),each=3))
K562_GSE311284_df=data.frame(K562_GSE311284_df)

HUDEP2_GSE314032_df=read.csv('../ref_data/bulk_RNAseq/erythroblast_cellline_data/GHUDEP2 GSE314032_sgControl_sgEGR1_gene_count_matrix.csv')
HUDEP2_GSE314032_df$gene=strsplit2(HUDEP2_GSE314032_df$gene_id,fixed = T,split = '|')[,1]
HUDEP2_GSE314032_df=HUDEP2_GSE314032_df[,-1]
HUDEP2_GSE314032_df=HUDEP2_GSE314032_df[rowSums(HUDEP2_GSE314032_df[,1:12]) >1,]
rownames(HUDEP2_GSE314032_df)=HUDEP2_GSE314032_df$gene
HUDEP2_GSE314032_df=HUDEP2_GSE314032_df[,1:12]
HUDEP2_GSE314032_df=log2(cpm(HUDEP2_GSE314032_df)+1)
HUDEP2_GSE314032_df=data.frame(HUDEP2_GSE314032_df)
c(fetal_markers,'HBD',adult_markers)[!c(fetal_markers,'HBD',adult_markers) %in% rownames(HUDEP2_GSE314032_df)]
HUDEP2_GSE314032_df[c("DLK1","MEG3","HMGA2","TNFSF13B"),1:dim(HUDEP2_GSE314032_df)[2]]=0




all_df=cbind(HEL_GSE203060_df[c(fetal_markers,'HBD',adult_markers),],K562_GSE311284_df[c(fetal_markers,'HBD',adult_markers),])
all_df=cbind(all_df,HUDEP2_GSE314032_df[c(fetal_markers,'HBD',adult_markers),])
colnames(all_df)=gsub(pattern ='HEL.' ,replacement = 'H',colnames(all_df))
colnames(all_df)=gsub(pattern ='SCR.' ,replacement = '',colnames(all_df))
colnames(all_df)=gsub(pattern ='SG2_' ,replacement = 'S2.',colnames(all_df))
sample_df=data.frame(row.names = colnames(all_df),celline=c(rep('HEL',dim(HEL_GSE203060_df)[2]),rep('K562',dim(K562_GSE311284_df)[2]),rep('HUDEP2',dim(HUDEP2_GSE314032_df)[2]) ),
                     source=c(rep('GSE203060',dim(HEL_GSE203060_df)[2]),rep('GSE311284',dim(K562_GSE311284_df)[2]),rep('GSE314032',dim(HUDEP2_GSE314032_df)[2])) )
sample_df$celline <- factor(sample_df$celline)
sample_df$source <- factor(sample_df$source)
annotation_colors <- list(
  celline = c(
    "HEL" = "#4DBBD5",
    "K562" = "#E64B35", 
    "HUDEP2" = "#3C5488"
  ),
  source = c(
    "GSE203060" = "#F39B7F",
    "GSE311284" = "#8491B4",
    "GSE314032" = "#00A087"
  ),
  type = c(
    "fetal" = "#F39B7F",
    "adult" = "#8491B4"
  )
)

library(ComplexHeatmap)
p=pheatmap(as.matrix(all_df),color=colorRampPalette(colors = c('navy','white','firebrick3'))(100),annotation_row =an_df,annotation_col =sample_df,annotation_colors = annotation_colors,cluster_cols = F,cluster_rows = F)
ggsave(as.ggplot(p),filename='../Protein_NRBC_marker/res_pic/main_figure6/cellline_fetal_marker_expression_heatmap.pdf',width = 8,height =8)


############################################################################################################################################################################
#-------------------------------------validated the DLK1 marker-------------------------#
############################################################################################################################################################################
####--------------------------sorted CD2351+DRAQ5+  NRBC----------------------------# 



smartseq_NRBC_exp_df=readRDS('../zx_lab_NRBC/experiment_data/normalized_NRBC_smartseq2_genecount.rds')
pheatmap(smartseq_NRBC_exp_df[c(fetal_markers,'HBD',adult_markers)[c(fetal_markers,'HBD',adult_markers) %in% rownames(smartseq_NRBC_exp_df) ],1:6],cluster_rows = F,cluster_cols = F)

pheatmap(smartseq_NRBC_exp_df[c(fetal_markers,'HBD',adult_markers)[c(fetal_markers,'HBD',adult_markers) %in% rownames(smartseq_NRBC_exp_df) ],2:5],cluster_rows = F,cluster_cols = F)

smartseq_NRBC_exp_df2=data.frame(smartseq_NRBC_exp_df[,'mPBMC2_NRBC']);colnames(smartseq_NRBC_exp_df2)='mPBMC1'
smartseq_NRBC_exp_df2$mPBMC2=rowMeans(smartseq_NRBC_exp_df[,c('mPBMC1_CD81H', 'mPBMC1_CD81L', 'mPBMC1_CD81M')])
smartseq_NRBC_exp_df2=smartseq_NRBC_exp_df2[rowSums(smartseq_NRBC_exp_df2) >0,]

markers=c(fetal_markers,'HBD',adult_markers)[c(fetal_markers,'HBD',adult_markers) %in% rownames(smartseq_NRBC_exp_df2) ]
p1=pheatmap(smartseq_NRBC_exp_df2[markers,],cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))
ggsave(as.ggplot(p1),filename='../Protein_NRBC_marker/res_pic/main_figure6/final_candidated_marker_mPBMC_NRBC_heatmap.pdf',width = 3,height =5)


Y_spefici_genes=c('SRY','ZFY','TSPY','DYS14','AMELY','RPS4Y1','TTTY15','DAZ1','RBMY','CDY','BPY2','HSFY','XKRY','PRY', 'USP9Y', 'UTY', 'KDM5D','ZFX')
p2=pheatmap(smartseq_NRBC_exp_df2[Y_spefici_genes[Y_spefici_genes %in% rownames(smartseq_NRBC_exp_df2) ],],cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))
ggsave(as.ggplot(p2),filename='../Protein_NRBC_marker/res_pic/main_figure6/Y_sex_gene_mPBMC_NRBC_heatmap.pdf',width = 3,height =5)


# new-------
samrtseq2_df=read.csv('zx_lab_NRBC/experiment_data/result_20260708_ycl_FACS_NRBC_smartseq2_genecount_new.csv',sep=',')
samrtseq2_df=samrtseq2_df[rowSums(samrtseq2_df[,3:dim(samrtseq2_df)[2]]) >0,]
du_samrtseq2_df=samrtseq2_df[duplicated(samrtseq2_df$gene_symbol),]

du_tmp_assay_tmp=data.frame()
for(symbol in unique(du_samrtseq2_df$gene_symbol)){
  temp=du_samrtseq2_df[du_samrtseq2_df$gene_symbol == symbol,3:7]
  temp=t(data.frame(colMaxs(as.matrix(temp))))
  rownames(temp)=symbol
  du_tmp_assay_tmp=rbind(du_tmp_assay_tmp,temp)
}



samrtseq2_df=samrtseq2_df[!duplicated(samrtseq2_df$gene_symbol),]
rownames(samrtseq2_df)=samrtseq2_df$gene_symbol
samrtseq2_df=samrtseq2_df[,3:7]
samrtseq2_df=rbind(samrtseq2_df,du_tmp_assay_tmp)
smaple_inf=c('Y1_AApos.2','Y1_AAneg','Y1_FDApos','Y1_FDApos.2','Y1_FDAneg')
colnames(samrtseq2_df)=smaple_inf

library(edgeR)
nr_samrtseq2_df=log2(edgeR::cpm(samrtseq2_df)+1)

cho_gene_samrtseq2_df=nr_samrtseq2_df[rownames(nr_samrtseq2_df) %in% c(top_fetal_marker,adult_markers,'HBD'),]
cho_gene_samrtseq2_df=cho_gene_samrtseq2_df[c('Y1_FDApos.2', 'Y1_FDAneg'),]
pheatmap(cho_gene_samrtseq2_df,cluster_cols = F,col=colorRampPalette(colors = c('navy','white','firebrick3'))(100))

DLK_pos_NRBC_samrtseq2_df=nr_samrtseq2_df[,c('Y1_FDApos','Y1_FDApos.2','Y1_FDAneg')]
DLK_pos_NRBC_samrtseq2_df=DLK_pos_NRBC_samrtseq2_df[rowSums(DLK_pos_NRBC_samrtseq2_df) >0,]
colnames(DLK_pos_NRBC_samrtseq2_df)

DLK1_pos_df=data.frame(cho_gene_samrtseq2_df[,c('Y1_FDApos.2','Y1_FDAneg')])
colnames(DLK1_pos_df)=c('DLK1_pos_rep1','DLK1_pos_rep2')
pheatmap(DLK1_pos_df[rowSums(DLK1_pos_df) >0,],cluster_cols = F,col=colorRampPalette(colors = c('navy','white','firebrick3'))(100))


Y_spefici_genes=c('SRY','ZFY','TSPY','DYS14','AMELY','RPS4Y1','TTTY15','DAZ1','RBMY','CDY','BPY2','HSFY','XKRY','PRY', 'USP9Y', 'UTY', 'KDM5D','ZFX')
rownames(nr_samrtseq2_df)[rownames(nr_samrtseq2_df) %in% Y_spefici_genes] # 测序数据量太低了，


