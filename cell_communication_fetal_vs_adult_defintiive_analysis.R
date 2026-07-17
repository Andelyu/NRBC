
################################################################################################################################################################
#---------------------------------------------------------prepare data---------------------------------------------------------------#
################################################################################################################################################################

library(Seurat)
library(CellChat)
library(nichenetr)
library(clusterProfiler)
library(reshape)
library(ggplot2)
library(patchwork )
library(pheatmap)
library(network)
library(ggnetwork)

library(RColorBrewer )
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]

setwd( "2021_NRBC_chlyu/NRBC_altas_CC" )
NRBC_subcelltype=c("BFUE/CFUE","ProE","Bas","Poly","Orth" )
 
NRBC_altas_LR_df=read.csv('res_data/filt_NRBC_altas_LR_df.csv',sep="\t")
Other2Ery_df=read.csv('res_data/filt_Other2Ery_df_new.csv',sep="\t")
all_NRBC_receptor_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']),'receptor.symbol'] # 取receptor gene
all_NRBC_ligand_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in%unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']),'ligand.symbol']# 取ligand gene
all_NRBC_receptor_genes=sort(unique(unlist(strsplit(all_NRBC_receptor_genes,split=','))))
all_NRBC_ligand_genes=sort(unique(unlist(strsplit(all_NRBC_ligand_genes,split=','))))
length(all_NRBC_receptor_genes);length(all_NRBC_ligand_genes)# 49， 80

filt_NBRC_altas_seu=readRDS('../20251125_filt_NBRC_altas_seu.rds')



# primary tissue data------------------#
BM_nRBC_MS_df=read.csv('../ref_data/Protein_NRBC/NC_2018_BM_hemo_celltype_LF_MS_data.tsv',header = T,sep="\t",skip = 1)
BM_nRBC_MS_df=BM_nRBC_MS_df[,c(colnames(BM_nRBC_MS_df)[1:4],'ERP.number.of.donors','ERP.normalized.LF.sum')]
BM_nRBC_MS_df=BM_nRBC_MS_df[!is.na(BM_nRBC_MS_df$ERP.normalized.LF.sum),]
BM_Ery_MS_protein_LF_genes=unique(unlist(strsplit(BM_nRBC_MS_df$gene.name[BM_nRBC_MS_df$ERP.normalized.LF.sum >0.1],split = ';')));length(BM_Ery_MS_protein_LF_genes)

adult_HSPC_Ery_MS_protein_df=read.table(file = '../Protein_NRBC_marker/recent_MS_Ery_protein_omics/dealt_adult_HSPC_Ery_MS_protein_df.tsv',sep="\t")
adult_HSPC_Ery_MS_protein_genes=unique(unlist(strsplit(adult_HSPC_Ery_MS_protein_df$Gene.names,split = ';')));length(adult_HSPC_Ery_MS_protein_genes)

############################################################################################################################################
#------------part1: analysis  LR in fetal vs adult definitive NRBC ---------------------------#
############################################################################################################################################
#------------part1.1: analysis conserved LR in fetal vs adult definitive NRBC ---------------------------#

conserved_LR_genes=readRDS('conserved_LR_genes.rds')
definitive_conserved_LR_genes=readRDS('definitive_conserved_LR_genes.rds')
pd_candidated_lr_genes=readRDS('res_data/pd_candidated_lr_genes.rds')
pd_hDEG_LR_genes=as.character(pd_candidated_lr_genes$gene[pd_candidated_lr_genes$cluster=='definitive'])
new_definitive_conserved_LR_genes=readRDS('res_data/new_definitive_conserved_LR_genes.rds')
fetal_adult_NRBC_whole_marker=readRDS('../Protein_NRBC_marker/res_data/fetal_adult_NRBC_whole_marker.rds')

temp_deg_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$avg_log2FC>0,]
temp_deg_df1=temp_deg_df[temp_deg_df$gene %in% conserved_LR_genes, ]
temp_deg_df1[conserved_LR_genes[!conserved_LR_genes %in% temp_deg_df1$gene] ,1:7]=0
temp_deg_df1[conserved_LR_genes[!conserved_LR_genes %in% temp_deg_df1$gene] ,'gene']=conserved_LR_genes[!conserved_LR_genes %in% temp_deg_df1$gene]  

temp_deg_df2=temp_deg_df[temp_deg_df$gene %in% pd_hDEG_LR_genes, ]
temp_deg_df2[pd_hDEG_LR_genes[!pd_hDEG_LR_genes %in% temp_deg_df2$gene] ,1:7]=0
temp_deg_df2[pd_hDEG_LR_genes[!pd_hDEG_LR_genes %in% temp_deg_df2$gene] ,'gene']=pd_hDEG_LR_genes[!pd_hDEG_LR_genes %in% temp_deg_df2$gene]  

temp_deg_df3=temp_deg_df[temp_deg_df$gene %in% new_definitive_conserved_LR_genes, ]
temp_deg_df3[new_definitive_conserved_LR_genes[!new_definitive_conserved_LR_genes %in% temp_deg_df3$gene] ,1:7]=0
temp_deg_df3[new_definitive_conserved_LR_genes[!new_definitive_conserved_LR_genes %in% temp_deg_df3$gene] ,'gene']=new_definitive_conserved_LR_genes[!new_definitive_conserved_LR_genes %in% temp_deg_df3$gene]  


temp_deg_df1=temp_deg_df1[order(temp_deg_df1$avg_log2FC,decreasing = T),]
temp_deg_df2=temp_deg_df2[order(temp_deg_df2$avg_log2FC,decreasing = T),]
temp_deg_df3=temp_deg_df3[order(temp_deg_df3$avg_log2FC,decreasing = T),]
temp_deg_df=rbind(temp_deg_df1,rbind(temp_deg_df2,temp_deg_df3))

temp_deg_df$gene=factor(temp_deg_df$gene,levels = unique(temp_deg_df$gene))
p1=ggplot(temp_deg_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+NoLegend()+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()

filt_NBRC_altas_seu$fa_type='fetal'
filt_NBRC_altas_seu$fa_type[filt_NBRC_altas_seu$tissue_stage=='ABM']='adult'
p3=VlnPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),features = temp_deg_df$gene,group.by = 'final_celltype',split.by = 'fa_type',stack = T)

p=p1+p3+plot_layout(ncol = 1,heights = c(0.6,1.2));p
ggsave(p,filename='res_pic/fa_definitive_conserved_LR_gene_expression.pdf',width = 30,height = 10)

sub_fetal_adult_all_Ery_tissue_markers=read.csv('../Protein_NRBC_marker/DE_marker/fetal_adult_all_Ery_RNA_markers.csv')
sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[,-1]

sub_fetal_adult_all_Ery_tissue_markers$celltype=factor(sub_fetal_adult_all_Ery_tissue_markers$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))
sub_temp_deg_df=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC>0 &  sub_fetal_adult_all_Ery_tissue_markers$gene %in% temp_deg_df$gene,]
sub_temp_deg_df$gene=factor(sub_temp_deg_df$gene,levels = temp_deg_df$gene)
sub_temp_deg_df$type='ligand'
sub_temp_deg_df$type[sub_temp_deg_df$gene %in% all_NRBC_receptor_genes]='receptor'
sub_temp_deg_df$type[sub_temp_deg_df$gene %in% intersect(all_NRBC_receptor_genes,all_NRBC_ligand_genes)]='Both'
p2=ggplot(sub_temp_deg_df,aes(x=gene,y=celltype,col=cluster,size=avg_log2FC,alpha=avg_log2FC,shape=type))+geom_point(stat = 'identity')+theme_classic()+
  theme(axis.text.x = element_blank(),legend.position = 'bottom')+ggtitle('defintive conserved LR genes')+xlab(label = '')
p2

# 
focused_fa_conserved_LR_genes=c('TGM2','EPOR','NR1H2','HLA-A','HLA-B','HLA-C','HLA-E','HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','CD74','LGALS9','PTPRC','IL1B')
p3=VlnPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),features = focused_fa_conserved_LR_genes,group.by = 'final_celltype',split.by = 'fa_type',stack = T)
ggsave(p3,filename='res_pic/focused_fa_definitive_conserved_LR_gene_expression.pdf',width = 8,height = 6)

key_ligand_fa_conserved_LR_genes=c('EPOR','IGF1R','IFNGR1','IFNGR2','TNFRSF1A','LTBR')
p=VlnPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),features = key_ligand_fa_conserved_LR_genes,group.by = 'final_celltype',split.by = 'fa_type',stack = T)


p=DotPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),features = key_ligand_fa_conserved_LR_genes,group.by = 'source_celltype',cols  = c('white','firebrick3'))+RotatedAxis()
p

#################-----------------------part2 heck the expression level in DEGs & canidated cmDEGs-------------------------------------------------------#

# 原条件筛选得到LR gene 太少了，由于adult基因中绝大多数是ealry阶段表达，整体比较时候，容易筛掉很多基因,放宽设置条件
cadidated_fa_markers=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$avg_log2FC >2 & fetal_adult_NRBC_whole_marker$pct.2<0.1 & fetal_adult_NRBC_whole_marker$pct.1>0.05,]
fa_candidated_lr_genes=cadidated_fa_markers[cadidated_fa_markers$gene %in% all_NRBC_LR_genes,]
fa_candidated_lr_genes=fa_candidated_lr_genes[!fa_candidated_lr_genes$gene %in% definitive_conserved_LR_genes,]
length(unique(fa_candidated_lr_genes$gene))
fa_candidated_lr_genes=fa_candidated_lr_genes[order(fa_candidated_lr_genes$cluster,fa_candidated_lr_genes$avg_log2FC,decreasing = T),]
fa_candidated_lr_genes$gene=factor(fa_candidated_lr_genes$gene,levels = fa_candidated_lr_genes$gene)
p1=ggplot(fa_candidated_lr_genes,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()
p1

temp_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% key_ligand_fa_conserved_LR_genes,]
temp_df=temp_df[temp_df$avg_log2FC >0,]
temp_df=temp_df[order(temp_df$avg_log2FC,decreasing = T),]
fa_candidated_lr_genes1=rbind(fa_candidated_lr_genes,temp_df)
fa_candidated_lr_genes1$gene=factor(fa_candidated_lr_genes1$gene,levels = c(levels(fa_candidated_lr_genes$gene),key_ligand_fa_conserved_LR_genes))
p1=ggplot(fa_candidated_lr_genes1,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()
p1



if(F){
  sub_fetal_adult_all_Ery_tissue_markers$celltype=factor(sub_fetal_adult_all_Ery_tissue_markers$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))
  
  sub_temp_deg_df=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC>0 &  sub_fetal_adult_all_Ery_tissue_markers$gene %in% unique(fa_candidated_lr_genes$gene),]
  sub_temp_deg_df$gene=factor(sub_temp_deg_df$gene,levels = fa_candidated_lr_genes$gene)
  sub_temp_deg_df$type='ligand'
  sub_temp_deg_df$type[sub_temp_deg_df$gene %in%unique(unlist(strsplit(CellChatDB.human$interaction$receptor.symbol,split = ', ')))  ]='receptor'
  both_genes=intersect(unique(unlist(strsplit(CellChatDB.human$interaction$ligand.symbol,split = ', '))),unique(unlist(strsplit(CellChatDB.human$interaction$receptor.symbol,split = ', '))))
  sub_temp_deg_df$type[sub_temp_deg_df$gene %in% both_genes  ]='Both'
  
  p2=ggplot(sub_temp_deg_df,aes(x=gene,y=celltype,col=cluster,size=avg_log2FC,alpha=avg_log2FC,shape=type))+geom_point(stat = 'identity')+theme_classic()+
    theme(legend.position = 'bottom')+ggtitle('fetal vs  ABM NRBC hDEGs-LR ')
  p2
  ggsave(p,filename='res_pic/main_figure5/fa_hDEG_LR_expression_DE_substage.pdf',width =6 ,height = 3)
  
}



p3=DotPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),group.by = 'source_celltype',scale = F,features =c(levels(fa_candidated_lr_genes$gene),key_ligand_fa_conserved_LR_genes))+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')

fa_candidated_lr_genes$MS='no_detected'
#fa_candidated_lr_genes$MS[fa_candidated_lr_genes$gene %in% adult_HSPC_Ery_MS_protein_genes]='aHSPC_nRBC_MS'
fa_candidated_lr_genes$MS[fa_candidated_lr_genes$gene %in% BM_nRBC_MS_df$Gene.Name]='ABM_nRBC_MS'
fa_candidated_lr_genes$MS[fa_candidated_lr_genes$gene %in% BM_Ery_MS_protein_LF_genes]='ABM_nRBC_MS'
fa_candidated_lr_genes$MS=factor(fa_candidated_lr_genes$MS,levels = c('ABM_nRBC_MS','no_detected'))

fa_candidated_lr_genes$HPA_IHC='.'
fa_candidated_lr_genes$HPA_IHC[fa_candidated_lr_genes$gene %in% c('PF4','ROBO2','DHCR7','TNFSF13B','ANXA1','PLXND1','HLA-DQB1','HLA-DRA','HLA-DRB1','HLA-DRB5','LGALS9') ]='positive' 
fa_candidated_lr_genes$HPA_IHC[fa_candidated_lr_genes$gene %in% c('APOA1','APOA2','DLK1','PTGR1','SEMA7A','C4A') ]='no_detected/negtive'
fa_candidated_lr_genes$HPA_IHC[fa_candidated_lr_genes$gene %in% c('WNT5B','RELN','ADORA2B','CXCR4') ]='not_available'
fa_candidated_lr_genes$HPA_IHC=factor(fa_candidated_lr_genes$HPA_IHC,levels = c('positive','no_detected/negtive','not_available','.'))
saveRDS(fa_candidated_lr_genes,file = 'res_data/fa_candidated_lr_genes.rds')

p4=ggplot(fa_candidated_lr_genes,aes(x=gene,y='MS',shape=MS))+geom_point(stat = 'identity')+theme_classic()+theme(axis.text.x = element_blank())+
  xlab(label = '')+ylab(label = '')+  scale_shape_manual(values = c(3, 1,4))

p5=ggplot(fa_candidated_lr_genes)+geom_point(aes(x=gene,y='HPA_IHC',shape=HPA_IHC),size=2)+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab(label = '')+scale_shape_manual(values = c(3, 1,4,20))

p=wrap_plots(p1,p3,heights = c(0.6,1.2));p
ggsave(p,filename='res_pic/main_figure5/fa_hDEG_LR_expression.pdf',width =14 ,height = 8)

#YS、FL：APOA1-/2-,PF4+-, FBM:C4A+
#  fetal:DLK1-,ROBO2+,WNT5B,RELN,PTGR1-,DHCR24-, ADORA2B,SEMA7A-,DHCR7+,PLXND1+
#adult:TNFSF13B+,ANXA1+,HLA-DQB1,
#candidated : CXCR4


#table(pd_candidated_lr_genes[,c('cluster','gene_type')])
if(F){
  fa_fetal_LR_gene_enrichGO_res=enrichGO(gene =fa_candidated_lr_genes$gene[fa_candidated_lr_genes$cluster=='fetal'],OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'MF' )
  saveRDS(fa_fetal_LR_gene_enrichGO_res,file = 'fa_fetal_LR_gene_enrichGO_res.rds')
  cnetplot(fa_fetal_LR_gene_enrichGO_res,30)
  
  fetal_cho_show_enrich_GO_df=fa_fetal_LR_gene_enrichGO_res@result # no DLK1
  cho_pathways=c('oxidoreductase activity, acting on the CH-CH group of donors, NAD or NADP as acceptor','peptide antigen binding','growth factor activity',
                 'sterol binding','virus receptor activity','heat shock protein binding','lipid transfer activity','intramolecular oxidoreductase activity')
  
  cho_show_enrich_GO_df=cho_show_enrich_GO_df[match(cho_pathways ,cho_show_enrich_GO_df$Description ),]
  cho_show_enrich_GO_df$Description=factor(cho_show_enrich_GO_df$Description,levels = cho_pathways)
  cho_show_enrich_GO_df$geneID=factor(cho_show_enrich_GO_df$geneID,levels =cho_show_enrich_GO_df$geneID)
  p1=ggplot(cho_show_enrich_GO_df,aes(y=Description,x=-log10(p.adjust) ))+geom_bar(stat = 'identity',fill='#00c6d9',color='black')+theme_classic()+
    geom_text(aes(y=Description, label=geneID),x =1.1,size=3 )
  
  
  fa_adult_LR_gene_enrichGO_res=enrichGO(gene =fa_candidated_lr_genes$gene[fa_candidated_lr_genes$cluster=='adult'],OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'MF' )
  cnetplot(fa_adult_LR_gene_enrichGO_res,30)
  
  saveRDS(fa_fetal_LR_gene_enrichGO_res,file = 'fa_adult_LR_gene_enrichGO_res.rds')
  
}

fa_all_NRBC_mDEGs_inLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$receptor.symbol %in% fa_candidated_lr_genes$gene[fa_candidated_lr_genes$gene %in% all_NRBC_receptor_genes ]])
fa_all_NRBC_mDEGs_outLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$ligand.symbol %in% fa_candidated_lr_genes$gene[fa_candidated_lr_genes$gene %in% all_NRBC_ligand_genes ]])
fa_all_NRBC_mDEGs_inLRs=fa_all_NRBC_mDEGs_inLRs[fa_all_NRBC_mDEGs_inLRs %in% NRBC_altas_LR_df$interaction_name]
fa_all_NRBC_mDEGs_outLRs=fa_all_NRBC_mDEGs_outLRs[fa_all_NRBC_mDEGs_outLRs %in% NRBC_altas_LR_df$interaction_name]
fa_candidated_lr_genes$gene[!fa_candidated_lr_genes$gene %in% c(CellChatDB.human$interaction$ligand.symbol ,CellChatDB.human$interaction$receptor.symbol )]# No factor

fa_candidated_lr_genes$gene_type='receptor'
fa_candidated_lr_genes$gene_type[fa_candidated_lr_genes$gene %in%all_NRBC_ligand_genes ]='ligand'

fa_NRBC_mDEGs_inLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$receptor.symbol %in% fa_candidated_lr_genes$gene])
fa_NRBC_mDEGs_inLRs=fa_NRBC_mDEGs_inLRs[fa_NRBC_mDEGs_inLRs %in% NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']]
fa_all_NRBC_mDEGs_LRs_df=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% fa_NRBC_mDEGs_inLRs & NRBC_altas_LR_df$target_type!='Ery2Other',]

fa_NRBC_mDEGs_outLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$ligand.symbol %in% fa_candidated_lr_genes$gene])
fa_NRBC_mDEGs_outLRs=fa_NRBC_mDEGs_outLRs[fa_NRBC_mDEGs_outLRs %in% NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Other2Ery']]
fa_all_NRBC_mDEGs_LRs_df=rbind(fa_all_NRBC_mDEGs_LRs_df,NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% fa_NRBC_mDEGs_outLRs & NRBC_altas_LR_df$target_type!='Other2Ery',])

# 获得new_definitive_conserved_LR 信息
new_definitive_conserved_LR_genes1=new_definitive_conserved_LR_genes[!new_definitive_conserved_LR_genes %in% fa_candidated_lr_genes$gene]
new_conserved_fa_NRBC_mDEGs_inLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$receptor.symbol %in% new_definitive_conserved_LR_genes1])
new_conserved_fa_NRBC_mDEGs_inLRs=new_conserved_fa_NRBC_mDEGs_inLRs[new_conserved_fa_NRBC_mDEGs_inLRs %in% NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']]
new_conserved_fa_all_NRBC_mDEGs_LRs_df1=NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% new_conserved_fa_NRBC_mDEGs_inLRs & NRBC_altas_LR_df$target_type!='Ery2Other',]

new_conserved_fa_NRBC_mDEGs_outLRs=unique(CellChatDB.human$interaction$interaction_name[ CellChatDB.human$interaction$ligand.symbol %in% new_definitive_conserved_LR_genes1])
new_conserved_fa_NRBC_mDEGs_outLRs=new_conserved_fa_NRBC_mDEGs_outLRs[new_conserved_fa_NRBC_mDEGs_outLRs %in% NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']]
new_conserved_fa_all_NRBC_mDEGs_LRs_df1=rbind(new_conserved_fa_all_NRBC_mDEGs_LRs_df1,NRBC_altas_LR_df[NRBC_altas_LR_df$interaction_name %in% new_conserved_fa_NRBC_mDEGs_outLRs & NRBC_altas_LR_df$target_type=='Ery2Other',])

defintive_specific_conserved_LR=readRDS('pd_defintive_specific_conserved_LR.rds')
new_conserved_fa_all_NRBC_mDEGs_LRs_df1=new_conserved_fa_all_NRBC_mDEGs_LRs_df1[!new_conserved_fa_all_NRBC_mDEGs_LRs_df1$interaction_name %in% defintive_specific_conserved_LR,]
new_conserved_fa_all_NRBC_mDEGs_LRs_df1=new_conserved_fa_all_NRBC_mDEGs_LRs_df1[!new_conserved_fa_all_NRBC_mDEGs_LRs_df1$interaction_name %in% fa_all_NRBC_mDEGs_LRs_df$interaction_name, ]
new_conserved_fa_all_NRBC_mDEGs_LRs_df1$type='01_new_conserved'


fa_all_NRBC_mDEGs_LRs_df$type='02_hDEG_LR'
fa_all_NRBC_mDEGs_LRs_df1=rbind(fa_all_NRBC_mDEGs_LRs_df,new_conserved_fa_all_NRBC_mDEGs_LRs_df1)

fa_all_NRBC_mDEGs_LRs_df1$target_type=factor(fa_all_NRBC_mDEGs_LRs_df1$target_type,levels = c('Other2Ery','Ery2Ery','Ery2Other'))
fa_all_NRBC_mDEGs_LRs_df1$cluster='fetal'
fa_all_NRBC_mDEGs_LRs_df1$cluster[fa_all_NRBC_mDEGs_LRs_df1$stage=='ABM']='adult'
fa_all_NRBC_mDEGs_LRs_df1=fa_all_NRBC_mDEGs_LRs_df1[fa_all_NRBC_mDEGs_LRs_df1$stage!='YS',]

fa_all_NRBC_mDEGs_LRs_df1=fa_all_NRBC_mDEGs_LRs_df1[order(fa_all_NRBC_mDEGs_LRs_df1$type,fa_all_NRBC_mDEGs_LRs_df1$cluster,fa_all_NRBC_mDEGs_LRs_df1$target_type,fa_all_NRBC_mDEGs_LRs_df1$ligand,fa_all_NRBC_mDEGs_LRs_df1$annotation),]

fa_all_NRBC_mDEGs_LRs_df1$interaction_name=factor(fa_all_NRBC_mDEGs_LRs_df1$interaction_name,levels =unique( fa_all_NRBC_mDEGs_LRs_df1$interaction_name))
p=ggplot(fa_all_NRBC_mDEGs_LRs_df1,aes(x=celltype ,y=interaction_name,size=prob ,color=annotation,shape=target_type))+geom_point(alpha=0.7)+facet_grid(~stage )+
  theme_classic()+scale_color_manual(values = cols[-2])+scale_shape_manual(values = c(0:2))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+facet_grid(~stage )+ggtitle('fetal vs adult LR based on hDEGs' )
p
ggsave(p,filename='res_pic/main_figure5/fa_fetal_adult_specific_LR_dotplot.pdf',width = 8,height = 10)

saveRDS(fa_all_NRBC_mDEGs_LRs_df1,file = 'res_data/fa_all_NRBC_mDEGs_LRs_df1.rds')



#cho_fetal_LR_interactions=fa_all_NRBC_mDEGs_LRs_df$interaction_name[grep('DLK1|WNT5B|RELN|PTGR1|ADORA2B|TNFSF13B|ANXA1|CXCR4|LGALS9|PLXND1|TBXAS1|ADGRG1|IL1B|SULT1A1',fa_all_NRBC_mDEGs_LRs_df$interaction_name)]
#ggplot(fa_all_NRBC_mDEGs_LRs_df[fa_all_NRBC_mDEGs_LRs_df$interaction_name %in% cho_fetal_LR_interactions, ],aes(x=celltype ,y=interaction_name,size=prob ,color=annotation,shape=target_type))+geom_point(alpha=0.7)+facet_grid(~stage )+theme_classic()+scale_color_manual(values = cols[-2])+
#  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+facet_grid(~stage )+ggtitle('primitive vs definitive LR based on FA_hDEGs-LR' )

netAnalysis_signalingRole_network(object.list2[['FL']], signaling = 'ncWNT', width = 36, height = 2.5, font.size = 10) #Endo & stroma cells
netAnalysis_signalingRole_network(object.list2[['FBM']], signaling = 'ncWNT', width = 40, height = 2.5, font.size = 10) # Endo & fibroblast

netAnalysis_signalingRole_network(object.list2[['ABM']], signaling = 'ANNEXIN', width = 40, height = 2.5, font.size = 10) # mono,myeloid cell ,Neu
netAnalysis_signalingRole_network(object.list2[['ABM']], signaling = 'BAFF', width = 40, height = 2.5, font.size = 10)# plasma
netAnalysis_signalingRole_network(object.list2[['ABM']], signaling = 'GALECTIN', width = 40, height = 2.5, font.size = 10)# sorts of immune cells



# check interaction mediated by MHC of NRBC
MHC_df3=NRBC_altas_LR_df[grep('HLA-',NRBC_altas_LR_df$interaction_name),]
MHC_df3$target=as.character(MHC_df3$target)
MHC_df3$source=factor(MHC_df3$source,levels = NRBC_subcelltype)
MHC_df3$target[MHC_df3$target=='PDC']='pDC'
#MHC_df3$target=gsub(pattern = 'CYCLING_',replacement ='' ,MHC_df3$target)
MHC_df3$target[MHC_df3$target=='CD8+ T-Cell']='CD8+T'
MHC_df3$target[MHC_df3$target=='CD4+ T-Cell']='CD4+T'
MHC_df3$target[MHC_df3$target=='TREG']='Treg'
MHC_df3$target[grep('MONOCYTE',  MHC_df3$target)]='MONOCYTE'
MHC_df3$target[grep('MACROPHAGE',MHC_df3$target)]='MACROPHAGE'
MHC_df3$target[grep('^DC|AS_DC', MHC_df3$target)]='DC'
MHC_df3$target[grep('CD8 T', MHC_df3$target)]='CD8+T'
MHC_df3$target[grep('CD4 T', MHC_df3$target)]='CD4+T'
MHC_df3$target[grep('NK', MHC_df3$target)]='NK'
MHC_df3$target[grep('ILC3', MHC_df3$target)]='ILC'
MHC_df3=MHC_df3[-grep('^CYCLING',MHC_df3$target),]

celltype_level=c("CD4+T","GZMK cytotoxic CD4 T", "Memory CD4 T" ,"CD8+T","Naive CD8 T",'ILC',"TYPE_1_INNATE_T"  , "GZMB CD8 T","GZMK CD8 T", "Treg",
                 "MONOCYTE","MACROPHAGE","NK","DC","pDC","OSTEOCLAST","THY1+ MSC","Fibro-MSC","Adipo-MSC","VSMC","ENDOTHELIUM_V"  )
unique(MHC_df3$target)[!unique(MHC_df3$target) %in% celltype_level] 

MHC_df3$target=factor(as.character(MHC_df3$target),levels = celltype_level)
ggplot(MHC_df3,aes(x=source,y=target,color=ligand,size=prob))+geom_point(alpha=0.6,stroke=1,shape=1)+theme_classic()+facet_grid(~stage)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))


############################################################################################################################################################
#-----------------------------part2:using the nichnet to analysis the key LR to target genes-------------------------#
############################################################################################################################################################
library(nichenetr)
library(ggplot2)
library(dplyr)

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")
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

draw_gonetcwork_pic_func=function(res=fa_fetal_LR_targetgene_enrichGO_res,showCategory=20,xlimits = c(-0.2, 2.5)){
  
  res <- as.data.frame(res)
  res <- head(res,showCategory)
  
  # 构建连接
  links <- do.call(rbind, lapply(1:nrow(res), function(i) {
    data.frame(gene = strsplit(res$geneID[i], "/")[[1]], 
               pathway = res$Description[i])}))
  
  # 计算位置
  links$y_gene <- match(links$gene, names(sort(table(links$gene), T))) / (length(unique(links$gene)) + 1)
  links$y_pathway <- match(links$pathway, res$Description[order(res$p.adjust)]) / (nrow(res) + 1)
  
  # 绘图
  p=ggplot() +  geom_segment(aes(x = 0, xend = 0.8, y = y_gene, yend = y_pathway), data = links,   color = "gray70", alpha = 0.4, size = 0.4) +
    geom_point(aes(x = 0, y = y_gene), data = distinct(links, gene, y_gene), size = 2, color = "#4878D0") +
    geom_text(aes(x = -0.05, y = y_gene, label = gene), data = distinct(links, gene, y_gene), hjust = 1, size = 3, color = "black") +
    geom_point(aes(x = 0.8, y = y_pathway), data = distinct(links, pathway, y_pathway), size = 2, color = "#EE854A") +
    geom_text(aes(x = 0.85, y = y_pathway, label = pathway), data = distinct(links, pathway, y_pathway),   hjust = 0, size = 2.5, color = "black") +
    scale_x_continuous(limits =xlimits, breaks = c(0, 1), labels = c("Genes", "Pathways")) +
    scale_y_continuous(breaks = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 10),
          plot.margin = margin(5, 15, 25, 5))
  print(p)
  return(p)
}


filt_NBRC_altas_seu$fa_celltype2=paste(as.character(filt_NBRC_altas_seu$fa_type),as.character(filt_NBRC_altas_seu$final_celltype),sep = "_")
fa_all_mexp_df=as.matrix(AverageExpression(subset(filt_NBRC_altas_seu,tissue_stage !='YS'),group.by = 'fa_celltype2',features =rownames(filt_NBRC_altas_seu) )$RNA)

known_key_regulator_erythropoiesis_genes=c('GATA1','GATA2','KLF1','TAL1','LMO2','FOG1','BCL11A','EPOR','HIF1A','HIF2A','STAT5A','ALAS2','SLC11A2','FOXO3')



# --------compare the LR on primitive and definitve NRBC -------------#
as.character(incoming_LR_type_df_list[['11_YS']])

LR_df=CellChatDB.human$interaction

fetal_NBRC_ligands=unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage %in% c('FL','FBM') & NRBC_altas_LR_df$target_type!='Ery2Other','interaction_name'])
fetal_NBRC_ligands=unique(LR_df$ligand.symbol[LR_df$interaction_name %in%fetal_NBRC_ligands ])
fetal_NBRC_ligands= unique(as.character(t(data.frame(strsplit(fetal_NBRC_ligands,split = ', ')))[,1]))
fetal_NBRC_ligands[!fetal_NBRC_ligands %in% rownames(filt_NBRC_altas_seu)]  

adult_NBRC_ligands=unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage %in% c('ABM') & NRBC_altas_LR_df$target_type!='Ery2Other','interaction_name'])
adult_NBRC_ligands=unique(LR_df$ligand.symbol[LR_df$interaction_name %in%adult_NBRC_ligands ])
adult_NBRC_ligands= unique(as.character(t(data.frame(strsplit(adult_NBRC_ligands,split = ', ')))[,1]))
adult_NBRC_ligands[!adult_NBRC_ligands %in% rownames(filt_NBRC_altas_seu)]  
adult_NBRC_ligands=c(adult_NBRC_ligands,'EPO')

# --------- HSPC-derived : fetal vs adult ------------------#
top_hspc_degs=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$avg_log2FC >1 & fetal_adult_NRBC_whole_marker$pct.1 >0.05 & fetal_adult_NRBC_whole_marker$pct.2 <0.3, ] %>% group_by(cluster) %>% top_n(wt =avg_log2FC,n = 500 )
top_hspc_sub_degs=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >1 & sub_fetal_adult_all_Ery_tissue_markers$pct.1>0.1,]  %>% group_by(cluster,celltype) %>% do(head(.,n=300))
table(top_hspc_sub_degs$cluster);table(top_hspc_degs$cluster);

fetal_target_genes=unique(c(top_hspc_degs$gene[top_hspc_degs$cluster=='fetal'],top_hspc_sub_degs$gene[top_hspc_sub_degs$cluster=='fetal']))
adult_target_genes=unique(c(top_hspc_degs$gene[top_hspc_degs$cluster=='adult'],top_hspc_sub_degs$gene[top_hspc_sub_degs$cluster=='adult']))
length(fetal_target_genes);length(adult_target_genes)


target_genes=fetal_target_genes;
target_genes=target_genes[rowMax(fa_all_mexp_df[target_genes,grep('fetal',colnames(fa_all_mexp_df))])>0.5]
length(target_genes)
potential_ligands=fetal_NBRC_ligands[fetal_NBRC_ligands %in% colnames(ligand_target_matrix)]

fetal_ligand_target_res_list=nichenet_predict_func(legend_title = 'fetal NRBC ligand to target gene score',potential_ligands =potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
fetal_ligand_target_res_list[[2]]
ggsave(fetal_ligand_target_res_list[[2]],filename = 'res_pic/main_figure5/fetal_recept_ligand_target_gene_heatmap_score.pdf',width = 10,height = 10)

saveRDS(fetal_ligand_target_res_list,file = 'res_data/fa_fetal_ligand_target_res_list.rds')

#
an_df=data.frame(fetal_ligand_target_res_list[[1]]$data[fetal_ligand_target_res_list[[1]]$data$y %in% c('EPO','IGF1','IFNG','IGF2','TNF','LTB','LTA'),])
an_df=data.frame(row.names = an_df$y,aupr_score=an_df$score)

p=pheatmap(t(fetal_ligand_target_res_list[[3]][,c('EPO','IGF1','IFNG','IGF2','TNF','LTB','LTA')]),annotation_row = an_df,border_color = 'white',cluster_rows = F,cluster_cols =F ,color =colorRampPalette(colors =  c("white" ,"#EF6548"))(100))
ggsave(as.ggplot(p),file='res_pic/main_figure5/fetal_key_ligand_target_genes_heatmap.pdf',width = 18,height = 3,dpi = 300)

fa_fetal_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(fetal_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
cnetplot(fa_fetal_LR_targetgene_enrichGO_res,showCategory =20,   color_category = "#E41A1C", color_gene = "#377EB8", cex_label_category = 0.6,cex_label_gene = 0.5, node_label = "all",  layout = "fr")
p=draw_gonetcwork_pic_func(res =fa_fetal_LR_targetgene_enrichGO_res,showCategory = 20 )
ggsave(p,filename='res_pic/main_figure5/fetal_target_gene_enrichGO_res.pdf',width =6,height = 6)
saveRDS(fa_fetal_LR_targetgene_enrichGO_res,file = 'res_data/fa_fetal_LR_targetgene_enrichGO_res.rds')


temp_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% rownames(fetal_ligand_target_res_list[[3]]) & fetal_adult_NRBC_whole_marker$avg_log2FC>1,]
#temp_df$gene=factor(temp_df$gene,levels = rownames(fetal_ligand_target_res_list[[3]]))
temp_df=temp_df[temp_df$pct.2<0.1 & temp_df$pct.1>0.1,]
temp_df=temp_df[order(temp_df$avg_log2FC,decreasing = T),]
temp_df$gene=factor(temp_df$gene,levels = temp_df$gene)
p1=ggplot(temp_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()
p2=DotPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),features =as.character(temp_df$gene),scale = F)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
p=p1+p2+plot_layout(heights = c(0.4,1.2),ncol = 1);p
ggsave(p,filename='res_pic/fa_fetal_hDEG_targetgene_expression.pdf',width = 6,height = 6)

cho_fetal_ligand_target_df=fetal_ligand_target_res_list[[3]][,c('EPO','IGF1','IFNG','IGF2','TNF','LTB','LTA')]
cho_fetal_ligand_target_list=list()
cho_fetal_ligand_target_list[['EPO']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'EPO'] >0.04]
cho_fetal_ligand_target_list[['IGF1']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IGF1'] >0.04]
cho_fetal_ligand_target_list[['IFNG']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IFNG'] >0.04]
cho_fetal_ligand_target_list[['TNF']]=names(cho_fetal_ligand_target_df[,'TNF'])[cho_fetal_ligand_target_df[,'TNF'] >0.04]
cho_fetal_ligand_target_list[['LTB']]=names(cho_fetal_ligand_target_df[,'EPO'])[cho_fetal_ligand_target_df[,'LTB'] >0.04]
cho_fetal_ligand_target_list[['LTB']]=names(cho_fetal_ligand_target_df[,'LTA'])[cho_fetal_ligand_target_df[,'LTA'] >0.04]

key_fa_fetal_LR_targetgene_enrichGO_res=enrichGO(gene =unique(unlist(cho_fetal_ligand_target_list)),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
p=draw_gonetcwork_pic_func(res =key_fa_fetal_LR_targetgene_enrichGO_res,showCategory = 20 )
ggsave(p,filename='res_pic/key_fa_fetal_hDEG_targetgene_expression.pdf',width = 6,height = 6)


#----------------------adult ------------------#
target_genes=adult_target_genes;
target_genes=target_genes[rowMax(fa_all_mexp_df[target_genes,grep('adult',colnames(fa_all_mexp_df))])>0.5]
length(target_genes)
potential_ligands=adult_NBRC_ligands[adult_NBRC_ligands %in% colnames(ligand_target_matrix)]
adult_ligand_target_res_list=nichenet_predict_func(legend_title = 'adult NRBC ligand to target gene score',potential_ligands =potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
p=adult_ligand_target_res_list[[2]]
p
ggsave(p,file='res_pic/main_figure5/adult_recept_ligand_target_gene_heatmap_score.pdf',height =10 ,width = 6)
saveRDS(adult_ligand_target_res_list,file = 'res_data/fa_adult_ligand_target_res_list.rds')

#
an_df=data.frame(adult_ligand_target_res_list[[1]]$data[adult_ligand_target_res_list[[1]]$data$y %in% c('EPO','CXCL12'),])
an_df=data.frame(row.names = an_df$y,aupr_score=an_df$score)
p=pheatmap(adult_ligand_target_res_list[[3]][,c('EPO','CXCL12')],annotation_col = an_df,border_color = 'white',cluster_rows = F,cluster_cols =F ,color =colorRampPalette(colors =  c("white" ,"#EF6548"))(100))
ggsave(as.ggplot(p),file='res_pic/main_figure5/adult_key_ligand_target_genes_heatmap.pdf',width = 4,height = 8,dpi = 300)


fa_adult_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(adult_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(fa_adult_LR_targetgene_enrichGO_res,file = 'res_data/fa_adult_LR_targetgene_enrichGO_res.rds')
cnetplot(fa_adult_LR_targetgene_enrichGO_res,20)
p=draw_gonetcwork_pic_func(res =fa_adult_LR_targetgene_enrichGO_res,showCategory = 20 ,xlimits = c(-0.5,2.5))
ggsave(p,filename='res_pic/main_figure5/adult_target_gene_enrichGO_res.pdf',width =6,height = 6)

temp_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% rownames(adult_ligand_target_res_list[[3]]) & fetal_adult_NRBC_whole_marker$avg_log2FC>1,]
#temp_df$gene=factor(temp_df$gene,levels = rownames(fetal_ligand_target_res_list[[3]]))
temp_df=temp_df[temp_df$pct.2<0.1 & temp_df$pct.1>0.1,]
temp_df=temp_df[order(temp_df$avg_log2FC,decreasing = T),]
temp_df$gene=factor(temp_df$gene,levels = temp_df$gene)
p1=ggplot(temp_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()
p2=DotPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),features =as.character(temp_df$gene),scale = F)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
p=p1+p2+plot_layout(heights = c(0.4,1.2),ncol = 1)

ggsave(p,filename='res_pic/main_figure5/fa_adult_hDEG_targetgene_expression.pdf',width = 8,height = 12)


cho_adult_ligand_target_df=adult_ligand_target_res_list[[3]][,c('EPO','CXCL12')]
cho_adult_ligand_target_list=list()
cho_adult_ligand_target_list[['EPO']]=names(cho_adult_ligand_target_df[,'EPO'])[cho_adult_ligand_target_df[,'EPO'] >0.04]
cho_adult_ligand_target_list[['CXCL12']]=names(cho_adult_ligand_target_df[,'CXCL12'])[cho_adult_ligand_target_df[,'CXCL12'] >0.04]
key_fa_adult_LR_targetgene_enrichGO_res=enrichGO(gene =unique(unlist(cho_adult_ligand_target_list)),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
p=draw_gonetcwork_pic_func(res =key_fa_adult_LR_targetgene_enrichGO_res,showCategory = 20 )
ggsave(p,filename='res_pic/main_figure5/key_fa_fetal_hDEG_targetgene_expression.pdf',width = 6,height = 6)




########################################################################################################################################################################
#----------------------part3: EPO、CXCL12、IGF1 、IGF2、IFNG、TNF(LTA\LTB) expression in niches ------------------------------------#
########################################################################################################################################################################


unique(NRBC_altas_LR_df[grep('IGF1',NRBC_altas_LR_df$ligand),'receptor']) # IGF1R
unique(NRBC_altas_LR_df[grep('IFNG',NRBC_altas_LR_df$ligand),'receptor']) # IFNGR1_IFNGR2
unique(NRBC_altas_LR_df[grep('TNF',NRBC_altas_LR_df$ligand),'receptor']) #  "TNFRSF1A"  "LTBR"      "TNFRSF21"  "TNFRSF13B" "TNFRSF13C" "TNFRSF17"
unique(NRBC_altas_LR_df[grep('LTA',NRBC_altas_LR_df$ligand),'receptor']) #  "TNFRSF1A"  "LTBR"   
unique(NRBC_altas_LR_df[grep('LTB',NRBC_altas_LR_df$ligand),'receptor']) #  LTB4R
unique(NRBC_altas_LR_df[grep('CXCL12',NRBC_altas_LR_df$ligand),'receptor']) #  CXCR4


p=DotPlot(subset(filt_NBRC_altas_seu,pd_celltype=='definitive'),scale=F,features = c('EPOR','CXCR4','IFNGR1','IFNGR2','IGF1R','TNFRSF1A',"LTBR",'LTB4R'))+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
ggsave(p,filename='res_pic/main_figure5/the_key_ligand_receptor_expression_in_definitive.pdf',width =6 ,height = 6)

#-----------------------------------------check the key secreted genes expression in FL/FBM/ABM---------------------#
FL_altas_seu=readRDS('../NRBC_FL_altas/tmp_FL_altas_seu.rds')

cho_celtype=c( "MONOCYTE","MACROPHAGE", "NK/T CELLS",'B CELLS','LMPP_MLP','DC' ,'ILC',"HEPATOCYTE","FIBROBLASTS","SMOOTH MUSCLE","SKELETAL MUSCLE","MESOTHELIUM","NEPHRON" )
cho_feature= c('EPO','CXCL12','IGF1','IGF2','IFNG','TNF','LTA','LTB')
p1=VlnPlot(subset(FL_altas_seu,subcelltype %in% cho_celtype),group.by = 'subcelltype',features =cho_feature,stack = T)+NoLegend()+ggtitle('FL ALTAS')
ggsave(p1,filename='res_pic/main_figure5/FL_key_ligand_expression_celltype_niche_vlnplot.pdf',width = 5,height = 5)


FL_altas_seu$age=factor(FL_altas_seu$age,levels = c("CS14_4PCW" ,"CS15_5PCW", "CS17_6PCW", "CS18","CS22","CS23",  "8PCW","8.1PCW","9.1PCW","9.7PCW",  "11PCW","11.4PCW" ,"12PCW","13.9PCW","14.4PCW","15PCW","16.3PCW","16PCW","17PCW"  ))
p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c('MACROPHAGE')),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FL  MACROPHAGE');p
ggsave(p,filename='res_pic/main_figure5/FL_Mac_key_ligand_expression_vlnplot.pdf',height = 6,width = 4)
rm(p);gc()

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c('MONOCYTE')),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FL  MONOCYTE ')
ggsave(p,filename='res_pic/main_figure5/FL_MONOCYTE_key_ligand_expression_vlnplot.pdf',height = 6,width = 4)
p=VlnPlot(subset(FL_altas_seu,anno_lvl_2_final_clean %in% c('MONOCYTE_I_CXCR4','MONOCYTE_II_CCR2','MONOCYTE_III_IL1B','PROMONOCYTE')),
        split.by = 'anno_lvl_2_final_clean',group.by = 'age',features =cho_feature[-1],stack = T)+ggtitle('FL  MONOCYTE ')
ggsave(p,filename='res_pic/main_figure5/FL_Mono_subcelltype_key_ligand_expression_vlnplot.pdf',height = 8,width = 6)

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("NK/T CELLS")),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FL  NK/T CELLS');p
ggsave(p,filename='res_pic/main_figure5/FL_NKT_key_ligand_expression_vlnplot.pdf',height = 6,width = 4)
rm(p);gc()

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("HEPATOCYTE")),group.by = 'age',features =cho_feature,stack = T)+NoLegend()+ggtitle('HEPATOCYTE');p
ggsave(p,filename='res_pic/main_figure5/FL_HEPATOCYTE_key_ligand_expression_vlnplot.pdf',height = 6,width =4)
rm(p);gc()

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("FIBROBLASTS")),group.by = 'age',features =c( "CXCL12", "IGF1","IGF2" ,  "TNF"  ,"LTB"   ),stack = T)+NoLegend()+ggtitle('FIBROBLASTS');p
ggsave(p,filename='res_pic/main_figure5/FL_FIBROBLASTS_key_ligand_expression_vlnplot.pdf',height = 6,width = 4)
rm(p);gc()

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("SMOOTH MUSCLE","SKELETAL MUSCLE")),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('MUSCLE');p
ggsave(p,filename='res_pic/main_figure5/FL_MUSCLE_key_ligand_expression_vlnplot.pdf',height = 6,width = 4)
rm(p);gc()


p=VlnPlot(subset(FL_altas_seu,anno_lvl_2_final_clean %in% unique(FL_altas_seu$anno_lvl_2_final_clean)[grep('MAC',unique(FL_altas_seu$anno_lvl_2_final_clean))]),group.by = 'anno_lvl_2_final_clean',
          features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FL MACROPHAGE');p

VlnPlot(subset(FL_altas_seu,anno_lvl_2_final_clean %in% c('MACROPHAGE_IRON_RECYCLING','MACROPHAGE_MHCII_HIGH')),group.by = 'age',split.by = 'anno_lvl_2_final_clean',
        features =c(cho_feature[-1],'HLA-DRA','FOLR2'),stack = T)+ggtitle('FL MACROPHAGE')
FL_altas_seu$MHCII='negtive'
FL_altas_seu@meta.data[colnames(subset(FL_altas_seu,subset =`HLA-DRA` >1)),'MHCII']='positve'
VlnPlot(subset(FL_altas_seu,anno_lvl_2_final_clean %in% unique(FL_altas_seu$anno_lvl_2_final_clean)[grep('MAC',unique(FL_altas_seu$anno_lvl_2_final_clean))]),group.by = 'anno_lvl_2_final_clean',split.by = 'MHCII',
        features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FL MACROPHAGE')


#MACROPHAGE_ERY MACROPHAGE_IRON_RECYCLING   MACROPHAGE_KUPFFER_LIKE     MACROPHAGE_LYVE1_HIGH     MACROPHAGE_MHCII_HIGH           MACROPHAGE_PERI  MACROPHAGE_PROLIFERATING          MACROPHAGE_TREM2 
#3072                     36851                      5091                      1260                     12156                         4                      4681                        12 

p=FeatureScatter(subset(FL_altas_seu,anno_lvl_2_final_clean %in% c('MACROPHAGE_IRON_RECYCLING','MACROPHAGE_KUPFFER_LIKE','MACROPHAGE_LYVE1_HIGH','MACROPHAGE_MHCII_HIGH','MACROPHAGE_PROLIFERATING')), 
               feature1 = "CXCL12", feature2 = "TNF",raster=FALSE,group.by ='anno_lvl_2_final_clean',cols = cols )
ggsave(p,filename='res_pic/main_figure5/FL_Mac_CXCL12_TNF_ligand_expression_vlnplot.pdf',height = 6,width = 6)

{
  FL_altas_seu$TNF_CXCL12='no'
  FL_altas_seu@meta.data[colnames(subset(FL_altas_seu,subset = CXCL12 < 0.001 & TNF >0.001)),'TNF_CXCL12']='TNF+'
  FL_altas_seu@meta.data[colnames(subset(FL_altas_seu,subset = CXCL12 > 0.001 & TNF <0.001)),'TNF_CXCL12']='CXCL12+'
  FL_altas_seu@meta.data[colnames(subset(FL_altas_seu,subset = CXCL12 > 0.001 & TNF >0.001)),'TNF_CXCL12']='CXCL12+TNF+'
  
  VlnPlot(subset(FL_altas_seu,anno_lvl_2_final_clean %in% unique(FL_altas_seu$anno_lvl_2_final_clean)[grep('MAC',unique(FL_altas_seu$anno_lvl_2_final_clean))]),group.by = 'anno_lvl_2_final_clean',split.by = 'TNF_CXCL12',
          features =cho_feature[-1],stack = T)+ggtitle('FL MACROPHAGE')
  
}

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("MESOTHELIUM")),group.by = 'age',features =c('CXCL12',"IGF1","IGF2"  , "LTA","LTB" ),stack = T)+NoLegend()+ggtitle('MESOTHELIUM');p
ggsave(p,filename='res_pic/main_figure5/FL_MESOTHELIUM_key_ligand_expression_vlnplot.pdf',height = 4,width = 6)
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
ggsave(as.ggplot(p),filename='res_pic/main_figure5/key_ligand_expression_inFL_niche_samples.pdf',width =6,height = 4 )

FL_cho_gene_aggregated_exp_df=t(FL_cho_gene_aggregated_exp)
FL_cho_gene_aggregated_exp_df=data.frame(FL_cho_gene_aggregated_exp_df/FL_cho_gene_aggregated_exp_df[,'IGF1'])
FL_cho_gene_aggregated_exp_df=FL_cho_gene_aggregated_exp_df[,c('IGF1','EPO','CXCL12','IFNG','TNF','IGF2','LTA','LTB')]
FL_cho_gene_aggregated_exp_df[FL_cho_gene_aggregated_exp_df==0]=NA
boxplot(FL_cho_gene_aggregated_exp_df,rm.na=T)
round(colMedians(as.matrix(FL_cho_gene_aggregated_exp_df),na.rm = T),digits = 1)
#IGF1    EPO CXCL12   IFNG    TNF   IGF2    LTA    LTB 
#1.0    0.7    1.4    1.2    1.2    1.2    0.9    1.4 
# FL_agrregated_expression_ref_IGF1_boxplot.pdf,6x6 

FeaturePlot(subset(FL_altas_seu,age %in% c('CS14_4PCW', 'CS15_5PCW', 'CS17_6PCW','CS18','CS22','CS23' )),features = 'CXCL12',split.by = 'age',ncol = 2)+NoLegend()
FeaturePlot(subset(FL_altas_seu,age %in% c('CS14_4PCW', 'CS15_5PCW', 'CS17_6PCW','CS18','CS22','CS23' )),features = 'TNF',split.by = 'age',ncol = 2)+NoLegend()

BM_altas_seu=readRDS('../NRBC_BM_altas/BM_altas_seu_v2.rds')#
VlnPlot(BM_altas_seu,group.by = 'age',features = cho_feature,stack = T)

table(BM_altas_seu$new_celltype[BM_altas_seu$stage=='EBM'])[table(BM_altas_seu$new_celltype[BM_altas_seu$stage=='EBM']) >0]
EBM_cho_celltype=c("MACROPHAGE", "Myoprogenitor","Osteoprogenitors","CHONDROBLAST",   "CHONDROCYTE","PMSC1","PMSC2","OCPs","BMSC1","BMSC2" )
p0=VlnPlot(subset(subset(BM_altas_seu,stage=='EBM'),new_celltype %in% EBM_cho_celltype ),features = cho_feature[-5],stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('EBM ALTAS')
ggsave(p0,filename='res_pic/main_figure5/EBM_celltype_key_ligand_expression_vlnplot.pdf',height = 4,width = 6)

levels(BM_altas_seu$new_celltype)

FBM_cho_cells=c("LMPP_MLP","DCs", "MACROPHAGE" ,"ILC","B CELLs" ,"NK/T CELLS", "CHONDROCYTE" , "OSTEOBLAST","FIBROBLASTS","SKELETAL MUSCLE","VSMC","Schwan cell"  )
p2=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% FBM_cho_cells ),features = cho_feature,stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('FBM ALTAS')

p=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% c('MACROPHAGE')),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FBM  MACROPHAGE');p
ggsave(p,filename='res_pic/main_figure5/FBM_Mac_key_ligand_expression_vlnplot.pdf',height = 6,width =4)
rm(p);gc()

p=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% c('MONOCYTE')),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FBM  MONOCYTE');p
ggsave(p,filename='res_pic/main_figure5/FBM_Mono_key_ligand_expression_vlnplot.pdf',height = 6,width =4)

VlnPlot(subset(BM_altas_seu,stage=='FBM'),features = c("IFNG",'TNF','LTA'),stack = T,group.by = 'anno_final_celltype2')+NoLegend()+ggtitle('FBM ALTAS')

# CD16 monocyte 和Proliferation T/NK 表达TNF 
VlnPlot(subset(BM_altas_seu,stage=='ABM'),features = c("IFNG",'TNF','LTA'),stack = T,group.by = 'ct')+NoLegend()+ggtitle('ABM ALTAS')


ABM_cho_cells=c("NEUTROPHIL","MACROPHAGE","Mac_Ery" ,"CLP","B CELLs","Plasma cell", "CD4 T","CD8 T","Treg","NK/T CELLS","OSTEOCLAST",
                "ENDOTHELIUM","Adipo-MSC","APOD+ MSC","Fibro-MSC","Osteo-MSC","THY1+ MSC","RNAlo MSC"  )
p3=VlnPlot(subset(subset(BM_altas_seu,stage=='ABM'),new_celltype %in% ABM_cho_cells),features = cho_feature,stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('ABM ALTAS')

p=p1+p2+p3;p
ggsave(p,filename='res_pic/main_figure5/key_ligand_expression_celltype_niche_vlnplot.pdf',width = 15,height = 5)


BM_altas_seu$donor[is.na(BM_altas_seu$donor)]=BM_altas_seu$sample[is.na(BM_altas_seu$donor)]
BM_altas_seu$id=paste(BM_altas_seu$donor,BM_altas_seu$age,sep="_")
table(BM_altas_seu$id)
BM_cho_gene_aggregated_exp=AggregateExpression(BM_altas_seu,group.by = 'id',features =cho_feature )$RNA
BM_cho_gene_aggregated_exp=log2(BM_cho_gene_aggregated_exp+1)
p=pheatmap(BM_cho_gene_aggregated_exp[,grep('H|F|CS',colnames(BM_cho_gene_aggregated_exp))],cluster_cols = F,cluster_rows = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure5/key_ligand_expression_in_BM_sample_heatmap.pdf',width =6,height = 4)

BM_cho_gene_aggregated_exp_df=t(BM_cho_gene_aggregated_exp[,grep('H|F',colnames(BM_cho_gene_aggregated_exp))])
BM_cho_gene_aggregated_exp_df=data.frame(BM_cho_gene_aggregated_exp_df/BM_cho_gene_aggregated_exp_df[,'IGF1'])
BM_cho_gene_aggregated_exp_df=BM_cho_gene_aggregated_exp_df[,c('IGF1','CXCL12','IFNG','TNF','IGF2','LTA','LTB')]
BM_cho_gene_aggregated_exp_df[BM_cho_gene_aggregated_exp_df==0]=NA
boxplot(BM_cho_gene_aggregated_exp_df[grep('CW',rownames(BM_cho_gene_aggregated_exp_df)),],rm.na=T)
round(colMedians(as.matrix(BM_cho_gene_aggregated_exp_df[grep('CW',rownames(BM_cho_gene_aggregated_exp_df)),]),na.rm = T),digits = 1)
#CXCL12 IGF1 IFNG  TNF                    LTA  LTB 
#1.4  1.0  0.7  1.3  1.3  1.1  1.9 
# FBM_agrregated_expression_ref_IGF1_boxplot.pdf,7.5 x 5.6

boxplot(BM_cho_gene_aggregated_exp_df[grep('y',rownames(BM_cho_gene_aggregated_exp_df)),],rm.na=T)
round(colMedians(as.matrix(BM_cho_gene_aggregated_exp_df[grep('y',rownames(BM_cho_gene_aggregated_exp_df)),]),na.rm = T),digits = 1)
# ABM_agrregated_expression_ref_IGF1_boxplot.pdf,7.5 x 5.6
#IGF1 CXCL12   IFNG    TNF   IGF2    LTA    LTB 
#1.0    1.2    0.6    0.6    0.8    0.4    1.0 

YS_altas_seu= readRDS('../NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' )
YS_altas_seu=NormalizeData(YS_altas_seu)
levels(YS_altas_seu$subcelltype)
YS_cho_celltype=c("EO/BASO/MAST" ,"MACROPHAGE","DEF_HSPC" ,"LMPP", "MOP","MONOCYTE","MOMO_MAC_DC","ELP","ILC","NK",  
                  "B_CELL", "MONOCYTE_MACROPHAGE" , "MESOTHELIUM","SMOOTH_MUSCLE" ,
                  "ENDODERM","ENDOTHELIUM","FIBROBLAST")
p4=VlnPlot(subset(YS_altas_seu,subcelltype %in% YS_cho_celltype),group.by = 'subcelltype',features =cho_feature,stack = T)+NoLegend()+ggtitle('YS ALTAS')
p4
ggsave(p4,filename='res_pic/main_figure5/key_ligand_expression_celltype_YS_vlnplot.pdf',width = 6,height = 6)

p4=VlnPlot(subset(YS_altas_seu,subcelltype %in% c('FIBROBLAST','SMOOTH_MUSCLE','MACROPHAGE')),cols = cols,group.by = 'stage',split.by = 'subcelltype',
           features =c('CXCL12','IGF1','IGF2','TNF'),stack = T)
p4
ggsave(p4,filename='res_pic/main_figure5/CXCL12_IGF2_expression_stagetime_YS_vlnplot.pdf',height = 6,width = 4)

YS_cho_gene_aggregated_exp=AggregateExpression(YS_altas_seu,features = cho_feature,group.by = 'id')$RNA
YS_cho_gene_aggregated_exp=log2(YS_cho_gene_aggregated_exp+1)
p=pheatmap(YS_cho_gene_aggregated_exp,cluster_rows = F,cluster_cols = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure5/key_ligand_expression_in_YS_sample_heatmap.pdf',height = 4,width = 6)

YS_cho_gene_aggregated_exp_df=t(YS_cho_gene_aggregated_exp)
YS_cho_gene_aggregated_exp_df=data.frame(YS_cho_gene_aggregated_exp_df/YS_cho_gene_aggregated_exp_df[,'IGF1'])
YS_cho_gene_aggregated_exp_df=YS_cho_gene_aggregated_exp_df[,c('IGF1','EPO','CXCL12','IFNG','TNF','IGF2','LTA','LTB')]
boxplot(as.matrix(YS_cho_gene_aggregated_exp_df))
round(colMedians(as.matrix(YS_cho_gene_aggregated_exp_df),na.rm = T),digits = 1)
#  IGF1    EPO CXCL12   IFNG    TNF   IGF2    LTA    LTB 
# 1.0    0.6    0.9    0.1    1.1    1.2    0.5    0.8 
# YS_agrregated_expression_ref_IGF1_boxplot.pdf,7.5 x 5.6

p=VlnPlot(subset(YS_altas_seu,subcelltype %in% c("MOP","MONOCYTE",'MONOCYTE_MACROPHAGE',"MOMO_MAC_DC")),cols = cols,group.by = 'stage',split.by = 'subcelltype',features =cho_feature[-1],stack = T)+ggtitle('YS MONOCYTE')
ggsave(as.ggplot(p),filename='res_pic/main_figure5/key_ligand_expression_celltype_YS_Mono_vlnplot.pdf',height = 6,width = 6)

#--------------------------------------IGF1+ Mac 在FL空间上的表达情况-------------------------------------------#
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

FL_H2A_sp_seu=Load10X_Spatial(data.dir ='../ref_data/ref_scRNAseq_data/E_MTAB_13067_FL_spational/H_2/H_2_A' ,filename = 'H_2_A-filtered_feature_bc_matrix.h5',image.name ="tissue_hires_image.png" )
DefaultAssay(FL_H2A_sp_seu)='Spatial'
FL_H2A_sp_seu=NormalizeData(FL_H2A_sp_seu)


boxplot(FL_H2A_sp_seu@assays$Spatial$data['IGF1',]) # 9 个样本
boxplot(FL_H5A_sp_seu@assays$Spatial$counts['IGF1',]) # only 1 

# 发现IGF1 表达为细胞几乎都是CD68+ Mac, 且都紧邻NRBC
p=SpatialFeaturePlot(FL_H2A_sp_seu, features = c("GYPA","IGF1", "CD68"),images = NULL,slot = 'data',pt.size.factor =5,image.scale='hires',crop = TRUE)
ggsave(p,filename='res_pic/main_figure5/FL_spational_IGF1_expression_raw.pdf',width = 10,height = 5)


# 提取坐标和表达数据
coords <- GetTissueCoordinates(FL_H2A_sp_seu, scale = NULL, which = "centroids",
                               cols = c("imagerow", "imagecol"))
expr <- GetAssayData(FL_H2A_sp_seu, layer = "data")[c("GYPA", "IGF1",'IGF2', "CD68",'CXCL12'), ]

# 构建数据框
df <- data.frame(
  x = coords[, 'x'], 
  y = coords[, 'y'],
  GYPA = as.numeric(expr["GYPA", rownames(coords)]),
  IGF1 = as.numeric(expr["IGF1", rownames(coords)]),
  CD68 = as.numeric(expr["CD68", rownames(coords)]),
  CXCL12 = as.numeric(expr["CXCL12", rownames(coords)]),
  IGF2 = as.numeric(expr["IGF2", rownames(coords)])
  
)

# 创建三个图层叠加
p_alpha <- ggplot() +
  geom_point(data = df[df$GYPA > 1, ], aes(x = x, y = y),  size = 2, alpha = 0.5,color='red' )+ #+scale_color_gradient(low = 'white',high =  'firebrick3')+
  geom_point(data = df[df$IGF1 > 0.1, ], aes(x = x, y = y), size = 4, alpha = 0.5,color='green')+ #+scale_color_gradient(low = 'white',high =  'green')+
  geom_point(data = df[df$CD68 > 0.1, ], aes(x = x, y = y), size = 3, alpha = 0.5,color='blue') +#+scale_color_gradient(low = 'white',high =  "blue")+
  geom_point(data = df[df$CXCL12 > 0.1, ], aes(x = x, y = y), size = 3, alpha = 0.5,color='gray') +#+scale_color_gradient(low = 'white',high =  "blue")+
  theme_void() +
  labs(title = "FL：GYPA (Red) | IGF1 (Green) | CD68 (Blue) |CXCL12(Purple)")

print(p_alpha)
ggsave(p_alpha,filename='res_pic/main_figure5/FL_spational_IGF1_expression.pdf',width = 5,height = 5)

p_alpha2 <- ggplot() +
  geom_point(data = df[df$GYPA > 1, ], aes(x = x, y = y),  size = 2, alpha = 0.5,color='red' )+ #+scale_color_gradient(low = 'white',high =  'firebrick3')+
  geom_point(data = df[df$CD68 > 0.1, ], aes(x = x, y = y), size = 3, alpha = 0.5,color='blue') +#+scale_color_gradient(low = 'white',high =  "blue")+
  geom_point(data = df[df$IGF2 > 0.1, ], aes(x = x, y = y), size = 3, alpha = 0.5,color='yellow') +#+scale_color_gradient(low = 'white',high =  "blue")+
  theme_void() +
  labs(title = "FL：GYPA (Red) | CD68 (Blue) |IGF2(Yellow)")

print(p_alpha2)






#FL_H5A_sp_seu=Load10X_Spatial(data.dir ='../ref_data/ref_scRNAseq_data/E_MTAB_13067_FL_spational/H_5/H_5_A' ,filename = 'H_5_A-filtered_feature_bc_matrix.h5')# 5 IGF1 :0
#DefaultAssay(FL_H5A_sp_seu)='Spatial'
#FL_H5A_sp_seu=NormalizeData(FL_H5A_sp_seu)

# ----------------IGF1 在 fetal Mac中的表达-----------------#

early_organogenesis_seu=readRDS('../ref_data/early_organogenesis_seu.rds')
# IGF2 广泛表达
VlnPlot(early_organogenesis_seu,group.by  = 'annotation',stack = T,
        features = c('IGF1','IGF2','CXCL12','TNF','LTA','LTB'))+ggtitle('CS12-CS16 embryo:blood systerm')+NoLegend()

p=VlnPlot(subset(early_organogenesis_seu,annotation %in% c('macrophage')),group.by  = 'dissection_part',stack = T,
          features = c('IGF1','IGF2','CXCL12','TNF','LTA','LTB'))+ggtitle('CS12-CS16 embryo:blood system')+NoLegend() # kupffer cell 不表达IGF1,TNF
ggsave(p,filename='res_pic/main_figure5/embryo_blood_Mac_IGF1_expression.pdf',width = 4,height = 6)

p=VlnPlot(subset(early_organogenesis_seu,annotation %in% c('macrophage')),group.by  = 'stage',stack = T,
          features = c('IGF1','IGF2','CXCL12','TNF','LTA','LTB'))+ggtitle('CS12-CS16 embryo:blood system')+NoLegend() # kupffer cell 不表达IGF1
p



# 在小孙部署的docker上完成的分析
# /mnt/data2/bio_projects_rawdata/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/fetal_macrophaege_IGF1.ipynb
if(F){
  human_Fetal_altas_seu=LoadH5Seurat('human_cell_atlas/human_immune_system_across_organs/PAN.A01.v01.raw_count.20210429.PFI.embedding.h5seurat', meta.data = FALSE,reduction='umap')
  human_Fetal_altas_seu_meta=read.csv('human_cell_atlas/human_immune_system_across_organs/PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.annotated.clean.csv',header =T )
  
  rownames(human_Fetal_altas_seu_meta)=human_Fetal_altas_seu_meta$X;human_Fetal_altas_seu_meta=human_Fetal_altas_seu_meta[,-1]
  dim(human_Fetal_altas_seu_meta)#  908178
  table(colnames(human_Fetal_altas_seu) %in%  rownames(human_Fetal_altas_seu_meta))# FALSE   TRUE ：3695 908178 
  human_Fetal_altas_seu=subset(human_Fetal_altas_seu,cells =rownames(human_Fetal_altas_seu_meta) )
  human_Fetal_altas_seu@meta.data=human_Fetal_altas_seu_meta[colnames(human_Fetal_altas_seu),]
  
  DimPlot(human_Fetal_altas_seu,group.by = c('organ','uniform_label_lvl0'),raster=FALSE)
  sort(unique(human_Fetal_altas_seu$anno_lvl_2_final_clean))
  cho_MAC=sort(unique(human_Fetal_altas_seu$anno_lvl_2_final_clean))[grep('MACROPHAGE',sort(unique(human_Fetal_altas_seu$anno_lvl_2_final_clean)))]
  cho_MAC
  human_Fetal_altas_seu=subset(human_Fetal_altas_seu,anno_lvl_2_final_clean %in% cho_MAC[-1:-2])
  human_Fetal_altas_seu=subset(human_Fetal_altas_seu,predicted_doublets=="False")
  DimPlot(human_Fetal_altas_seu,group.by = c('organ'),raster=FALSE)

  VlnPlot(human_Fetal_altas_seu,feature=c('IGF1',"CXCL12","IFNG","TNF"),stack=T,group.by='organ',split.by='anno_lvl_2_final_clean',drop = FALSE)
  
  VlnPlot(human_Fetal_altas_seu,feature=c('IGF1',"CXCL12","IFNG","TNF"),group.by='organ',stack=T,layer='data')+xlim(c(0,15))+ggtitle('the Macrophage from 4PWC-17PWC fetus tissues')
}else{
  human_Fetal_altas_Mac=readRDS('../ref_data/ref_scRNAseq_data/human_cell_atlas/human_immune_system_across_organs/human_Fetal_altas_Mac_seu.rds')
  cho_org=unique(human_Fetal_altas_Mac$organ)
  cho_org=cho_org[!cho_org %in% c('YS','LI','BM')]
  human_Fetal_altas_Mac_YS_FL_FBM=subset(human_Fetal_altas_Mac,organ %in% c('YS','LI','BM'))
  
  p1=VlnPlot(subset(human_Fetal_altas_Mac_YS_FL_FBM,organ %in% 'YS'),feature=c('IGF1','IGF2',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='age',stack=T,layer='data',cols = cols)+xlim(c(0,15))+ggtitle('YS  Macrophage')+NoLegend()
  p1 # 样本太少了
  
  YS_altas_seu=readRDS('../NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds')
  YS_altas_seu=NormalizeData(YS_altas_seu)
  p1=VlnPlot(subset(YS_altas_seu,subcelltype=='MACROPHAGE'),feature=c('IGF1','IGF2',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='stage',stack=T,layer='data',cols = cols)+xlim(c(0,15))+ggtitle('YS  Macrophage')+NoLegend()
  p1 # 样本太少了
  
  p2=VlnPlot(subset(human_Fetal_altas_Mac_YS_FL_FBM,organ %in% 'LI'),feature=c('IGF1','IGF2',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='age',stack=T,layer='data',cols = cols)+xlim(c(0,15))+ggtitle('FL Macrophage')+NoLegend()
  p2
  # 采用FL altas 数据
  FL_altas_seu=readRDS('../NRBC_FL_altas/tmp_FL_altas_seu.rds')
  p2=VlnPlot(subset(FL_altas_seu,subcelltype %in% 'MACROPHAGE'),feature=c('IGF1',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='age',stack=T,layer='data',cols = cols)+xlim(c(0,15))+ggtitle('FL Macrophage')+NoLegend()
  p2
  
  p3=VlnPlot(subset(human_Fetal_altas_Mac_YS_FL_FBM,organ %in% 'BM'),feature=c('IGF1',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='age',stack=T,layer='data',cols = cols)+xlim(c(0,15))+ggtitle('FBM Macrophage')+NoLegend()
  p3
  
  p1+p2+p3
  
  human_Fetal_altas_Mac=subset(human_Fetal_altas_Mac,organ %in% cho_org)
  human_Fetal_altas_Mac=NormalizeData(human_Fetal_altas_Mac)
  human_Fetal_altas_Mac=FindVariableFeatures(human_Fetal_altas_Mac);human_Fetal_altas_Mac=ScaleData(human_Fetal_altas_Mac)
  p=VlnPlot(human_Fetal_altas_Mac,feature=c('IGF1','IGF2',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='organ',stack=T,layer='data',cols = cols)+ggtitle('the Macrophage from fetus tissues')+NoLegend()
  p
  ggsave(p,filename='res_pic/main_figure5/fetus_other_organ_tissue_Mac_IGF1_lvnplot.pdf',width = 4,height = 6)
 
  p=VlnPlot(human_Fetal_altas_Mac,feature=c('IGF1','IGF2',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='age',stack=T,layer='data',cols = cols)+ggtitle('the Macrophage from fetus tissues')+NoLegend()
  p
  ggsave(p,filename='res_pic/main_figure5/fetus_other_age_Mac_IGF1_lvnplot.pdf',width = 4,height = 6)
  
  p=VlnPlot(human_Fetal_altas_Mac,feature=c('IGF1','IGF2',"CXCL12","IFNG","TNF",'LTA','LTB'),group.by='age',split.by = 'organ',stack=T,layer='data',cols = cols)+
    ggtitle('the Macrophage from 4PWC-17PWC fetus tissues')
  p
  ggsave(p,filename='res_pic/main_figure5/fetus_other_age_organ_Mac_IGF1_lvnplot.pdf',width = 4,height = 6)
  table(human_Fetal_altas_Mac@meta.data[,c('organ','age')])
  
}

########################################################################################################################################################################
###################----------------------------------------------the signaling pathway -------------------------------------###################
########################################################################################################################################################################
library(nichenetr)
library(network)
library(ggnetwork)

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


# =--------------------测试：CXCL12->CXCR4->信号传导->TF->ANXA1-------------------------#
# TF : GATA1
ligands_oi <- "CXCL12" # this can be a list of multiple ligands if required
targets_oi <- c( "CXCR4","ANXA1")#"IFI16","CITED2","PIM1","SGK1","CXCR4","MT1E","PFKFB4","SLC38A2", "JUN" ,
CXCL12_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )

cho_fetal_ligand_target_df=fetal_ligand_target_res_list[[3]][,c('IGF1','IGF2','IFNG','TNF','LTA','LTB','EPO')]
cho_fetal_ligand_target_list=list()
cho_fetal_ligand_target_list[['IGF1']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IGF1'] >0.04]
cho_fetal_ligand_target_list[['IGF2']]=names(cho_fetal_ligand_target_df[,'IGF2'])[cho_fetal_ligand_target_df[,'IGF2'] >0.04]
cho_fetal_ligand_target_list[['IFNG']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IFNG'] >0.04]
cho_fetal_ligand_target_list[['TNF']]=names(cho_fetal_ligand_target_df[,'TNF'])[cho_fetal_ligand_target_df[,'TNF'] >0.04]
cho_fetal_ligand_target_list[['LTA']]=names(cho_fetal_ligand_target_df[,'LTA'])[cho_fetal_ligand_target_df[,'LTA'] >0.04]
cho_fetal_ligand_target_list[['LTB']]=names(cho_fetal_ligand_target_df[,'LTB'])[cho_fetal_ligand_target_df[,'LTB'] >0.04]
cho_fetal_ligand_target_list[['EPO']]=names(cho_fetal_ligand_target_df[,'EPO'])[cho_fetal_ligand_target_df[,'EPO'] >0.04]


ligands_oi <- c("IGF1") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("IGF1")]))))
IGF1_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# IGF1_signaling_network.pdf, 6x6


ligands_oi <- c("IGF2") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("IGF2")]))))
IGF2_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# IGF2_signaling_network.pdf, 5x5


ligands_oi <- c("IFNG") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("IFNG")]))))
IFNG_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# IFNG_signaling_network.pdf, 6x6

ligands_oi <- c("TNF") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("TNF")]))))
TNF_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# TNF_signaling_network.pdf, 6x6


ligands_oi <- c("LTA") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("LTA")]))))
LTA_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# LTA_signaling_network.pdf, 5x5

ligands_oi <- c("LTB") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("LTB")]))))
LTB_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# LTB_signaling_network.pdf, 5x5


BM_altas_seu=readRDS('NRBC_BM_altas/BM_altas_seu_v2.rds')#
VlnPlot(subset(BM_altas_seu,stage=='FBM'),group.by = 'new_celltype',features = c('EPO','CXCL12','IGF1','IGF2','IFNG','TNF'),stack = T)+NoLegend()+ggtitle('FBM ALTAS')
VlnPlot(subset(BM_altas_seu,stage=='ABM'),group.by = 'new_celltype',features = c('EPO','CXCL12','IGF1','IGF2','IFNG','TNF'),stack = T)+NoLegend()+ggtitle('ABM ALTAS')

rm(list = c('BM_altas_seu','FL_altas_seu'));gc()


cho_adult_ligand_target_df=adult_ligand_target_res_list[[3]][,c('EPO','CXCL12')]
cho_adult_ligand_target_list=list()
cho_adult_ligand_target_list[['EPO']]=names(cho_adult_ligand_target_df[,'EPO'])[cho_adult_ligand_target_df[,'EPO'] >0.04]
cho_adult_ligand_target_list[['CXCL12']]=names(cho_adult_ligand_target_df[,'CXCL12'])[cho_adult_ligand_target_df[,'CXCL12'] >0.04]

ligands_oi <- 'EPO' # this can be a list of multiple ligands if required
targets_oi <- cho_adult_ligand_target_list[['EPO']]
adult_EPO_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# EPO_signaling_adult_network.pdf,5 X 5 

targets_oi <- cho_fetal_ligand_target_list[['EPO']]
fetal_EPO_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# EPO_signaling_fetal_network.pdf,5 X 5 

targets_oi <-unique(c( cho_adult_ligand_target_list[['EPO']], cho_fetal_ligand_target_list[['EPO']]))
fetal_adulkt_EPO_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# EPO_signaling_fetal_adult_network.pdf

ligands_oi <- 'CXCL12' # this can be a list of multiple ligands if required
targets_oi <- cho_adult_ligand_target_list[['CXCL12']]
CXCL12_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# CXCL12_signaling_network.pdf ,5 x5 

# analysis the expression of Signaling mediators based on ligand signals
Signaling_mediator_genelist=list()
Signaling_mediator_genelist[['CXCL12']]=unique(CXCL12_screated_network$data$label)
Signaling_mediator_genelist[['CXCL12']]=Signaling_mediator_genelist[['CXCL12']][ !Signaling_mediator_genelist[['CXCL12']] %in% cho_adult_ligand_target_list[['CXCL12']] ]
Signaling_mediator_genelist[['EPO']]=unique(fetal_adulkt_EPO_screated_network$data$label)
Signaling_mediator_genelist[['EPO']]=Signaling_mediator_genelist[['EPO']][ !Signaling_mediator_genelist[['EPO']] %in% unique(c( cho_adult_ligand_target_list[['EPO']], cho_fetal_ligand_target_list[['EPO']]))]

ligand='IGF1'
Signaling_mediator_genelist[[ligand]]=unique(IGF1_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% cho_fetal_ligand_target_list[[ligand]]]
ligand='IGF2'
Signaling_mediator_genelist[[ligand]]=unique(IGF1_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% cho_fetal_ligand_target_list[[ligand]]]
ligand='IFNG'
Signaling_mediator_genelist[[ligand]]=unique(IFNG_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% cho_fetal_ligand_target_list[[ligand]]]
ligand='TNF'
Signaling_mediator_genelist[[ligand]]=unique(TNF_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% cho_fetal_ligand_target_list[[ligand]]]
ligand='LTA'
Signaling_mediator_genelist[[ligand]]=unique(LTA_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% cho_fetal_ligand_target_list[[ligand]]]
ligand='LTB'
Signaling_mediator_genelist[[ligand]]=unique(LTB_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% cho_fetal_ligand_target_list[[ligand]]]
length(unique(unlist(Signaling_mediator_genelist)))
all_genes=unique(unlist(Signaling_mediator_genelist))
all_genes=all_genes[!all_genes %in% c('CXCL12','IGF1','IGF2','EPO','LTA','LTB','IFNG','TNF')]
p=DotPlot(object = filt_NBRC_altas_seu,features =unique(unlist(Signaling_mediator_genelist)),group.by = 'source_celltype',scale = F )
temp_df=p$data
temp_df=temp_df[-grep('YS',temp_df$id),]
rownames(temp_df)=NULL


temp_df=temp_df[temp_df$avg.exp >0.1 & temp_df$pct.exp>10,]
length(unique(temp_df$features.plot))
ggplot(temp_df,aes(y=id,x=features.plot,size=pct.exp,color=avg.exp))+geom_point()+theme_classic()+scale_color_gradient(low = 'white',high = 'firebrick3')+RotatedAxis()

order_Signaling_mediator_genes=c(unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['EPO']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in%    Signaling_mediator_genelist[['IGF1']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['IFNG']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['TNF']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['IGF2']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['CXCL12']]]
                                 )

order_Signaling_mediator_genes=as.character(unique(order_Signaling_mediator_genes));length(order_Signaling_mediator_genes)
p1=DotPlot(object = filt_NBRC_altas_seu,features =unique(temp_df$features.plot),group.by = 'source_celltype',scale = F,cols = c('gray','firebrick3') )+RotatedAxis()
