#############################################################
#Seurat deal
dt=readRDS('merge_CAF.RDS')

##CellChat
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(ggpubr)

run_cellchat_single <- function(seuratObject, cellinfo, species="human", population.size=TRUE) {
    # seuratObject: seurat object
    # cellinfo: names in seurat meta.data
    # species: human or mouse
    # population.size: sorting-enriched single cells-->FALSE, unsorted single-cell transcriptomes-->TRUE
    library(Seurat)
    library(CellChat)
    # seuratObject@active.ident = factor(seuratObject@meta.data[[cellinfo]])
    # names(seuratObject@active.ident) = colnames(seuratObject)
    # data.input <- GetAssayData(seuratObject, assay = "RNA", slot = "data") # normalized data matrix
    # labels <- Idents(seuratObject)
    # identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
    cellchat <- createCellChat(object = seuratObject, group.by = cellinfo, assay = "RNA")
    #cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
    #cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    #levels(cellchat@idents) # show factor levels of the cell labels

    if (species == "human") {
        CellChatDB <- CellChatDB.human
        ppi = PPI.human
    } else {
        CellChatDB <- CellChatDB.mouse
        ppi = PPI.mouse
    }

    #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
    cellchat@DB <- CellChatDB

    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    future::plan("multiprocess", workers = 6) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, ppi)

    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
    return(cellchat)
}

all=readRDS('All.RDS')

Cancer=subset(all,subset=SampleType=="Tumor")
Paracancer=subset(all,subset=SampleType=="Normal")

Cancer_chat<-run_cellchat_single(Cancer,"CellType","human",population.size=TRUE)
saveRDS(Cancer_chat,"Cancer_chat.RDS")

Paracancer_chat<-run_cellchat_single(Paracancer,"CellType","human",population.size=TRUE)
saveRDS(Paracancer_chat,"Paracancer_chat.RDS")

cco.list <- list(Cancer = Cancer_chat,Paracancer=Paracancer_chat)
cellchat <- mergeCellChat(cco.list,add.names = names(cco.list),cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf('Tumor_Normal_numbers_strenth.pdf',width=7,height=5)
cowplot::plot_grid(gg1,gg2,nrow=1)
dev.off()

p1=netVisual_diffInteraction(cellchat, weight.scale = T)
p2=netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
pdf('Tumor_Normal_netVisual.pdf',width=9,height=5)
cowplot::plot_grid(p1,p2,nrow=1)
dev.off()

weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
pdf('Tumor_Normal_netVisual.pdf',width=9,height=5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
dev.off()

weight.max <- getMaxWeight(cco.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
pdf('Tumor_Normal_diffInteraction.pdf',width=9,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
dev.off()

#cibersortx
exp=readRDS('../Bulk_T_vs_N/TCGA/TCGA_HCC_Expr.RDS')
info=readRDS('../Bulk_T_vs_N/TCGA/TCGA_HCC_Info.RDS')

HCC_samples=rownames(info)[info$TYPE=="HCC"]

ciber_rs=read.table("CIBERSORTx_Job115_Adjusted.txt",header=T,row.names=1,sep="\t")
colnames(ciber_rs)=gsub("\\.","_",colnames(ciber_rs))

ciber_rs=ciber_rs[intersect(rownames(ciber_rs),HCC_samples),1:25]

library(corrplot)
matrix_corr <- cor(as.matrix(ciber_rs),method="spearman")
rownames(matrix_corr)=gsub('_','+',rownames(matrix_corr))
colnames(matrix_corr)=gsub('_','+',colnames(matrix_corr))
orders=c(levels(factor(CAF$CellType)),levels(Mye$CellType))
matrix_corr=matrix_corr[orders,orders]

res1 <- cor.mtest(matrix_corr,conf.level = .95)
pdf('Mye_CAF_Ciber_corr.pdf',width=20,height=20)
corrplot(matrix_corr,order="original", p.mat = res1$p,col=rev(COL2('RdYlBu')), method = "circle", diag = T, type = "lower",sig.level = c(0.001,.01, .05), insig = "label_sig",
tl.cex = 1, tl.col = "black",tl.srt = 1,number.cex= .8)
dev.off()

#################################################
#Nichenet
library(nichenetr)
library(Seurat)
library(tidyverse)
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")

CAF=readRDS('merge_CAF.RDS')

mg=merge(Mye,CAF)

mg@active.ident<-as.factor(mg$CellType)

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

##
min_pct=0.1

## receiver
receiver = unique(Mye$CellType)
expressed_genes_receiver = get_expressed_genes(receiver, mg, pct = min_pct)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = unique(CAF$CellType)
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, mg, min_pct) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

seurat_obj_receiver= subset(mg, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["SampleType"]])

condition_oi = "Tumor"
condition_reference = "Normal" 

future::plan("multicore", workers = 5)
options(future.globals.maxSize=3221225472)
DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = min_pct) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

pdf('mg_CellType_top30_ligand_expr_dotplot.pdf',width=10,height=8)
DotPlot(mg, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
dev.off()

#Active target gene inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.8)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
pdf('mg_CellType_ligand_target_network.pdf',width=20,height=8)
p_ligand_target_network
dev.off()

#Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network[,order_ligands] %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
pdf('mg_CellType_ligand_receptor_network.pdf',width=20,height=8)
p_ligand_receptor_network
dev.off()

#Summary visualizations of the NicheNet analysis
ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
pdf('mg_top30_ligand_activity.pdf',width=5,height=8)
p_ligand_aupr
dev.off()

# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
rotated_dotplot = DotPlot(mg %>% subset(CellType %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
pdf('mg_top30_ligand_dotplot.pdf',width=5,height=6)
rotated_dotplot
dev.off()

figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  #p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(0.5,1,3))

legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
    #ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h", rel_widths = c(0.5,  1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
pdf('mg_CAF_Mye_Nichenet_result.pdf',width=20,height=10)
cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
dev.off()

#IL6 TGFB1 receptors expression
pdf('Mye_IL6_TGFB1_receptors_expr_dotplot.pdf',width=10,height=6)
DotPlot(Mye,group.by='CellType', features = c('HRH1','IL6R','IL6ST','ITGAV','ITGB5','TGFBR2','ENG','VDR','TLR2','TGFBR1','F11R','FCN1','VTN') , cols = "RdYlBu") + RotatedAxis()
dev.off()


