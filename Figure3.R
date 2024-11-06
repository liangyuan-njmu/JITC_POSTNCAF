#############################################################
#Seurat deal
dt=readRDS('merge_CAF.RDS')

##################
#Monocle2
##################
devtools::load_all("~/R/library/4.2.0/monocle")

dat <- Seurat::as.CellDataSet(dt)
dat <- estimateSizeFactors(dat)
dat <- detectGenes(dat, min_expr = 10)
fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(dat),num_cells_expressed >= 10))
dat <- reduceDimension(dat,
                          max_components = 2,
                          norm_method = 'log',
                          num_dim = 20,
                          reduction_method = 'tSNE',
                          verbose = T,
                          check_duplicates=F)
dat <- clusterCells(dat,verbose = F)
clustering_DEG_genes <- differentialGeneTest(dat[expressed_genes,],
           fullModelFormulaStr = '~Cluster',
           cores = 10)

#dat <- setOrderingFilter(dat, ordering_genes = clustering_DEG_genes$gene_short_name[clustering_DEG_genes$use_for_ordering=="TRUE"])
dat <- setOrderingFilter(dat, ordering_genes = unique(top50$gene))
#dat <- setOrderingFilter(dat, ordering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000])
dat <- reduceDimension(dat, method = 'DDRTree')
dat <- orderCells(dat)
saveRDS(dat,"CAF_monocle_top50_2.RDS")

p<-plot_cell_trajectory(dat, cell_size=1, show_backbone=FALSE,cell_link_size = 1,color_by="CellType",theta=0.8)+
	facet_wrap(~CellType, nrow = 3)+scale_color_manual(values=CAFcol)
ggsave('./CAF_CellType_wrap_top50.pdf',p)
p<-plot_cell_trajectory(dat, cell_size=0.5,show_tree=T, show_backbone=FALSE,cell_link_size = 1,color_by="CellType")+
	scale_color_manual(values=CAFcol)
ggsave('./CAF_CellType_top50.pdf',p,width=6,height=6)

p<-plot_cell_trajectory(dat, cell_size=1, show_backbone=FALSE,cell_link_size = 1,color_by="CellType")+
	facet_wrap(~SampleType2, nrow = 3)+scale_color_manual(values=Myecol)
ggsave('./CAF_SampleType_wrap_top50.pdf',p,useDingbat=F)

p<-plot_cell_trajectory(dat, cell_size=1, show_backbone=FALSE,cell_link_size = 1,color_by="CellType")+
	facet_wrap(~sampletype2, nrow = 1)+scale_color_manual(values=Myecol)
ggsave('./CAF_SampleType2_wrap_top50.pdf',p,width=14,height=7,useDingbat=F)

dat$Pseudotime<-max(dat$Pseudotime)-dat$Pseudotime
p<-plot_cell_trajectory(dat, cell_size=1, show_backbone=FALSE,cell_link_size = 1,color_by="Pseudotime")+viridis::scale_color_viridis()
ggsave('./CAF_Pseudotime_top50.pdf',p,useDingbat=F)

dat <- detectGenes(dat)
dat_subset <- dat[fData(dat)$num_cells_expressed >=10, ]
BEAM_res <- BEAM(dat_subset, branch_point=2, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

p<-plot_genes_branched_heatmap(dat_subset[rownames(BEAM_res)[BEAM_res$qval<0.001],],
                branch_point = 2,
                num_clusters = 4,
                cores = 1,
                return_heatmap=T,
                show_rownames = F)
p_clu5=plot_genes_branched_heatmap(dat_subset[rownames(BEAM_res)[BEAM_res$qval<0.001],],
                branch_point = 2,
                num_clusters = 5,
                cores = 1,
                return_heatmap=T,
                show_rownames = F)
pdf('CAF_top50_Branch2_clu5_pheatmap.pdf')
p_clu5$ph_res
dev.off()

p<-plot_cell_trajectory(dat, cell_size=1,show_tree=T, show_backbone=FALSE,cell_link_size = 1,color_by="State",theta=0.8)+
	scale_color_manual(values=jjAnno::useMyCol('calm',n=6))
ggsave('./CAF_State_top50.pdf',p,width=6,height=6)

tmp=sapply(1:5,function(x){
  sapply(levels(factor(dat$CellType)),function(y){
    length(which(dat$CellType==y & dat$State==x))
  })
})
colnames(tmp)=paste0('State',1:5)

plot_lists=lapply(colnames(tmp),function(x){
  df=data.frame(
    value=as.numeric(as.vector(tmp[,x])),
    group=factor(rownames(tmp),levels=levels(factor(dat$CellType)))
  )
  ggplot(df, aes(x="",y=value,fill=group))+
  geom_bar(stat="identity",width = 1)+
  scale_fill_manual(values=CAFcol)+labs(title=x)+
  coord_polar("y")+theme_minimal()+theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))+Seurat::NoLegend()

})
pdf('Monocle2_State_CellType_pie.pdf')
cowplot::plot_grid(plotlist=plot_lists,nrow=1)
dev.off()

#############################
#Scenic
dir.create("./scenic/")
setwd("scenic/")

badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB",  "-", ":", "\\.", '^KIAA',
  "^DNAJ"
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(dt),perl=T,value=T))

ct=dt@assays$RNA@counts[-which(rownames(dt) %in% badgenelt),]

mat<-t(as.matrix(dt@assays$RNA@data[rownames(ct),]))
write.table(mat,"./CAF.csv",quote=F,row.names=T,col.names=T,sep=",")

#After pyscenic deal
read.csv("./auc_mtx.csv",header=T,row.names=1,sep=",")->auc
colnames(auc)<-gsub("\\.","",colnames(auc))
t(auc)->auc
as.matrix(auc)->auc

sapply(levels(factor(dt$CellType)),function(x){
  tmp<-auc[,which(colnames(auc) %in% colnames(dt)[which(dt$CellType==x)])]
  apply(tmp,1,mean)->mean
  mean
})->a
library(pheatmap)
pdf("./dt_TF_pheatmap.pdf",height=50,width=8)
pheatmap(a,scale="row",cellheight = 10,cellwidth=10,angle_col=45,show_rownames=T,show_colnames=T,cluster_rows = T,cluster_cols=T,color = rev(heatmaply::RdBu(100)))
dev.off()





