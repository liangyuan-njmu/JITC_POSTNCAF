#############################################################
#Seurat deal
dt=readRDS('merge_CAF.RDS')

#GSEA
future::plan("multicore", workers = 5)

library(org.Hs.eg.db)
library(clusterProfiler)

diff_rs=FindMarkers(dt,ident.1="POSTN+CAF",group.by="CellType",logfc.threshold=0.2,min.pct=0.2,only.pos=F)

diff_rs$type=ifelse(diff_rs$avg_log2FC>=0.1 & diff_rs$p_val_adj<=0.001,"High",
             ifelse(diff_rs$avg_log2FC<=(-0.1) & diff_rs$p_val_adj<=0.001,"Low","Nosig"))
diff_rs$gene=rownames(diff_rs)
diff_rs=diff_rs[order(diff_rs$avg_log2FC,decreasing=T),]

diff_rs$logP= -log10(diff_rs$p_val_adj)
diff_rs$logP[which(diff_rs$logP=="Inf")]=max(diff_rs$logP[diff_rs$logP!="Inf"])+10
diff_rs$gene=rownames(diff_rs)
labels=c(head(diff_rs,20)$gene,tail(diff_rs,20)$gene)
diff_rs$label=ifelse(diff_rs$gene %in% labels,diff_rs$gene,"")

p=ggplot(diff_rs,aes(x=avg_log2FC,y=logP,color=type))+
    ggrastr::rasterise(geom_point(aes(color=type),size=1.5),dpi=150)+
    geom_vline(xintercept = c(-0.25,0.25),size=1,linetype="dashed")+
    geom_hline(yintercept = c(-log10(0.001)),size=1,linetype="dashed")+
    #xlim(-1,1)+
    #ggrepel::geom_text_repel(aes(label=label),size=3,max.overlaps=Inf)+
    scale_color_manual(values=c('High'="#89288F",'Low'="lightblue",'Nosig'="grey87"))+
    ylab("-log10(p_val_adj)")+ggtitle("POSTN+CAF vs Other Cells")+xlab("average log2FC")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),legend.position = "none",
      axis.title.x.bottom = element_text(size = 16,color = "black"),axis.text.y.left = element_text(size = 14,color = "black"),
      axis.title.y.left = element_text(size = 16),plot.title = element_text(size = 16,hjust = 0.5)
    )
ggsave("POSTN+CAF_vs_OtherCells_diffgenes.pdf",p,width=8,height=4)

genelist <- bitr(diff_rs$gene, fromType="SYMBOL",toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- dplyr::distinct(genelist,SYMBOL,.keep_all=TRUE) 
up_df=diff_rs
up_df$SYMBOL<-up_df$gene
tmp <- up_df %>% inner_join(genelist,by="SYMBOL")
ENTREZID <- tmp$avg_log2FC
names(ENTREZID) <- tmp$ENTREZID

ego<- gseGO(geneList= ENTREZID,OrgDb= org.Hs.eg.db,ont= "BP",
        minGSSize= 0.05,pvalueCutoff = 1,verbose= FALSE)

ekegg<-gseKEGG(ENTREZID,organism="hsa",pvalueCutoff=1,minGSSize=0.05)

library("msigdbr")
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(ont = gs_name, gene = entrez_gene)
egsea=GSEA(ENTREZID, TERM2GENE = msig_h, pvalueCutoff=1)

grep("hypo",egsea@result$Description[egsea@result$pvalue<0.05],value=T,ignore.case=T)
ekegg@result$Description[ekegg@result$pvalue<0.05]
ego@result$Description[ego@result$pvalue<0.05][1:30]

select_kegg_pathways=c('ECM-receptor interaction','Focal adhesion','TGF-beta signaling pathway','PI3K-Akt signaling pathway')

library(GseaVis)
keggterms <- ekegg@result$ID[ekegg@result$Description %in% select_kegg_pathways]
lapply(keggterms, function(x){
  gseaNb(object = ekegg,geneSetID = x,addPval = T,pvalX = 0.5,pvalY = 0.65,pCol = 'black',pHjust = 0)
}) -> gseaList
pdf('POSTN+CAF_KEGG_GSEA.pdf',width=24,height=6)
cowplot::plot_grid(plotlist = gseaList,ncol = 4,align = 'hv')
dev.off()

#################################################
#T cells exclusion
#TCGA LIHC
exp=readRDS('~/TCGARNASeq/LIHC/expr.RDS')
tide_rs=read.table('TCGA.LIHC.RNASeq.norm_subtract.OS_base',header=T,row.names=1,sep="\t")
rownames(tide_rs)=paste0(rownames(tide_rs),"-01")

inte_samples=intersect(rownames(tide_rs),colnames(exp))

exp=exp[,inte_samples]
tide_rs=tide_rs[inte_samples,]

diffgenes=readRDS("CAF_CellType_diffgenes.RDS")
top=diffgenes %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

genelt=lapply(levels(top$cluster),function(x){
  top$gene[top$cluster==x]
})
names(genelt)=levels(top$cluster)

ssgsea_rs=GSVA::gsva(as.matrix(exp),gset.idx.list=genelt,method="ssgsea",min.sz=2, max.sz=100, verbose=T,parallel.sz=1)

for(i in 1:6){
  x=names(genelt)[i]
  df<-data.frame(
      CAF=as.numeric(as.vector(ssgsea_rs[x,])),
      Tex=as.numeric(as.vector(tide_rs[,'Exclusion']))
  )
  p<-ggplot(df,aes(x=CAF,y=Tex))+ 
      ggpubr::stat_cor(method = "pearson",label.sep = "\n", size = 10)+
      geom_point(color=CAFcol[i],size=2)+
      geom_smooth(color="darkred",method="lm",level = 0.99,size=4)+
      xlab(x)+ylab("T Cells Exclusion")+
      theme_bw()+theme(text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor=element_blank(), legend.position="right")
  assign(paste0("p",i),p)
}
pdf(paste0("TCGALIHC_","6CAF_Exclusion_corr",".pdf"),width=18,height=12)
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2)
dev.off()



