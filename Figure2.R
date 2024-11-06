#############################################################
#Seurat deal
dt=readRDS('merge_CAF.RDS')

#Sankey
sapply(levels(factor(dt$SampleType)),function(x){
    sapply(levels(factor(dt$CellType)),function(y){
        round(length(which(dt$SampleType==x & dt$CellType==y))/length(which(dt$SampleType==x)),4)
    })
})->tmp
df<-data.frame(value=as.numeric(as.vector(unlist(tmp))),
			  type=factor(rep(colnames(tmp),rep(nrow(tmp),ncol(tmp))),levels=c("Normal","Tumor")),
			  cell=factor(rep(rownames(tmp),ncol(tmp)),levels=levels(factor(unique(dt$CellType)))))

library(ggalluvial)
library(reshape)
dat<-data.frame(group=levels(factor(unique(dt$CellType))),Normal=df$value[df$type=="Normal"],Tumor=df$value[df$type=="Tumor"])
dat=melt(dat)
dat$order = c(1:6,1:6)
dat$group=factor(dat$group,levels=levels(factor(unique(dt$CellType))))
pdf("Sankey.pdf",width=6,height=8)
ggplot(dat,
       aes(x = variable, stratum =group, 
        alluvium = order,y = value,
        fill = group, label = group)) +
  scale_x_discrete(expand = c(.1, .1)) +
  scale_fill_manual(values=CAFcol)+
  geom_flow() +
  geom_stratum(alpha = 1) +
  ylab('Percent')+
  #geom_text(stat = "stratum", size = 2.5) +
  theme(panel.border=element_rect(color='black', fill=NA), 
       panel.grid.major =element_blank(), 
       panel.grid.minor = element_blank(),
       panel.background = element_blank(), 
       axis.line = element_line(colour = "black"),
       axis.text=element_text(colour="black")) 
dev.off()

#Cluster sampletype barplot
samplecol=c('#EC478B','#3F94D7')
names(samplecol)=c('Tumor','Normal')

tmp<-sapply(unique(dt$CellType),function(x){
    sapply(unique(dt$SampleType),function(y){
        round(length(which(dt$CellType==x & dt$SampleType==y))/length(which(dt$CellType==x)),6)
    })
})
df<-data.frame(
    value=as.numeric(as.vector(unlist(tmp))),
    cluster=factor(rep(colnames(tmp),rep(nrow(tmp),ncol(tmp)))),
    sampletype=factor(rep(rownames(tmp),ncol(tmp)))
)
df2<-df %>% group_by(cluster) %>% filter(sampletype=="Tumor") %>% arrange(desc(value))
df$cluster<-factor(df$cluster,levels=rev(df2$cluster))
df$sampletype=factor(df$sampletype,levels=c("Tumor","Normal"))
p<-ggplot(df,aes(x=cluster,y=value,fill=sampletype))+
    geom_bar(stat = 'identity',width=0.4)+
    coord_flip() + 
    scale_fill_manual(values =samplecol)+
    theme_bw()+theme(text = element_text(size = 15),axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
    labs(x="", y = "Fraction")
ggsave("CAF_CellType_sampletype_fraction.pdf",p,width=3.5,height=3.5)

#########
future::plan("multicore", workers = 10)
library(dplyr)
dt@active.ident<-factor(dt$CellType)
dt.markers <- FindAllMarkers(dt,features=rownames(dt)[-which(rownames(dt) %in% badgenelt)], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dt.markers$pct_diff<-dt.markers$pct.1-dt.markers$pct.2
dt.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) ->top50
saveRDS(dt.markers,"CAF_CellType_diffgenes.RDS")

##DotPlot
dt.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) ->top3
diff3genes=c(
  'CD36','PDZD2','PDE1C',
  'IGFBP3','CXCL12','CXCL14',
  'MYH11','MUSTN1','RGS6',
  'LUM','POSTN','VCAN',
  'SEPT7','NEAT1','ALDOA',
  'STMN1','TOP2A','CENPF'
)
p<-DotPlot(dt,cluster.idents=F,group.by="CellType",features=diff3genes,cols=c("darkblue","yellow"))+scale_color_gradientn(colours=rev(heatmaply::RdBu(100)))+
    xlab("")+ylab("")+ggpubr::rotate_x_text(angle=45)
ggsave("./CAF_CellType_top3_diffgenes_dotplot.pdf",p,width=8,height=4)

############################
#CAF signatures Heatmap
############################
VascularDevelopment=c('MYH11','DSTN','MUSTN1','ADIRF')
ECM=c('POSTN','COL6A3','FN1','MMP14','MMP11','DCN','LUM')
LipidProcess=c('APOC3','APOC1','APOA2','APOC2','FABP1')
AP=c('CD74','HLA-DRB1','CXCL12','CCL21','HLA-DRA')
Proliferation=c('STMN1','TOP2A','MKI67','CCNB1','TYMS','HMGB2','PCNA','DEK')
Inflammation=c("DPT",'CXCL9','CXCL10','CXCL11',"EFEMP1","PDPN","PDGFRA")

#p<-DotPlot(dt,cluster.idents=F,group.by="CellType",features=Inflammation,cols=c("darkblue","yellow"))+scale_color_gradientn(colours=rev(heatmaply::RdBu(100)))+
#    xlab("")+ylab("")+ggpubr::rotate_x_text(angle=45)
#ggsave("./CAF_CellType_Inflammation_dotplot.pdf",p,width=8,height=4)

filter_genes=c(Inflammation,VascularDevelopment,ECM,LipidProcess,AP,Proliferation)
dt$type<-paste0(dt$CellType,"_",dt$SampleType)
order_type<-c()
for(i in levels(factor(dt$CellType))){
    order_type<-c(order_type,paste0(i,"_",c("Tumor","Normal")))
}
mat<-as.matrix(dt@assays$RNA@data)
a<-sapply(filter_genes,function(x){
    sapply(order_type,function(y){
        mean(as.numeric(as.vector(mat[x,which(dt$type==y)])))
    })
})
t(a)->a

annotation_col = data.frame(
  Sample = factor(rep(c("Tumor","Normal"), 6),levels=c("Tumor","Normal")), 
  CellType=factor(rep(levels(factor(dt$CellType)),rep(2,6)),levels=levels(factor(dt$CellType)))
)
rownames(annotation_col) = colnames(a)

samplecol<-c('#99CCFF','#CC0066')
names(samplecol)<-c('Normal','Tumor')

ann_colors = list(
  Sample = samplecol,
  CellType=CAFcol
)

library(pheatmap)
pdf("dt_Cluster_sampletype_6signature_pheatmap.pdf",width=10,height=12)
pheatmap(a,scale="row",cellwidth=10,cellheight=10,border_color=NA,annotation_col=annotation_col,annotation_colors = ann_colors,angle_col=315,cluster_rows=F,cluster_cols=F,show_rownames=T,show_colnames=T,color=rev(heatmaply::RdBu(100)))
dev.off()

##Enrichment
dt.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) ->top100

dir.create('./Enrich')
setwd('./Enrich')
source("~/Biotools/Enrich_Radarchart.R")
for(i in 1:length(levels(top100$cluster))){
  genes=top100$gene[top100$cluster==levels(top100$cluster)[i] & top100$p_val_adj<0.001]
  OneGroup_Radarchart(genelt=genes,enrichFunc=c('enrich'),db="human",DataBase=c('KEGG'),return_enrich_rs=F,
            col=CAFcol[i],savedir="./",savename=paste0(levels(top100$cluster)[i],"_enrichKEGG_Radarchart.pdf"))
}
for(i in 1:length(levels(top100$cluster))){
  genes=top100$gene[top100$cluster==levels(top100$cluster)[i] & top100$p_val_adj<0.001]
  OneGroup_Radarchart(genelt=genes,enrichFunc=c('enrich'),db="human",DataBase=c('GO'),return_enrich_rs=F,
            col=CAFcol[i],savedir="./",savename=paste0(levels(top100$cluster)[i],"_enrichGO_Radarchart.pdf"))
}

