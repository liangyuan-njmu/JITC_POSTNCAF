#############################################################
#Seurat deal example one ST
library(SpaGene)
#P1T
P1T=readRDS('JH_TIB/P1T_Spatial.rds')

count=P1T@assays$Spatial@counts
location=GetTissueCoordinates(P1T,image='image_P11_T')

sg<-SpaGene(count,location)
P1T_pattern<-FindPattern(sg,nPattern=10)

pdf('P1T_pattern.pdf',width=24,height=20)
PlotPattern(P1T_pattern,location)
dev.off()

#P1T
diet.seurat = Seurat::DietSeurat(P1T, graphs = "pca") #slim down Seurat obj prior to conversion
sce = as.SingleCellExperiment(diet.seurat) #convert seurat to SCE
#P1T@images$image@coordinates[,c(4,5)]=GetTissueCoordinates(P1T)
colData(sce) = cbind(colData(sce), P1T@images$image_P11_T@coordinates) #add spatial info to SCE

set.seed(12345)
sce <- spatialPreprocess(sce, platform="ST", 
                              n.PCs=30, n.HVGs=2000, log.normalize=TRUE)

sce <- qTune(sce, qs=seq(2, 20), platform="ST", d=7)
pdf('P1T_qTune.pdf')
qPlot(sce)
dev.off()

sce <- spatialCluster(sce, q=10, platform="ST", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)

sce.enhanced <- spatialEnhance(sce, q=10, platform="ST", d=7,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

markers <- c('COL1A1',"POSTN",'CD68','SPP1', 'CD3D','CD8B')
sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                     feature_names=markers,
                                     nrounds=0)

plot_list= purrr::map(markers, function(x) featurePlot(sce.enhanced, x,is.enhanced=T))
pdf('P1T_6genes_enhanced.pdf',width=30,height=8)
patchwork::wrap_plots(plot_list, ncol=6)
dev.off()

markers <- c('CD4',"CD19",'GZMA','GNLY')
sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                     feature_names=markers,
                                     nrounds=0)

plot_list= purrr::map(markers, function(x) featurePlot(sce.enhanced, x,is.enhanced=T))
pdf('P1T_4genes_enhanced.pdf',width=20,height=8)
patchwork::wrap_plots(plot_list, ncol=4)
dev.off()
