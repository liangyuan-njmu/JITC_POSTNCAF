#############################################################
#scanpy deal
`python`
import anndata as ad
import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

adata = sc.read_h5ad('./All.h5ad')
adata.raw = adata
adata.layers["counts"] = adata.raw.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.highly_variable_genes(adata, min_mean=0.0125,n_top_genes=4000, max_mean=3, min_disp=0.25)
sc.tl.pca(adata, svd_solver='arpack')
sce.pp.harmony_integrate(adata, 'orig.ident')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

sc.pl.umap(adata,color="leiden",legend_loc='on data',title='Harmony Leiden',save='.pdf')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False,save='.pdf')

Tcells<-c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B")
BCells<-c("CD79A",'MS4A1','MZB1')
Plasma<-c("MZB1","CD38","PRDM1","IGHG1","IGHG2","IGHG3","JCHAIN","IGLC3")
DC<-c("CD1C","IDO1","CLEC9A","CLEC10A",'ITGAX','CD83','CD86')
cDC1<-c("CLEC9A")
cDC2<-c('CD1C','CLEC10A')
Monocyte<-c('CD14',"FCN1","CEBPD",'CLEC4A')
pDC<-c('LILRA4',"JCHAIN")
mDC<-c('CCR7',"LAMP3")
Endothelial=c('PECAM1','CLDN5','VWF')
Epithelial=c('KRT19','CK18','EPCAM','KRT8','CLDN4')
Cholangiocytes=c('SOX9','LCN2','ELF3','CXCL1','FXYD2')
Hepatocytes=c('CYP2A7',"ALB",'AFP','APOA2')
Macrophage=c('CD163','MRC1','CSF1R','CD68')
NK=c('KLRD1','NKG7','PRF1','KLRC1')
Fibroblast=c('COL1A1','DCN','COL1A2','COL3A1')
Mast=c('CPA3','KIT','TPSAB1','TPSB2','MS4A2')
gdT=c('TRDC')
Megakaryocytes=c('GNG11')
Proliferation=c("TOP2A","MKI67","STMN1")
Neutrophils=c('CSF3R','S100A8','S100A9')

new_cluster_names = {
  '9':'Fibroblast',
  '0':'T/NK',
  '1':'T/NK',
  '7':'T/NK',
  '26':'T/NK',
  '5':'T/NK',
  '8':'T/NK',
  '24':'T/NK',
  '18':'T/NK',
  '10':'T/NK',
  '20':'T/NK',
  '22':'T/NK',
  '12':'BCells',
  '13':'PlasmaCells',
  '27':'PlasmaCells',
  '6':'Myeloid',
  '14':'Myeloid',
  '2':'Myeloid',
  '21':'Myeloid',
  '23':'Myeloid',
  '17':'Myeloid',
  '32':'Myeloid',
  '4':'Endothelial',
  '15':'Endothelial',
  '11':'Epithelial',
  '29':'Epithelial',
  '19':'Epithelial',
  '16':'Epithelial',
  '3':'Epithelial',
  '30':'Epithelial',
  '31':'Epithelial',
  '28':'RBC',
  '25':'Myeloid'
}
adata.obs['CellType'] = adata.obs['leiden'].map(new_cluster_names).astype('category')

adata = adata[-adata.obs['CellType'].isin(['RBC'])]

sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')
import matplotlib as plt
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": (300)}):sc.pl.umap(adata,color="CellType",legend_loc=None,title='CellType',save='CellType.pdf')
adata.write("./All_HCC.h5ad")

adata_subset = adata[adata.obs['CellType'].isin(['Endothelial'])]
adata_subset.write("./HCC_Endothelial.h5ad")

adata_subset = adata[adata.obs['CellType'].isin(['BCells','PlasmaCells'])]
adata_subset.write("./HCC_B.h5ad")

adata_subset = adata[adata.obs['CellType'].isin(['Myeloid','MastCells'])]
adata_subset.write("./HCC_Myeloid.h5ad")

adata_subset = adata[adata.obs['CellType'].isin(['Epithlial'])]
adata_subset.write("./HCC_Epithelial.h5ad")

marker_genes=[
  'CD79A','CD79B',
  'PECAM1','VWF',
  'ALB','APOA2',
  'COL1A1','COL1A2',
  'LYZ','C1QB','CD1C','S100A9',
  'MZB1','IGHG1',
  'CD3D','CD3E'
]
sc.pl.stacked_violin(adata, marker_genes, groupby='CellType', rotation=90,save="marker_violin.pdf")

sc.pl.dotplot(adata, marker_genes, groupby='CellType',save="Dotplot.pdf")

#
df=data.frame(
  value=Cell_Numbers,
  celltype=c('BCells','Endothelial','Epithelial','Fibroblast','Myeloid','PlasmaCells','TNK')
)
df2=df %>% arrange(desc(value))
df$celltype=factor(df$celltype,levels=rev(df2$celltype))
p=ggplot(df,aes(x=celltype,y=value))+
    geom_bar(stat = 'identity',fill='#80E6EC',width=0.5)+
    coord_flip() + 
    theme_bw()+theme(text = element_text(size = 15),axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
    labs(x="", y = "Cell Nums")
ggsave("CellType_CellNumbers.pdf",p,width=3,height=3)

df=data.frame(
  value=Cell_Fractions,
  sampletype=factor(rep(c('Tumor','Normal'),7),levels=c('Normal','Tumor')),
  celltype=factor(rep(c('BCells','Endothelial','Epithelial','Fibroblast','Myeloid','PlasmaCells','TNK'),rep(2,7)),levels=)
)
df$celltype=factor(df$celltype,levels=rev(df2$celltype))
p=ggplot(df,aes(x=celltype,y=value,fill=sampletype))+
    geom_bar(stat = 'identity',width=0.5)+
    coord_flip() + 
    scale_fill_manual(values=c('Tumor'='#FF0200','Normal'='#1600FF'))+
    theme_bw()+theme(text = element_text(size = 15),axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
    labs(x="", y = "Fraction")
ggsave("CellType_sampletype_fraction.pdf",p,width=3,height=4)

