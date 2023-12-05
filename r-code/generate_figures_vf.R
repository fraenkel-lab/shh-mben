####################################################################################
## Max Gold
## Code to generate figures from data
####################################################################################

#############################################
## Import packages
#############################################
library(ggplot2)
library(gridExtra)
library(monocle3)
library(RColorBrewer)
library(Seurat)
library(SeuratWrappers)
library(stringr)

#############################################
## Load Data
#############################################

mben.big <- readRDS('fdata/mben_big.rds')
mben.tums <- readRDS('fdata/mben_tums.rds')
cds <- readRDS('fdata/mben_bigtums_ptime.rds')
rmip <- readRDS('fdata/mben_tums_pseudotime.rds')
mtf <- subset(mben.tums, TumCellType=='Ribosomal', invert=TRUE)
bt9n <- readRDS('fdata/BT9_big.rds')
btg <- readRDS('fdata/BT9_tum.rds')
p14g <- readRDS('fdata/p14_gn.rds')

#############################################
## Set parameter variables
#############################################

big.point.size = 2
point.size <- 0.5
small.point.size <- 0.1

hs <- 20
bs <- 16
ts <- 12
ss <- 10

ftheme <- theme(plot.title=element_text(size=ts, hjust=0.5,  face="bold"), axis.title.x =element_text(size=ss), axis.title.y =element_text(size=ss), legend.text=element_text(size=ss), axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5), axis.text.x = element_text(size=ss-4,  colour = "black"), axis.text.y = element_text(size=ss-4,colour = "black")) 
bftheme <- theme(plot.title=element_text(size=hs, hjust=0.5,  face="bold"), axis.title.x =element_text(size=bs), axis.title.y =element_text(size=bs), legend.text=element_text(size=bs), axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5), axis.text.x = element_text(size=bs-4,  colour = "black"), axis.text.y = element_text(size=bs-4,colour = "black")) 


#############################################
## Cells to be removed for mtf plots
#############################################

# remove very small outlier cells that cluster with ribosomal set in UMAP but assigned to other groups
oo <- as.data.frame(mtf@reductions$umap@cell.embeddings)
oof <- oo[((oo$UMAP_1 >-6.5) & (oo$UMAP_1 <0) & (oo$UMAP_2 <0) ),]
gsamp <- setdiff(rownames(oo), rownames(oof))

#############################################
## Plotting functions
#############################################

reset_featplot <- function(plot, gene){
	plot$data <- plot$data[order(plot$data[gene]),]
	return(plot)
}

reset_featplot_gene <- function(plot, gene){
	plot$data <- plot$data[order(plot$data[gene]),]
	plot <- plot + ggtitle(gene)
	return(plot)
}

get_sampplot <- function(data, sample, colors){
	cc <- rownames(data@meta.data[data@meta.data$sample == sample, ])
	dp <- DimPlot(data, label=TRUE, label.size=5, cells=cc) + NoLegend() + ggtitle(sample) + theme(plot.title=element_text(size=bts, hjust = 0.5), axis.title.x =element_text(size=ss), axis.title.y =element_text(size=ss), legend.text=element_text(size=ss))
	return(dp)
}

####################
## Figure 1
####################

### Figure 1A  (Cell Types)

bcol = brewer.pal(length(unique(mben.big@meta.data$CellType)), 'Dark2')
bcol[7] = '6BAED6' ## change last color to blue
mben.big.dimplot.annot <- DimPlot(mben.big, group.by = 'CellType', cols=bcol)+ ggtitle(element_blank()) + ftheme+ theme(legend.text=element_text(size=ss-4), legend.key.size = unit(0.15, 'cm'), legend.key.width = unit(0.05, 'cm'), legend.key.height = unit(0.5, 'cm'))
ggsave(mben.big.dimplot.annot, file = 'png_images/mben_allct.png', width = 3.75, height = 3, dpi = 300, units = "in", device='png')
ggsave(mben.big.dimplot.annot, file = 'pdf_images/mben_allct.pdf', width = 3.75, height = 3, dpi = 300, units = "in", device='pdf')


### Figure 1B  (Cell Types and Pseudotime)

trajectory.plot <- plot_cells(cds, label_branch_points = FALSE, label_leaves = FALSE, label_roots=FALSE, color_cells_by='pseudotime', cell_size = 0.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + NoLegend() +ggtitle("Pseudotime")  +xlab('UMAP 1') + ylab('UMAP 2') + ftheme + theme(axis.title.x =element_text(size=ss-4), axis.title.y =element_text(size=ss-4), plot.title=element_text(size=ss-2, hjust=0.5,  face="bold"))

neun.plot <- FeaturePlot(rmip, c('RBFOX3'), min.cutoff='q5')+ NoLegend()+ggtitle("NeuN") +xlab('UMAP 1') + ylab('UMAP 2')+ NoLegend() +ftheme + theme(axis.title.x =element_text(size=ss-2), axis.title.y =element_text(size=ss-2), plot.title=element_text(size=ss, hjust=0.5,  face="bold"))
gli2.plot <- FeaturePlot(rmip, c('GLI2'), min.cutoff='q5') + NoLegend() +ggtitle("GLI2") +xlab('UMAP 1') + ylab('UMAP 2')+ NoLegend()+ ftheme + theme(axis.title.x =element_text(size=ss-2), axis.title.y =element_text(size=ss-2), plot.title=element_text(size=ss, hjust=0.5,  face="bold"))
top2a.plot <- FeaturePlot(rmip, c('TOP2A'), min.cutoff='q5')+ NoLegend() +ggtitle("TOP2A") +xlab('UMAP 1') + ylab('UMAP 2')+ NoLegend() + ftheme + theme(axis.title.x =element_text(size=ss-2), axis.title.y =element_text(size=ss-2), plot.title=element_text(size=ss, hjust=0.5,  face="bold"))
tums.plot <- grid.arrange(top2a.plot, gli2.plot,neun.plot,trajectory.plot, nrow=2)
ggsave(tums.plot, file = 'png_images/mben_tum_mark_ptime.png', width = 3, height = 3, dpi = 300, units = "in", device='png')
ggsave(tums.plot, file = 'pdf_images/mben_tum_mark_ptime.pdf', width = 3, height = 3, dpi = 300, units = "in", device='pdf')


### Figure 1D (Marker Genes for Each Stage)

bus = 1
lorder <- c('GCP-Cycling', 'GCP-SHH', 'GN-Premigratory', 'GN-Migrating', 'GN-Postmigratory')
corder <- c('red2',  'gold3', 'green4', 'deepskyblue4', 'mediumorchid')
mben.dimplot <- DimPlot(mtf,group.by='TumCellType',cells=gsamp, pt.size=point.size ) + NoLegend() + ggtitle("MBEN Tumor Nuclei") +ftheme + theme(axis.title.x =element_text(size=ss+bus), axis.title.y =element_text(size=ss+bus), plot.title=element_text(size=ss+bus+4, hjust=0.5,  face="bold"))
mben.dimplot$data$TumCellType <- factor(x=mben.dimplot$data$TumCellType, levels = lorder)
mben.dimplot <- mben.dimplot + scale_color_manual(labels= lorder, values =corder)

marker.genes <- c( 'TOP2A', 'GLI2', 'SEMA6A', 'GRIN2B', 'GRIN2C')
pl <- lapply(marker.genes, function(x) reset_featplot_gene(FeaturePlot(mtf, c(x), cells=gsamp,pt.size=point.size, min.cutoff='q80' ) , x) + NoLegend() +ftheme +  theme(axis.title.x =element_text(size=ss+bus), axis.title.y =element_text(size=ss+bus), plot.title=element_text(size=ss+bus+4, hjust=0.5,  face="bold"))) 
rpl <- c(list('celltype'=mben.dimplot), pl)
tum.featplot <- grid.arrange(grobs=rpl, nrow=2)

ggsave(tum.featplot, file = 'png_images/tum_granule_dev.png', width = 6, height = 4, dpi = 300, units = "in", device='png')
ggsave(tum.featplot, file = 'pdf_images/tum_granule_dev.pdf', width = 6, height = 4, dpi = 300, units = "in", device='pdf')


####################
## Figure 4
####################

### Figure 4B (SHHb signature plot)
plot_shhb <- FeaturePlot(mtf, features=c('SHHb1'), min.cutoff=0) + ggtitle("SHHb signature") + ftheme
ggsave(plot_shhb, file = 'png_images/sc_shhb_sig.png', width = 4, height = 3, dpi = 300, units = "in", device='png')
ggsave(plot_shhb, file = 'png_images/sc_shhb_sig.pdf', width = 4, height = 3, dpi = 300, units = "in", device='pdf')

## Figure 4C (plot individual samples to organize later)
samples <- unique(mtf@meta.data$sample)

for (samp in samples){
	stum <- subset(mtf, sample == samp)
	lorder <- c('GCP-Cycling',  'GCP-SHH', 'GN-Premigratory', 'GN-Migrating', 'GN-Postmigratory')
	corder <- c('red2',  'gold3', 'green4', 'deepskyblue4', 'mediumorchid')
	mben.dimplot <- DimPlot(stum,group.by='TumCellType') + NoLegend() + ggtitle(toupper(samp)) +ftheme + xlim(-10,9)
	mben.dimplot$data$TumCellType <- factor(x=mben.dimplot$data$TumCellType, levels = lorder)
	mben.dimplot <- mben.dimplot + scale_color_manual(labels= lorder, values =corder)
	ggsave(mben.dimplot, file = paste0('png_images/' ,samp, '_gndev.png'), width = 3, height = 3, dpi = 300, units = "in", device='png')
	ggsave(mben.dimplot, file = paste0('pdf_images/' ,samp, '_gndev.pdf'), width = 3, height = 3, dpi = 300, units = "in", device='pdf')
}


####################
## Supp Figures
####################
             
             
### Supp Figure 1A  (Cell Types)

cl <- length(unique(mben.big@meta.data$RNA_snn_res.0.25))
scol = colorRampPalette(brewer.pal(cl, "Accent"))

marker.genes <- c('ALDH1L1', 'CD163', 'COL3A1', 'OLIG1', 'PTPRC', 'VWF')
pl <- lapply(marker.genes, function(x) reset_featplot(FeaturePlot(mben.big, c(x), pt.size=point.size, min.cutoff='q70'), x) + NoLegend() + ftheme
mben.big.featplot <- grid.arrange(grobs=pl, nrow=2)
ggsave(mben.big.featplot, file = 'png_images/tumor_stroma_markers.png', width = 12, height = 6, dpi = 300, units = "in", device='png')
ggsave(mben.big.featplot, file = 'pdf_images/tumor_stroma_markers.pdf', width = 12, height = 6, dpi = 300, units = "in", device='pdf')


### Supp Figure 1B (Clustering by samples)
samples <- unique(mben.big@meta.data$sample)
si <- lapply(samples, function(s) get_sampplot(mben.big, s, scol(cl)))
mben.big.dimplot.split <- grid.arrange(grobs=si, nrow=3)
ggsave(mben.big.dimplot.split, file = 'real_images/big_samp_splitplot.png', width = 12, height = 6, dpi = 300, units = "in", device='png')
ggsave(mben.big.dimplot.split, file = 'real_images/big_samp_splitplot.pdf', width = 12, height = 6, dpi = 300, units = "in", device='pdf')


### Supp Figure 2A (Tumor Clusters)
mbt.dimplot <- DimPlot(mben.tums, label=TRUE, group.by='seurat_clusters') +NoLegend() + ggtitle("Clusters") + bftheme
ggsave(mbt.dimplot, file = 'png_images/mben_tum_clusters.png', width = 6, height = 6, dpi = 300, units = "in", device='png')
ggsave(mbt.dimplot, file = 'pdf_images/mben_tum_clusters.pdf', width = 6, height = 6, dpi = 300, units = "in", device='pdf')

### Supp Figure 2B (Marker Genes For GN Dev)
marker.genes.supp <- c( 'ETV1', 'PRKCB','SLC17A6', 'SLC17A7')
qpl <- lapply(marker.genes.supp , function(x) reset_featplot_gene(FeaturePlot(mtf, c(x), pt.size=point.size, min.cutoff='q80'), x) + NoLegend() +bftheme )
tum.featplot.extra <- grid.arrange(grobs=qpl, nrow=1)
ggsave(tum.featplot.extra, file = 'png_images/tum_granule_dev_extra.png', width = 16, height = 4, dpi = 300, units = "in", device='png')
ggsave(tum.featplot.extra, file = 'pdf_images/tum_granule_dev_extra.pdf', width = 16, height = 4, dpi = 300, units = "in", device='pdf')

### Supp Figure 3A (BT2019 annotation plot)
bt9_dimplot <- DimPlot(bt9n, group.by='CellType')
ggsave(bt9_dimplot, file = 'png_images/bt9n_annot.png', width = 6, height = 6, dpi = 300, units = "in", device='png')
ggsave(bt9_dimplot, file = 'pdf_images/bt9n_annot.png', width = 6, height = 6, dpi = 300, units = "in", device='pdf')

### Supp Figure 3B (Marker genes for BT2019)
marker.genes <- c('TOP2A', 'SFRP1',  'STMN2', 'PRKCB', 'GABRD', 'VSNL1')
pl <- lapply(marker.genes, function(x) reset_featplot(FeaturePlot(btg, c(x), pt.size=point.size, min.cutoff='q40'), x) + NoLegend() + ftheme)
b9.featplot <- grid.arrange(grobs=pl, nrow=2)
ggsave(b9.featplot, file = 'png_images/bt2019_featplot.png', width = 9, height = 6, dpi = 300, units = "in", device='png')
ggsave(b9.featplot, file = 'pdf_images/bt2019_featplot.png', width = 9, height = 6, dpi = 300, units = "in", device='pdf')



### Supp Figure 5A (Mouse Dimplot and features)
p14.dimplot <- DimPlot(p14g, label=FALSE, group.by='GnCellType') + NoLegend()  + ggtitle('P14 Cell Type') + ftheme
p14.dimplot <- LabelClusters(p14.dimplot, id='GnCellType',fontface='bold', color='black',size=3)

marker.genes <- c('Top2a', 'Sfrp1', 'Stmn2', 'Grin2b', 'Gabra6')
pl <- lapply(marker.genes, function(x) reset_featplot(FeaturePlot(p14g, c(x), pt.size=point.size, min.cutoff='q60'), x) + NoLegend() + ftheme)
rpl <- c(list('celltype'=p14.dimplot), pl)
p14g.featplot <- grid.arrange(grobs=rpl, nrow=2)

ggsave(p14g.featplot, file = 'png_images/p14_mouse_markers.png', width = 9, height = 6, dpi = 300, units = "in", device='png')
ggsave(p14g.featplot, file = 'pdf_images/p14_mouse_markers.pdf', width = 9, height = 6, dpi = 300, units = "in", device='pdf')


### Supp Figure 5B (Mouse featureplot for ribosomal genes)
marker.genes <- c('Rps24', 'Rpl10', 'Rpl34')
pl <- lapply(marker.genes, function(x) reset_featplot(FeaturePlot(p14g, c(x), pt.size=point.size, min.cutoff='q80'), x) + NoLegend() + ftheme)
p14g.ribo <- grid.arrange(grobs=pl, nrow=1)
ggsave(p14g.ribo, file = 'png_images/p14_ribo.png', width = 9, height = 3, dpi = 300, units = "in", device='png')
ggsave(p14g.ribo, file = 'pdf_images/p14_ribo.pdf', width = 9, height = 3, dpi = 300, units = "in", device='pdf')
             
             
## Supplementary Figure 17
fmrp_plot <- FeaturePlot(mts, 'fmrp_targets') + ftheme
ggsave(fmrp_plot, file = 'png_images/sc_fmrp_sig.png', width = 3, height = 3, dpi = 300, units = "in", device='png')
ggsave(fmrp_plot, file = 'pdf_images/sc_fmrp_sig.pdf', width = 3, height = 3, dpi = 300, units = "in", device='pdf')             
             
## Supplementary Figure 18
marker.genes <- c('CNTN1','MAP2', 'MKI67',  'VSNL1')
pl <- lapply(marker.genes, function(x) reset_featplot(FeaturePlot(mtf, c(x), pt.size=point.size, min.cutoff='q80'), x) + NoLegend() + ftheme)
cycif.featplot <- grid.arrange(grobs=pl, nrow=1)
ggsave(vsnl1.featplot, file = 'png_images/cycif_markers.png', width = 12, height = 3, dpi = 300, units = "in", device='png')
ggsave(cycif.featplot, file = 'pdf_images/cycif_markers.pdf', width = 12, height = 3, dpi = 300, units = "in", device='png')
