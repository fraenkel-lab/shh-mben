####################################################################################
## Max Gold
## Code to process snRNA-seq for validation cohort
####################################################################################

#############################################
## Import packages
#############################################
library(data.table)
library(ggplot2)
library(gridExtra)
library(harmony)
library(Seurat)
library(stringr)

#######################################
## general methods
#######################################

import_data <- function(folder, data_file, assay, species, project){
	if (project == 'taylor'){
		counts <- Read10X_h5(paste0(folder,data_file), use.names = TRUE, unique.features = TRUE)
		data <- CreateSeuratObject(counts, assay='RNA', min.cells = 5, min.features = 300)
		samp <- strsplit(strsplit(data_file, '_')[[1]][5], '\\.')[[1]][1]
	} else if (project == 'mouse'){
		counts <- Read10X(paste0(folder,data_file))
		data <- CreateSeuratObject(counts, assay='RNA', min.cells = 5, min.features = 300)
		samp <- data_file
	}

	if (assay == 'cell') {mt.cutoff <- 25} else {mt.cutoff <- 5}
	## add scrublet features
	scrub <- read.csv(paste0('full_scrublet_files/', samp, '.csv'), row.names=1)

	### add metadata
	data <- AddMetaData(data, samp, col.name = 'sample')
	data <- AddMetaData(data, scrub)
	data <- AddMetaData(data, assay, col.name = 'assay')
	data <- AddMetaData(data, species, col.name = 'species')
	data <- AddMetaData(data, paste0(project,'-',assay), col.name='project')

	if (species == 'mouse'){
		mtv <- PercentageFeatureSet(data, pattern = "^mt-")
	} else if (species == 'human') {
		mtv <- PercentageFeatureSet(data, pattern = '^MT-')
	}
	data <- AddMetaData(data, mtv, col.name='percent.mt')

	# ## subset the data
	data <- get_quantile_cutoffs(data)
	cdata <- subset_data(data, mt.cutoff)
	return(cdata)
}

get_quantile_cutoffs <- function(data){
	cq5 <- quantile(data@meta.data$nCount_RNA, 0.05)[[1]]
	cq95 <- quantile(data@meta.data$nCount_RNA, 0.95)[[1]]
	fq5 <- quantile(data@meta.data$nFeature_RNA, 0.05)[[1]]
	fq95 <- quantile(data@meta.data$nFeature_RNA, 0.95)[[1]]

	cs <- sapply(data@meta.data$nCount_RNA, function(x) between(x, cq5, cq95))
	fs <- sapply(data@meta.data$nCount_RNA, function(x) between(x, fq5, fq95))
	data <- AddMetaData(data, cs, col.name = 'ncount_call')
	data <- AddMetaData(data, cs, col.name = 'nfeat_call')
	return(data)
}

# ### get high quality samples
subset_data <- function(data, mt.cutoff){
	cdata <- subset(data, subset = (nfeat_call == TRUE) & (ncount_call == TRUE) & (Call90 == 'Good') & (percent.mt < mt.cutoff))
	return(cdata)
}
## normalize, standardize, UMAP, and clustering
process_data <- function(data, res,  nfeat, har = FALSE, hcat='sample', regress=TRUE){
	DefaultAssay(data) <- 'RNA'
	data <- NormalizeData(data)
	data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeat, verbose = FALSE)
	if (regress == TRUE){
		data <- ScaleData(data, verbose=FALSE, vars.to.regress='nCount_RNA')
	} else {
		data <- ScaleData(data, verbose=FALSE)
	}
	data <- RunPCA(data, verbose = FALSE)

	if (har == TRUE){
		data <- RunHarmony(data, hcat)
		data <- RunUMAP(data, reduction = 'harmony',dims = 1:50)
		data <- FindNeighbors(data, reduction = "harmony", dims = 1:50)
	} else {
		data <- RunUMAP(data, reduction = 'pca',dims = 1:40)
		data <- FindNeighbors(data, reduction = "pca", dims = 1:40)
	}

	data <- FindClusters(data, resolution = res)
	return(data)
}

#############################################
## snRNA-seq loading and processing
#############################################
 

## collect old MBEN data and process
main.folder <- "../data/taylor_data/"
main.files <-c("filtered_feature_bc_matrix_mb0595.h5",
			"filtered_feature_bc_matrix_mb0601.h5",  
			"filtered_feature_bc_matrix_mb2112.h5",  
			"filtered_feature_bc_matrix_mb3201.h5",  
			"filtered_feature_bc_matrix_mb3400.h5",  
			"filtered_feature_bc_matrix_mb4113.h5",
			"filtered_feature_bc_matrix_mb3612.h5")
main.list <- lapply(main.files, function(file) import_data(main.folder, file, 'nuc', 'human', 'taylor'))  

## collect new validation data and process
val.folder <- "../data/val_data/"
val.files <-c("filtered_feature_bc_matrix_mb002.h5",
             "filtered_feature_bc_matrix_mb019.h5",
             "filtered_feature_bc_matrix_mb084.h5",
             "filtered_feature_bc_matrix_mb005.h5",
             "filtered_feature_bc_matrix_mb015.h5",
             "filtered_feature_bc_matrix_mb009.h5")

val.list <- lapply(val.files, function(file) import_data(val.folder, file, 'nuc', 'human', 'taylor')) 
      
## merge datasets and put into one big object
mben.list <- c(main.list,v al.list)

mben.big <- merge(mben.list[[1]], y = sapply(2:length(mben.list), function(x) mben.list[x]), project='BigMBEN')
rm(mben.list)
rm(olv.list)
rm(rev.list)

res = 0.25
nfeat = 2500
mben.big <- process_data(mben.big, res, nfeat, har=TRUE)

saveRDS(mben.big, file='fdata/mben_validation_combined.rds')

#############################################
## Plotting
#############################################


## Load MBEN Object
mben.big <- readRDS(file='fdata/mben_validation_combined.rds')


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
## Plotting functions
#############################################

reset_featplot_gene <- function(plot, gene){
	plot$data <- plot$data[order(plot$data[gene]),]
	plot <- plot + ggtitle(gene)
	return(plot)
}

get_sampplot <- function(data, sample, colors){
	cc <- rownames(data@meta.data[data@meta.data$sample == sample, ])
	dp <- DimPlot(data, label=TRUE, label.size=5, cells=cc) + NoLegend() + ggtitle(sample) + theme(plot.title=element_text(size=hs, hjust = 0.5), axis.title.x =element_text(size=ss), axis.title.y =element_text(size=ss), legend.text=element_text(size=ss))
	return(dp)
}


#############################################
## Validation Plots
#############################################

### Sample Plots for Two Representatives From MBEN Cohort
samples <- c('mb4113', 'mb0601')
si <- lapply(samples, function(s) get_sampplot(mben.big, s, scol(cl)))
mben.orig.dimplot.split <- grid.arrange(grobs=si, nrow=2)
ggsave(mben.orig.dimplot.split, file = 'png_images/mben_splitplot.png', width = 3, height = 6, dpi = 300, units = "in", device='png')
ggsave(mben.orig.dimplot.split, file = 'pdf_images/mben_splitplot.pdf', width = 3, height = 6, dpi = 300, units = "in", device='pdf')


### SHHa Plots
samples <- c('mb002', 'mb009', 'mb019')
si <- lapply(samples, function(s) get_sampplot(mben.big, s, scol(cl)))
shha.dimplot.split <- grid.arrange(grobs=si, nrow=1)
ggsave(shha.dimplot.split, file = 'png_images/shha_splitplot.png', width = 9, height = 3, dpi = 300, units = "in", device='png')
ggsave(shha.dimplot.split, file = 'pdf_images/shha_splitplot.pdf', width = 9, height = 3, dpi = 300, units = "in", device='pdf')


### SHHb Plots
samples <- c('mb005', 'mb015', 'mb084')
si <- lapply(samples, function(s) get_sampplot(mben.big, s, scol(cl)))
shhb.dimplot.split <- grid.arrange(grobs=si, nrow=1)
ggsave(shhb.dimplot.split, file = 'png_images/shhb_splitplot.png', width = 9, height = 3, dpi = 300, units = "in", device='png')
ggsave(shhb.dimplot.split, file = 'pdf_images/shhb_splitplot.pdf', width = 9, height = 3, dpi = 300, units = "in", device='pdf')


### Marker Genes
marker.genes.supp <- c('GLI2','SEMA6A', 'GABRD')
qpl <- lapply(marker.genes.supp , function(x) reset_featplot_gene(FeaturePlot(mben.big, c(x), pt.size=point.size, min.cutoff='q80'), x) + bftheme )
tum.featplot.extra <- grid.arrange(grobs=qpl, nrow=1)
ggsave(tum.featplot.extra, file = 'png_images/integrated_marker.png', width = 12, height = 4, dpi = 300, units = "in", device='png')
ggsave(tum.featplot.extra, file = 'pdf_images/integrated_marker.pdf', width = 12, height = 4, dpi = 300, units = "in", device='pdf')

