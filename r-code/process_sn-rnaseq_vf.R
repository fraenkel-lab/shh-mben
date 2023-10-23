####################################################################################
## Max Gold
## Code to run through sn/sc-RNA-Seq clustering and signature analysis
####################################################################################

# Import packages
library(biomaRt)
library(data.table)
library(ggplot2)
library(gridExtra)
library(harmony)
library(monocle3)
library(RColorBrewer)
library(readxl)
library(Seurat)
library(SeuratWrappers)
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
	scrub <- read.csv(paste0('scrublet_files/', samp, '.csv'), row.names=1)

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

######################################
## Process MBEN-Single-NUC
######################################

taylor.folder <- "../data/taylor_data/"
sn.files <-c("filtered_feature_bc_matrix_mb0595.h5",
			"filtered_feature_bc_matrix_mb0601.h5",  
			"filtered_feature_bc_matrix_mb2112.h5",  
			"filtered_feature_bc_matrix_mb3201.h5",  
			"filtered_feature_bc_matrix_mb3400.h5",  
			"filtered_feature_bc_matrix_mb4113.h5",
			"filtered_feature_bc_matrix_mb3612.h5")

## collect MBEN data and merge
mben.list <- lapply(sn.files, function(file) import_data(taylor.folder, file, 'nuc', 'human', 'taylor'))

mben.big <- merge(mben.list[[1]], y = sapply(2:length(mben.list), function(x) mben.list[x]), project='BigMBEN')
rm(mben.list)

res = 0.25
nfeat = 2500
mben.big <- process_data(mben.big, res, nfeat, har=TRUE)


annot.list <- list('7'= 'Fibroblast', '8'= 'Macrophage', '9' = 'Astrocyte', '10' = 'Oligodendrocyte', '11'= 'Microglia', '13'='Immune-Other')
mben.idents <- as.character(Idents(mben.big))
annots <- as.character(sapply(mben.idents, function(x) if (x %in% names(annot.list)){annot.list[[x]]} else {'Tumor'}))
mben.big <- AddMetaData(mben.big, annots, col.name = 'CellType')

saveRDS(mben.big, file='fdata/mben_big.rds')


######################################
## Analyze the Tumor Cells
######################################

## ignore 6 (only in MB4113) and 12 (only in MB2112)
mben.tums <- subset(mben.big, seurat_clusters %in% c('0', '1', '2', '3', '4', '5'))
mben.tums <- FindNeighbors(mben.tums, reduction = 'harmony', dims = 1:50)
mben.tums <- RunUMAP(mben.tums, reduction = 'harmony',dims = 1:50)
mben.tums <- FindClusters(mben.tums, resolution = 0.25)


tum.annot.list <- list('0'= 'GN-Postmigratory', '1'= 'GCP-SHH', '2' = 'GN-Migrating', '3' = 'GN-Premigratory', '4'= 'Ribosomal', '5'='GCP-Cycling', '6'= 'GN-Postmigratory')
tum.idents <- as.character(Idents(mben.tums))
tum.annots <- as.character(sapply(tum.idents, function(x) tum.annot.list[[x]] ))
mben.tums <- AddMetaData(mben.tums, tum.annots, col.name = 'TumCellType')

saveRDS(mben.tums, file='fdata/mben_tums.rds')


## Create clean dataset for pseudotime
## Remove few "tumor" cells that cluster in other areas

mip <- subset(mben.big, seurat_clusters %in% c('0', '1', '2', '3', '4', '5'))
gt <- c()
u2 <- as.double(mip@reductions$umap@cell.embeddings[,'UMAP_2'])
u1 <- as.double(mip@reductions$umap@cell.embeddings[,'UMAP_1'])

for (i in 1:length(u1)){
	u1v = u1[i]
	u2v = u2[i]
	if (u1v < -8){
		val <- FALSE
	} else if (u2v < -2){
		if (u1v > 0){
			val <- FALSE
		} else {
			val <- TRUE
		}
	} else {
		val <- TRUE
	}
	gt <- c(gt, val)
}

mip <- AddMetaData(mip, gt, 'TumCell')
rmip <- subset(mip, TumCell == TRUE)

######################################
## Pseudotime Analysis
######################################

## Perform Pseudotime Analysis on clean dataset using monocle
cds <- as.cell_data_set(DietSeurat(rmip, graphs = "umap"))
cds <- cluster_cells(cds)
cds <- learn_graph(cds, learn_graph_control=list(minimal_branch_len=15))
cyc_cells <- subset(rmip, seurat_clusters == '5')
uemb = as.data.frame(cyc_cells@reductions$umap@cell.embeddings)
min.avp <- which.max(uemb$UMAP_1)
min.avp <- rownames(uemb)[min.avp]
cds <- order_cells(cds, root_cells = min.avp)
ptime <- as.numeric(cds@principal_graph_aux$UMAP$pseudotime)
rmip <- AddMetaData(rmip, ptime, 'Pseudotime')
saveRDS(rmip, file='fdata/mben_tums_pseudotime.rds')
saveRDS(cds, file="fdata/mben_bigtums_ptime.rds")

######################################
## Process bt2019110 (scRNA-Seq)
######################################

bt9_file <- 'filtered_feature_bc_matrix_bt2019110.h5'
bt9 <- import_data(taylor.folder,bt9_file, 'cell', 'human', 'taylor')
bt9n <- process_data(bt9, 0.5, 2500, har=FALSE, regress=FALSE)
saveRDS(bt9n, file='fdata/BT9_big.rds')


annot.list.bt9 <- list('2'='Myeloid', '3'='Myeloid', '9'='Astrocyte', '10'='Oligodendrocyte', '11'='Immune-Other',  '12'='Fibroblast', '13'='Microglia')
bt9.idents <- as.character(Idents(bt9n))
bt9.annots <- as.character(sapply(bt9.idents, function(x) if (x %in% names(annot.list.bt9)){annot.list.bt9[[x]]} else {'Tumor'}))
bt9n <- AddMetaData(bt9n, bt9.annots, col.name = 'CellType')


btg<- subset(bt9n, seurat_clusters %in% c('0', '1', '4', '5', '6', '7','8'))
btg <- FindNeighbors(btg, reduction = 'pca', dims = 1:40)
btg <- RunUMAP(btg, reduction = 'pca',dims = 1:40)
btg <- FindClusters(btg, resolution = 0.2)

saveRDS(btg, file='fdata/BT9_tum.rds')


######################################
## Process P14 Mouse From Vladoiu 2019
######################################

mouse.folder <- 'taylor_mouse/'
p14 <- import_data(mouse.folder,'P14', 'cell', 'mouse', 'mouse')
p14n <- process_data(p14, 0.5, 2500, har=FALSE, regress=FALSE)
saveRDS(p14n, 'fdata/p14_all.rds')


p14g <- subset(p14n, seurat_clusters %in% c('0','2', '3', '4', '5', '9','11'))
p14g <- FindNeighbors(p14g, reduction = 'pca', dims = 1:40)
p14g <- RunUMAP(p14g, reduction = 'pca',dims = 1:40)
p14g <- FindClusters(p14g, resolution = 0.3)


mouse.annot.list <- list('1'= 'GCP', '2'= 'GCP-Cycling', '4'='GCP-Cycling', '6'='GCP-Cycling',  '0' = 'GN-Premigratory', '3'= 'GN-Migrating', '5'= 'GN-Postmigratory')
mouse.idents <- as.character(Idents(p14g))
mouse.annots <- as.character(sapply(mouse.idents, function(x) mouse.annot.list[[x]] ))
p14g <- AddMetaData(p14g, mouse.annots, col.name = 'GnCellType')

saveRDS(p14g, 'fdata/p14_gn.rds')

