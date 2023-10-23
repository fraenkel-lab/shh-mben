####################################################################################
## Max Gold
## Code to look at gene sets and markers from published data
####################################################################################

#############################################
## Functions
#############################################

library(Seurat)

## data already comes log2 normalized so start with variable features and then scaling, etc
process_hov <- function(file){
	data <- read.table(file, row.names=1, sep='\t')
	data <- CreateSeuratObject(data)
	data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2500, verbose = FALSE)
	data <- ScaleData(data, verbose=FALSE)
	data <- RunPCA(data, verbose = FALSE)
	data <- RunUMAP(data, reduction = 'pca',dims = 1:40) 
	data <- FindNeighbors(data, reduction = "pca", dims = 1:40)
	data <- FindClusters(data, resolution = 0.3)	
	return(data)
}

process_counts_rie <- function(pbmc){
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nCount_RNA > 500 & nCount_RNA < 40000 & percent.mt < 30) ### criteria used in the paper
	return(pbmc)
}

normalize_samp_rie <- function(mben.big, res,  nn, har = FALSE, hcat='sample'){
	DefaultAssay(mben.big) <- 'RNA'
	mben.big <- NormalizeData(mben.big)
	mben.big <- FindVariableFeatures(mben.big, selection.method = "vst", nfeatures = nn, verbose = FALSE) ## choose 4000 to match them
	mben.big <- ScaleData(mben.big,verbose = FALSE)
	mben.big <- RunPCA(mben.big, verbose = FALSE)

	if (har == TRUE){
		mben.big <- RunHarmony(mben.big, hcat, theta=1.5)
		mben.big <- RunUMAP(mben.big, reduction = 'harmony',dims = 1:50)
		mben.big <- FindNeighbors(mben.big, reduction = "harmony", dims = 1:50)
	} else {
		mben.big <- RunUMAP(mben.big, reduction = 'pca',dims = 1:30)
		mben.big <- FindNeighbors(mben.big, reduction = "pca", dims = 1:30)
	}

	mben.big <- FindClusters(mben.big, resolution = res)
	mben.big <- AddMetaData(mben.big, paste0(mben.big$sample, '-', mben.big$seurat_clusters), 'orig_cluster')
	mben.big <- AddMetaData(mben.big, mben.big$seurat_clusters, 'ooo')
	return(mben.big)
}


make_plot <- function(marker.genes, data, output){
	point.size <- 0.5
	bftheme <- theme(plot.title=element_text(size=hs, hjust=0.5,  face="bold"), axis.title.x =element_text(size=bs), axis.title.y =element_text(size=bs), legend.text=element_text(size=bs), axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5), axis.text.x = element_text(size=bs-4,  colour = "black"), axis.text.y = element_text(size=bs-4,colour = "black")) 	
	plot.list <- list()
	lgray <- CustomPalette(low = "lightgrey", high = "lightgrey", mid = NULL, k = 5)

	## if there are all 0's, then just make everyhting light gray
	for (x in marker.genes){
		mx = max(data[x,]@assays$RNA@counts)
		if (mx != 0){
			tplot = reset_featplot_gene(FeaturePlot(data, c(x), pt.size=point.size, min.cutoff='q20') , x) +  bftheme
		} else {
			tplot = reset_featplot_gene(FeaturePlot(data, c(x), pt.size=point.size, min.cutoff='q20', cols=lgray) , x) +  bftheme
		}	
	}
	plot.list[[x]] <- tplot
	tum.featplot <- grid.arrange(grobs=pl, nrow=2)
	ggsave(tum.featplot, file = paste0('png_images/', output, '_tum-markers.png'), width = 21, height = 12, dpi = 300, units = "in", device='png')
	ggsave(tum.featplot, file = paste0('pdf_images/', output, '_tum-markers.pdf'), width = 21, height = 12, dpi = 300, units = "in", device='pdf')
}



#############################################
## Load Data
#############################################

### Two SHH Tumors from Vladoiu 2019 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6675628/)

load("vlad_data/SM4217_SHH_res.0.1.Robj")
sm4217 <- UpdateSeuratObject(scRNAseq)
rm(scRNAseq)

load("vlad_data/SHH_BT2017017_res0.1.Robj")
bt2017 <-  UpdateSeuratObject(scRNAseq)
rm(scRNAseq)


### Three SHH tumors from Hovestadt 2019 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6754173/)
muv41 <- process_hov('hov_data/GSM3905417_MUV41.txt.gz')
sj577 <- process_hov('hov_data/GSM3905425_SJ577.txt.gz')
sj454 <- process_hov('hov_data/GSM3905423_SJ454.txt.gz')


### 9 SHH Tumors from Riemondy 2022 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8804892/)

rcounts <- fread('riemondy/GSE155446_human_raw_counts.csv.gz', sep=',',header=TRUE)
rcounts <- as.matrix(rcounts, rownames=1)
meta <- read.csv('riemondy/GSE155446_human_cell_metadata.csv.gz', row.names=1, header=TRUE, sep=',')
rdata <- CreateSeuratObject(rcounts, assay='RNA', project='Riemondy', meta.data=meta)
rm(rcounts)
riemondy <- subset(rdata, subgroup == 'SHH')
rm(rdata)
riemondy <- AddMetaData(riemondy, 'riemondy', col.name = 'dataset')
rt <- subset(riemondy, coarse_cell_type == 'malignant')
rm(riemondy)

rt <- process_counts_rie(rt)
rt <- AddMetaData(rt, rt@meta.data$geo_sample_id, 'sample')
rt <- normalize_samp_rie(rt, 0.3, 4000, har=TRUE, hcat='sample')

### Load BT2019
btg <- readRDS("fdata/BT9_tum.rds")


#############################################
## Make Plots
#############################################

marker.genes <- c('TOP2A', 'SFRP1', 'STMN2', 'SEMA6A', 'GABRD', 'VSNL1')

## Hovestadt (Supp Figures 9, 10, & 11)
make_plot(marker.genes, muv41, 'hov-muv41')
make_plot(marker.genes, sj577, 'hov-sj577')
make_plot(marker.genes, sj454, 'hov-sj454')

## Vladoiu (Supp Figures 12, & 13)
make_plot(marker.genes, sm4217, 'vlad-sm4217')
make_plot(marker.genes, bt2017, 'vlad-bt2017')

## Riemondy (Supp Figure 14)
make_plot(marker.genes, rt, 'riemondy-all')

## BT2019 MBEN (Supp Figure 15)
make_plot(marker.genes, btg, 'MBEN-bt2019')
