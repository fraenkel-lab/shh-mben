####################################################################################
## Max Gold
## Code to perform GSVA and consensus clustering for Figure 2
####################################################################################

library(biomaRt)
library(ConsensusClusterPlus)
library(GSVA)
library(readxl)
library(Seurat)


get_differential <- function(obj, annot,name){
	clust <- as.character(unique(obj@meta.data[,annot]))
	dl <- sapply(clust, function(c) FindMarkers(obj, ident.1=c, only.pos=TRUE, group.by=annot), simplify = FALSE,USE.NAMES = TRUE)
	names(dl) <- as.character(sapply(names(dl), function(x) paste0(name, '_', x)))
	return(dl)
}

#############################################
## Load Data
#############################################

p14g <- readRDS('fdata/p14_gn.rds')
btg <- readRDS('fdata/BT9_tum.rds')
mben.tums <- readRDS('fdata/mben_tums.rds')

#############################################
## Downsample to minimum cell
#############################################

gcmin <- min(table(mben.tums@meta.data$TumCellType))
bmin <- min(table(btg@meta.data$seurat_clusters))
mmin <- min(table(p14g@meta.data$GnCellType))

Idents(p14g) <- 'CellType'
msmall <- subset(p14g, downsample=mmin)

Idents(btg) <- 'seurat_clusters'
bsmall <- subset(btg, downsample=bmin)

Idents(mben.tums) <- 'TumCellType'
gsmall <- subset(mben.tums, downsample=gcmin)

#############################################
## Calculate differential features
#############################################

mdiff_s <- get_differential(msmall,'GnCellType', 'Mouse')
bdiff_s <- get_differential(bsmall,'seurat_clusters', 'MBEN-Cell')
gcdiff_s <- get_differential(gsmall,'TumCellType', 'MBEN-Nuc')

#############################################
## Find Mouse/Human Corresponding Genes
#############################################

### get mouse to human
library("biomaRt")
# closeAllConnections()
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

x <- rownames(p14g)
genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
             values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

## check with the microarray dataset
cav_df <- read.csv("bulk_data/shh_cav_microarray_mean.csv.gz", row.names=1)
cav_genes <- rownames(cav_df)


mouse.human.list <- list()
for (row in 1:nrow(genesV2)){
	mouse.gene <- genesV2[row, "MGI.symbol"]
	human.gene <- genesV2[row, "HGNC.symbol"]
	if (!is.null(human.gene)){
		mouse.human.list[[mouse.gene]] <- human.gene
	}
}

### Apply the genes rules (either in sheet or capitalized version in Cavalli)
for (n in names(mdiff_s)){
	ll <- as.vector(sapply(rownames(mdiff_s[[n]]), function(x) x %in% names(mouse.human.list)))
	fll <- c()
	for (name in rownames(mdiff_s[[n]])){
		if (name %in% names(mouse.human.list)){
			val <- mouse.human.list[[name]]
		} else if (toupper(name) %in% rownames(cav_df)){
			val <- toupper(name)
		} else {
			val <- FALSE
		}
		fll <- c(fll, val)
	}

	mdiff_s[[n]][['HumanGene']] <- fll
}

#############################################
## Load other gene sets
#############################################

## proteomics
shha <- rownames(read.csv('cellsigs/shhab_comp_shha_200.csv', row.names=1))
shhb <- rownames(read.csv('cellsigs/shhab_comp_shhb_200.csv', row.names=1))

## Cavalli
alpha <- rownames(read.csv('cellsigs/Alpha_cavalli_sigs_200.csv', row.names=1))
beta <- rownames(read.csv('cellsigs/Beta_cavalli_sigs_200.csv', row.names=1))
delta <- rownames(read.csv('cellsigs/Delta_cavalli_sigs_200.csv', row.names=1))
gamma <- rownames(read.csv('cellsigs/Gamma_cavalli_sigs_200.csv', row.names=1))

## Korshunov
tcl1 <- rownames(read.table('cellsigs/mben_1_signature.txt', row.names=1,  sep=' ', header=TRUE))
tcl2 <- rownames(read.table('cellsigs/mben_2_signature.txt', row.names=1,  sep=' ', header=TRUE))

## Riemondy
a1 <- as.data.frame(read_excel("cellsigs/riemondy_cell_sigs.xlsx", sheet = "SHH A1"))$feature
a2 <- as.data.frame(read_excel("cellsigs/riemondy_cell_sigs.xlsx", sheet = "SHH A2"))$feature
b1 <- as.data.frame(read_excel("cellsigs/riemondy_cell_sigs.xlsx", sheet = "SHH B1"))$feature
b2 <- as.data.frame(read_excel("cellsigs/riemondy_cell_sigs.xlsx", sheet = "SHH B2"))$feature
c1 <- as.data.frame(read_excel("cellsigs/riemondy_cell_sigs.xlsx", sheet = "SHH C1"))$feature
c2 <- as.data.frame(read_excel("cellsigs/riemondy_cell_sigs.xlsx", sheet = "SHH C2"))$feature

hgc1 <- as.data.frame(read_excel("cellsigs/Okonechnikov_signatures.xlsx", sheet = "gc_diff_1"))$Feature
hgc2 <- as.data.frame(read_excel("cellsigs/Okonechnikov_signatures.xlsx", sheet = "gc_diff_2"))$Feature
hgcd <- as.data.frame(read_excel("cellsigs/Okonechnikov_signatures.xlsx", sheet = "GC_defined"))$Feature
hgcp <- as.data.frame(read_excel("cellsigs/Okonechnikov_signatures.xlsx", sheet = "GCP"))$Feature

other.gs <- list('Archer_SHHa'= shha, 'Archer_SHHb'= shhb, 'WHO_SHH-3(α)'=alpha, 'WHO_SHH-1(β)'=beta, 'WHO_SHH-4(δ)'=delta, 'WHO_SHH-2(γ)'=gamma, 
	'Korshunov_TCL1'=tcl1, 'Korshunov_TCL2'=tcl2, 'Riemondy_A1'=a1,'Riemondy_A2'=a2,'Riemondy_B1'=b1,'Riemondy_B2'=b2,'Riemondy_C1'=c1,'Riemondy_C2'=c2, 
	'Human_GCP'=hgcp, 'Human_GC-early'=hgc1, 'Human_GC-middle'=hgc2, 'Human_GC-late'=hgcd)


#############################################
## Create gene lists for each size
#############################################

gene.nums <- c(50,100,150,200)
big.gnl <- list()
for (gn in gene.nums){
	small.gnl <- list()
	## process outside gs
	for (name in names(other.gs)){
		small.gnl[[name]] <- head(other.gs[[name]], gn)
	}
	## process bt9 genesets
	for (name in names(bdiff_s)){
		small.gnl[[name]] <- rownames(head(bdiff_s[[name]], gn))
	}
	## process sn-mben genesets
	for (name in names(gcdiff_s)){
		small.gnl[[name]] <- rownames(head(gcdiff_s[[name]], gn))
	}
	## process P14 mouse genesets
	for (name in names(mdiff_s)){
		o <- mdiff_s[[name]]
		clean.mouse <- o[o[,'HumanGene']!=FALSE,]
		small.gnl[[name]] <- head(clean.mouse,gn)$HumanGene
	}	
	big.gnl[[as.character(gn)]] <- small.gnl
}

saveRDS(big.gnl, file='big_genelist.rds')

#############################################
## Run GSVA For Each Size
#############################################

for (name in names(big.gnl)){
	print(name)
	g.gs <- gsva(data.matrix(cav_df), gset.idx.list=big.gnl[[name]], method='gsva', verbose=FALSE)
	write.csv(g.gs, file = gzfile(paste0('gsva_data/FINAL_GSVA_genesigs', name, '.csv.gz')))
}


#############################################
## Run Consensus Clustering For all
#############################################

perl <- c(0.3, 0.5, 0.7, 0.9)
gs <-c(50,100,150,200)

pacl <- list()

for (perc in perl){
	for (gn in gs){
		file <- paste0('gsva_data/FINAL_GSVA_genesigs',gn,'.csv.gz')
		data <- read.csv(file, row.names=1, header=TRUE)
		d <- as.matrix(t(data))
		rownames(d) <- colnames(data)
		results = ConsensusClusterPlus(d,maxK=10,reps=1000,pItem=1,pFeature=perc,clusterAlg="km",distance="euclidean",seed=500,plot="png")
		rname <- paste0('cc_results/ccres_gn-', gn, '_per-', perc, '.csv')
		write.csv(results[[5]]$consensusMatrix, file=rname)
	}
}

