
# (Maria) script to identify Sox2+/Sox10+ and its variations from mice SS2 data
# Input is an RDS with a seurat object
# QC performed following Oscars pipeline:
# samples processed independently in age batches: 1) E17  2) PN or together, depending on args2
# QC: keep cells with < 10% mitochondrial genes, 200 < number_of_genes < 10000
# exclude genes expressed in fewer than 3 cells
# PN_n_cells aprox 6274
# E17 n_cells aprox 4166
# regress total read counts and mitochondrial content
# cell cycle regressed using Tirosh et. al 2016 gene signature
# Model HAS TO HAVE some vars.to.regress or batch variables
# 10 nearest neighbors
# 40 PC
# clustering method Louvain
# normalized expression: log(counts_per_10000 +1) used to get DEGs
# DEGs significance Welch's BH correction with FDR threshold = 0-01

#-------------------------------------------------------------------------#
rm(list=ls())
set.seed(1234)

library("dplyr")
library("ggplot2")
library("Seurat")
library("SeuratDisk") # library(anndata)
library("patchwork")
library("ComplexHeatmap")
library("RaceID")
library("RColorBrewer")
library("SeuratWrappers")
#-----------------------------------#
cat("READING PASSED VARIABLES\n")
args <-commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
        assign(paste("args",i,sep=""),eval(parse(text=args[i])))
        cat(paste("args",i,sep=""),":\n")
        str(eval(parse(text=paste("args",i,sep=""))))
}
cat("Done",Sys.time(),"-----------------\n")

cat("input data\n")
query=readRDS(args1)
#-----------------------------------------#

# Housekeep

# Do we need to filter out any cells? Provide meta.data field and values to keep
## otherwise "no"
## it only accepts one field. Look at my other script for an expanded version 

if(toupper(args2)!="NO"){
	idx=query@meta.data %>% filter(!!as.name(names(args2))%in%args2[[paste(names(args2))]]) %>% rownames()
	query=query[,idx]	 
}

# Remove spike ins
idx=grep(x=rownames(query),pattern="^ERCC")
if(length(idx)>0) query=query[-idx,]
#-----------------------------------------#
# QC
# nFeatures
idx= query@meta.data %>% filter(nFeature_RNA >10000 | nFeature_RNA < 200) %>% rownames()

# % mito
mito.genes = query %>% grep(pattern = "^MT-", x = rownames(.), value = TRUE, ignore.case=TRUE)
percent.mito = Matrix::colSums(query[["RNA"]]@counts[mito.genes, ])/Matrix::colSums(query[["RNA"]]@counts)
query = AddMetaData(object = query, metadata = percent.mito, col.name = "percent.mito")
idx = c(idx, which((percent.mito*100)>10) %>% names())

if(length(idx)>0) { query = query %>% .[,!colnames(.)%in%idx] }

# features expressed
idx=rownames(query)[which(rowSums(query[["RNA"]]@counts)<3)]

if(length(idx)>0) { query=query[!rownames(query)%in%idx,] }

# Remove ERCC


# Cell cycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015. There is an updated list -- cc.genes.updated.2019 -- but Oscar used the old one
data(cc.genes) 

# Basic function to convert human to mouse gene names
# This solution comes inspired from: https://github.com/satijalab/seurat/issues/2493 and biomaRt::getLDS() 
# it needs an online server
# I set a different mirror and older host due to connection issues on their side - will most likely not impact results
# Another option would be to just get the list convertion online and then read in text file
convertHumanGeneList <- function(x){
	require("biomaRt")
	human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror="useast", host = "dec2021.archive.ensembl.org")
	mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror="useast", host = "dec2021.archive.ensembl.org")
	genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
	humanx <- unique(genesV2[, 2])
	return(humanx)
}

m.s.genes <- convertHumanGeneList(cc.genes$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)

query = CellCycleScoring(query, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)

#-----------------------------------------#
# Fit model
if(toupper(args4) == "NO") {
	query = NormalizeData(query, normalization.method = "LogNormalize", scale.factor = 10000) %>%
		FindVariableFeatures(., selection.method = "vst", nfeatures = 2500) %>%
		ScaleData(., features = rownames(.),vars.to.regress = args5) %>%
		RunPCA(., features = VariableFeatures(object = .)) %>% 
		FindNeighbors(., dims = 1:40,k.param=20) %>%
		FindClusters(., resolution = 0.6, n.iter=100) %>%
		RunUMAP(.,dims=1:40)
}
if(toupper(args4) == "YES") {
	query = NormalizeData(query, normalization.method = "LogNormalize", scale.factor = 10000) %>%
		FindVariableFeatures(., selection.method = "vst", nfeatures = 2500) %>%
		.@assays$RNA@scale.data <- as.matrix(0) %>%
		RunFastMNN(object.list = SplitObject(., split.by = args5))
		RunPCA(., features = VariableFeatures(object = .)) %>% 
		FindNeighbors(., dims = 1:40,k.param=20, reduction="mnn") %>%
		FindClusters(., resolution = 0.6, n.iter=100) %>%
		RunUMAP(.,dims=1:40, reduction="mnn")	
}

# Format for cellxgene
i <- query@meta.data %>% sapply(., is.factor)
query@meta.data[i] <- lapply(query@meta.data[i], as.character)

#-----------------------------------#
#outnm="/crex2/proj/sllstore2017016/Paper2024_mmWT/output-d/"
save.image(paste(args3,".RData",sep=""))
saveRDS(object=query,file=paste(args3,".RDS",sep=""))

SaveH5Seurat(query, filename = paste(args3,".h5Seurat",sep=""))
Convert(paste(args3,".h5Seurat",sep=""), dest = "h5ad")
#sceasy::convertFormat(query, from="seurat", to="anndata", outFile=paste(outnm,"01_main_",args3,".h5ad",sep=""))
file.remove(paste(args3,".h5Seurat",sep=""))


cat("DONE",Sys.time(),"-----------------------------------------------------\n")

#-----------------------------------------#
#------------Annotate Clusters------------#
#---------------Jiacheng Zhu--------------#
cat("Add More Meta Information and Annotate the Clusters\n")
Mouse_Embryo <- readRDS("/crex/proj/sllstore2017016/joeyzhu/Project/Spatial_NB/00_Data/scMouseData/allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_oAllOthrMnk_scnpy.depNmtchndrlRgrs_ccRgrssn.pca.data.louvain_n_leiden.Seurat.rds")
Mouse_Adult <- readRDS("/crex/proj/sllstore2017016/joeyzhu/Project/Spatial_NB/00_Data/scMouseData/ISS_Bin_Deconv_Refs2/allMouse_expnd_Kif1bW2r_nWyLate_eP5nP14nP90_sAllOthrMnk_scnpy.depNmtchndrlRgrs_ccRgrssn.pca.data.louvain_n_leiden.Seurat.rds")
Idents(Mouse_Embryo) <- "louvain"
new.cluster.ids <- c("ADR Chromaffin","NOR Chromaffin","Glial","Neuroblast","ADR Chromaffin WT","Cycling Transitory","Cycling Neuroblast","Bridge-Neuroblast","Cycling Hybrid","Cycling Hybrid","Bridge","Unknown","Unknown")
names(new.cluster.ids) <- levels(Mouse_Embryo)
Mouse_Embryo <- RenameIdents(Mouse_Embryo, new.cluster.ids)
Mouse_Embryo$Celltype <- Idents(Mouse_Embryo)

Idents(Mouse_Adult) <- "leiden"
new.cluster.ids <- c("Macrophage", "Macrophage", "Endothelial", "Macrophage", "Zona Glomerulosa", "Neoplastic-Pheo", "Endothelial", "Zona Glomerulosa", "Neoplastic-Hybrid",
                     "Endothelial", "Dendritic", "NK", "Endothelial", "Glial", "Pre-Neoplastic", "T", "Zona Fasciculata", "NOR Chromaffin",
                     "ADR Chromaffin", "Mesenchymal", "Neoplastic-Neuroblast like", "Glial", "Macrophage", "Dendritic", "Mesenchymal", "CAF", "Capsule",
                     "Macrophage", "Dendritic", "Zona Fasciculata", "Macrophage", "Unknown", "Endothelial", "B", "Mast", "Unknown")
names(new.cluster.ids) <- levels(Mouse_Adult)
Mouse_Adult <- RenameIdents(Mouse_Adult, new.cluster.ids)
Mouse_Adult$Celltype <- Idents(Mouse_Adult)

MetaInfor <- rbind(Mouse_Adult@meta.data[,c("samples", "nCounts_RNA", "nFeaturess_RNA", "stage", "outcome", "percent_mito", "louvain", "leiden", "wenyu_names", "Celltype")],
                   Mouse_Embryo@meta.data[,c("samples", "nCounts_RNA", "nFeaturess_RNA", "stage", "outcome", "percent_mito", "louvain", "leiden", "wenyu_names", "Celltype")])
MetaInfor$Dev_Stage <- c(rep("PN", times = dim(Mouse_Adult)[2]), rep("E", times = dim(Mouse_Embryo)[2]))
MetaInfor$Celltype_Stage <- paste0(MetaInfor$Dev_Stage, "_", MetaInfor$Celltype)

rm(Mouse_Embryo, Mouse_Adult)
# Add label to the new data
MyCells <- colnames(query)[which(colnames(query) %in% rownames(MetaInfor))]

query$OldCelltype <- paste0(query$wenyu_names, "_New")
query$OldCelltype[MyCells] <- as.character(MetaInfor[MyCells,"Celltype"])

query$OldCelltype_Stage <- paste0(query$wenyu_names, "_New")
query$OldCelltype_Stage[MyCells] <- as.character(MetaInfor[MyCells,"Celltype_Stage"])

query$NewCelltype <- ""
query$NewCelltype[which(query$seurat_clusters == 0)] <- "0-E_ADR Chromaffin"
query$NewCelltype[which(query$seurat_clusters == 1)] <- "1-PN_NOR Chromaffin"
query$NewCelltype[which(query$seurat_clusters == 2)] <- "2-PN_ADR Chromaffin"
query$NewCelltype[which(query$seurat_clusters == 3)] <- "3-PN_Endothelial"
query$NewCelltype[which(query$seurat_clusters == 4)] <- "4-PN_Macrophage"
query$NewCelltype[which(query$seurat_clusters == 5)] <- "5-E_NOR Chromaffin"
query$NewCelltype[which(query$seurat_clusters == 6)] <- "6-E_Glial/SCP"
query$NewCelltype[which(query$seurat_clusters == 7)] <- "7-PN_Cortex ZG"
query$NewCelltype[which(query$seurat_clusters == 8)] <- "8-E_Cycling cell"
query$NewCelltype[which(query$seurat_clusters == 9)] <- "9-PN_Cortex X-Zone"
query$NewCelltype[which(query$seurat_clusters == 10)] <- "10-E_Neuroblast"
query$NewCelltype[which(query$seurat_clusters == 11)] <- "11-PN_Glial/Sustentacular"
query$NewCelltype[which(query$seurat_clusters == 12)] <- "12-PN_Cortex ZF"
query$NewCelltype[which(query$seurat_clusters == 13)] <- "13-PN_Macrophage"
query$NewCelltype[which(query$seurat_clusters == 14)] <- "14-PN_Cortex ZG"
query$NewCelltype[which(query$seurat_clusters == 15)] <- "15-PN_Dendritic cell"
query$NewCelltype[which(query$seurat_clusters == 16)] <- "16-PN_T cell"
query$NewCelltype[which(query$seurat_clusters == 17)] <- "17-PN_NK cell"
query$NewCelltype[which(query$seurat_clusters == 18)] <- "18-PN_Mesenchymal"
query$NewCelltype[which(query$seurat_clusters == 19)] <- "19-PN_Dendritic cell"

NewCelltypes <- sort(unique(query$NewCelltype))
NewCelltypesOrdered <- NewCelltypes[c(1,2,13:20,3:12)]

query$NewCelltype <- factor(query$NewCelltype, levels = NewCelltypesOrdered)

#Rnm <- "01.2_main_E17andPN_WT_AddMeta_Annotated"
save.image(paste(args3,"_AddMeta_Annotated.RData",sep=""))
saveRDS(object=query,file=paste(args3,"_AddMeta_Annotated.RDS",sep=""))

query_copy <- query
query_copy$NewCelltype <- as.character(query_copy$NewCelltype)
query_copy@assays$RNA@scale.data <- as_matrix(query@assays$RNA@data)

SaveH5Seurat(query, filename = paste(args3,"_AddMeta_Annotated.h5Seurat",sep=""))
Convert(paste(args3,"_AddMeta_Annotated.h5Seurat",sep=""), dest = "h5ad")

cat("Done",Sys.time(),"-----------------\n")

sessionInfo()

