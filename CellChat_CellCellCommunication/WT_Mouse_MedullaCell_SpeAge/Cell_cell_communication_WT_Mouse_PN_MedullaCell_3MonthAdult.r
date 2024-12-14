set.seed(1234)
setwd("./WT_Mouse_MedullaCell_SpeAge")

library("dplyr")
library("ggplot2")
library("Seurat")
library("SeuratDisk") 
library("patchwork")
library("ComplexHeatmap")
library("RaceID")
library("cowplot")
library("RColorBrewer")
library("SeuratWrappers")
library("CellChat")
options(stringsAsFactors = FALSE)

query <- readRDS("../01.2_main_E17andPN_WT_AddMeta_Annotated.RDS")

## Subset Specific Age Group
SelectedAge <- "3MonthAdult"
query_subset <- subset(query, subset = wenyu_names %in% c("WT_AG_3m") & seurat_clusters %in% c(1, 2, 11))
Cluster_Number <- table(query_subset$seurat_clusters)
query_subset <- subset(query_subset, subset = seurat_clusters %in% names(Cluster_Number)[which(Cluster_Number > 10)])
query_subset$NewCelltype <- as.character(query_subset$NewCelltype)

## Get Cluster Colors
a <- DimPlot(query, reduction = "umap", group.by = "NewCelltype")
pbuild <- ggplot2::ggplot_build(a)
pdata <- pbuild$data[[1]]
pdata <-  pdata[order(pdata$group), ]
ucols <- unique(pdata$colour) 
names(ucols) <- levels(query$NewCelltype)

# Create a CellChat object
cellChat <- createCellChat(object = query_subset, group.by = "NewCelltype", assay = "RNA")
ucols_subset <- ucols[levels(cellChat@idents)]

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

cellChat@DB <- CellChatDB

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) 

cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

# Compute the communication probability and infer cellular communication network
ptm = Sys.time()
cellChat <- computeCommunProb(cellChat, type = "triMean") 

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

cellChat <- filterCommunication(cellChat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
#  Infer the cell-cell communication at a signaling pathway level
cellChat <- computeCommunProbPathway(cellChat)

#  Calculate the aggregated cell-cell communication network
ptm = Sys.time()
cellChat <- aggregateNet(cellChat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

saveRDS(cellChat, paste0("WT_Mouse_",SelectedAge,"_MedullaCell_CellChat.RDS"))

# Compute network centrality scores
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")

# Save Data
saveRDS(cellChat, paste0("WT_Mouse_",SelectedAge,"_MedullaCell_CCC.rds"))
rm(query, query_subset)
save.image(paste0("WT_Mouse_",SelectedAge,"_MedullaCell_CCC.RData"))
