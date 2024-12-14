library(dplyr)
library(CellChat)
library(patchwork)
library(network)
library(igraph)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(extrafont)
library(grid)
library(gridExtra)
library(gridGraphics)
library(ComplexHeatmap)

set.seed(1234)

## Load CellChat Object of Each Dataset
cellChat_AgedAdult <- readRDS("./WT_Mouse_AgedAdult_MedullaCell_CCC.rds")
## Add Meta Infors
cellChat_AgedAdult@meta$ReNewCelltype <- ""
cellChat_AgedAdult@meta$ReNewCelltype[which(cellChat_AgedAdult@meta$NewCelltype == "1-PN_NOR Chromaffin")] <- "NOR"
cellChat_AgedAdult@meta$ReNewCelltype[which(cellChat_AgedAdult@meta$NewCelltype == "2-PN_ADR Chromaffin")] <- "ADR"
cellChat_AgedAdult@meta$ReNewCelltype[which(cellChat_AgedAdult@meta$NewCelltype == "11-PN_Glia/Sustentacular")] <- "Glia"

## Load CellChat Object of Each Dataset
cellChat_3MonthAdult <- readRDS("./WT_Mouse_3MonthAdult_MedullaCell_CCC.rds")
## Add Meta Infors
cellChat_3MonthAdult@meta$ReNewCelltype <- ""
cellChat_3MonthAdult@meta$ReNewCelltype[which(cellChat_3MonthAdult@meta$NewCelltype == "1-PN_NOR Chromaffin")] <- "NOR"
cellChat_3MonthAdult@meta$ReNewCelltype[which(cellChat_3MonthAdult@meta$NewCelltype == "2-PN_ADR Chromaffin")] <- "ADR"
cellChat_3MonthAdult@meta$ReNewCelltype[which(cellChat_3MonthAdult@meta$NewCelltype == "11-PN_Glia/Sustentacular")] <- "Glia"

## Load CellChat Object of Each Dataset
cellChat_P14 <- readRDS("./WT_Mouse_PN_14D_MedullaCell_CCC.rds")
## Add Meta Infors
cellChat_P14@meta$ReNewCelltype <- ""
cellChat_P14@meta$ReNewCelltype[which(cellChat_P14@meta$NewCelltype == "1-PN_NOR Chromaffin")] <- "NOR"
cellChat_P14@meta$ReNewCelltype[which(cellChat_P14@meta$NewCelltype == "2-PN_ADR Chromaffin")] <- "ADR"
cellChat_P14@meta$ReNewCelltype[which(cellChat_P14@meta$NewCelltype == "11-PN_Glia/Sustentacular")] <- "Glia"

## Load CellChat Object of Each Dataset
cellChat_P5 <- readRDS("./WT_Mouse_PN_5D_MedullaCell_CCC.rds")
## Add Meta Infors
cellChat_P5@meta$ReNewCelltype <- ""
cellChat_P5@meta$ReNewCelltype[which(cellChat_P5@meta$NewCelltype == "1-PN_NOR Chromaffin")] <- "NOR"
cellChat_P5@meta$ReNewCelltype[which(cellChat_P5@meta$NewCelltype == "2-PN_ADR Chromaffin")] <- "ADR"
cellChat_P5@meta$ReNewCelltype[which(cellChat_P5@meta$NewCelltype == "11-PN_Glia/Sustentacular")] <- "Glia"

## Load CellChat Object of Each Dataset
cellChat_E17 <- readRDS("./WT_Mouse_E_17_MedullaCell_CCC.rds")
## Add Meta Infors
cellChat_E17@meta$ReNewCelltype <- ""
cellChat_E17@meta$ReNewCelltype[which(cellChat_E17@meta$NewCelltype == "0-E_ADR Chromaffin")] <- "ADR"
cellChat_E17@meta$ReNewCelltype[which(cellChat_E17@meta$NewCelltype == "5-E_NOR Chromaffin")] <- "NOR"
cellChat_E17@meta$ReNewCelltype[which(cellChat_E17@meta$NewCelltype == "6-E_Glia/SCP")] <- "Glia"

## Re-label the cells
cellChat_AgedAdult <- setIdent(cellChat_AgedAdult, ident.use = "ReNewCelltype")
cellChat_3MonthAdult <- setIdent(cellChat_3MonthAdult, ident.use = "ReNewCelltype")
cellChat_P14 <- setIdent(cellChat_P14, ident.use = "ReNewCelltype")
cellChat_P5 <- setIdent(cellChat_P5, ident.use = "ReNewCelltype")
cellChat_E17 <- setIdent(cellChat_E17, ident.use = "ReNewCelltype")

## Re-do the analysis on the new labels
cellChat_AgedAdult <- computeCommunProb(cellChat_AgedAdult, type = "triMean") 
cellChat_AgedAdult <- filterCommunication(cellChat_AgedAdult, min.cells = 10)
cellChat_AgedAdult <- computeCommunProbPathway(cellChat_AgedAdult)
cellChat_AgedAdult <- aggregateNet(cellChat_AgedAdult)
cellChat_AgedAdult <- netAnalysis_computeCentrality(cellChat_AgedAdult, slot.name = "netP")
saveRDS(cellChat_AgedAdult, file = "WT_Mouse_AgedAdult_MedullaCell_CCC_ReNamed.rds")

cellChat_3MonthAdult <- computeCommunProb(cellChat_3MonthAdult, type = "triMean") 
cellChat_3MonthAdult <- filterCommunication(cellChat_3MonthAdult, min.cells = 10)
cellChat_3MonthAdult <- computeCommunProbPathway(cellChat_3MonthAdult)
cellChat_3MonthAdult <- aggregateNet(cellChat_3MonthAdult)
cellChat_3MonthAdult <- netAnalysis_computeCentrality(cellChat_3MonthAdult, slot.name = "netP")
saveRDS(cellChat_3MonthAdult, file = "WT_Mouse_3MonthAdult_MedullaCell_CCC_ReNamed.rds")

cellChat_P14 <- computeCommunProb(cellChat_P14, type = "triMean") 
cellChat_P14 <- filterCommunication(cellChat_P14, min.cells = 10)
cellChat_P14 <- computeCommunProbPathway(cellChat_P14)
cellChat_P14 <- aggregateNet(cellChat_P14)
cellChat_P14 <- netAnalysis_computeCentrality(cellChat_P14, slot.name = "netP")
saveRDS(cellChat_P14, file = "WT_Mouse_PN_14D_MedullaCell_CCC_ReNamed.rds")

cellChat_P5 <- computeCommunProb(cellChat_P5, type = "triMean") 
cellChat_P5 <- filterCommunication(cellChat_P5, min.cells = 10)
cellChat_P5 <- computeCommunProbPathway(cellChat_P5)
cellChat_P5 <- aggregateNet(cellChat_P5)
cellChat_P5 <- netAnalysis_computeCentrality(cellChat_P5, slot.name = "netP")
saveRDS(cellChat_P5, file = "WT_Mouse_PN_5D_MedullaCell_CCC_ReNamed.rds")

cellChat_E17 <- computeCommunProb(cellChat_E17, type = "triMean") 
cellChat_E17 <- filterCommunication(cellChat_E17, min.cells = 10)
cellChat_E17 <- computeCommunProbPathway(cellChat_E17)
cellChat_E17 <- aggregateNet(cellChat_E17)
cellChat_E17 <- netAnalysis_computeCentrality(cellChat_E17, slot.name = "netP")
saveRDS(cellChat_E17, file = "WT_Mouse_E_17_MedullaCell_CCC_ReNamed.rds")

object.list <- list(E17 = cellChat_E17, P5 = cellChat_P5, P14 = cellChat_P14, P90 = cellChat_3MonthAdult, Aged = cellChat_AgedAdult) 
cellChat_Merged <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

saveRDS(cellChat_Merged, "WT_Mouse_Merged_AllAgeGroup_CCC.Rds")
