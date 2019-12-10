#!/usr/bin/env Rscript

# Aaron M. Allen, 2019.06.06

library(Seurat)
library(dplyr)
library(Matrix)



NumberOfCCAs <- 45
nGeneLowCutOff <- 200
nGeneHighCutOff <- Inf
nUMILowCutOff <- 1200
nUMIHighCutOff <- 10000
ClusterResolution <- 12  
MitoLowCutOff <- -Inf
MitoHighCutOff <- 0.15
Perplex <- 30
Theta <- 0.05

AnalysisName <- paste0(
    "noDF_mt(0-",MitoHighCutOff,
    ")_G(",nGeneLowCutOff,"-",nGeneHighCutOff,
    ")_U(",nUMILowCutOff,"-",nUMIHighCutOff,
    ")_CCA(",NumberOfCCAs,
    ")_res(",ClusterResolution,
    ")_Perplex(",Perplex,
    ")_theta(",Theta#,
    ")_TopVarGenes"
    )

DateAndTime <- format(Sys.time(), "%Y%m%dT%H%M")
AnalysisDir <- paste0(DateAndTime,"_",AnalysisName)
dir.create(file.path(AnalysisDir))

setwd(file.path(AnalysisDir))
WorkingDirectory <- getwd()


dir.create(file.path("R_Objects"))

# Load Data
###########
MaleVNCRep1.data <- read.table("../InputFiles/12_MaleVNCRep1_merged_gene_exon_tagged_top30000cells.dge.txt.gz", header = TRUE, sep = "\t", row.names = 1)
FemaleVNCRep1.data <- read.table("../InputFiles/12_FemaleVNCRep1_merged_gene_exon_tagged_top30000cells.dge.txt.gz", header = TRUE, sep = "\t", row.names = 1)
MaleVNCRep2.data <- read.table("../InputFiles/12_MaleVNCRep2_merged_gene_exon_tagged_top30000cells.dge.txt.gz", header = TRUE, sep = "\t", row.names = 1)
FemaleVNCRep2.data <- read.table("../InputFiles/12_FemaleVNCRep2_merged_gene_exon_tagged_top30000cells.dge.txt.gz", header = TRUE, sep = "\t", row.names = 1)

MaleVNCRep1.data <- Matrix(as.matrix(MaleVNCRep1.data), sparse= TRUE)
FemaleVNCRep1.data <- Matrix(as.matrix(FemaleVNCRep1.data), sparse= TRUE)
MaleVNCRep2.data <- Matrix(as.matrix(MaleVNCRep2.data), sparse= TRUE)
FemaleVNCRep2.data <- Matrix(as.matrix(FemaleVNCRep2.data), sparse= TRUE)

MaleVNCRep1 <- CreateSeuratObject(raw.data = MaleVNCRep1.data, min.cells = 3, min.genes = nGeneLowCutOff, project = "VNC")
FemaleVNCRep1 <- CreateSeuratObject(raw.data = FemaleVNCRep1.data, min.cells = 3, min.genes = nGeneLowCutOff, project = "VNC")
MaleVNCRep2 <- CreateSeuratObject(raw.data = MaleVNCRep2.data, min.cells = 3, min.genes = nGeneLowCutOff, project = "VNC")
FemaleVNCRep2 <- CreateSeuratObject(raw.data = FemaleVNCRep2.data, min.cells = 3, min.genes = nGeneLowCutOff, project = "VNC")

MaleVNCRep1@meta.data$Replicate <- "MaleRep1"
MaleVNCRep2@meta.data$Replicate <- "MaleRep2"
FemaleVNCRep1@meta.data$Replicate <- "FemaleRep1"
FemaleVNCRep2@meta.data$Replicate <- "FemaleRep2"



length(WhichCells(MaleVNCRep1, subset.name = "nUMI", accept.low = nUMILowCutOff))
length(WhichCells(MaleVNCRep2, subset.name = "nUMI", accept.low = nUMILowCutOff))
length(WhichCells(FemaleVNCRep1, subset.name = "nUMI", accept.low = nUMILowCutOff))
length(WhichCells(FemaleVNCRep2, subset.name = "nUMI", accept.low = nUMILowCutOff))




#### control for mitochondrial genes
####################################
mito.genes <- grep(pattern = "^mt:", x = rownames(x = MaleVNCRep1@data), value = TRUE)
prop.mito <- Matrix::colSums(MaleVNCRep1@raw.data[mito.genes, ]) / Matrix::colSums(MaleVNCRep1@raw.data)
MaleVNCRep1 <- AddMetaData(object = MaleVNCRep1, metadata = prop.mito, col.name = "prop.mito")

mito.genes <- grep(pattern = "^mt:", x = rownames(x = MaleVNCRep2@data), value = TRUE)
prop.mito <- Matrix::colSums(MaleVNCRep2@raw.data[mito.genes, ]) / Matrix::colSums(MaleVNCRep2@raw.data)
MaleVNCRep2 <- AddMetaData(object = MaleVNCRep2, metadata = prop.mito, col.name = "prop.mito")

mito.genes <- grep(pattern = "^mt:", x = rownames(x = FemaleVNCRep1@data), value = TRUE)
prop.mito <- Matrix::colSums(FemaleVNCRep1@raw.data[mito.genes, ]) / Matrix::colSums(FemaleVNCRep1@raw.data)
FemaleVNCRep1 <- AddMetaData(object = FemaleVNCRep1, metadata = prop.mito, col.name = "prop.mito")

mito.genes <- grep(pattern = "^mt:", x = rownames(x = FemaleVNCRep2@data), value = TRUE)
prop.mito <- Matrix::colSums(FemaleVNCRep2@raw.data[mito.genes, ]) / Matrix::colSums(FemaleVNCRep2@raw.data)
FemaleVNCRep2 <- AddMetaData(object = FemaleVNCRep2, metadata = prop.mito, col.name = "prop.mito")



####################################
MaleVNC.combined <- MergeSeurat(object1 = MaleVNCRep1, 
                                object2 = MaleVNCRep2, 
                                do.normalize = FALSE,
                                add.cell.id1 = "Rep1", 
                                add.cell.id2 = "Rep2", 
                                project = "VNC")
FemaleVNC.combined <- MergeSeurat(object1 = FemaleVNCRep1, 
                                  object2 = FemaleVNCRep2,
                                  do.normalize = FALSE, 
                                  add.cell.id1 = "Rep1", 
                                  add.cell.id2 = "Rep2", 
                                  project = "VNC")
MaleVNC.combined@meta.data$Sex <- "Male"
FemaleVNC.combined@meta.data$Sex <- "Female"

rm(MaleVNCRep1.data)
rm(FemaleVNCRep1.data)
rm(MaleVNCRep2.data)
rm(FemaleVNCRep2.data)
rm(MaleVNCRep1)
rm(FemaleVNCRep1)
rm(MaleVNCRep2)
rm(FemaleVNCRep2)
gc()




#### filter, normalize, and scale data
######################################
MaleVNC.combined <- FilterCells(object = MaleVNC.combined, 
                                subset.names = c("nGene","nUMI","prop.mito"), 
                                low.thresholds = c(nGeneLowCutOff,nUMILowCutOff,MitoLowCutOff), 
                                high.thresholds = c(nGeneHighCutOff,nUMIHighCutOff,MitoHighCutOff)
)
MaleVNC.combined <- NormalizeData(MaleVNC.combined)
MaleVNC.combined <- ScaleData(object = MaleVNC.combined, vars.to.regress = c("Replicate", "nUMI","prop.mito"), do.par = TRUE, num.cores = 12)


FemaleVNC.combined <- FilterCells(object = FemaleVNC.combined, 
                                  subset.names = c("nGene","nUMI","prop.mito"), 
                                  low.thresholds = c(nGeneLowCutOff,nUMILowCutOff,MitoLowCutOff), 
                                  high.thresholds = c(nGeneHighCutOff,nUMIHighCutOff,MitoHighCutOff)
)
FemaleVNC.combined <- NormalizeData(FemaleVNC.combined)
FemaleVNC.combined <- ScaleData(object = FemaleVNC.combined, vars.to.regress = c("Replicate", "nUMI","prop.mito"), do.par = TRUE, num.cores = 12)
######################################



###### use variable genes
MaleVNC.combined <- FindVariableGenes(MaleVNC.combined, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001)
FemaleVNC.combined <- FindVariableGenes(FemaleVNC.combined, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001)

saveRDS(MaleVNC.combined, file = "R_Objects/MaleVNC.combined.rds")
saveRDS(FemaleVNC.combined, file = "R_Objects/FemaleVNC.combined.rds")

g.1 <- (MaleVNC.combined@var.genes)
g.2 <- (FemaleVNC.combined@var.genes)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(MaleVNC.combined@scale.data))
genes.use <- intersect(genes.use, rownames(FemaleVNC.combined@scale.data))
length(genes.use)


NumberOfCCAs <- 100
VNC.combined <- RunCCA(MaleVNC.combined, FemaleVNC.combined, genes.use = genes.use, num.cc = NumberOfCCAs, add.cell.id1 = "Male", add.cell.id2 = "Female")

VNC.combined <- AlignSubspace(VNC.combined, reduction.type = "cca", grouping.var = "Sex", dims.align = 1:NumberOfCCAs, do.par = TRUE, num.cores = 12)
saveRDS(VNC.combined, file = "R_Objects/VNC.combined.rds")

NumberOfCCAs <- 45
VNC.combined <- RunTSNE(VNC.combined,
                            num_threads = 12,
                            verbose = T,
                            reduction.use = "cca.aligned",
                            dims.use = 1:NumberOfCCAs,
                            do.fast = FALSE,
                            check_duplicates = F,
                            perplexity = Perplex,
                            max_iter = 20000,
                            theta = Theta)
VNC.combined <- FindClusters(VNC.combined,
                                 reduction.type = "cca.aligned",
                                 resolution = ClusterResolution,
                                 dims.use = 1:NumberOfCCAs,
                                 save.SNN = TRUE)
saveRDS(VNC.combined, file = "R_Objects/VNC.combined.rds")








# Cluster Markers
#################
setwd(WorkingDirectory)

dir.create(file.path("ClusterMarkers"))
setwd(file.path("ClusterMarkers"))

AllClusterMarkers <- FindAllMarkers(object = VNC.combined,
                                    test.use = "negbinom",
                                    latent.vars = c("Replicate"),
                                    print.bar = TRUE,
                                    only.pos = TRUE)

AllClusterMarkers <- AllClusterMarkers %>% 
    filter(avg_logFC>0.5) %>% 
    filter(p_val_adj<0.05)

AllClusterMarkersFileName <- paste0("AllCluster_Markers_",AnalysisName,".csv")
write.csv(AllClusterMarkers, AllClusterMarkersFileName)

##################




