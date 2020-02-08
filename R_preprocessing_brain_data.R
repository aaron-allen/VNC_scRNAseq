library("Seurat")
library("tidyverse")
library("loomR")




# waddell mid brain data

waddell_loom <- connect(filename = "waddell_CentralBrain_10k.loom", skip.validate = TRUE, mode = "r+")
cell_id <- waddell_loom$col.attrs$CellID[]
waddell_raw_data <- waddell_loom[["matrix"]][,]
gene <- waddell_loom$row.attrs$Gene[]
waddell_raw_data <- as.data.frame(waddell_raw_data)
rownames(waddell_raw_data) <- cell_id
colnames(waddell_raw_data) <- gene

write_csv(waddell_raw_data, "waddell_loom_raw_data.csv")

raw_data_t <- as.data.frame(t(raw_data))
rownames(raw_data_t) <- colnames(raw_data)
colnames(raw_data_t) <- rownames(raw_data)
write_csv(raw_data_t, "waddell_loom_raw_data_t.csv")

write_csv(as.data.frame(genotype), "waddell_loom_genotype.csv")
write_csv(as.data.frame(cell_id), "waddell_loom_cellid.csv")
write_csv(as.data.frame(gene), "waddell_loom_gene.csv")
waddell_loom$close.all()

waddell_raw_data <- read_csv("waddell_loom_raw_data.csv")
waddell_raw_data <- t(waddell_raw_data)
waddell_raw_data <- as.data.frame(waddell_raw_data)
waddell_cellid <- read_csv("waddell_loom_cellid.csv")
waddell_gene <- read_csv("waddell_loom_gene.csv")
colnames(waddell_raw_data) <- waddell_cellid$cell_id
rownames(waddell_raw_data) <- waddell_gene$gene

waddell_raw_data <- Matrix(as.matrix(waddell_raw_data), sparse= TRUE)
waddell <- CreateSeuratObject(raw.data = waddell_raw_data, project = "waddell")

rm(waddell_raw_data)
rm(waddell_gene)
rm(waddell_cellid)
gc()


waddell@meta.data$Replicate <- waddell@ident
waddell_genotype <- read_csv("waddell_loom_genotype.csv")
waddell@meta.data$Genotype <- waddell_genotype$genotype
rm(waddell_genotype)
gc()

mito.genes <- grep(pattern = "^mt:", x = rownames(x = waddell@data), value = TRUE)
prop.mito <- Matrix::colSums(waddell@raw.data[mito.genes, ]) / Matrix::colSums(waddell@raw.data)
waddell <- AddMetaData(object = waddell, metadata = prop.mito, col.name = "prop.mito")
rm(mito.genes)
rm(prop.mito)
gc()

waddell <- NormalizeData(object = waddell, normalization.method = "LogNormalize", scale.factor = 10000)
waddell <- ScaleData(object = waddell, vars.to.regress = c("nUMI", "prop.mito"), do.par = TRUE, num.cores = 8)
waddell <- FindVariableGenes(waddell, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001)
length(waddell@var.genes)

saveRDS(waddell,"waddell_brain.rds")
gc()

waddell <- RunPCA(object = waddell, pcs.compute = 50, pc.genes = waddell@var.genes, do.print = TRUE)
waddell <- RunTSNE(waddell,
                  num_threads = 7,
                  verbose = T,
                  dims.use = 1:50,
                  do.fast = FALSE,
                  check_duplicates = F,
                  max_iter = 20000,
                  theta = 0.1)
gc()
waddell <- FindClusters(waddell,
                       resolution = 2.5,
                       dims.use = 1:50,
                       save.SNN = TRUE)

saveRDS(waddell,"waddell_brain.rds")
gc()

TSNEPlot(object = waddell,coord.fixed = T, no.axes = T, do.label = T, no.legend = F, pt.size = 1)
waddell_markers <- FindAllMarkers(object = waddell,
                                 logfc.threshold = 0.5,
                                 test.use = "negbinom",
                                 latent.vars = c("Replicate"),
                                 print.bar = TRUE,
                                 only.pos = TRUE)
write_csv(waddell_markers, "waddell_brain_markers.csv")
gc()

waddell_markers_filtered <- waddell_markers %>% 
    filter(p_val_adj<0.05)









# aerts brain data

aerts <- connect(filename = "aerts_Fly_AdultBrain_Filtered_57k.loom", skip.validate = TRUE, mode = "r+")
cell_id <- aerts$col.attrs$CellID[]
raw_data <- aerts[["matrix"]][,]
Replicate <- aerts$col.attrs$Replicate[]
genotype <- aerts$col.attrs$Genotype[]
age <- aerts$col.attrs$Age[]
gene <- aerts$row.attrs$Gene[]
raw_data <- as.data.frame(raw_data)
rownames(raw_data) <- cell_id
colnames(raw_data) <- gene

write_csv(raw_data, "aerts_raw_data.csv")

raw_data_t <- as.data.frame(t(raw_data))
rownames(raw_data_t) <- colnames(raw_data)
colnames(raw_data_t) <- rownames(raw_data)
write_csv(raw_data_t, "aerts_raw_data_t.csv")

write_csv(as.data.frame(Replicate), "aerts_replicate.csv")
write_csv(as.data.frame(genotype), "aerts_genotype.csv")
write_csv(as.data.frame(cell_id), "aerts_cellid.csv")
write_csv(as.data.frame(gene), "aerts_gene.csv")

aerts$close_all()


aerts_raw_data <- read_csv("aerts_raw_data.csv")
aerts_raw_data <- t(aerts_raw_data)
aerts_raw_data <- as.data.frame(aerts_raw_data)
aerts_cellid <- read_csv("~/Desktop/temp/old_seurat_objects/v3/aerts_cellid.csv")
aerts_gene <- read_csv("~/Desktop/temp/old_seurat_objects/v3/aerts_gene.csv")
colnames(aerts_raw_data) <- aerts_cellid$cell_id
rownames(aerts_raw_data) <- aerts_gene$gene

aerts_raw_data <- Matrix(as.matrix(aerts_raw_data), sparse= TRUE)
aerts <- CreateSeuratObject(raw.data = aerts_raw_data, project = "aerts")

rm(aerts_cellid)
rm(aerts_gene)
rm(aerts_raw_data)
gc()

aerts_replicate <- read_csv("aerts_replicate.csv")
aerts@meta.data$Replicate <- aerts_replicate$Replicate
aerts_genotype <- read_csv("aerts_genotype.csv")
aerts@meta.data$Genotype <- aerts_genotype$genotype
rm(aerts_replicate)
rm(aerts_genotype)
gc()

mito.genes <- grep(pattern = "^mt:", x = rownames(x = aerts@data), value = TRUE)
prop.mito <- Matrix::colSums(aerts@raw.data[mito.genes, ]) / Matrix::colSums(aerts@raw.data)
aerts <- AddMetaData(object = aerts, metadata = prop.mito, col.name = "prop.mito")
rm(mito.genes)
rm(prop.mito)
gc()

aerts <- NormalizeData(object = aerts, normalization.method = "LogNormalize", scale.factor = 10000)
aerts <- ScaleData(object = aerts, vars.to.regress = c("nUMI", "prop.mito"), do.par = TRUE, num.cores = 7)
aerts <- FindVariableGenes(aerts, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001)
length(aerts@var.genes)

saveRDS(aerts,"aerts_brain.rds")
gc()

aerts <- RunPCA(object = aerts, pcs.compute = 82, pc.genes = aerts@var.genes, do.print = TRUE)
aerts <- RunTSNE(aerts,
                 num_threads = 7,
                 verbose = T,
                 dims.use = 1:82,
                 do.fast = FALSE,
                 check_duplicates = F,
                 max_iter = 20000,
                 theta = 0.1)
gc()
aerts <- FindClusters(aerts,
                      resolution = 2,
                      dims.use = 1:82,
                      save.SNN = TRUE)

saveRDS(aerts,"aerts_brain.rds")
gc()

aerts_markers <- FindAllMarkers(object = aerts,
                                logfc.threshold = 0.5,
                                test.use = "negbinom",
                                latent.vars = c("Replicate"),
                                print.bar = TRUE,
                                only.pos = TRUE)
write_csv(aerts_markers, "aerts_brain_markers.csv")
gc()

aerts_markers_filtered <- aerts_markers %>% 
    filter(p_val_adj<0.05)


write_csv(aerts_markers, "aerts_brain_markers_filtered.csv")