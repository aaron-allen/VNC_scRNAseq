library("Seurat")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("clustree")
library("colorspace")
library("ggcorrplot")
library("ggpubr")
library("circlize")
library("eulerr")
library("corrr")
library("ggpubr")



VNC.combined <- readRDS("R_Objects/VNC.combined.rds")
cluster_markers <- read_csv(file = "AllCluster_Markers_noDF_mt(0-0.15)_G(200-Inf)_U(1200-10000)_CCA(45)_res(12)_Perplex(30)_theta(0.05).csv",col_names = T)





##################################################
##    figure 1    ################################
##################################################


# Figure 1B
TSNEPlot(VNC.combined,
               do.return = T, 
               pt.size = 0.01,
               no.axes = T,
               do.label = F, 
               coord.fixed=T)


# Figure 1C


top5 <- cluster_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
ClusMarExpByCluster <- AverageExpression(object = VNC.combined,genes.use = unique(top5$gene), use.scale = T)
pheatmap(ClusMarExpByCluster,
         breaks=seq(-3,3,0.05),
         color = diverging_hcl(120, "Vik"),
         border_color = NA,
         cellwidth = 3, 
         cellheight = 3,
         labels_row="Clusters",
         labels_col="Genes",
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 16,
         angle_col = 90)





# Figure 1D

length(WhichCells(object = VNC.combined))

exp_genes <- data.frame(mean.exp = rowMeans(expm1(x = VNC.combined@data)))
exp_genes <- exp_genes[order(exp_genes$mean.exp, decreasing = TRUE),, drop = FALSE] 
exp_genes <- rownames_to_column(exp_genes,var="gene")
exp_genes <- exp_genes %>% filter(mean.exp>0)
length(exp_genes$gene)


non.neuron <- c(18,23,24,55,70,98,99,105,106,116)
non.neuron.markers <- cluster_markers %>% filter(cluster %in% non.neuron)
neuron.markers <- cluster_markers %>% filter(!(cluster %in% non.neuron))
length(unique(non.neuron.markers$gene))
length(unique(neuron.markers$gene))




# Figure 1-Figure supplement 1A
vnc_raw <- readRDS("R_objects/vnc_raw.rds")
VlnPlot(object = vnc_raw,features.plot = c("nGene","nUMI","prop.mito"),point.size.use = 0,group.by = "Replicate",nCol = 3,x.lab.rot = T)

# Figure 1-Figure supplement 1B
cellxnUMI <- data.frame(nUMI=vnc_raw@meta.data[["nUMI"]],
                        nGene=vnc_raw@meta.data[["nGene"]])
cellxnUMI %>%
    filter(nUMI>0) %>% 
    ggplot(aes(x=nUMI))+
    geom_histogram(color="darkblue", fill="lightblue", bins = 50) + 
    xlim(0,5000) +
    ylim(0,4000) +
    theme_gray()

# Figure 1-Figure supplement 1C
VlnPlot(object = VNC.combined, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)


# Figure 1-Figure supplement 1D
genefpkm <- read_csv("flyatlas2/genefpkm.csv",col_names = FALSE)
colnames(genefpkm) <- c("submitted_id","TissueID","FPKM","Replicate1","Replicate2","Replicate3","SD","Status")
FBgn <- unique(genefpkm$submitted_id)
write.csv(FBgn,"FBgnList.csv")

FBgn_to_symbol <- read_delim("flyatlas2/FBgn_to_symbol.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
FBgn_to_symbol <- FBgn_to_symbol %>% select(-c("current_id","converted_id"))

flyatlas.vns <- genefpkm %>% filter(TissueID == 7 | TissueID ==8)
fa.vns.sym <- left_join(x = flyatlas.vns[,1:3],y = FBgn_to_symbol,by="submitted_id")
fa.vns.sym2 <- fa.vns.sym %>% group_by(current_symbol) %>% summarize(FPKM=mean(FPKM))

scGene_to_newGene <- read_delim("geneLists/scGene_to_newGene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)           # read the file in
scGene_to_newGene <- scGene_to_newGene[- grep("D*{3}[\\]", scGene_to_newGene$current_symbol),]                  # flybase spits out some genes from other species
scGene_to_newGene <- scGene_to_newGene[! scGene_to_newGene$current_id=="unknown ID", ]                          # remove the rows that flybase couldn't find current symbols for
scGene_to_newGene <- drop_na(scGene_to_newGene)                                                                 # ...
scGene_to_newGene <- scGene_to_newGene[! (scGene_to_newGene$current_symbol!=scGene_to_newGene$submitted_id 
                                          & duplicated(scGene_to_newGene$current_symbol)),]                     # remove rows where the input id had two ouput ids and the input and output don't match

norm.sc.data <- data.frame(submitted_id=VNC.combined@data@Dimnames[[1]], norm_exp=rowSums(VNC.combined@data))
norm.sc.data <- left_join(x = norm.sc.data,y = scGene_to_newGene,by = "submitted_id")
norm.sc.data <- drop_na(norm.sc.data)

sc.fa2.join <- full_join(x = fa.vns.sym2,y = norm.sc.data,by = "current_symbol")
sc.fa2.join <- drop_na(sc.fa2.join)
sc.fa2.join <- sc.fa2.join %>% filter(FPKM>0) %>% filter(norm_exp>0)
sc.fa2.join$FPKM <- log(sc.fa2.join$FPKM)
sc.fa2.join$norm_exp <- log(sc.fa2.join$norm_exp)

ggplot(sc.fa2.join, aes(x = FPKM,y=norm_exp)) + 
    geom_point() +
    geom_smooth(method=lm) +
    stat_cor(label.y = 14, size =5) +
    stat_regline_equation(label.y = 13, size =5) +
    coord_fixed() +
    theme_gray() +
    theme(text = element_text(size=20))


# Figure 1-Figure supplement 2
rep1 <- SubsetData(object = VNC.combined,subset.name = "Replicate",accept.value = "FemaleRep1",subset.raw = T)
rep1_avgexp_byclust <- AverageExpression(object = rep1,use.raw = F,show.progress = T)
rep1_avgexp_byclust_gathered <- rep1_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep1")
rep2 <- SubsetData(object = VNC.combined,subset.name = "Replicate",accept.value = "FemaleRep2",subset.raw = T)
rep2_avgexp_byclust <- AverageExpression(object = rep2,use.raw = F,show.progress = T)
rep2_avgexp_byclust_gathered <- rep2_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep2")
rep3 <- SubsetData(object = VNC.combined,subset.name = "Replicate",accept.value = "MaleRep1",subset.raw = T)
rep3_avgexp_byclust <- AverageExpression(object = rep3,use.raw = F,show.progress = T)
rep3_avgexp_byclust_gathered <- rep3_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep3")
rep4 <- SubsetData(object = VNC.combined,subset.name = "Replicate",accept.value = "MaleRep2",subset.raw = T)
rep4_avgexp_byclust <- AverageExpression(object = rep4,use.raw = F,show.progress = T)
rep4_avgexp_byclust_gathered <- rep4_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep4")

rep1_bulk <- tibble(genes = rep1@data@Dimnames[[1]], rep1_exp = rowSums(rep1@data))
rep2_bulk <- tibble(genes = rep2@data@Dimnames[[1]], rep2_exp = rowSums(rep2@data))
rep3_bulk <- tibble(genes = rep3@data@Dimnames[[1]], rep3_exp = rowSums(rep3@data))
rep4_bulk <- tibble(genes = rep4@data@Dimnames[[1]], rep4_exp = rowSums(rep4@data))

bulk_by_rep <- full_join(full_join(full_join(rep1_bulk,rep2_bulk),rep3_bulk),rep4_bulk)

bulk_by_rep %>% 
    column_to_rownames("genes") %>% 
    log1p() %>% 
    # arrange(rep4_exp) %>% 
    t() %>% 
    pheatmap(breaks=seq(0,10,0.1),
             labels_col = "Genes",
             color = sequential_hcl(100, "inferno"),
             border_color = NA,
             cluster_rows = F,
             cluster_cols = F,
             fontsize = 16,
             angle_col = 90)


rep1_v_rep2 <- bulk_by_rep %>% 
    ggplot(aes(x = log1p(rep1_exp),y=log1p(rep2_exp))) + 
    geom_point() +
    coord_fixed() +
    geom_smooth(method=lm) +
    stat_cor(label.y = 14, size =5) +
    stat_regline_equation(label.y = 13, size =5) +
    ylim(0,15) +
    xlim(0,15) +
    theme_gray() +
    theme(text = element_text(size=20))
rep1_v_rep3 <- bulk_by_rep %>% 
    ggplot(aes(x = log1p(rep1_exp),y=log1p(rep3_exp))) + 
    geom_point() +
    coord_fixed() +
    geom_smooth(method=lm) +
    stat_cor(label.y = 14, size =5) +
    stat_regline_equation(label.y = 13, size =5) +
    ylim(0,15) +
    xlim(0,15) +
    theme_gray() +
    theme(text = element_text(size=20))
rep1_v_rep4 <- bulk_by_rep %>% 
    ggplot(aes(x = log1p(rep1_exp),y=log1p(rep4_exp))) + 
    geom_point() +
    coord_fixed() +
    geom_smooth(method=lm) +
    stat_cor(label.y = 14, size =5) +
    stat_regline_equation(label.y = 13, size =5) +
    ylim(0,15) +
    xlim(0,15) +
    theme_gray() +
    theme(text = element_text(size=20))
rep2_v_rep3 <- bulk_by_rep %>% 
    ggplot(aes(x = log1p(rep2_exp),y=log1p(rep3_exp))) + 
    geom_point() +
    coord_fixed() +
    geom_smooth(method=lm) +
    stat_cor(label.y = 14, size =5) +
    stat_regline_equation(label.y = 13, size =5) +
    ylim(0,15) +
    xlim(0,15) +
    theme_gray() +
    theme(text = element_text(size=20))
rep2_v_rep4 <- bulk_by_rep %>% 
    ggplot(aes(x = log1p(rep2_exp),y=log1p(rep4_exp))) + 
    geom_point() +
    coord_fixed() +
    geom_smooth(method=lm) +
    stat_cor(label.y = 14, size =5) +
    stat_regline_equation(label.y = 13, size =5) +
    ylim(0,15) +
    xlim(0,15) +
    theme_gray() +
    theme(text = element_text(size=20))
rep3_v_rep4 <- bulk_by_rep %>% 
    ggplot(aes(x = log1p(rep3_exp),y=log1p(rep4_exp))) + 
    geom_point() +
    coord_fixed() +
    geom_smooth(method=lm) +
    stat_cor(label.y = 14, size = 5) +
    stat_regline_equation(label.y = 13, size =5) +
    ylim(0,15) +
    xlim(0,15) +
    theme_gray() +
    theme(text = element_text(size=20))

ggarrange(plotlist = list(rep1_v_rep2,
                          rep1_v_rep3,
                          rep1_v_rep4,
                          rep2_v_rep3,
                          rep2_v_rep4,
                          rep3_v_rep4),
          labels = c("rep1 v rep2",
                     "rep1 v rep3",
                     "rep1 v rep4",
                     "rep2 v rep3",
                     "rep2 v rep4",
                     "rep3 v rep4"),
          hjust = -0.7,
          ncol = 3,
          nrow = 2)





# Figure 1-Figure supplement 3A
VNC.combined.res1 <- FindClusters(VNC.combined,reduction.type = "cca.aligned",resolution = 1,dims.use = 1:45,save.SNN = TRUE)
TSNEPlot(VNC.combined.res1,no.legend = T,no.axes = T,pt.size = 0.1,do.label = TRUE,coord.fixed=T)

# Figure 1-Figure supplement 3B
VNC.combined.res3 <- FindClusters(VNC.combined,reduction.type = "cca.aligned",resolution = 6,dims.use = 1:45,save.SNN = TRUE)
TSNEPlot(VNC.combined.res3,no.legend = T,no.axes = T,do.return = T,pt.size = 0.1,do.label = TRUE,coord.fixed=T)

# Figure 1-Figure supplement 3C
VNC.combined.res12 <- FindClusters(VNC.combined,reduction.type = "cca.aligned",resolution = 12,dims.use = 1:45,save.SNN = TRUE)
TSNEPlot(VNC.combined.res12,no.legend = T,no.axes = T,do.return = T,pt.size = 0.1,do.label = TRUE,coord.fixed=T)

# Figure 1-Figure supplement 4D
clustTreePlot <- clustree(VNC.combined, prefix = "res.") + theme(legend.position = "bottom")

# Figure 1-Figure supplement 4
TSNEPlot(object = VNC.combined,
         no.legend = F,
         no.axes = T, 
         do.label = T,
         pt.size = 0.2,
         coord.fixed=T
)


# Figure 1-Figure supplement 5A
RepAlignScore <- CalcAlignmentMetric(object = VNC.combined, reduction.use = "cca.aligned", dims.use = 50, grouping.var = "Replicate")
RepAlignScore
VNC.combined@meta.data$res.12 <- factor(VNC.combined@meta.data$res.12, levels=c(seq(0,119)))
ggplot(VNC.combined@meta.data, aes(x=res.12, fill=Replicate)) + geom_bar(position = "fill")

# Figure 1-Figure supplement 5B
p1 <- TSNEPlot(object = VNC.combined,
               plot.order = c("FemaleRep1"),
               do.return = T,
               group.by = "Replicate",
               colors.use = c("lightgrey","lightgrey","lightgrey","#F8766D"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)
p2 <- TSNEPlot(object = VNC.combined,
               plot.order = c("FemaleRep2"),
               do.return = T,
               group.by = "Replicate",
               colors.use = c("lightgrey","lightgrey","lightgrey","#7CAE00"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)
p3 <- TSNEPlot(object = VNC.combined,
               plot.order = c("MaleRep1"),
               do.return = T,
               group.by = "Replicate",
               colors.use = c("lightgrey","lightgrey","lightgrey","#00BFC4"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)
p4 <- TSNEPlot(object = VNC.combined,
               plot.order = c("MaleRep2"),
               do.return = T,
               group.by = "Replicate",
               colors.use = c("lightgrey","lightgrey","lightgrey","#C77CFF"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)

ggarrange(plotlist = list(p1,p2,p3,p4),
          labels = c("Rep1","Rep2","Rep3","Rep4"),
          hjust = -0.7,
          ncol = 4,
          nrow = 1)


# Figure 1-Figure supplement 5C
all_avgexp_byclust <- bind_rows(rep1_avgexp_byclust_gathered,
                                rep2_avgexp_byclust_gathered,
                                rep3_avgexp_byclust_gathered,
                                rep4_avgexp_byclust_gathered)
mydata <- tibble(comparison = character(),
                 r = numeric(),
                 cluster = character())
for (clust_num in unique(all_avgexp_byclust$cluster)) {
    sub_data <- all_avgexp_byclust %>%
        filter(cluster == clust_num) %>%
        spread(replicate, expression) %>%
        column_to_rownames("gene") %>%
        select(-cluster) %>%
        correlate() %>%
        shave() %>%
        stretch() %>%
        unite(col = "comparison",x:y,sep = ":") %>%
        drop_na() %>%
        add_column(cluster = clust_num)
    mydata <- bind_rows(mydata,sub_data)
}
mydata$cluster <- as.numeric(mydata$cluster)
mydata_spread <- mydata %>%
    spread(comparison,r) %>%
    column_to_rownames("cluster") %>%
    t()

pheatmap(mydata_spread,
         breaks=seq(-1,1,0.05),
         color = diverging_hcl(40, "Blue-Red 3"),
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 16,
         angle_col = 90)







# Figure 1-Figure supplement 6A
FeaturePlot(object = VNC.combined,
            min.cutoff = 0,
            pt.size = 0.01,
            cols.use = c("lightgrey","black"),
            nCol = 3,
            no.axes = T,
            no.legend = F,
            coord.fixed = T,
            features.plot = c("nGene","nUMI","prop.mito")
)


# Figure 1-Figure supplement 6B
cluster_by_num_markers <- cluster_markers %>% group_by(cluster) %>% summarise(num_mar = length(gene))
temp2 <- tibble(cluster = as.numeric(VNC.combined@meta.data$res.12))
temp2 <- temp2 %>% left_join(cluster_by_num_markers)

VNC.combined@meta.data$num_markers <- temp2$num_mar
FeaturePlot(object = VNC.combined,
                         features.plot = "num_markers",
                         pt.size = 1,
                         cols.use = c("lightgrey","black"),
                         coord.fixed = T,
                         no.axes = T,
                         no.legend = F,
                         do.return = TRUE,
                         max.cutoff = 100)

custer_by_max_logfc <- cluster_markers %>% group_by(cluster) %>% summarise(max_logFC = max(avg_logFC))
temp2 <- temp2 %>% left_join(cluster_by_max_logfc)

VNC.combined@meta.data$max_logFC <- temp2$max_logFC
FeaturePlot(object = VNC.combined,
                         features.plot = "max_logFC",
                         pt.size = 1,
                         cols.use = c("lightgrey","black"),
                         coord.fixed = T,
                         no.axes = T,
                         no.legend = F,
                         do.return = TRUE,
                         max.cutoff = 4)


# Figure 1-Figure supplement 7A
FeaturePlot(object = VNC.combined,
            no.axes = T,
            nCol = 5,
            coord.fixed = T,
            pt.size = 0.1,
            no.legend = T,
            cols.use = c("lightgrey","black"),
            features.plot = c("elav","nSyb","noe","repo","MRE16"))


# Figure 1-Figure supplement 7B
NeuronExpByCluster <- AverageExpression(object = vnc,
										genes.use = c("noe",
													"elav",
													"nSyb",
													"para",
													"VAChT",
													"ChAT",
													"Gad1",
													"VGAT",
													"VGlut",
													"MRE16",
													"repo"), 
										use.scale = T)
pheatmap(NeuronExpByCluster,
         breaks=seq(-2,2,0.05),
         color = diverging_hcl(80, "Vik"),
         border_color = NA,
         treeheight_row = 50,
         treeheight_col = 50,
         cluster_rows = F,
         cluster_cols = T,
         fontsize = 16,
         angle_col = 90)




##################################################
##    figure 2    ################################
##################################################


tale.markers <- cluster_markers %>% 
    filter(gene %in% homeo.list$gene) %>% 
    select(gene,cluster) %>% 
    left_join(y = homeo.list,by = "gene") %>% 
    filter(group=="tale") %>% 
    arrange(gene)
tale.markers
circos.par(start.degree = 90, clock.wise = FALSE)
grid.col = c(qualitative_hcl(length(unique(tale.markers$gene)), "Dark 3"), rep("grey",length(unique(tale.markers$cluster))))
chordDiagram(x = tale.markers[1:2],annotationTrack = "grid",preAllocateTracks = 1,grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
circos.clear()



lim.markers <- cluster_markers %>% 
    filter(gene %in% homeo.list$gene) %>% 
    select(gene,cluster) %>% 
    left_join(y = homeo.list,by = "gene") %>% 
    filter(group=="lim") %>% 
    arrange(gene)
circos.par(start.degree = 90, clock.wise = FALSE)
grid.col = c(qualitative_hcl(length(unique(lim.markers$gene)), "Dark 3"), rep("grey",length(unique(lim.markers$cluster))))
chordDiagram(x = lim.markers[1:2],annotationTrack = "grid",preAllocateTracks = 1,grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
circos.clear()





# Figure 2-Figure supplement 1A
waddell <- readRDS("R_Objects/waddell_brain.rds")
TSNEPlot(object = waddell,coord.fixed = T,no.axes = T,do.label = T,pt.size = 1)

# Figure 2-Figure supplement 1B
aerts <- readRDS("R_Objects/aerts_brain.rds")
TSNEPlot(object = aerts,coord.fixed = T,no.axes = T,do.label = T,pt.size = 0.2)



# Figure 2-Figure supplement 1C
aerts_markers <- read_csv("aerts_brain_markers.csv",col_names = T)
waddell_markers <- read_csv("waddell_brain_markers.csv",col_names = T)
goodwin_markers <- cluster_markers

aerts_markers <- aerts_markers %>% filter(avg_logFC > 0.5 & p_val_adj < 0.05)
waddell_markers <- waddell_markers %>% filter(avg_logFC > 0.5 & p_val_adj < 0.05)
goodwin_markers <- goodwin_markers %>% filter(avg_logFC > 0.5 & p_val_adj < 0.05)

length(unique(aerts_markers$gene))
length(unique(waddell_markers$gene))
length(unique(goodwin_markers$gene))

goodwin_markers <- goodwin_markers %>% 
    filter(cluster != 116) %>% 
    filter(cluster != 99)
length(unique(goodwin_markers$gene))

fit <- euler(c("Aerts" = length(setdiff(setdiff(unique(aerts_markers$gene),unique(waddell_markers$gene)),unique(goodwin_markers$gene))),
               "Waddell" = length(setdiff(setdiff(unique(waddell_markers$gene),unique(aerts_markers$gene)),unique(goodwin_markers$gene))),
               "Goodwin" = length(setdiff(setdiff(unique(goodwin_markers$gene),unique(aerts_markers$gene)),unique(waddell_markers$gene))),
               "Aerts&Waddell" = length(setdiff(intersect(unique(aerts_markers$gene),unique(waddell_markers$gene)),unique(goodwin_markers$gene))),
               "Aerts&Goodwin" = length(setdiff(intersect(unique(aerts_markers$gene),unique(goodwin_markers$gene)),unique(waddell_markers$gene))),
               "Waddell&Goodwin" = length(setdiff(intersect(unique(waddell_markers$gene),unique(goodwin_markers$gene)),unique(aerts_markers$gene))),
               "Aerts&Waddell&Goodwin" = length(intersect(unique(aerts_markers$gene),intersect(unique(waddell_markers$gene),unique(goodwin_markers$gene))))))
plot(fit,
     quantities = T,
     fills = list(fill = c("#F8766D", "#00BA38", "#619CFF"), alpha = 1)
)






# Figure 2-Figure supplement 2A
igsf_list <- read_tsv("geneLists/igsf_converted.tsv",col_names = T)
exp_igsfs <- intersect(igsf_list$current_symbol,exp_genes$gene)
nonexp_igsfs <- setdiff(igsf_list$current_symbol,exp_genes$gene)
igsfs_clusMarkers <- intersect(igsf_list$current_symbol,unique(cluster_markers$gene))

myDotPlot(object = VNC.combined,
          x.lab.rot = TRUE,
          col.min = -2,
          col.max = 2,
          dot.min = 0.01,
          plot.legend = TRUE,
          cols.use = c("#5A90BB", "#F1F1F1", "#A9854D"),
          do.return = TRUE,
          genes.plot = igsfs_clusMarkers) + 
    coord_flip()



# Figure 2-Figure supplement 1B
DIPdpr.list <- read_tsv(file = "geneLists/DIPdpr.tsv",col_names = F)
dipdpr.markers <- cluster_markers %>% filter(gene %in% DIPdpr.list$X1) %>% select(gene,cluster) %>% arrange(gene)
new_order <- c("dpr20","dpr18","dpr15","dpr11","dpr17","dpr3","dpr2","dpr1","dpr10","dpr6","dpr8","dpr21","dpr9","dpr13","CG45781","DIP-zeta","DIP-epsilon","DIP-delta","CG31814","DIP-eta","DIP-theta","DIP-gamma","DIP-beta","DIP-alpha")

dipdpr.markers <- dipdpr.markers %>%
    mutate(gene =  factor(gene, levels = new_order)) %>%
    arrange(gene)

circos.par(start.degree = 90, clock.wise = FALSE)
grid.col = c(rep("#E16A86",14),rep("#00AD9A",10), rep("grey",length(unique(dipdpr.markers$cluster))))
chordDiagram(x = dipdpr.markers,annotationTrack = "grid",preAllocateTracks = 1,grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)



# Figure 2-Figure supplement 1C
beatside.list <- read_tsv(file = "geneLists/beatslide.tsv",col_names = F)
beat.markers <- cluster_markers %>% filter(gene %in% beatside.list$X1) %>% select(gene,cluster) %>% arrange(gene)
new_order <- c("beat-VII","beat-VI","beat-IIa","beat-Ia","beat-Ic","beat-Ib","beat-Vc","beat-Va","beat-Vb","beat-IV","beat-IIIc","beat-IIIb","beat-IIIa","CG34371","CG34114","CG12484","side","CG34113","CG14372","CG42313")
beat.markers <- beat.markers %>%
    mutate(gene =  factor(gene, levels = new_order)) %>%
    arrange(gene)

circos.par(start.degree = 90, clock.wise = FALSE)
grid.col = c(rep("#E16A86",13),rep("#00AD9A",7), rep("grey",length(unique(beat.markers$cluster))))
chordDiagram(x = beat.markers,annotationTrack = "grid",preAllocateTracks = 1,grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
circos.clear()



# Figure 2-Figure supplement 3A
gpcr_list <- read_tsv("geneLists/GPCRs.tsv",col_names = F)
exp_gpcrs <- intersect(gpcr_list$X1,exp_genes$gene)
nonexp_gpcrs <- setdiff(gpcr_list$X1,exp_genes$gene)
gpcrs_clusMarkers <- intersect(gpcr_list$X1,unique(cluster_markers$gene))

myDotPlot(object = VNC.combined,
          x.lab.rot = TRUE,
          col.min = -2,
          col.max = 2,
          dot.min = 0.01,
          plot.legend = TRUE,
          cols.use = c("#5A90BB", "#F1F1F1", "#A9854D"),
          do.return = TRUE,
          genes.plot = gpcrs_clusMarkers) + 
    coord_flip()


# Figure 2-Figure supplement 4A
ion_list <- read_tsv("geneLists/ionchannels.tsv",col_names = F)
exp_ions <- intersect(ion_list$X1,exp_genes$gene)
nonexp_ions <- setdiff(ion_list$X1,exp_genes$gene)
ions_clusMarkers <- intersect(ion_list$X1,unique(cluster_markers$gene))

myDotPlot(object = VNC.combined,
        x.lab.rot = TRUE,
        col.min = -2,
        col.max = 2,
        dot.min = 0.01,
        plot.legend = TRUE,
        cols.use = c("#5A90BB", "#F1F1F1", "#A9854D"),
        do.return = TRUE,
        genes.plot = ions_clusMarkers) + 
    coord_flip()




# Figure 2-Figure supplement 5A
hdtf_list <- read_tsv("geneLists/homeodomainTF.tsv",col_names = TRUE)
exp_hdtfs <- intersect(hdtf_list$gene,exp_genes$gene)
nonexp_hdtfs <- setdiff(hdtf_list$gene,exp_genes$gene)
hdtf_clusMarkers <- intersect(hdtf_list$gene,unique(cluster_markers$gene))

myDotPlot(object = VNC.combined,
          x.lab.rot = TRUE,
          col.min = -2,
          col.max = 2,
          dot.min = 0.01,
          plot.legend = TRUE,
          cols.use = c("#5A90BB", "#F1F1F1", "#A9854D"),
          do.return = TRUE,
          genes.plot = hdtf_clusMarkers) + 
    coord_flip()



# Figure 2-Figure supplement 6A
number_of_class_in_group <- function(object,
                                     gene_list,
                                     accept_low = 0,
                                     cells_use = NULL,
                                     cluster = FALSE,
                                     use_raw = TRUE,
                                     use_scaled = FALSE,
                                     thresh_min = 0,
                                     det_thresh = 0) {
    
    if (length(setdiff(gene_list, object@raw.data@Dimnames[[1]]))>0) {
        cat("The following genes were not found in the dataset:\n",
            setdiff(gene_list, object@raw.data@Dimnames[[1]]),
            "\n")
    }
    exp_gene_list <- intersect(gene_list, object@raw.data@Dimnames[[1]])
    if (cluster) {
        det_gene_data <- dplyr::as_tibble(t(AverageDetectionRate(object = object,
                                                                 thresh.min = thresh_min))) %>%
            select(one_of(exp_gene_list))
        exp_gene_data <- dplyr::as_tibble(t(AverageExpression(object = object,
                                                              genes.use = exp_gene_list,
                                                              use.raw = use_raw,
                                                              use.scale = use_scaled)))
        exp_gene_data[det_gene_data < det_thresh] = 0
        
        counts <- dplyr::tibble(cluster_id = rownames(exp_gene_data)-1,
                                count = rowSums(exp_gene_data>accept_low))
        return(counts)
    } else {
        exp_gene_data <- FetchData(object = object,
                                   vars.all = exp_gene_list,
                                   cells.use = cells_use,
                                   use.raw = use_raw,
                                   use.scaled = use_scaled)
        counts <- dplyr::tibble(cell_id = rownames(exp_gene_data),
                                count = rowSums(exp_gene_data>accept_low))
        return(counts)
    }
}



igsf_AvgDet <- dplyr::as_tibble(t(AverageDetectionRate(object = VNC.combined,
                                                       thresh.min = 0))) %>%
    select(one_of(exp_igsf))
igsf_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = VNC.combined,
                                                    genes.use = exp_igsf,
                                                    use.raw = T,
                                                    use.scale = F)))

igsf_big_counts <- tibble(cluster = character(),
                          percent_cells = numeric(),
                          min_UMI = numeric(),
                          counts = numeric())

temp_seu <- vnc
for (umi_thresh in seq(0,8,2)) {
    temp_seu@raw.data[temp_seu@raw.data < umi_thresh] = 0
    igsf_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = temp_seu,
                                                        genes.use = exp_igsf,
                                                        use.raw = T,
                                                        use.scale = F)))
    
    for (det_thresh in seq(0,0.6,0.1)) {
        temp <- igsf_AvgExp
        temp[igsf_AvgDet < det_thresh] = 0
        
        temp_counts <- tibble(cluster = rownames(temp),
                              percent_cells = det_thresh,
                              min_UMI = umi_thresh,
                              counts = rowSums(temp>0))
        
        igsf_big_counts <- bind_rows(igsf_big_counts,temp_counts)
    }
}

igsf_big_counts_sub <- igsf_big_counts %>% filter(percent_cells < 0.6)
igsf_big_counts_sub$percent_cells <-as.character(igsf_big_counts_sub$percent_cells) 

for (umi_thresh in seq(0,8,2)) {
    igsf_per_cell <- number_of_class_in_group(object = VNC.combined,gene_list = igsf_list$current_symbol,accept_low = umi_thresh,use_raw = TRUE)
    temp_counts <- tibble(cluster = igsf_per_cell$cell_id,
                          percent_cells = "single_cell",
                          min_UMI = umi_thresh,
                          counts = igsf_per_cell$count)
    
    igsf_big_counts_sub <- bind_rows(igsf_big_counts_sub,temp_counts)
    
}

igsf_big_counts_sub$min_UMI <-as.character(igsf_big_counts_sub$min_UMI) 
ggplot(igsf_big_counts_sub, aes(x=percent_cells, y=counts, fill=min_UMI)) + 
    geom_boxplot() +
    theme(text = element_text(size=20))



# Figure 2-Figure supplement 6B
gpcr_AvgDet <- dplyr::as_tibble(t(AverageDetectionRate(object = VNC.combined,
                                                       thresh.min = 0))) %>%
    select(one_of(exp_gpcrs))
gpcr_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = VNC.combined,
                                                    genes.use = exp_gpcrs,
                                                    use.raw = T,
                                                    use.scale = F)))


gpcr_big_counts <- tibble(cluster = character(),
                          percent_cells = numeric(),
                          min_UMI = numeric(),
                          counts = numeric())

temp_seu <- vnc
for (umi_thresh in seq(0,8,2)) {
    temp_seu@raw.data[temp_seu@raw.data < umi_thresh] = 0
    gpcr_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = temp_seu,
                                                        genes.use = exp_gpcrs,
                                                        use.raw = T,
                                                        use.scale = F)))
    
    for (det_thresh in seq(0,0.6,0.1)) {
        temp <- gpcr_AvgExp
        temp[gpcr_AvgDet < det_thresh] = 0
        
        temp_counts <- tibble(cluster = rownames(temp),
                              percent_cells = det_thresh,
                              min_UMI = umi_thresh,
                              counts = rowSums(temp>0))
        
        gpcr_big_counts <- bind_rows(gpcr_big_counts,temp_counts)
    }
}


gpcr_big_counts_sub <- gpcr_big_counts %>% filter(percent_cells < 0.6)
gpcr_big_counts_sub$percent_cells <-as.character(gpcr_big_counts_sub$percent_cells) 

for (umi_thresh in seq(0,8,2)) {
    gpcrs_per_cell <- number_of_class_in_group(object = VNC.combined,gene_list = gpcr_list$X1,accept_low = umi_thresh,use_raw = TRUE)
    temp_counts <- tibble(cluster = gpcrs_per_cell$cell_id,
                          percent_cells = "single_cell",
                          min_UMI = umi_thresh,
                          counts = gpcrs_per_cell$count)
    
    gpcr_big_counts_sub <- bind_rows(gpcr_big_counts_sub,temp_counts)
    
}
gpcr_big_counts_sub$min_UMI <-as.character(gpcr_big_counts_sub$min_UMI) 

ggplot(gpcr_big_counts_sub, aes(x=percent_cells, y=counts, fill=min_UMI)) + 
    geom_boxplot() +
    theme(text = element_text(size=20))






# Figure 2-Figure supplement 6C
ionCh_AvgDet <- dplyr::as_tibble(t(AverageDetectionRate(object = VNC.combined,
                                                        thresh.min = 0))) %>%
    select(one_of(exp_ions))
ionCh_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = VNC.combined,
                                                     genes.use = exp_ions,
                                                     use.raw = T,
                                                     use.scale = F)))



ionCh_big_counts <- tibble(cluster = character(),
                           percent_cells = numeric(),
                           min_UMI = numeric(),
                           counts = numeric())


temp_seu <- vnc
for (umi_thresh in seq(0,8,2)) {
    temp_seu@raw.data[temp_seu@raw.data < umi_thresh] = 0
    ionCh_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = temp_seu,
                                                         genes.use = exp_ions,
                                                         use.raw = T,
                                                         use.scale = F)))
    
    for (det_thresh in seq(0,0.6,0.1)) {
        temp <- ionCh_AvgExp
        temp[ionCh_AvgDet < det_thresh] = 0
        
        temp_counts <- tibble(cluster = rownames(temp),
                              percent_cells = det_thresh,
                              min_UMI = umi_thresh,
                              counts = rowSums(temp>0))
        
        ionCh_big_counts <- bind_rows(ionCh_big_counts,temp_counts)
    }
}

ionCh_big_counts_sub <- ionCh_big_counts %>% filter(percent_cells < 0.6)
ionCh_big_counts_sub$percent_cells <-as.character(ionCh_big_counts_sub$percent_cells) 

for (umi_thresh in seq(0,8,2)) {
    ionCh_per_cell <- number_of_class_in_group(object = VNC.combined,gene_list = ionCh_list$X1,accept_low = umi_thresh,use_raw = TRUE)
    temp_counts <- tibble(cluster = ionCh_per_cell$cell_id,
                          percent_cells = "single_cell",
                          min_UMI = umi_thresh,
                          counts = ionCh_per_cell$count)
    
    ionCh_big_counts_sub <- bind_rows(ionCh_big_counts_sub,temp_counts)
    
}

ionCh_big_counts_sub$min_UMI <-as.character(ionCh_big_counts_sub$min_UMI) 
ggplot(ionCh_big_counts_sub, aes(x=percent_cells, y=counts, fill=min_UMI)) + 
    geom_boxplot() +
    theme(text = element_text(size=20))






# Figure 2-Figure supplement 6D
hdtf_AvgDet <- dplyr::as_tibble(t(AverageDetectionRate(object = VNC.combined,
                                                       thresh.min = 0))) %>%
    select(one_of(exp_hdtfs))
hdtf_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = VNC.combined,
                                                    genes.use = exp_hdtfs,
                                                    use.raw = T,
                                                    use.scale = F)))

hdtf_big_counts <- tibble(cluster = character(),
                          percent_cells = numeric(),
                          min_UMI = numeric(),
                          counts = numeric())


temp_seu <- vnc
for (umi_thresh in seq(0,8,2)) {
    temp_seu@raw.data[temp_seu@raw.data < umi_thresh] = 0
    hdtf_AvgExp <- dplyr::as_tibble(t(AverageExpression(object = temp_seu,
                                                        genes.use = exp_hdtfs,
                                                        use.raw = T,
                                                        use.scale = F)))
    
    for (det_thresh in seq(0,0.6,0.1)) {
        temp <- hdtf_AvgExp
        temp[hdtf_AvgDet < det_thresh] = 0
        
        temp_counts <- tibble(cluster = rownames(temp),
                              percent_cells = det_thresh,
                              min_UMI = umi_thresh,
                              counts = rowSums(temp>0))
        
        hdtf_big_counts <- bind_rows(hdtf_big_counts,temp_counts)
    }
}

hdtf_big_counts_sub <- hdtf_big_counts %>% filter(percent_cells < 0.6)
hdtf_big_counts_sub$percent_cells <-as.character(hdtf_big_counts_sub$percent_cells) 

for (umi_thresh in seq(0,8,2)) {
    hdtfs_per_cell <- number_of_class_in_group(object = VNC.combined,gene_list = hdtfs_list$gene,accept_low = umi_thresh,use_raw = TRUE)
    temp_counts <- tibble(cluster = hdtfs_per_cell$cell_id,
                          percent_cells = "single_cell",
                          min_UMI = umi_thresh,
                          counts = hdtfs_per_cell$count)
    
    hdtf_big_counts_sub <- bind_rows(hdtf_big_counts_sub,temp_counts)
    
}

hdtf_big_counts_sub$min_UMI <-as.character(hdtf_big_counts_sub$min_UMI) 
ggplot(hdtf_big_counts_sub, aes(x=percent_cells, y=counts, fill=min_UMI)) + 
    geom_boxplot() +
    theme(text = element_text(size=20))





##################################################
##    figure 3    ################################
##################################################


# Figure 3A
homeo.list <- read_tsv(file = "geneLists/homeodomainTF.tsv",col_names = T)
hox.markers <- cluster_markers %>% 
    filter(gene %in% homeo.list$gene) %>% 
    select(gene,cluster) %>% 
    left_join(y = homeo.list,by = "gene") %>% 
    filter(group=="hox-like") %>% 
    arrange(gene)
hox.markers
hox.markers <- hox.markers %>% filter(gene != "eve" & gene != "pb")
circos.par(start.degree = 90, clock.wise = FALSE)
grid.col = c(qualitative_hcl(length(unique(hox.markers$gene)), "Dark 3"), rep("grey",length(unique(hox.markers$cluster))))
chordDiagram(x = hox.markers[1:2],annotationTrack = "grid",preAllocateTracks = 1,grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
circos.clear()



# Figure 3B
hoxGenes <- FetchData(object = VNC.combined,vars.all = c("Antp","Ubx","abd-A","Abd-B"))
hoxGenes.Corr <- round(cor(hoxGenes,method = "pearson"), 3)
hoxGenes.Corr.p <- cor_pmat(hoxGenes)
hoxGenesPlot <- ggcorrplot(hoxGenes.Corr, show.diag = FALSE, p.mat = hoxGenes.Corr.p, insig = "blank", outline.color = "white")
hoxGenesPlot + scale_fill_gradient2(limit = c(-0.4,0.4), low = "blue", high =  "red", mid = "white", midpoint = 0)


# Figure 3C
FeaturePlot(object = VNC.combined,
            nCol = 4,
            features.plot = c("Antp","Ubx","abd-A","Abd-B"),
            pt.size = 0.4,
            cols.use = c("lightgrey","black"),
            coord.fixed = T,
            no.axes = T,
            no.legend = F)



# Figure 3E
hoxGenes <- hoxGenes %>% mutate(Cells = seuObject@data@Dimnames[[2]])
hoxGenes <- hoxGenes %>% mutate(
    hox = if_else((Antp>Ubx & Antp>`abd-A` & Antp>`Abd-B` & Antp>0), "MesoNm",
                  if_else((Ubx>Antp & Ubx>`abd-A` & Ubx>`Abd-B` & Ubx>0), "MetaNm",
                          if_else(((`abd-A`>Antp & `abd-A`>Ubx & `abd-A`>0) | (`Abd-B`>Antp & `Abd-B`>Ubx & `Abd-B`>0)), "ANm",
                                  if_else((Antp==0 & Ubx==0 & `abd-A`==0 & `Abd-B`==0),"ProNm","unknown")))))

VNC.combined@meta.data$Hox <- SubGenes$hox
seuObject <- VNC.combined
new.cluster.ids <- as_factor(seuObject@meta.data$Hox)
seuObject <- SetIdent(object = seuObject, ident.use = new.cluster.ids)

p1 <- TSNEPlot(object = VNC.combined,
               plot.order = c("ProNm"),
               do.return = T,
               group.by = "Hox",
               colors.use = c("lightgrey","lightgrey","lightgrey","lightgrey","#F8766D"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)
p2 <- TSNEPlot(object = VNC.combined,
               plot.order = c("MesoNm"),
               do.return = T,
               group.by = "Hox",
               colors.use = c("lightgrey","lightgrey","lightgrey","lightgrey","#7CAE00"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)
p3 <- TSNEPlot(object = VNC.combined,
               plot.order = c("MetaNm"),
               do.return = T,
               group.by = "Hox",
               colors.use = c("lightgrey","lightgrey","lightgrey","lightgrey","#00BFC4"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)
p4 <- TSNEPlot(object = VNC.combined,
               plot.order = c("ANm"),
               do.return = T,
               group.by = "Hox",
               colors.use = c("lightgrey","lightgrey","lightgrey","lightgrey","#C77CFF"),
               pt.size = 0.1,
               coord.fixed = TRUE,
               no.axes = TRUE,
               do.label = FALSE,
               no.legend = TRUE)

ggarrange(plotlist = list(p1,p2,p3,p4),
          labels = c("ProNm","MesoNm","MetaNm","ANm"),
          ncol = 2,
          nrow = 2)

rm(seuObject)
gc()











##################################################
##    figure 4    ################################
##################################################



# Figure 4A
lin_by_cluster <- tibble(hemilineage = c("lin0A","lin0A","lin0A",
                                         "lin11B",
                                         "lin12B","lin12B","lin12B","lin12B","lin12B",
                                         "lin1B","lin1B",
                                         "lin9A","lin9A","lin9A","lin9A",
                                         "lin4B","lin4B",
                                         "lin21A","lin21A",
                                         "lin14A","lin14A","lin14A",
                                         "lin8A & 24B","lin8A & 24B","lin8A & 24B",
                                         "lin23B","lin23B","lin23B","lin23B",
                                         "lin13B","lin13B",
                                         "lin5B & 6B","lin5B & 6B","lin5B & 6B","lin5B & 6B","lin5B & 6B",
                                         "lin13A & 19A","lin13A & 19A","lin13A & 19A","lin13A & 19A","lin13A & 19A","lin13A & 19A","lin13A & 19A","lin13A & 19A","lin13A & 19A",
                                         "lin6A","lin6A","lin6A",
                                         "lin9B",
                                         "lin8B","lin8B","lin8B","lin8B","lin8B",
                                         "lin10B","lin10B","lin10B",
                                         "lin1A",
                                         "lin16B","lin16B",
                                         "lin20A & 22A","lin20A & 22A","lin20A & 22A","lin20A & 22A","lin20A & 22A",
                                         "lin15B","lin15B","lin15B"),
                         cluster = as.factor(c(22,88,112,
                                               38,
                                               83,94,73,30,81,
                                               12,47,
                                               50,56,57,31,
                                               100,0,
                                               1,107,
                                               41,74,13,
                                               69,6,110,
                                               35,51,67,93,
                                               25,17,
                                               89,3,20,97,87,
                                               79,75,48,19,82,59,96,102,65,
                                               28,9,66,
                                               54,
                                               117,53,8,49,76,
                                               39,91,68,
                                               16,
                                               5,46,
                                               108,33,34,78,14,
                                               36,52,80)))


cellid_by_cluster <- as_tibble(rownames_to_column(as.data.frame(VNC.combined@ident)))
cellid_by_cluster <- cellid_by_cluster %>% rename(cell_id = "rowname", cluster = "VNC.combined@ident")
cell_cluster_lineage <- full_join(cellid_by_cluster,lin_by_cluster)   
cell_cluster_lineage$hemilineage <- cell_cluster_lineage$hemilineage %>% replace_na("unknown")    

VNC.combined@meta.data$predicted_lineage <- cell_cluster_lineage$hemilineage

lin_markers <- read_tsv("geneLists/lineageMarkers.tsv",col_names = FALSE)
lin_markers <- c(lin_markers$X1, "ap","VAChT","VGlut","Gad1")

temp_vnc <- VNC.combined
class(temp_vnc@ident)
class(temp_VNC.combined@meta.data$predicted_lineage)
temp_vnc <- SetIdent(object = temp_vnc,
                     ident.use = as.factor(temp_VNC.combined@meta.data$predicted_lineage)
)

linExpByCluster <- AverageExpression(object = temp_vnc,genes.use = lin_markers, use.scale = T)
pheatmap(t(linExpByCluster),
         breaks=seq(-2,2,0.1),
         color = diverging_hcl(40, "Vik"),
         border_color = NA,
         cellwidth = 20,
         cellheight = 20,
         treeheight_row = 100,
         treeheight_col = 100,
         cluster_cols = T,
         cluster_rows = F,
         fontsize = 16,
         angle_col = 90)




# Figure 4B
TSNEPlot(object = SubsetData(object = VNC.combined,
                             ident.remove = c(23,24,70,98,99,105,106,116)),
         coord.fixed = T,
         no.axes = T,
         pt.size = 0.4,
         do.label = T,
         group.by = "predicted_lineage"
)



# Figure 4C
hemilineage_markers <- FindAllMarkers(object = temp_vnc,
                                      logfc.threshold = 0.5,
                                      test.use = "negbinom",
                                      only.pos = TRUE,
                                      do.print = TRUE,
                                      latent.vars = "Replicate")

hemilineage_markers_filtered <- hemilineage_markers %>% 
    filter(p_val_adj < 0.05)
write.csv(hemilineage_markers,"hemilineage_markers.csv")
write.csv(hemilineage_markers_filtered,"hemilineage_markers_filtered.csv")

top10 <- hemilineage_markers_filtered %>% group_by(cluster) %>% top_n(10, avg_logFC)
lin_markers <- c(lin_markers,"ChAT","VGAT")
top_filtered <- hemilineage_markers_filtered %>%
    filter(cluster != "unknown") %>% 
    filter(pct.1 > 0.2) %>% 
    filter(pct.2 < 0.1) %>% 
    filter(!gene %in% lin_markers) %>% 
    group_by(cluster) %>% 
    top_n(10, avg_logFC)
genes_to_plot <- unique(top_filtered$gene)

ExpByCluster <- AverageExpression(object = temp_vnc,genes.use = genes_to_plot, use.scale = T)
pheatmap(t(ExpByCluster),
         breaks=seq(-2,2,0.1),
         color = diverging_hcl(40, "Vik"),
         border_color = NA,
         cellwidth = 20,
         cellheight = 20,
         treeheight_row = 100,
         treeheight_col = 100,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 16,
         angle_col = 90)





# Figure 4-Figure supplement 1A
lin_markers <- setdiff(lin_markers,c("VGAT","ChAT"))
linExpByCluster <- AverageExpression(object = SubsetData(object = VNC.combined,cells.use = WhichCells(vnc,ident.remove = c(18,23,24,55,70,98,99,105,106,116))),genes.use = lin_markers, use.scale = T)
pheatmap(linExpByCluster,
         breaks=seq(-2,2,0.1),
         color = diverging_hcl(40, "Vik"),
         border_color = NA,
         cellwidth = 14,
         cellheight = 14,
         cluster_cols = T,
         cluster_rows = T,
         fontsize = 16,
         angle_col = 90)



# Figure 4-Figure supplement 1B
lin_markers <- sort(setdiff(lin_markers,c("VAChT","VGlut","Gad1")))
FeaturePlot(object = VNC.combined,
            nCol = 5,
            features.plot = lin_markers,
            pt.size = 0.4,
            cols.use = c("lightgrey","black"),
            coord.fixed = T,
            no.axes = T,
            no.legend = T)





# Figure 4-Figure supplement 2
abdb.cells <- WhichCells(object = VNC.combined,subset.name = "Abd-B",accept.low = 0)
abda.cells <- WhichCells(object = VNC.combined,subset.name = "abd-A",accept.low = 0)
anm.cells <- union(abda.cells,abdb.cells)
thorax.cells <- setdiff(WhichCells(object = VNC.combined),anm.cells)

pros.pos.cells <- WhichCells(object = VNC.combined,subset.name = "pros",accept.low = 0)
pros.neg.cells <- setdiff(WhichCells(object = VNC.combined), pros.pos.cells)
imp.pos.cells <- WhichCells(object = VNC.combined,subset.name = "Imp",accept.low = 0)
imp.neg.cells <- setdiff(WhichCells(object = VNC.combined), imp.pos.cells)
embryonic.cells <- intersect(imp.pos.cells,pros.neg.cells)

postemb.thor.cells <- intersect(thorax.cells,pros.pos.cells)

#0A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(22,88,112)))
#11B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(38)))
#12B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(83,94,73,30,81)))
#1B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(12,47)))
#9A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(50,56,57,31)))
#4B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(100,0)))
#21A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(1,107)))
#14A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(41,74,13)))
#8A 24B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(69,6,110)))
#23B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(35,51,67,93)))
#13B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(25,17)))
#5B 6B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(89,3,20,97,87)))
#13A 19A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(79,75,48,19,82,59,96,102,65)))
#6A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(28,9,66)))
#9B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(54)))
#8B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(117,53,8,49,76)))
#10B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(39,91,68)))
#1A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(16)))
#16B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(5,46)))
#20A/22A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(108,33,34,78,14)))
#15B
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(36,52,80)))
#12A
length(WhichCells(object = VNC.combined, cells.use = postemb.thor.cells,ident = c(40)))


lin.cell.counts <- read_csv("lineage_cell_counts.csv")
ggplot(lin.cell.counts, aes(x=lacin_counts, y=my_counts)) +
    geom_point() + 
    theme_grey() + 
    geom_smooth(method='lm',formula=y~x, col="black")





# Figure 4-Figure supplement 2B
FeaturePlot(object = VNC.combined,
            nCol = 1,
            coord.fixed = T,
            no.axes = T,
            no.legend = F,
            pt.size = 0.4,
            cols.use = c("lightgrey","black"),
            features.plot = "br")

# Figure 4-Figure supplement 2C
FeaturePlot(object = VNC.combined,
            nCol = 1,
            coord.fixed = T,
            no.axes = T,
            no.legend = F,
            pt.size = 0.4,
            cols.use = c("lightgrey","black"),
            features.plot = c("pros","dati","mamo","Imp"))



# Figure 4-Figure supplement 2D
devGenes <- FetchData(object = VNC.combined,vars.all = c("pros","dati","Imp","mamo"))
devGenes.Corr <- round(cor(devGenes,method = "pearson"), 3)
devGenes.Corr.p <- cor_pmat(devGenes)
devGenesPlot <- ggcorrplot(devGenes.Corr, show.diag = TRUE, p.mat = devGenes.Corr.p, insig = "blank", outline.color = "white",lab = TRUE)
devGenesPlot + scale_fill_gradient2(limit = c(-0.5,0.5), low = "blue", high =  "red", mid = "white", midpoint = 0)

# Figure 4-Figure supplement 2E
devGenes <- FetchData(object = VNC.combined,vars.all = c("pros","Imp","Antp","Ubx","abd-A","Abd-B"))
devGenes.Corr <- round(cor(devGenes,method = "pearson"), 3)
devGenes.Corr.p <- cor_pmat(devGenes)
devGenesPlot <- ggcorrplot(devGenes.Corr, show.diag = TRUE, p.mat = devGenes.Corr.p, insig = "blank", outline.color = "white",lab = TRUE)
devGenesPlot + scale_fill_gradient2(limit = c(-0.5,0.5), low = "blue", high =  "red", mid = "white", midpoint = 0)

# Figure 4-Figure supplement 2F
devGenes <- FetchData(object = VNC.combined,vars.all = c("pros","Imp","VAChT","VGlut","Gad1"))
devGenes.Corr <- round(cor(devGenes,method = "pearson"), 3)
devGenes.Corr.p <- cor_pmat(devGenes)
devGenesPlot <- ggcorrplot(devGenes.Corr, show.diag = TRUE, p.mat = devGenes.Corr.p, insig = "blank", outline.color = "white",lab = TRUE)
devGenesPlot + scale_fill_gradient2(limit = c(-0.5,0.5), low = "blue", high =  "red", mid = "white", midpoint = 0)







##################################################
##    figure 5    ################################
##################################################


# Figure 5B
FeaturePlot(object = VNC.combined,
            nCol = 1,
            coord.fixed = T,
            no.axes = T,
            no.legend = F,
            pt.size = 0.4,
            cols.use = c("lightgrey","black"),
            features.plot = "acj6")


# Figure 5C
FeaturePlot(object = SubsetData(object = VNC.combined,ident.use = c(35,51,67,93)),
            nCol = 5,
            coord.fixed = T,
            no.axes = T,
            no.legend = F,
            pt.size = 0.4,
            cols.use = c("lightgrey","black"),
            features.plot = c("acj6","unc-4","VAChT","kn","twz"))



# Figure 5-Figure supplement 2A
FeaturePlot(object = VNC.combined,
            nCol = 3,
            coord.fixed = T,
            no.axes = T,
            no.legend = F,
            pt.size = 0.4,
            cols.use = c("lightgrey","black"),
            features.plot = c("Lim3","kn","twz"))







##################################################
##    figure 6    ################################
##################################################




# Figure 6A
FeaturePlot(object = VNC.combined,
            pt.size = 0.01,
            cols.use = c("lightgrey","black"),
            nCol = 5,
            no.axes = T,
            no.legend = F,
            coord.fixed = T,
            features.plot = c("VAChT","ChAT","VGlut","VGAT","Gad1")
            )


# Figure 6B
fan.list <- tibble(gene = c("VAChT","VGlut","Gad1"))
fan.markers <- cluster_markers %>% filter(gene %in% fan.list$gene) %>% select(gene,cluster) %>% arrange(gene)
circos.par(start.degree = 90, clock.wise = FALSE)
grid.col = c(qualitative_hcl(length(unique(fan.markers$gene)), "Dark 3"), rep("grey",length(unique(fan.markers$cluster))))
chordDiagram(x = fan.markers,annotationTrack = "grid",preAllocateTracks = 1,grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
circos.clear()



# Figure 6C
fanGenes <- FetchData(object = VNC.combined,vars.all = c("VAChT","ChAT","VGlut","Gad1"))
ClusterIdent <- data.frame(seuObject@ident)
fanGenes <- fanGenes %>% mutate(
    Identity = ClusterIdent$seuObject.ident,
    Cells = seuObject@data@Dimnames[[2]])

Gad1Clusters <- c(64,66,109,45,81,83,64,12,47,3,112,31,75,87,20,9,28,17,88,96,56,57,38,25,26,22,90,19,48,30,79,102,73,97,50,82,104,119)
VGlutClusters <- c(80,36,52,113,46,54,5,118,41,15,69,110,6,13,1,86,74,107)
VAChTClusters <- c(67,101,78,8,108,76,49,63,33,35,39,51,71,0,7,37,77,68,114,11,10,85,34,14,42,58,2,16,95,117,105,91,103,40,100,62,44,53,111,92,21,93)
mixedClusters <- c(27,4,61,43,29,32,55,18,115,72,84)
intersect(Gad1Clusters,VGlutClusters)
intersect(Gad1Clusters,VAChTClusters)
intersect(VAChTClusters,VGlutClusters)
intersect(mixedClusters,VGlutClusters)
intersect(mixedClusters,VAChTClusters)
intersect(mixedClusters,Gad1Clusters)

Gad1Cells <- WhichCells(object = VNC.combined, ident = Gad1Clusters)
VGlutCells <- WhichCells(object = VNC.combined, ident = VGlutClusters)
VAChTCells <- WhichCells(object = VNC.combined, ident = VAChTClusters)
AllCells <- WhichCells(VNC.combined,ident.remove = c(23,24,70,98,99,106,116))
unknownCells <- setdiff(setdiff(setdiff(AllCells,Gad1Cells),VAChTCells),VGlutCells)

fanGenes <- fanGenes %>% mutate(
    FAN = if_else(is.element(Identity,unlist(VAChTClusters)), "Cholinergic",
                  if_else(is.element(Identity,unlist(VGlutClusters)), "Glutamatergic",
                          if_else(is.element(Identity,unlist(Gad1Clusters)), "GABAergic",
                                  if_else((VAChT>0 & ChAT>0), "Cholinergic",
                                          if_else((VAChT>2 | ChAT>2), "Cholinergic",
                                                  if_else(VGlut>2 & VGlut>Gad1 & VGlut>VAChT, "Glutamatergic",
                                                          if_else(Gad1>2 & Gad1>VAChT & Gad1>VGlut, "GABAergic","unknown"
                                                                  )
                                                          )
                                                  )
                                          )
                                  )
                          )
                  )
)

VNC.combined@meta.data$FAN <- fanGenes$FAN
TSNEPlot(object = SubsetData(object = VNC.combined,
                             cells.use = WhichCells(VNC.combined,
                                                    ident.remove = c(18,23,24,55,70,98,99,105,106,116))),
         group.by = "FAN",
         do.label = F,
         pt.size = 1,
         coord.fixed=T)


# Figure 6D
NeuronTypes <- data.frame(CellType=c("GABAergic","Glutamatergic","Cholinergic","none"),
                          Count=c(length(WhichCells(object = VNC.neurons, subset.name = "FAN",accept.value = "GABAergic")),
                                  length(WhichCells(object = VNC.neurons, subset.name = "FAN",accept.value = "Glutamatergic")),
                                  length(WhichCells(object = VNC.neurons, subset.name = "FAN",accept.value = "Cholinergic")),
                                  length(WhichCells(object = VNC.neurons, subset.name = "FAN",accept.value = "unknown"))
                          )
)
NeuronTypes %>% 
    mutate(prop=round(x = 100*Count/sum(Count),digits = 0)) %>%
    mutate(lab.ypos = cumsum(prop) - 0.5*prop) %>%
    ggplot(aes(x=CellType, y=prop, fill=fct_relevel(CellType,"none","Cholinergic","Glutamatergic","GABAergic"))) +
    geom_bar(stat = "identity", color = "white",size = 4) +
    theme(legend.title=element_blank(),
        legend.text=element_text(size=20),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        plot.subtitle = element_text(size = 12, face = "bold",hjust = 0.5)
    ) +
    ylim(0, 50)






# Figure 6-Figure supplement 1A
fanGenes <- FetchData(object = VNC.combined,vars.all = c("VAChT","ChAT","VGlut","VGAT","Gad1"))
fanGenes.Corr <- round(cor(fanGenes,method = "pearson"), 3)
fanGenes.Corr.p <- cor_pmat(fanGenes)
fanGenesPlot <- ggcorrplot(fanGenes.Corr, show.diag = TRUE, p.mat = fanGenes.Corr.p, insig = "blank", outline.color = "white",lab = TRUE)
fanGenesPlot + scale_fill_gradient2(limit = c(-0.5,0.5), low = "blue", high =  "red", mid = "white", midpoint = 0)



# Figure 6-Figure supplement 1B
Gad1Cells <- WhichCells(object = VNC.combined,subset.name = "Gad1",accept.low = 0)
VGlutCells <- WhichCells(object = VNC.combined,subset.name = "VGlut",accept.low = 0)
VAChTCells <- WhichCells(object = VNC.combined,subset.name = "VAChT",accept.low = 0)


fit <- euler(c("VAChT" = round(100*length(setdiff(setdiff(VAChTCells,VGlutCells),Gad1Cells))/length(WhichCells(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116)))),
               "VGlut" = round(100*length(setdiff(setdiff(VGlutCells,VAChTCells),Gad1Cells))/length(WhichCells(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116)))),
               "Gad1" = round(100*length(setdiff(setdiff(Gad1Cells,VAChTCells),VGlutCells))/length(WhichCells(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116)))),
               "VAChT&VGlut" = round(100*length(setdiff(intersect(VAChTCells,VGlutCells),Gad1Cells))/length(WhichCells(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116)))),
               "VAChT&Gad1" = round(100*length(setdiff(intersect(VAChTCells,Gad1Cells),VGlutCells))/length(WhichCells(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116)))),
               "VGlut&Gad1" = round(100*length(setdiff(intersect(VGlutCells,Gad1Cells),VAChTCells))/length(WhichCells(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116)))),
               "VAChT&VGlut&Gad1" = round(100*length(intersect(VAChTCells,intersect(VGlutCells,Gad1Cells)))/length(WhichCells(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116))))))
plot(fit,
     quantities = T,
     fills = list(fill = c("#F8766D", "#00BA38", "#619CFF"), alpha = 1)
)






# Figure 6-Figure supplement 1C
FANExpByCluster <- AverageExpression(object = SubsetData(object = VNC.combined,
                                                         ident.remove = c(18,23,24,55,70,98,99,105,106,116)
                                                         ),
                                     genes.use = c("VAChT","ChAT","VGlut","Gad1","VGAT"), 
                                     use.scale = T)
pheatmap(FANExpByCluster,
         breaks=seq(-2,2,0.05),
         color = diverging_hcl(80, "Vik"),
         border_color = NA,
         treeheight_row = 100,
         treeheight_col = 100,
         cluster_cols = T,
         cluster_rows = F,
         fontsize = 16,
         angle_col = 90)




# Figure 6-Figure supplement 1D
seuObject <- VNC.combined
new.cluster.ids <- as_factor(seuObject@meta.data$FAN)
seuObject <- SetIdent(object = seuObject, ident.use = new.cluster.ids)
Gad1Cells <- WhichCells(seuObject,subset.name = "Gad1",accept.low = 0)
VAChTCells <- WhichCells(seuObject,subset.name = "VAChT",accept.low = 0)
VGlutCells <- WhichCells(seuObject,subset.name = "VGlut",accept.low = 0)
GenePlot(object = seuObject,gene1 = "VGlut",gene2 = "VAChT",cell.ids = intersect(VAChTCells,VGlutCells))
GenePlot(object = seuObject,gene1 = "Gad1",gene2 = "VAChT",cell.ids = intersect(VAChTCells,Gad1Cells))
GenePlot(object = seuObject,gene1 = "Gad1",gene2 = "VGlut",cell.ids = intersect(Gad1Cells,VGlutCells))










##################################################
##    figure 7    ################################
##################################################



# Figure 7B
FeaturePlot(object = VNC.combined,
            coord.fixed = T,
            pt.size = 0.1,
            cols.use = c("lightgrey","black"),
            no.axes = T,
            no.legend = F,
            min.cutoff = 0,
            features.plot = "Vmat"
)



# Figure 7C
MaleVNC.combined <- readRDS("R_Objects/MaleVNC.combined.rds")
FemaleVNC.combined <- readRDS("R_Objects/FemaleVNC.combined.rds")
Vmat.cells <- WhichCells(object = SubsetData(object = VNC.combined,ident.use = c(72,84)),subset.name = "Vmat",accept.low = 0)
simple.Vmat.cells <- gsub(pattern = ".*ale_", replacement = "", x = Vmat.cells)
male.Vmat.cells <- intersect(simple.Vmat.cells, WhichCells(object = MaleVNC.combined,subset.name = "Vmat",accept.low = 0))
female.Vmat.cells <- intersect(simple.Vmat.cells, WhichCells(object = FemaleVNC.combined,subset.name = "Vmat",accept.low = 0))

Male.Vmat <- SubsetData(object = MaleVNC.combined,cells.use = male.Vmat.cells)
Female.Vmat <- SubsetData(object = FemaleVNC.combined,cells.use = female.Vmat.cells)
Male.Vmat <- FindVariableGenes(Male.Vmat, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001,do.plot = F)
Female.Vmat <- FindVariableGenes(Female.Vmat, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001,do.plot = F)

g.1 <- (Male.Vmat@var.genes)
g.2 <- (Female.Vmat@var.genes)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(Male.Vmat@scale.data))
genes.use <- intersect(genes.use, rownames(Female.Vmat@scale.data))
length(genes.use)

NumberOfCCAs <- 40
Vmat.combined <- RunCCA(Male.Vmat, Female.Vmat, genes.use = genes.use, num.cc = NumberOfCCAs, add.cell.id1 = "Male", add.cell.id2 = "Female")
Vmat.combined <- AlignSubspace(Vmat.combined, reduction.type = "cca", grouping.var = "Sex", dims.align = 1:NumberOfCCAs, do.par = TRUE, num.cores = 8)
saveRDS(Vmat.combined, file = "R_Objects/Vmat.combined.rds")

NumberOfCCAs <- 7
ClusterResolution <- 1.2  
Vmat.combined <- RunTSNE(Vmat.combined,
                     num_threads = 0,
                     verbose = T,
                     reduction.use = "cca.aligned",
                     dims.use = 1:NumberOfCCAs,
                     do.fast = FALSE,
                     check_duplicates = F,
                     max_iter = 10000,
                     theta = 0)
Vmat.combined <- FindClusters(Vmat.combined,
                          reduction.type = "cca.aligned",
                          resolution = ClusterResolution,
                          dims.use = 1:NumberOfCCAs,
                          save.SNN = TRUE)
TSNEPlot(Vmat.combined,
               no.axes = T,
               pt.size = 8,
               do.label = F,
               coord.fixed=T)




# Figure 7D
FeaturePlot(object = Vmat.combined,
            nCol = 2,
            no.axes = T,
            coord.fixed = T,
            pt.size = 4,
            no.legend = T,
            cols.use = c("lightgrey","black"),
            features.plot = c("DAT","SerT","Tbh","Hdc"))




# Figure 7E
Vmat.markers <- FindAllMarkers(object = Vmat.combined,
                                     test.use = "negbinom",
                                     only.pos = T,
                                     print.bar = T,
                                     latent.vars = c("Replicate"),
                                     logfc.threshold = 0.5)
Vmat.markers <- Vmat.markers %>% filter(p_val_adj<0.05)
write.csv(Vmat.markers,"Vmat_markers.csv")

top10 <- Vmat.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
myDotPlot(object = Vmat.combined,genes.plot = top10$gene,x.lab.rot = T,plot.legend = T,col.min = -2,col.max = 2,cols.use = c("lightgrey","white","#F8766D"))
myDotPlot(object = Vmat.combined,genes.plot = top10$gene,x.lab.rot = T,plot.legend = T,col.min = -2,col.max = 2,cols.use = c("lightgrey","white","#7CAE00"))
myDotPlot(object = Vmat.combined,genes.plot = top10$gene,x.lab.rot = T,plot.legend = T,col.min = -2,col.max = 2,cols.use = c("lightgrey","white","#00BFC4"))
myDotPlot(object = Vmat.combined,genes.plot = top10$gene,x.lab.rot = T,plot.legend = T,col.min = -2,col.max = 2,cols.use = c("lightgrey","white","#C77CFF"))





# Figure 7-Figure supplement 1A
TSNEPlot(Vmat.combined,
         group.by ="Replicate",
         no.axes = T,
         pt.size = 8,
         do.label = F,
         coord.fixed=T)



# Figure 7-Figure supplement 1B
rep1 <- SubsetData(object = Vmat.combined,subset.name = "Replicate",accept.value = "FemaleRep1",subset.raw = T)
rep1_avgexp_byclust <- AverageExpression(object = rep1,use.raw = F,show.progress = T)
rep1_avgexp_byclust_gathered <- rep1_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep1")
rep2 <- SubsetData(object = Vmat.combined,subset.name = "Replicate",accept.value = "FemaleRep2",subset.raw = T)
rep2_avgexp_byclust <- AverageExpression(object = rep2,use.raw = F,show.progress = T)
rep2_avgexp_byclust_gathered <- rep2_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep2")
rep3 <- SubsetData(object = Vmat.combined,subset.name = "Replicate",accept.value = "MaleRep1",subset.raw = T)
rep3_avgexp_byclust <- AverageExpression(object = rep3,use.raw = F,show.progress = T)
rep3_avgexp_byclust_gathered <- rep3_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep3")
rep4 <- SubsetData(object = Vmat.combined,subset.name = "Replicate",accept.value = "MaleRep2",subset.raw = T)
rep4_avgexp_byclust <- AverageExpression(object = rep4,use.raw = F,show.progress = T)
rep4_avgexp_byclust_gathered <- rep4_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep4")


all_avgexp_byclust <- bind_rows(rep1_avgexp_byclust_gathered,
                                rep2_avgexp_byclust_gathered,
                                rep3_avgexp_byclust_gathered,
                                rep4_avgexp_byclust_gathered)

mydata <- tibble(comparison = character(),
                 r = numeric(),
                 cluster = character())
for (clust_num in unique(all_avgexp_byclust$cluster)) {
    sub_data <- all_avgexp_byclust %>%
        filter(cluster == clust_num) %>%
        spread(replicate, expression) %>%
        column_to_rownames("gene") %>%
        select(-cluster) %>%
        correlate() %>%
        shave() %>%
        stretch() %>%
        unite(col = "comparison",x:y,sep = ":") %>%
        drop_na() %>%
        add_column(cluster = clust_num)
    mydata <- bind_rows(mydata,sub_data)
}
mydata$cluster <- as.numeric(mydata$cluster)
mydata_spread <- mydata %>%
    spread(comparison,r) %>%
    column_to_rownames("cluster") %>%
    t()

pheatmap(mydata_spread,
         breaks=seq(-1,1,0.05),
         color = diverging_hcl(40, "Blue-Red 3"),
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 16,
         angle_col = 90)
















##################################################
##    figure 8    ################################
##################################################




# Figure 8A
All.NPs <- read.table(file = "geneLists/neuropeptides.tsv",header = F,stringsAsFactors = FALSE)

max.umi <- apply(VNC.combined@raw.data, 1, function(x) max(x))
max.umi <- as.data.frame(max.umi)
max.umi <- rownames_to_column(max.umi)
exp.genes.raw.EgtX <- max.umi$rowname[max.umi$max.umi>10]
exp.NPs.raw.EgtX <- intersect(exp.genes.raw.EgtX,All.NPs$V1)

NPExpByCluster <- AverageExpression(object = VNC.combined,genes.use = exp.NPs.raw.EgtX, use.scale = F,use.raw = F)

colMax <- function(data) sapply(data, max, na.rm = TRUE)
subNPExpByCluster <- NPExpByCluster[, -which(colMax(NPExpByCluster) < 10)]
subNPExpByCluster.log <- log1p(subNPExpByCluster)
highlimit <- 4
pheatmap(subNPExpByCluster.log,
         breaks=seq(0,highlimit,0.1),
         color = colorRampPalette(rev(brewer.pal(11,"Spectral")))(10*highlimit),
         border_color = NA,
         cluster_rows = T,
         cluster_cols = T,
         treeheight_row = 100,
         treeheight_col = 100,
         fontsize = 16,
         angle_col = 90)





# Figure 8B
seuTest <- VNC.combined
SubGenes <- tibble(cell = VNC.combined@data@Dimnames[[2]],
                   identity = VNC.combined@meta.data[["res.12"]])
AllCellsandID <- SubGenes %>% 
    mutate(c84 = if_else(identity == 84,"Yes","no"))
seuTest@meta.data$c84 <- AllCellsandID$c84
TSNEPlot(object = seuTest, group.by = "c84",no.legend = T,do.label = F, pt.size = 0.1,coord.fixed=T,no.axes = T,colors.use = c("lightgrey","black"))
rm(seuTest)





# Figure 8D
cluster84 <- SubsetData(object = VNC.inf.combined,ident.use = 88,subset.raw = T)
FeaturePlot(object = cluster84,
            nCol = 4,
            coord.fixed = T,
            no.axes = T,
            no.legend = T,
            pt.size = 0.4,
            cols.use = c("lightgrey","red"),
            features.plot = c("amon","svr","Pal2","Phm"))



# Figure 8E
highExpNP <- rownames(subNPExpByCluster.log)[subNPExpByCluster.log[18] > 3]
FeaturePlot(object = cluster84,
            min.cutoff = 3,
            nCol = length(highExpNP),
            coord.fixed = T,
            no.axes = T,
            no.legend = T,
            pt.size = 0.4,
            cols.use = c("lightgrey","red"),
            features.plot = highExpNP)




# Figure 8-Figure supplement 1A
median(VNC.combined@meta.data[["nUMI"]])
np.seu <- SubsetData(object = VNC.combined,ident.use = 84)
median(np.seu@meta.data[["nUMI"]])

VlnPlot(object = VNC.combined,group.by = "orig.ident",point.size.use = 0,x.lab.rot = T,y.log = F,use.raw = T,features.plot = "nUMI")
VlnPlot(object = VNC.combined,ident.include = 84,group.by = "orig.ident",point.size.use = 0,x.lab.rot = T,y.log = F,use.raw = T,features.plot = "nUMI")



# Figure 8-Figure supplement 1B
All.NPs <- read.table(file = "geneLists/neuropeptides.tsv",header = F,stringsAsFactors = FALSE)
exp.NPs <- intersect(exp_genes$gene,All.NPs$V1)

npGenes <- FetchData(object = VNC.combined,vars.all = exp.NPs)
npGenes <- gather(data = npGenes,key = "gene",value = "expression")
ggplot(data = SubGenes, aes(x = gene,y = expression)) +
    geom_boxplot(aes(fill = gene)) +
     scale_y_continuous(trans = 'log10') +
    theme(legend.title=element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle = 90,hjust = 0))



# Figure 8-Figure supplement 1C
neuro.10 <- CreateSeuratObject(raw.data = neuro.Seu@raw.data,min.cells = 10)
neuro.absexp <- data.frame(as.matrix(neuro.10@raw.data))
neuro.absexp <- rownames_to_column(neuro.absexp,"gene")
neuro.absexp <- gather(data = neuro.absexp,key = "cell",value = "expression",-gene)

neuro.absexp.avg.max.10 <- neuro.absexp %>%
    group_by(gene) %>%
    top_n(10, expression) %>%
    summarise(avg.top.exp = mean(expression),stdev = sd(expression))

neuro.absexp.avg.max.10 <- neuro.absexp.avg.max.10 %>% 
    filter(avg.top.exp>0) %>% 
    mutate(NP = ifelse(
        neuro.absexp.avg.max.10$gene %in% All.NPs$V1,
        1,
        0)
    )

avgtopabsexp.barplot <- ggplot(neuro.absexp.avg.max.10, aes(x=reorder(gene,-avg.top.exp), y=avg.top.exp, fill=NP))+
    geom_bar(width = 1, stat = "identity") +
    scale_y_continuous(trans = 'log1p') +
    theme(axis.text.x = element_blank())
avgtopabsexp.barplot



# Figure 8-Figure supplement 1D
neuro.Seu <- SubsetData(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116),subset.raw = T)
neuro.freqs <- scale(neuro.Seu@raw.data, center = FALSE, scale = colSums(neuro.Seu@raw.data))
neuro.freqs <- data.frame(neuro.freqs)
neuro.freqs <- rownames_to_column(neuro.freqs,"gene")
neuro.freqs <- gather(data = neuro.freqs,key = "cell",value = "frequency",-gene)

neuro.freqs.avg.max.10 <- neuro.freqs %>%
    group_by(gene) %>%
    top_n(10, frequency) %>%
    summarise(avg.top.freq = mean(frequency),stdev = sd(frequency))

neuro.freqs.avg.max.10 <- neuro.freqs.avg.max.10 %>% 
    mutate(NP = ifelse(
        neuro.freqs.avg.max.10$gene %in% All.NPs$V1,
        1,
        0)
    )

avgtopfreqs.barplot <- ggplot(neuro.freqs.avg.max.10, aes(x=reorder(gene,-avg.top.freq), y=avg.top.freq*100, fill=NP))+
    geom_bar(width = 1, stat = "identity") +
    scale_y_continuous(trans = 'log1p') +
    theme(axis.text.x = element_blank())
avgtopfreqs.barplot





# Table 1

neuro.Seu <- SubsetData(object = VNC.combined,ident.remove = c(18,23,24,55,70,98,99,105,106,116),subset.raw = T)
neuro.Seu2 <- CreateSeuratObject(raw.data = neuro.Seu@raw.data,min.cells = 10,is.expr = 1)
neuro.Seu2 <- NormalizeData(object = neuro.Seu2)
neuro.Seu2 <- ScaleData(object = neuro.Seu2)

gene_max_exp <- tibble(gene = neuro.Seu2@raw.data@Dimnames[[1]],
                       raw_exp = apply(neuro.Seu2@raw.data, 1, function(x) max(x)),
                       norm_exp = apply(neuro.Seu2@data, 1, function(x) max(x))
)
gene_max_exp <- gene_max_exp %>%
    filter(raw_exp > 0) %>% 
    mutate(raw_percentile = percent_rank(raw_exp) * 100,
           norm_percentile = percent_rank(norm_exp) * 100) %>% 
    arrange(desc(norm_percentile))
gene_max_exp <- rowid_to_column(gene_max_exp, "ranked_exp")

exp_genes <- gene_max_exp$gene
np_genes <- read_tsv("geneLists/neuropeptides.tsv",col_names = F)
exp_np <- intersect(np_genes$X1,exp_genes)

np_max_exp <- gene_max_exp %>%
    filter(gene %in% exp_np) %>% 
    filter(raw_exp > 5) %>% 
    arrange(desc(norm_percentile))

mean(gene_max_exp$raw_exp)
median(gene_max_exp$raw_exp)

write_csv(np_max_exp, "neuropeptide_expression_percentile.csv")

genes_99_percentile <- gene_max_exp %>%
    filter(norm_percentile > 99) %>% 
    arrange(desc(norm_percentile))










# Figure 8-Figure supplement 2
# Same as previous analysis, but with the high UMI cutoff set to Inf



# Figure 8-Figure supplement 3
VlnPlot(object = VNC.combined,features.plot = "Orcokinin",x.lab.rot = T,use.raw=T,y.log = T,point.size.use = 0)








##################################################
##    figure 9    ################################
##################################################



# Figure 9A
seuTest <- VNC.combined
SubGenes <- NULL
SubGenes <- tibble(cell = seuTest@data@Dimnames[[2]],
                   identity = seuTest@meta.data[["res.12"]])
AllCellsandID <- SubGenes %>% 
    mutate(glia = if_else((identity == 23 | identity == 24 | identity == 70 | identity == 98 | identity == 106),"Yes","no"))
seuTest@meta.data$glia <- AllCellsandID$glia
TSNEPlot(object = seuTest, group.by = "glia",no.legend = T,do.label = F, pt.size = 0.1,coord.fixed=T,no.axes = T,colors.use = c("lightgrey","black"))
rm(seuTest)



# Figure 9B
male.glia.cells <- str_remove(string = WhichCells(object = SubsetData(object = VNC.combined,subset.name = "Sex",accept.value = "Male"),ident = c(23,24,70,98,106)),pattern = "Male_")
female.glia.cells <- str_remove(string = WhichCells(object = SubsetData(object = VNC.combined,subset.name = "Sex",accept.value = "Female"),ident = c(23,24,70,98,106)),pattern = "Female_")

elav.cells <- gsub(pattern = ".*ale_",replacement = "",x = WhichCells(object = SubsetData(object = VNC.combined,ident.use = c(23,24,70,98,106)),subset.name = "elav",accept.low = 0))
nSyb.cells <- gsub(pattern = ".*ale_",replacement = "",x = WhichCells(object = SubsetData(object = VNC.combined,ident.use = c(23,24,70,98,106)),subset.name = "nSyb",accept.low = 0))
VAChT.cells <- gsub(pattern = ".*ale_",replacement = "",x = WhichCells(object = SubsetData(object = VNC.combined,ident.use = c(23,24,70,98,106)),subset.name = "VAChT",accept.low = 0))
VGlut.cells <- gsub(pattern = ".*ale_",replacement = "",x = WhichCells(object = SubsetData(object = VNC.combined,ident.use = c(23,24,70,98,106)),subset.name = "VGlut",accept.low = 0))
Gad1.cells <- gsub(pattern = ".*ale_",replacement = "",x = WhichCells(object = SubsetData(object = VNC.combined,ident.use = c(23,24,70,98,106)),subset.name = "Gad1",accept.low = 0))


neuronal.cells <- union(VAChT.cells,union(VGlut.cells,union(Gad1.cells,union(elav.cells,nSyb.cells))))
male.glia.cells <- setdiff(male.glia.cells,neuronal.cells)
female.glia.cells <- setdiff(female.glia.cells,neuronal.cells)

Male.glia <- SubsetData(object = MaleVNC.combined,cells.use = male.glia.cells)
Female.glia <- SubsetData(object = FemaleVNC.combined,cells.use = female.glia.cells)
Male.glia <- FindVariableGenes(Male.glia, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001,do.plot = F)
Female.glia <- FindVariableGenes(Female.glia, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001,do.plot = F)

g.1 <- (Male.glia@var.genes)
g.2 <- (Female.glia@var.genes)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(Male.glia@scale.data))
genes.use <- intersect(genes.use, rownames(Female.glia@scale.data))
length(genes.use)

NumberOfCCAs <- 40
Glia.combined <- RunCCA(Male.glia, Female.glia, genes.use = genes.use, num.cc = NumberOfCCAs, add.cell.id1 = "Male", add.cell.id2 = "Female")
Glia.combined <- AlignSubspace(Glia.combined, reduction.type = "cca", grouping.var = "Sex", dims.align = 1:NumberOfCCAs, do.par = TRUE, num.cores = 8)

NumberOfCCAs <- 6
Glia.combined <- RunTSNE(Glia.combined,
                         num_threads = 0,
                         verbose = T,
                         reduction.use = "cca.aligned",
                        dims.use = 1:NumberOfCCAs,
                        do.fast = FALSE,
                        check_duplicates = F,
                        max_iter = 10000,
                        theta = 0.01)
Glia.combined <- FindClusters(Glia.combined,
                             reduction.type = "cca.aligned",
                             resolution = 0.9,
                             dims.use = 1:NumberOfCCAs,
                             save.SNN = TRUE)
TSNEPlot(Glia.combined,
               do.return = T, 
               pt.size = 4,
               no.axes = T,
               do.label = TRUE, 
               coord.fixed=T)

# Figure 9C
FeaturePlot(object = Glia.combined,
            no.axes = T,
            nCol = 2,
            coord.fixed = T,
            pt.size = 4,
            no.legend = T,
            cols.use = c("lightgrey","black"),
            features.plot = c("alrm","Indy","wrapper","Eaat2"))



# Figure 9D
glial.markers <- FindAllMarkers(object = Glia.combined,
                                test.use = "negbinom",
                                logfc.threshold = 0.5,
                                only.pos = T,
                                latent.vars = c("Replicate"),
                                do.print = T,
                                print.bar = T)
glial.markers <- glial.markers %>% filter(p_val_adj<0.05)
top10 <- glial.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DotPlot(object = Glia.combined,genes.plot = top10$gene,x.lab.rot = T,plot.legend = T,col.min = 0,col.max = 2,cols.use = c("white","blue"))





# Figure 9-Figure supplement 1A
TSNEPlot(Glia.combined,
         group.by ="Replicate",
         no.axes = T,
         pt.size = 8,
         do.label = F,
         coord.fixed=T)



# Figure 9-Figure supplement 1B
rep1 <- SubsetData(object = Glia.combined,subset.name = "Replicate",accept.value = "FemaleRep1",subset.raw = T)
rep1_avgexp_byclust <- AverageExpression(object = rep1,use.raw = F,show.progress = T)
rep1_avgexp_byclust_gathered <- rep1_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep1")
rep2 <- SubsetData(object = Glia.combined,subset.name = "Replicate",accept.value = "FemaleRep2",subset.raw = T)
rep2_avgexp_byclust <- AverageExpression(object = rep2,use.raw = F,show.progress = T)
rep2_avgexp_byclust_gathered <- rep2_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep2")
rep3 <- SubsetData(object = Glia.combined,subset.name = "Replicate",accept.value = "MaleRep1",subset.raw = T)
rep3_avgexp_byclust <- AverageExpression(object = rep3,use.raw = F,show.progress = T)
rep3_avgexp_byclust_gathered <- rep3_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep3")
rep4 <- SubsetData(object = Glia.combined,subset.name = "Replicate",accept.value = "MaleRep2",subset.raw = T)
rep4_avgexp_byclust <- AverageExpression(object = rep4,use.raw = F,show.progress = T)
rep4_avgexp_byclust_gathered <- rep4_avgexp_byclust %>%
    rownames_to_column("gene") %>%
    gather("cluster", "expression", -gene) %>%
    add_column(replicate = "rep4")


all_avgexp_byclust <- bind_rows(rep1_avgexp_byclust_gathered,
                                rep2_avgexp_byclust_gathered,
                                rep3_avgexp_byclust_gathered,
                                rep4_avgexp_byclust_gathered)

mydata <- tibble(comparison = character(),
                 r = numeric(),
                 cluster = character())
for (clust_num in unique(all_avgexp_byclust$cluster)) {
    sub_data <- all_avgexp_byclust %>%
        filter(cluster == clust_num) %>%
        spread(replicate, expression) %>%
        column_to_rownames("gene") %>%
        select(-cluster) %>%
        correlate() %>%
        shave() %>%
        stretch() %>%
        unite(col = "comparison",x:y,sep = ":") %>%
        drop_na() %>%
        add_column(cluster = clust_num)
    mydata <- bind_rows(mydata,sub_data)
}
mydata$cluster <- as.numeric(mydata$cluster)
mydata_spread <- mydata %>%
    spread(comparison,r) %>%
    column_to_rownames("cluster") %>%
    t()

pheatmap(mydata_spread,
         breaks=seq(-1,1,0.05),
         color = diverging_hcl(40, "Blue-Red 3"),
         border_color = NA,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 16,
         angle_col = 90)


