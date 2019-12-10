#!/usr/bin/env Rscript

# Aaron M. Allen, 2019.08.23


scRNA_bulkFlyAtlas2_comparison <- function(seuObject = VNS.combined,tissue1 = 7, tissue2 = 8, use.raw = TRUE) {

    flyatlas.tissue <- genefpkm %>% filter(TissueID == tissue1 | TissueID == tissue2)
    fa.tissue.sym <- left_join(x = flyatlas.tissue[,1:3],y = FBgn_to_symbol,by="submitted_id")
    fa.tissue.sym2 <- fa.tissue.sym %>% group_by(current_symbol) %>% summarize(FPKM=mean(FPKM))
    
    scGene_to_newGene <- read_delim("geneLists/scGene_to_newGene.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)           # read the file in
    scGene_to_newGene <- scGene_to_newGene[- grep("D*{3}[\\]", scGene_to_newGene$current_symbol),]                  # flybase spits out some genes from other species
    scGene_to_newGene <- scGene_to_newGene[! scGene_to_newGene$current_id=="unknown ID", ]                          # remove the rows that flybase couldn't find current symbols for
    scGene_to_newGene <- drop_na(scGene_to_newGene)                                                                 # ...
    scGene_to_newGene <- scGene_to_newGene[! (scGene_to_newGene$current_symbol!=scGene_to_newGene$submitted_id 
                                              & duplicated(scGene_to_newGene$current_symbol)),]                     # remove rows where the input id had two ouput ids and the input and output don't match
    
    if (use.raw == TRUE) {
           sc.data <- data.frame(submitted_id=seuObject@raw.data@Dimnames[[1]],UMI=rowSums(seuObject@raw.data))
           sc.data <- left_join(x = sc.data,y = scGene_to_newGene,by = "submitted_id")
           sc.data <- drop_na(sc.data)
           sc.fa2.join <- full_join(x = fa.tissue.sym2,y = sc.data,by = "current_symbol")
           sc.fa2.join <- drop_na(sc.fa2.join)
           sc.fa2.join <- sc.fa2.join %>% filter(FPKM>0) %>% filter(UMI>0)
           sc.fa2.join$FPKM <- log(sc.fa2.join$FPKM)
           sc.fa2.join$UMI <- log(sc.fa2.join$UMI)
           ggscatter(sc.fa2.join, x = "FPKM", y = "UMI", 
                     add = "reg.line", conf.int = TRUE, #conf.int.level = 95,
                     cor.coef = TRUE, cor.method = "spearman",
                     xlab = "log(FPKM)", ylab = "log(nUMI)",
                     alpha = I(1 / 20),
                     ggtheme = theme_gray())
           
    } else {
           sc.data <- data.frame(submitted_id=VNS.combined@data@Dimnames[[1]],norm=rowSums(VNS.combined@data))
           sc.data <- left_join(x = sc.data,y = scGene_to_newGene,by = "submitted_id")
           sc.data <- drop_na(sc.data)
           sc.fa2.join <- full_join(x = fa.tissue.sym2,y = sc.data,by = "current_symbol")
           sc.fa2.join <- drop_na(sc.fa2.join)
           sc.fa2.join <- sc.fa2.join %>% filter(FPKM>0) %>% filter(norm>0)
           sc.fa2.join$FPKM <- log(sc.fa2.join$FPKM)
           sc.fa2.join$norm <- log(sc.fa2.join$norm)
           ggscatter(sc.fa2.join, x = "FPKM", y = "norm", 
                     add = "reg.line", conf.int = TRUE, #conf.int.level = 95,
                     cor.coef = TRUE, cor.method = "spearman",
                     xlab = "log(FPKM)", ylab = "log(norm.expression)",
                     alpha = I(1 / 20),
                     ggtheme = theme_gray())
           
    }
}


# 
# scRNA_bulkFlyAtlas2_comparison()
# scRNA_bulkFlyAtlas2_comparison(use.raw = F)
# scRNA_bulkFlyAtlas2_comparison(tissue1 = 5, tissue2 = 6,use.raw = FALSE)
# scRNA_bulkFlyAtlas2_comparison(tissue1 = 19, tissue2 = 20)
# 


