# modified DotPlot function from Seurat v2.3.4


myDotPlot <- function (object, genes.plot, cols.use = c("blue","white","red"), 
          col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
          scale.by = "radius", scale.min = NA, scale.max = NA, group.by, 
          plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE) 
{
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
                         radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    if (!missing(x = group.by)) {
        object <- SetAllIdent(object = object, id = group.by)
    }
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    colnames(x = data.to.plot) <- genes.plot
    data.to.plot$cell <- rownames(x = data.to.plot)
    data.to.plot$id <- object@ident
    data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
                                            value = expression, -c(cell, id))
    data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
        summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, 
                                                                                threshold = 0))
    data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
        mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, 
                                                                                     max = col.max, min = col.min))
    data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                      levels = rev(x = genes.plot))
    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
    p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                   y = id)) + geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) + 
        scale.func(range = c(0, dot.scale), limits = c(scale.min, 
                                                       scale.max)) + theme(axis.title.x = element_blank(), 
                                                                           axis.title.y = element_blank())
    if (length(x = cols.use) == 1) {
        p <- p + scale_color_distiller(palette = cols.use)
    }
    else {

        # modified this next line by adding the mid arguement
        p <- p + scale_color_gradient2(low = cols.use[1], mid = cols.use[2], high = cols.use[3])
        #

    }
    if (!plot.legend) {
        p <- p + theme(legend.position = "none")
    }
    if (x.lab.rot) {
        p <- p + theme(axis.text.x = element_text(angle = 90, 
                                                  vjust = 0.5))
    }
    suppressWarnings(print(p))
    if (do.return) {
        return(p)
    }
}

    
    
PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
}