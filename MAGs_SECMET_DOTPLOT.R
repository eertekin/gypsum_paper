## MAGs SECMET DOTPLOT ##

dev.new(width=12, height=10, unit="in")
pdf("secmet_dotplot_3.pdf")

secmet_dotplot = ggplot(secmet.sum.lineage  , aes( x = lineage ,  y = number.mbp , fill = site)) + 
        geom_dotplot( binaxis= 'y', binwidth = 0.15 , binpositions = "all" , stackratio= 1 ,
                      stackdir = "center", stackgroups = TRUE) + 
        theme_bw() +
        scale_fill_manual(values = alpha(c("gray" , "#00B0F0" , "#92D050"))) + 
        theme( text = element_text(size = 15) , 
               plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in") , legend.title = element_blank() , axis.text.x = element_text(angle = 90 , face = "italic" , hjust = 1) ) + 
        facet_grid( . ~ lineage , scales = "free" , space = "free") +
        theme(strip.text.x = element_blank()) + labs(y = "Number of clusters per assembled genome" , x = "") + 
        ylab(expression(paste("Number of clusters per \n assembled genome")))


ggsave("secmet_dotplot.pdf" , width = 5 , height = 5)
