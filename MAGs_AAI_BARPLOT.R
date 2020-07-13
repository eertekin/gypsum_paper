## MAGs AAI BARPLOT##
library(ggplot2)
pdf("MAGs_AAI.pdf")
ggplot(MAGs_phylogenetic_information , aes(x = MAG , y = AAI)) +
        geom_bar(stat = "identity" , color = "gray21" , fill = "gray41" ) + theme_bw() + 
        theme(aspect.ratio = 0.5 , text = element_text(size = 15) , axis.text.x = element_text(angle = 90) , 
              legend.text = element_text(face = "italic") , legend.position = "top" ,
              legend.title = element_blank()) + scale_y_continuous( limits = c(0, 100) ) +
        labs(y = "AAI(%)" , x = "")

dev.off()
