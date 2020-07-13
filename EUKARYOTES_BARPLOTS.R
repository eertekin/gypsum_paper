## EUKARYOTES BARPLOTS ##
        #FUNGI BARPLOT#
library(stringr) ; library(ggplot2) ; library(tidyr) ; library(RColorBrewer)

gypsum_fungi = gypsum.phylum.summary[c(15,3),]
colnames(gypsum_fungi) = c("phylum" , "CL1" ,"CL2","CL3","KM1","KM2","KM3","MTQ1","MTQ2","MTQ3")
gypsum_fungi_tidy = gather( gypsum_fungi , key = "site" , value = "ra" , 2:10 )
gypsum_fungi_tidy$site = factor(gypsum_fungi_tidy$site , levels = c("CL1" ,"CL2","CL3", "MTQ1","MTQ2","MTQ3" , "KM1","KM2","KM3"))

fungi_colors = c("indianred1" ,"orange" )

pdf("fungi_barplot.pdf")
p1 = ggplot(gypsum_fungi_tidy , aes(x = site , y = ra , fill = phylum)) +
        geom_bar(stat = "identity") + theme_bw() + 
        theme(aspect.ratio = 1 , text = element_text(size = 15) ,
              axis.title.y = element_text(size = 15) , axis.text.x = element_text(angle = 90) ,
              legend.text = element_text(face = "italic") , legend.position = "top" ,
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in") ,
              legend.title = element_blank()) + scale_fill_manual(values = fungi_colors) + 
        labs(y = "Relative abundance(%)" , x = "")
dev.off()
        #ALGAE BARPLOT#
        
gypsum_algae = chlorophyta_genus[1:2,]
gypsum_algae = gypsum_algae[order(gypsum_algae$KM1 , decreasing = TRUE),]

gypsum_algae_tidy = gather(gypsum_algae , key = site , value = ra , 2:10)
 
gypsum_algae_tidy$site = factor(gypsum_algae_tidy$site , levels = c("CL1" ,"CL2","CL3", "MTQ1","MTQ2","MTQ3" , "KM1","KM2","KM3"))

algae_colors = c("seagreen" , "turquoise" , "cornflowerblue" , "deepskyblue" , "cadetblue1" , "gray")

set.cwd = "~/Desktop"
ggplot(gypsum_algae_tidy , aes(x = site , y = ra , fill = taxa)) +
        geom_bar(stat = "identity") + theme_bw() +
        theme(aspect.ratio = 1 , text = element_text(size = 15) , 
              axis.title.y = element_text(size = 15) , axis.text.x = element_text(angle = 90) , 
              legend.text = element_text(face = "italic") , legend.position = "top" , 
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in") ,
              legend.title = element_blank()) + scale_fill_manual(values = algae_colors) + 
        labs(y = "Relative abundance(%)" , x = "")
ggsave("algae_barplot.pdf" , width = 5 , height = 5)

graphics.off()


p1 ; ggsave("fungi_barplot.pdf" , width = 5 , height = 5)
 












