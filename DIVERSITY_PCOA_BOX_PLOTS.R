## OTU PCOA PLOT ##

OTU_pcoa = pcoa(bray_curtis_CSS_normalized_otu_table)

otu_vectors = as.data.frame(OTU_pcoa$vectors[,1:2])
otu_vectors$site = c(rep("CL" , 17) , rep("MTQ" , 8) , rep("KM" , 9))

site_colors = c("#FFC000" , "#00B0F0" , "#92D050")

pdf("~/Desktop/otu_pcoa.pdf")
ggplot(otu_vectors, aes( x = Axis.1 , y = Axis.2 , fill = site)) +
        geom_point(size = 5 , shape = 21 , color = "gray41") + scale_fill_manual(values = site_colors) + theme_bw() +
        theme( axis.title=element_text(size = 30) , axis.text = element_text(size = 15)  , aspect.ratio = 1 , legend.position = "none" )  +
        xlab("Variance explained = 10%") + ylab("Variance explained = 47%")
dev.off()

## OTU OBSERVED RICHNESS BOX PLOT ##
pdf("~/Desktop/otu_boxplot.pdf")
ggplot(Obs_richness_site , aes(x = Site , y = Richness , fill = Site)) +
        geom_boxplot(color = "gray41") + scale_fill_manual(values = site_colors) + 
        scale_x_discrete(limits=c("CL" , "MTQ" , "KM")) + theme_bw() +
        theme( axis.title=element_text(size = 30) , axis.text = element_text(size = 15)  , aspect.ratio = 1  , legend.position = "none") +
        xlab("") + ylab("Observed richness")
dev.off()
## METAGENOME PCOA PLOT ##

gypsum_dist = vegdist(t(gypsum.phylum.summary[,2:ncol(gypsum.phylum.summary)]) , method = "bray")
gypsum.phylum.pcoa = pcoa(gypsum_dist)
phylum_vectors =  as.data.frame(gypsum.phylum.pcoa$vectors[,1:2])
phylum_vectors$site = c(rep("CL" , 3) , rep("KM" , 3) , rep("MTQ" , 3))

pdf("~/Desktop/kaiju_phylum_pcoa.pdf")

ggplot(phylum_vectors, aes( x = Axis.1 , y = Axis.2 , fill = site)) +
        geom_point(size = 5 , shape = 21 , color = "gray41") + scale_fill_manual(values = site_colors) + theme_bw() +
        theme( axis.title=element_text(size = 30) , axis.text = element_text(size = 15) , aspect.ratio = 1 , legend.position = "none" )  +
        xlab("Variance explained = 59%") + ylab("Variance explained = 30%")
dev.off()

## FUNCTIONS OBSERVED RICHNESS BOXPLOT##

KO_richness = data.frame(community = c("CL1" , "CL2" , "CL3" , "KM1" , "KM2" , "KM3" , "MTQ1" , "MTQ2" , "MTQ3") ,
                         site = c(rep("CL" , 3) , rep("KM" , 3) , rep("MTQ" , 3)) , 
                                richness = c(4505,4566,4689,4318,4210,4504,4495,4508,4695) )
pdf("~/Desktop/functions_boxplot.pdf")


ggplot(KO_richness , aes(x = site , y = richness , fill = site)) + 
        geom_boxplot(color = "gray41") + scale_fill_manual(values = site_colors) + 
        scale_x_discrete(limits=c("CL" , "MTQ" , "KM")) + theme_bw() +
        theme( axis.title=element_text(size=30) , axis.text = element_text(size = 15) , aspect.ratio = 1 , legend.position = "none") +
        xlab("") + ylab("Observed richness")

dev.off()
## FUNCTIONS PCOA PLOT ##

functions_dist = vegdist(t(gypsum_kegg_sums[,2:ncol(gypsum_kegg_sums)]) , method = "bray")
functions_pcoa = pcoa(functions_dist)
functions_vectors =  as.data.frame(functions_pcoa$vectors[,1:2])
functions_vectors$site = c(rep("CL" , 3) , rep("KM" , 3) , rep("MTQ" , 3))

pdf("~/Desktop/functions_pcoa.pdf")
ggplot(functions_vectors, aes( x = Axis.1 , y = Axis.2 , fill = site)) +
        geom_point(size = 5 , shape = 21 , color = "gray41") + scale_fill_manual(values = site_colors) + theme_bw() +
        theme( axis.title=element_text(size=30) , axis.text = element_text(size = 15) ,aspect.ratio = 1 , legend.position = "none" )  +
        xlab("Variance explained = 67%") + ylab("Variance explained = 21%")

dev.off()



