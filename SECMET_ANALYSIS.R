## SECONDARY METABOLITE CLUSTERS ANALYSIS ##
library(pheatmap) ; library(dplyr) ; library(tidyr) ; library(stringr) ; library(RColorBrewer)
       

colnames(all_secmet_clusters) = c("contig_img" , "cluster")
all_secmet_clusters = merge(all_secmet_clusters , gene_quant[,c(2,4:12)] , by = "contig_img")
all_secmet_clusters = unique(all_secmet_clusters)


secmet_collapse = as.data.frame(
        apply(all_secmet_clusters[,3:11] , 2 , function(x){
                as.data.frame(aggregate(x ~ cluster , data = all_secmet_clusters , FUN = sum))
                })) 
secmet_collapse = secmet_collapse[,c(1,cols2select)]
colnames(secmet_collapse) = c("cluster" , "CL1" , "CL2" , "CL3" , "KM1" , "KM2" , "KM3" , "MTQ1" , "MTQ2" , "MTQ3" )

secmet = gather(secmet_collapse , key = "sample" , value = "abundance" , 2:10 )
colnames(secmet) = c("cluster" , "sample" , "abundance")


for( i in 1:nrow(secmet) ) {
        if( str_detect(secmet[,2][i] ,  "CL")) {
                secmet$community[i] = "CL"} else if (str_detect(secmet[,2][i] ,  "KM")) {
                        secmet$community[i] = "KM"} else if (str_detect(secmet[,2][i] ,  "MTQ")) {secmet$community[i] = "MTQ"}          
}

secmet = gather(secmet_collapse , key = "sample" , value = "abundance" , 2:10 )
colnames(secmet) = c("cluster" , "sample" , "abundance")

# GENERATING THE HEATMAP 
annotation_column = data.frame( site = factor(sites)) ; row.names(annotation_column) = c("CL1" , "CL2" , "CL3" , "KM1" , "KM2" , "KM3" , "MTQ1" , "MTQ2" , "MTQ3" ) 

secmet_colors = rev(colorRampPalette(brewer.pal(n = 10 , name = "Spectral"))(256))
secmet_data = log10(secmet_collapse[,2:10] + 1) ; row.names(secmet_data) = secmet_collapse$cluster

par(mar = c(0,0,0,0) , oma = c(0,0,0,0))

pheatmap(secmet_data , color = secmet_colors , border_color = NA , 
         annotation_col = annotation_column , annotation_colors = mat_colors ,
         clustering_method = "average")

# BINS SECMET ANALYSIS BARPLOT #

secmet.sum.bin$lineage = str_to_title(secmet.sum.bin$lineage , locale = "en" )
secmet.sum.bin.simple = secmet.sum.bin
secmet_categories = c("terpene" , "NRPS" , "NRPS-like" , "T1PKS" , "T3PKS" , "bacteriocin" , "betalactone" , "hserlactone" , "lassopeptide" , "siderophore" , "other")

for (i in 1:nrow(secmet.sum.bin.simple) ) { 
        if ( ! (secmet.sum.bin.simple[i,2] %in% secmet_categories )) { 
                secmet.sum.bin.simple[i,2] = "other" }
        }
secmet.sum.bin.simple$cluster = factor(secmet.sum.bin.simple$cluster , levels = secmet_categories)

secmet.sum.bin.simple$lineage = factor(secmet.sum.bin.simple$lineage , levels = c("Cyanobacteria" , "Actinobacteria" , "Proteobacteria" , "Chloroflexi" , "Deinococcus" , "Gemmatimonidetes" )) 
secmet_colors = c("forestgreen" ,
                  rev(brewer.pal( n = 4 , name = "Reds")) ,
                  rev(brewer.pal(n = 4 , name = "Blues")) ,
                  "mediumorchid1" , "gray41" )
pdf("secmet_clusters_bins.pdf")

ggplot(secmet.sum.bin.simple , aes(x = lineage , y = number , fill = cluster)) +
        geom_bar(stat = "identity" , width = 0.5 ) + coord_flip() + theme_bw() + scale_fill_manual( values = secmet_colors ) +
        theme(axis.text = element_text(size = 15) , axis.text.y = element_text(face = "italic"),
              aspect.ratio = 0.75 , legend.position = "bottom" , legend.text = element_text(size = 12))
ggsave("secmet_clusters_bins.pdf" , height = 5 , width = 8)







