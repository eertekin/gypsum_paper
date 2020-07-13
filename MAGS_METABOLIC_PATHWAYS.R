library(dplyr) ; library(ggplot2) ; library(RColorBrewer) ; library(stringr) ; library(wesanderson) ; library(tidyr) ; library(cowplot) ; library(pheatmap)
library(viridis)
## MAGs HEATMAP ##

abundance_data = as.data.frame(log10(abundance_table[,2:10] + 1))
rownames(abundance_data) = abundance_table$MAG

annotation_row_mags = data.frame(
        MAG = factor( abundance_table$lineage) )
annotation_column = data.frame(
        site = factor( c("CL","CL","CL","KM","KM","KM","MTQ","MTQ","MTQ")) )
annotation_column$site = factor(annotation_column$site , levels = c("CL" , "KM" , "MTQ"))

row.names(annotation_column) = c("CL1","CL2","CL3","KM1","KM2","KM3","MTQ1","MTQ2","MTQ3")

row.names(annotation_row_mags) = row.names(abundance_data)
rownames(annotation_column) = c("CL1","CL2","CL3","KM1","KM2","KM3","MTQ1","MTQ2","MTQ3")
colnames(abundance_data) = c("KM1" , "CL3"  , "KM3" ,  "CL1" , "KM2" , "CL2" ,"MTQ2","MTQ3" ,"MTQ1")

phylum_colors = list( 
        MAG = c(cyanobacteria = "#92D050" , actinobacteria = "#FFC100" , deinococcus = "lightseagreen" ,
                gemmatimonidetes = "indianred1" , proteobacteria = "#7DBEFF" ,
                chloroflexi = "orchid" ) ,
        site = c(CL = "gold" , KM = "skyblue" , MTQ = "forestgreen")
        )

 pheatmap( abundance_data , border_color = NA , 
          color = colorRampPalette(viridis(n = 8 , option = 1))(250)  , 
          annotation_col = annotation_column , annotation_row = annotation_row_mags  , 
          annotation_colors = phylum_colors )

## MAGs METABOLIC PATHWAYS ##

 for (i in 1:nrow(bins_lineage)) {
         if(bins_lineage$bin[i] %in% chloroflexi) { bins_lineage$lineage[i] = "chloroflexi"} else { bins_lineage$lineage[i] = bins_lineage$lineage[i]}
 }  
 
 
nitrogen_met = c( c("K00370" , "K00371" , "K00372" , "K00374" , "K02567" , "K02568" , "K00367" , "K10534" , "K00360", "K00362" , "K00363") , c("K17877" , "K00366") , c("K00368" , "K15864" , "K04561" , "K02305") )
sulfur_met = c( "K00392" , "K00381" ,"K00380" , "K00956" , "K00955" , "K00390" , "K00860" , "K00958" , "K00957" , "K17227" , "K17226" , "K17218" , "K00387" )
anox_ps = c( "K11333" ,  "K04641" , "K11334" , "K11335" , "K08926" , "K08928" , "K08929" )
ox_ps = c("K02689" , "K02690" , "K02703" , "K02704" )
c_fix = c("K01601" , "K01602" , "K00855" , "K15231" , "K01965" , "K01966")
        
#pathways = c( c("K01601" , "K01602" , "K00855" , "K15231" , "K01965" , "K01966") , c("K02689" , "K02690" , "K02703" , "K02704" ) , c("K00362" , "K00363") , c( "K11333" , "K11334" , "K11335" , "K08926" , "K08928" , "K08929" ) ,
#              c("K00370" , "K00371" , "K00374" , "K02567" , "K02568" , "K00367" , "K10534" , "K00372" , "K00360") , c("K17877" , "K00366") , c("K00368" , "K15864" , "K04561" , "K02305") ,
#              c("K13811" , "K00958" , "K00955" , "K00956" , "K00957") , c("K13811" , "00955" , "K00860" , "K00390" ) , c("K00380" , "K00381" , "K00392" ) , c("K17226" , "K17227") , c("K00387" , "K17218" ) ,
#              c( "K11333" ,  "K04641" , "K11334" , "K11335" , "K08926" , "K08928" , "K08929" ) )

pathways = c(nitrogen_met , sulfur_met , anox_ps , ox_ps , c_fix)

MAGs_energy_carbon = subset(bins_kegg_count , bins_kegg_count$kegg_no %in% pathways)
MAGs_energy_carbon = merge(MAGs_energy_carbon  , bins_lineage , by = "bin" , all.x = TRUE)
level3 = data.frame( level3 = "text" , stringsAsFactors = FALSE)
d = MAGs_energy_carbon
d = cbind(d , level3)

for (i in 1:nrow(d) ) {
        if( d$kegg_no[i] %in% anox_ps) {
                d$level3[i] = "anoxygenic photosynthesis"
        } else if ( d$kegg_no[i] %in% ox_ps ) { 
                d$level3[i] = "oxygenic photosynthesis"
        } else if ( d$kegg_no[i] %in% c_fix ) { 
                d$level3[i] = "carbon fixation"
        } else if (d$kegg_no[i] %in% nitrogen_met) {
                d$level3[i] = "nitrogen metabolism"
        } else if(d$kegg_no[i] %in% sulfur_met) {
                d$level3[i] = "sulfur metabolism"
        }
}

a = subset(d, level3 %in% c("carbon fixation" , "anoxygenic photosynthesis" , "oxygenic photosynthesis" , "sulfur metabolism"  , "nitrogen metabolism" ))      
a$kegg_no = factor(a$kegg_no)
a$level3 = factor(a$level3)
a$bin = str_replace(a$bin , pattern = "bin." , replacement = "MAG")
a$level3 = factor(a$level3 , levels = c("oxygenic photosynthesis" , "anoxygenic photosynthesis" , "carbon fixation" , "sulfur metabolism" , "nitrogen metabolism"))
a[is.na(a)] = 0
a$enzyme = as.character(a$enzyme)

a$kegg_no = factor( a$kegg_no , levels = c(
        c("K00370" , "K00371" , "K00372" , "K00374" , "K00360" , "K02567" , "K02568" , "K00367" , "K10534"   , "K00362" , "K00363") , c("K17877" , "K00366") , c("K00368" , "K15864" , "K04561" , "K02305")  ,
         c( "K00392" , "K00381" ,"K00380" , "K00956" , "K00955" , "K00390" , "K00860" , "K00958" , "K00957" , "K17227" , "K17226" , "K17218" , "K00387" ) ,
       c( "K11333" ,  "K04641" , "K11334" , "K11335" , "K08926" , "K08928" , "K08929" ) ,
       c("K02689" , "K02690" , "K02703" , "K02704" ) ,
       c("K01601" , "K01602" , "K00855" , "K15231" , "K01965" , "K01966") )
        
        
        
)

a$name = paste( a$kegg_no, a$enzyme , sep = " | ")
a$name = factor(a$name , levels = c(rev(paste(kegg_gene_levels$kegg_no , kegg_gene_levels$gene , sep = " | " ) )))
        


a$bin = factor(a$bin ,  levels = c(paste("MAG" , c(35,17,16,13,4,28,40,29)  , sep = "") , paste("MAG" , c(1,25,19,33) , sep = "") , paste("MAG" , c(23,41,30,20,9,43,31,36,37,26,5,12,8,27,3,24,34,42,11,44) , sep = "") ,
                                   paste("MAG" , c(7,22,39) , sep = "") , paste("MAG" , c(10,14,2,21,32,38,18) , sep = "") , paste("MAG" , c(6,15) , sep = "") ))
a$number = as.numeric(a$number)
for ( i in 1:nrow(a) ) { if( a$number[i] > 0 ) { a$number[i] = 1} }
for ( i in 1:nrow(a) ) { if( a$enzyme[i] == "nirA" ) { a$kegg_no[i] = "K00366" } }

a$number = as.factor(a$number)
svg("mags_metabolic_pathways.svg")

ggplot(a, aes( x = bin , y = name, fill = number  )) +
        geom_tile(color = "gray41") + facet_grid(level3~lineage , scales = "free" , space = "free") + 
        theme_bw() + theme(axis.text.x = element_text(angle = 90) , text = element_text(size = 10)) + scale_fill_manual(values = c("gray" , "black") )

dev.off()








