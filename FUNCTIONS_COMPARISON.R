library(pheatmap) ; library(viridis) ; library(tidyr) ; library(dplyr)
library(RColorBrewer) ; library(viridis) ; library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(9, "Set1"))

## DIFFERENTIALLY ABUNDANT FUNCTIONS ##

        # ANOVA ANALYSIS #
gypsum_kegg_pathways = gene_quant[,4:14]
gypsum_kegg_pathways = as.data.frame(apply( gene_quant[,4:12] , 2 , function(x) {
        aggregate( x ~ gene_quant$kegg_no , data = gene_quant , FUN = sum )} )   )

gypsum_kegg_pathways = gypsum_kegg_pathways[,c(1,cols2select)]

gypsum_kegg_pathways_t = as.data.frame(t(gypsum_kegg_pathways[,-1]))
colnames(gypsum_kegg_pathways_t) = gypsum_kegg_pathways[,1]
gypsum_kegg_pathways_t$community = c("CL","CL","CL","KM","KM","KM" ,"MTQ","MTQ","MTQ")

ncol(gypsum_kegg_pathways_t)
keggs.anova = apply(gypsum_kegg_pathways_t[1:4791] , 2 , function(x){ 
        TukeyHSD(aov( x ~ gypsum_kegg_pathways_t$community , data = gypsum_kegg_pathways_t))
} )
anova.res = t(as.data.frame(lapply(keggs.anova , function(x){ as.data.frame(x[[1]][10:12]) })))
colnames(anova.res) = c("CL-KM" , "CL-MTQ" , "KM-MTQ") ;

p.vals = gather(anova.res[,c(1:3,7) ] , key = "comparison" , value = "pval" , 1:3 )
p.vals$p.adj = p.adjust(p.vals$pval , method = "fdr") 

CL.KM = p.vals %>% filter(comparison == "CL-KM" )  
CL.MTQ = p.vals %>% filter( comparison == "CL-MTQ" )
KM.MTQ = p.vals %>% filter( comparison == "KM-MTQ")
kegg.padj = cbind(CL.KM[,1:4] , CL.MTQ$p.adj , KM.MTQ$p.adj)
colnames(kegg.padj) = c("kegg_no" , "comparison" , "pval" , "CL.KM" , "CL.MTQ" , "KM.MTQ")

kegg_sums = unique(gypsum_img_kegg_sums[,c(1,4:12)])

kegg.padj.sf = kegg.padj %>% filter(kegg.padj$CL.KM < 0.05 & kegg.padj$KM.MTQ < 0.05 & kegg.padj$CL.MTQ > 0.05 )
kegg.padj.sf = merge(kegg.padj.sf , kegg_enzyme_desc , by = "kegg_no" )
kegg.padj.sf = merge(kegg.padj.sf , kegg_sums , by = "kegg_no")

## ENRICHED KOS HEATMAP ##

z.scores.kegg.sf = t(apply( keggs.data.signif[,2:10] , 1 , function(x) {
        (x - mean(x)) / sd(x)
} ))
community = c("CL1" , "CL2" ,"CL3", "KM1","KM2","KM3","MTQ1","MTQ2","MTQ3") 

kos_annotation_sites = data.frame( row.names = c("CL1" , "CL2" ,"CL3", "KM1","KM2","KM3","MTQ1","MTQ2","MTQ3") , 
                                   site = c(rep("CL",3) , rep("KM",3) , rep("MTQ",3)))
                                   
                                  
annotation_colors = list( 
        site = c(CL = "#FFC000" , KM = "#00B0F0", MTQ = "#92D050" ) )

par(oma = c(5,5,5,5))
svg("enriched_kos_heatmap.svg")
pheatmap( z.scores.kegg.sf , color = colorRampPalette(viridis(n = 10 , option = 2 ))(256) ,
          annotation_col = kos_annotation_sites , annotation_colors = annotation_colors)
dev.off()

## ENRICHED KOS PATHWAYS HEATMAP ##

kegg.sf.pathways = merge(kegg.padj.sf , kegg2brite , by = "kegg_no" , all.x = TRUE )
kegg.sf.pathway.sums = as.data.frame(apply( kegg.sf.pathways[,9:17] , 2 , function(x) {
        aggregate( x ~ kegg.sf.pathways$level3 , data = kegg.sf.pathways , FUN = sum )} )   )

kegg.sf.pathway.sums = kegg.sf.pathway.sums[,c(1,cols2select)]
colnames(kegg.sf.pathway.sums) = c("level3","CL1","CL2","CL3","KM1","KM2","KM3" , "MTQ1","MTQ2","MTQ3")

row.names(kegg.sf.pathway.sums) = 1:nrow(kegg.sf.pathway.sums)
kegg.sf.pathway.sums = kegg.sf.pathway.sums[-c(31,72),]

kegg.sf.pathway.sums = unique(merge(kegg.sf.pathway.sums , kegg2brite[,c(1,5)] , by = "level3" , all.x = TRUE))
kegg.sf.pathway.sums$SUMS = rowSums(kegg.sf.pathway.sums[,2:10])

kegg.sf.pathway.sums = kegg.sf.pathway.sums[order(kegg.sf.pathway.sums$SUMS , decreasing = TRUE),]
kegg.sf.path.sums.top30 = kegg.sf.pathway.sums[1:30,]

# GENERATING THE HEATMAP #

z.scores.path.sums = t(apply(kegg.sf.path.sums.top30[,2:10] , 1 , function(x) {
        (x - mean(x)) / sd(x)
} ) )

row.names(z.scores.path.sums) = kegg.sf.path.sums.top30$level3
annotation_row_pathsums = data.frame(
        level1 = factor(kegg.sf.path.sums.top30$level1) )  

row.names(annotation_row_pathsums) = row.names(z.scores.path.sums)

sites = c(rep("CL" , times = 3) , rep("KM" , times = 3) , rep("MTQ" , times = 3))

annotation_column_pathsums = data.frame( site = factor(sites) ) 
row.names(annotation_column_pathsums ) = c( "CL1","CL2","CL3","KM1","KM2","KM3","MTQ1","MTQ2","MTQ3")
l1_c = brewer.pal(3,"Set1")
annotation_colors = list(
        site = c(CL = "#FFC000" , KM = "#00B0F0", MTQ = "#92D050" ) ,
        level1 = c( "Cellular Processes" =  l1_c[1] , "Environmental Information Processing" =  l1_c[3] , "Metabolism" =  l1_c[2])
)


pheatmap( z.scores.path.sums , border_color = "NA" , 
          color = colorRampPalette(viridis(n=10 , option=1))(256) , fontsize = 10 ,
          annotation_row = annotation_row_pathsums  ,
          annotation_col = annotation_column_pathsums  ,
 annotation_colors = annotation_colors, width = 10 , height = 6 ,
 filename = "~/Desktop/enriched_pathways_heatmap.pdf" )










