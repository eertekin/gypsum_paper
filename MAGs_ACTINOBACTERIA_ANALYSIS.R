library(ggplot2) ; library(tidyr) ; library(dplyr) ; library(RColorBrewer)

## MAGs ACTINOBACTERIA ANALYSIS ##
        #GENES ENRICHED IN KM
MAG.actino.sf.KM = MAG.actino.pvals.sf %>% filter(KM > MTQ)
hyd.ox.kos = c("K04651","K04652","K04653","K04656","K06281","K06282","K07321")
MAG.actino.sf.KM.hyd = MAG.actino.sf.KM %>% filter(kegg_no %in% hyd.ox.kos )
MAG.actino.sf.KM.hyd.tidy = gather(MAG.actino.sf.KM.hyd , key = "MAG" , value = "number" , 5:22 )
MAG.actino.sf.KM.hyd.tidy = merge(MAG.actino.sf.KM.hyd.tidy , abundance_table[,c(1,14)], by = "MAG")
MAG.actino.sf.KM.hyd.tidy = MAG.actino.sf.KM.hyd.tidy[order(MAG.actino.sf.KM.hyd.tidy$site),]
setwd("~/Desktop")
pdf("MAG_actino_enriched_KM.pdf")
ggplot(MAG.actino.sf.KM.hyd.tidy , aes( x = MAG , y = enzyme , fill = number )) + facet_grid( .~site , scales = "free" , space = "free") + 
        geom_tile() + theme_bw() + theme(axis.text.x = element_text(angle = 90 , size = 10 ) , aspect.ratio = 1  , text = element_text( size = 15)) +
        scale_fill_manual(values = c("gray20","darkred","firebrick3", "firebrick1" , "orangered" , "tomato" , "orange"  , "gold"))

dev.off()
        #GENES ENRICHED IN MTQ
MAG.actino.sf.MTQ = MAG.actino.pvals.sf %>% filter(KM < MTQ)
signaling.kos = c("K02480","K00945","K03100","K02519")


aa.uptake.degradation.kos = c("K11069" , "K11070" , "K11072" , 
                              "K00663" , "K00274" , "K00130" ,
                              "K01619" , "K01811") 

MAG.actino.sf.MTQ.aa.tidy$enzyme = factor(MAG.actino.sf.MTQ.aa.tidy$enzyme , levels = rev(c("potC","potD","potA",
                                                                                           "MAO","betB","aacA","deoC","xylS")))
 

MAG.actino.sf.MTQ.aa = MAG.actino.sf.MTQ %>% filter(kegg_no %in% aa.uptake.degradation.kos)
MAG.actino.sf.MTQ.aa.tidy = gather(MAG.actino.sf.MTQ.aa , key = "MAG" , value = "number" , 5:22 )
MAG.actino.sf.MTQ.aa.tidy = merge(MAG.actino.sf.MTQ.aa.tidy , abundance_table[,c(1,14)], by = "MAG")
MAG.actino.sf.MTQ.aa.tidy = MAG.actino.sf.MTQ.aa.tidy[order(MAG.actino.sf.MTQ.aa.tidy$site),]
MAG.actino.sf.MTQ.aa.tidy$number = as.factor(MAG.actino.sf.MTQ.aa.tidy$number)

pdf("MAG_actino_enriched_MTQ.pdf")
ggplot(MAG.actino.sf.MTQ.aa.tidy , aes( x = MAG , y = enzyme , fill = number )) + facet_grid( .~site , scales = "free" , space = "free") + 
        geom_tile() + theme_bw() + theme(axis.text.x = element_text(angle = 90 , size = 10 ) , aspect.ratio = 1  , text = element_text( size = 15)) +
        scale_fill_manual(values = c("gray20","darkred","firebrick3", "firebrick1" , "orangered" , "tomato" , "orange"  , "gold"))
dev.off()

















