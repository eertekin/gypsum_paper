library(dplyr)
MAG1 = gene_quant_bin %>% select(gene , kegg_no , bin ) %>% filter(bin == "bin.1")
MAG1 = MAG1[complete.cases(MAG1$kegg_no),]
write.table(MAG1[,1:2] , file = "~/Desktop/MAG1.keggs.txt" ,row.names = FALSE , sep = "\t")


MAG4 = gene_quant_bin %>% select(gene , kegg_no , bin ) %>% filter(bin == "bin.4")
MAG4 = MAG4[complete.cases(MAG4$kegg_no),]
write.table(MAG4[,1:2] , file = "~/Desktop/MAG4.keggs.txt" ,row.names = FALSE , sep = "\t")

MAG11 = gene_quant_bin %>% select(gene , kegg_no , bin ) %>% filter(bin == "bin.11")
MAG11 = MAG11[complete.cases(MAG11$kegg_no),]
write.table(MAG11[,1:2] , file = "~/Desktop/MAG11.keggs.txt" ,row.names = FALSE , sep = "\t")

MAG20 = gene_quant_bin %>% select(gene , kegg_no , bin ) %>% filter(bin == "bin.20")
MAG20 = MAG20[complete.cases(MAG20$kegg_no),]
write.table(MAG20[,1:2] , file = "~/Desktop/MAG20.keggs.txt" ,row.names = FALSE , sep = "\t")



kegg_binned = gene_quant_bin[complete.cases(gene_quant_bin$gene),]
library(plyr)
count(kegg_binned$kegg_no)
collapse_MAGs = function(x) {
        for ( i in levels(x$bin)) {
             assign( i , subset(kegg_binned , x$bin == "i" , c(13,15,16)) )
         
        }
}

collapse_MAGs(kegg_binned)
bins_list = mget(ls(pattern = "kegg_bin.")) ; bins_kegg_quant = lapply(bins_list , FUN = function(x) {count(x$kegg_no)})
bins_kegg_count.2 = Reduce(function(x,y) {merge(x,y , by = "x", all = TRUE )} , bins_kegg_quant ) 
colnames(bins_kegg_count.2) = c("kegg_no",paste("bin." , c(1,10:19,2,20:29,3,30:39,4,40:44,5,6,7,8,9) , sep = ""))
bins_kegg_count_spread = merge( kegg_enzyme_desc , bins_kegg_count.2  , by = "kegg_no" , all.y = TRUE)
bins_kegg_count_spread[is.na(bins_kegg_count_spread)] = 0
chloroflexi = c("bin.1" , "bin.25" , "bin.19" , "bin.33")
c.fix = c("K01601" , "K01602" , "K01965" , "K01966" , "K00174" )
ox.ps = c("K02689" , "K02690" , "K02703" , "K02704" ) 

dnra = c("K00362" , "K00363")
no3.no2 = paste(c("K00370" , "K00371" , "K00374" , "K02567" , "K02568" , "K00367" , "K10534" , "K00372" , "K00360") , collapse = "|")  
no2.nh3 = paste(c("K17877" , "K00366") , collapse = "|")
no2.n2o = paste(c("K00368" , "K15864" , "K04561" , "K02305") , collapse = "|")


so4.paps = c("K13811" , "K00958" , "K00955" , "K00956" , "K00957")
paps.so3 = paste(c("K13811" , "00955" , "K00860" , "K00390" ) , collapse = "|" )
so3.h2s = c("K00380" , "K00381" , "K00392" )
thio.ox = c("K17226" , "K17227")
so3.ox = paste(c("K00387" , "K17218" ) , collapse = "|" )
anox.ps = c( "K11333" , "K11334" , "K11335" , "K08926" , "K08928" , "K08929" )

MAG_pathways = data.frame( MAG = character() , c.fix = numeric() , ox.ps = character() , anox.ps = character() , dnra = character() , no3.no2 = character() , no2.nh3 = character() , no2.n2o = character() ,
                           so4.paps = character() , paps.so3 = character() , so3.h2s = character() , thio.ox = character() , so3.ox = character() ,
                           lineage = character() , stringsAsFactors = FALSE)

        d = bins_kegg_count_taxa ; d$bin = factor(d$bin) ; d$lineage = factor(d$lineage)
for ( i in 1:length( unique(levels(d$bin) ) ) ) {
        b = d %>% filter(bin == levels(d$bin)[i]) ; 
        b = b[complete.cases(b$number),]
        MAG_pathways[i,1] = levels(d$bin)[i]
        MAG_pathways[i,14] = as.character(b$lineage[1])
        ## carbon fixation
        if((c.fix[1] %in% b$kegg_no & c.fix[2] %in% b$kegg_no )|(c.fix[3] %in% b$kegg_no & c.fix[4] %in% b$kegg_no)|(c.fix[5] %in% b$kegg_no) ) {
                MAG_pathways[i,2] = 1 } else{ MAG_pathways[i,2] = 0 }
        ##oxygenic ps##
        if((ox.ps[1] %in% b$kegg_no & ox.ps[2] %in% b$kegg_no) | (ox.ps[3] %in% b$kegg_no & ox.ps[4] %in% b$kegg_no)) {
                MAG_pathways[i,3] = 1 } else{ MAG_pathways[i,3] = 0 }
        ## dnra
        if(dnra[1] %in% b$kegg_no & dnra[2] %in% b$kegg_no) { MAG_pathways[i,5] = 1} else{ MAG_pathways[i,5] = 0 }
        ## nitrate reduction
        if( any(grepl(no3.no2 , b$kegg_no)) ) { MAG_pathways[i,6] = 1 } else{ MAG_pathways[i,6] = 0 }
        ## nitrite reduction
        if( any(grepl(no2.nh3 , b$kegg_no) ) ) { MAG_pathways[i,7] = 1 } else{ MAG_pathways[i,7] = 0 }
        ##dnit ##
        if( any(grepl(no2.n2o , b$kegg_no) ) ) { MAG_pathways[i,8] = 1 } else{ MAG_pathways[i,8] = 0 }
        ##so4.paps##
        if( so4.paps[1] %in% b$kegg_no | so4.paps[2] %in% b$kegg_no | ( (so4.paps[3] %in% b$kegg_no) & (so4.paps[4] %in% b$kegg_no) & (so4.paps[5] %in% b$kegg_no) ) )  { 
                MAG_pathways[i,9] = 1 } else{ MAG_pathways[i,9] = 0 }
        ##
        if( any(grepl(paps.so3 , b$kegg_no))){ MAG_pathways[i,10] = 1 } else{ MAG_pathways[i,10] = 0 }
        ##so3.h2s##
        if((so3.h2s[1] %in% b$kegg_no & so3.h2s[2] %in% b$kegg_no) | ( so3.h2s[3] %in% b$kegg_no )) { MAG_pathways[i,11] = 1 } else{ MAG_pathways[i,11] = 0 }
        ## thio.ox ##
        if(thio.ox[1] %in% b$kegg_no & thio.ox[2] %in% b$kegg_no)  { MAG_pathways[i,12] = 1 } else{ MAG_pathways[i,12] = 0}
        ##so3.ox##
        if( any(grepl(so3.ox, b$kegg_no)))  { MAG_pathways[i,13] = 1 } else{ MAG_pathways[i,13] = 0}
        ## anox.ps ##
        if( ( anox.ps[1] %in% b$kegg_no & anox.ps[2] %in% b$kegg_no & 
              anox.ps[3] %in% b$kegg_no & anox.ps[4] %in% b$kegg_no & anox.ps[4] %in% b$kegg_no & anox.ps[5] %in% b$kegg_no & 
              anox.ps[6] %in% b$kegg_no & anox.ps[4] %in% b$kegg_no )  ) { MAG_pathways[i,4] = 1 } else { MAG_pathways[i,4] = 0}
        
}

#MAG_taxa = bins_kegg_count_brite[,c(6,8)] 
#MAG_taxa = unique(MAG_taxa) ; MAG_taxa$bin = str_replace(MAG_taxa$bin , pattern = "bin." , replacement = "MAG")
#colnames(MAG_taxa) = c("MAG" , "lineage")
MAG_pathways_tidy = gather(MAG_pathways , key = "pathway" , value = "presence" , 2:13)
#MAG_pathways_tidy = merge(MAG_pathways_tidy , MAG_taxa , by = "bin" , all.x = TRUE)

for (i in 1:nrow(MAG_pathways_tidy)) {
        if(MAG_pathways_tidy$MAG[i] %in% chloroflexi) { MAG_pathways_tidy$lineage[i] = "chloroflexi"} else { MAG_pathways_tidy$lineage[i] = MAG_pathways_tidy$lineage[i]}
} 

for (i in 1:nrow(MAG_pathways_tidy)) {
        if(MAG_pathways_tidy$pathway[i] == "c.fix") { MAG_pathways_tidy$metabolism[i] = "autotrophy" } else if(MAG_pathways_tidy$pathway[i] == "dnra" | MAG_pathways_tidy$pathway[i] == "no3.no2" | MAG_pathways_tidy$pathway[i] == "no2.nh3" | MAG_pathways_tidy$pathway[i] == "no2.n2o" ) {
                MAG_pathways_tidy$metabolism[i] = "nitrogen" } else if (MAG_pathways_tidy$pathway[i] == "anox.ps") { MAG_pathways_tidy$metabolism[i] = "anoxygenic photosynthesis" 
                } else if(MAG_pathways_tidy$pathway[i] == "ox.ps") { MAG_pathways_tidy$metabolism[i] = "oxygenic photosynthesis" } else { MAG_pathways_tidy$metabolism[i] = "sulfur" }
}

MAG_pathways_tidy$MAG = str_replace(MAG_pathways_tidy$MAG , pattern = 'bin.' , replacement = "MAG" )
MAG_pathways_tidy$pathway = factor( MAG_pathways_tidy$pathway , levels = c("anox.ps" , "ox.ps" , "c.fix" , rev(c("no3.no2" , "no2.nh3" , "dnra", "no2.n2o")) , rev(c( "so4.paps" , "paps.so3" ,  "so3.h2s" , "so3.ox" , "thio.ox"))  ) )
MAG_pathways_tidy$MAG = factor(MAG_pathways_tidy$MAG ,  levels = c(paste("MAG" , c(35,17,16,13,4,28,40,29)  , sep = "") , paste("MAG" , c(1,25,19,33) , sep = "") , paste("MAG" , c(23,41,30,20,9,43,31,36,37,26,5,12,8,27,3,24,34,42,11,44) , sep = "") ,
                                                                   paste("MAG" , c(7,22,39) , sep = "") , paste("MAG" , c(10,14,2,21,32,38,18) , sep = "") , paste("MAG" , c(6,15) , sep = "") ))
ggplot(MAG_pathways_tidy , aes(MAG , pathway , fill = presence)) +
        geom_tile(color = "gray21") + facet_grid(metabolism ~ lineage , space = "free" , scales = "free") + scale_fill_manual(values = c("gray91" , "royalblue")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) 











