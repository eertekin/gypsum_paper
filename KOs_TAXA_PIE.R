## ENRICHED KOS TAXA PIECHART ##
library(dplyr) ; library(plyr)
ko.sf.taxa = kegg.padj.sf.phylo[,c(1,9:17,21)]
ko.sf.taxa = separate(ko.sf.taxa , 11 , into = c("domain" , "phylum" , "class" , "order", "family" , "genus" , "species") , sep = ";")
ko.sf.taxa.KM = ko.sf.taxa %>% filter(kegg_no %in% kegg.padj.sf_KM$kegg_no)
ko.sf.taxa.KM.count = count(ko.sf.taxa.KM$class)
colnames(ko.sf.taxa.KM.count) = c("class" , "KM")


ko.sf.taxa.MTQ = ko.sf.taxa %>% filter(kegg_no %in% kegg.padj.sf_MTQ$kegg_no)
ko.sf.taxa.MTQ.count = count(ko.sf.taxa.MTQ$class)
colnames(ko.sf.taxa.MTQ.count) = c("class" , "MTQ")

ko.sf.taxa.counts = merge(ko.sf.taxa.MTQ.count , ko.sf.taxa.KM.count , by = "class" , all = TRUE)
ko.sf.taxa.counts[is.na(ko.sf.taxa.counts)] = 0
ko.sf.taxa.counts = ko.sf.taxa.counts[(order(ko.sf.taxa.counts$MTQ , decreasing = TRUE )),]

        
ko.sf.taxa.counts$class = factor(ko.sf.taxa.counts$class , levels = ko.sf.taxa.counts$class )
generate_others = function(d) {
        d[,1] = factor(d[,1] , levels = c(levels(d[,1]) , "others")) 
        d = d[order(d[,2], decreasing = TRUE),]
        sums = colSums(d[1:10,-1])
        others = colSums(d[,-1]) - sums
        d_others = data.frame(rbind(d[1:10,] , c("others" , others ) ) )
        
}

ko.sf.taxa.counts_others = generate_others(ko.sf.taxa.counts)
ko.sf.taxa.counts_others$class = factor(ko.sf.taxa.counts_others$class , levels = as.character(ko.sf.taxa.counts_others$class))

ko.sf.taxa.counts_others$KM = as.numeric(ko.sf.taxa.counts_others$KM)
pie_colors = c("tomato" , "palegreen2" , "skyblue1" , "olivedrab2" , "aquamarine1" ,"lightpink" , "coral" ,   "orchid1" , "orange1" , "chocolate1" , "gray91" )

pdf("taxa_pie_MTQ.pdf")
p1 = ggplot(ko.sf.taxa.counts_others[,1:2] , aes(x = "" , y = MTQ , fill = class)) +
        geom_bar(stat = "identity" , width=0.25, color="gray60" ) + 
        coord_polar("y" , start=0) + theme_void() +
        theme(legend.text = element_text(face = "italic" , size = 15) , 
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in")) +
        scale_fill_manual( values = pie_colors)
dev.off()
pdf("taxa_pie_KM.pdf")
p2 = ggplot(ko.sf.taxa.counts_others[,c(1,3)] , aes(x = "" , y = KM , fill = class)) +
        geom_bar(stat = "identity" , width=0.25, color="gray60" ) + 
        coord_polar("y" , start=0) + theme_void() +
        theme(legend.text = element_text(face = "italic" , size = 15) , 
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in")) +
        scale_fill_manual( values = pie_colors)
dev.off()
ggsave(p1,p2, filename = "pies.pdf")



