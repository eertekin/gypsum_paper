library(stringr) ; library(RColorBrewer)

## PHEATMAPS ##


actino_data = CSS_actino_otu_notax
chloro_data  = CSS_otu_chloro_notax
proteo_data = CSS_normalized_otu_all_proteo_notax
cyano_data = CSS_otu_cyano_notax
annotation_sites = data.frame( row.names = colnames(actino_data[,2:35]))
annotation_sites$community = c(rep("CL" , 17) , rep("MTQ" , 8 ) , rep("KM" , 9))
annotation_colors = list(
        community = c(CL = "#FFC000" , KM = "#00B0F0", MTQ = "#92D050" ))


pheatmap(actino_data[,2:35] , annotation_col = annotation_sites , border_color = NA ,
         annotation_colors = annotation_colors ,show_colnames = FALSE , show_rownames = FALSE ,
         color = colorRampPalette(c("white" , "blue" ,"yellow" , "red"))(250) )
pheatmap(proteo_data[,2:35] , annotation_col = annotation_sites , border_color = NA ,
         annotation_colors = annotation_colors ,show_colnames = FALSE , show_rownames = FALSE ,
         color = colorRampPalette(c("white" , "blue" ,"yellow" , "red"))(250) )
pheatmap(chloro_data[,2:35] , annotation_col = annotation_sites , border_color = NA ,
         annotation_colors = annotation_colors ,show_colnames = FALSE , show_rownames = FALSE ,
         color = colorRampPalette(c("white" , "blue" ,"yellow" , "red"))(250) )
pheatmap(cyano_data[,2:35] , annotation_col = annotation_sites , border_color = NA ,
         annotation_colors = annotation_colors ,show_colnames = FALSE , show_rownames = FALSE ,
         color = colorRampPalette(c("white" , "blue" ,"yellow" , "red"))(250) )
