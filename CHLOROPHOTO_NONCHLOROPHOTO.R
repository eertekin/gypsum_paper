## CHLOROPHOTOTROPHS HETEROTROPHS BARPLOT ##

photo_nonphoto_tidy = gather(photo_nonphoto , key = "site" , value = "abundance" , 2:10 )
photo_nonphoto_tidy$site = factor(photo_nonphoto_tidy$site , levels = c("CL1","CL2","CL3","MTQ1","MTQ2","MTQ3","KM1","KM2","KM3"))
ggplot(photo_nonphoto_tidy , aes(x = site , y = abundance , fill = group )) +
        geom_bar(stat= "identity") +
         theme_bw() + 
        theme(aspect.ratio = 1 , text = element_text(size = 15) , 
              axis.title.y = element_text(size = 15) , axis.text.x = element_text(angle = 90) , 
              legend.position = "top" ,
              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "in") ,
              legend.title = element_blank()) + scale_fill_manual(values = c("tomato" , "yellowgreen") ,
                                                                  labels = c("non-chlorophotrophs" , "chlorophototrophs")) +
        labs(y = "Relative abundance(%)" , x = "")
ggsave("photo_nonphoto.pdf" , height = 5 , width = 5)
