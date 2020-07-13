library(ggplot2) ; library(stringr) ; library(dplyr) ; library(tidyr)
x.names = c( 10, 50 , expression(10^2) , bquote(paste("5 x " , 10^2)) , 
                 expression(10^3) , bquote(paste("5 x " , 10^3)) , expression(10^4) , expression(10^5) , expression(10^6) , expression(10^7) )

ct.freq = CT.FREQ
ct.freq = gather(ct.freq , key = site , value = freq , 2:3 )
ct.freq$VOLUME = factor(ct.freq$VOLUME)

ggplot(ct.freq , aes(x = VOLUME , y = freq , fill = site)) + 
        geom_bar(stat  = "identity" , position = "dodge") +
        theme_bw() + theme( axis.title=element_text(size=20) , axis.text = element_text(size = 15) ,
        aspect.ratio = 0.5 , legend.position = "none" , plot.margin = unit(c(0.5,0.5,0.5,0.5) , "in") ) +
        scale_x_discrete(labels = x.names) + 
        theme(axis.text.x = element_text(angle = 90) , legend.position = "top") + 
        xlab( bquote(paste("Volume categories " , (µm^3))) ) + ylab("Frequence (%)") +
        scale_fill_manual( labels=c("CL-MTQ", "KM") , values = c("hotpink", "dodgerblue") )
                          
ggsave("ct_scan_feq.pdf" , height = 5 , width = 8)

