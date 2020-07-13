## WEATHER FIGURE ##

weather.data = MTQ_KM_humidity_data 
weather.data$range = factor(weather.data$range , levels = unique(weather.data$range))

pdf("climate_data.pdf")

ggplot(weather.data, aes(x = weather.data$range)) +
        geom_bar( aes(y = weather.data$KM) , stat = "identity" , alpha=.5,fill='royalblue',color='royalblue') +
        geom_bar( aes(y = weather.data$MTQ) , stat = "identity"  , alpha=.3,fill='hotpink',color='hotpink')   +
        theme_bw() + theme(aspect.ratio = 1 , axis.title = element_text(size = 15) , axis.text = element_text(size = 15) ,
                           axis.text.x = element_text(angle = 90) ) + theme(legend.position = "top") +
        labs(x = " Relative Humidity (%)" , y = "Total number of hours in a year")                           
dev.off()
