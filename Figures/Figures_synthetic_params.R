library(ggplot2)
library(latex2exp)

library(extrafont)
library(patchwork)
highlight <- function(x, value, col.value, col=NA, ...){
  hst <- hist(x, ...)
  idx <- findInterval(value, hst$breaks)
  cols <- rep(col, length(hst$counts))
  cols[idx] <- col.value
  hist(x, col=cols, ...)
}
colnames<-c('rho','beta','mean','var')
li <- read.table('../Synthetic_data/li.txt',col.names=colnames)
lr <- read.table('../Synthetic_data/lr.txt',col.names=colnames)
li_lr <- read.table('../Synthetic_data/li_lr.txt',col.names=colnames)

len = length(head(li$rho,1000))

 rho <- data.frame(Likelihood = c(rep("value_1",len), rep("value_2",len), rep("value_3",len) ) ,
                   value = c(head(li$rho,1000), head(lr$rho,1000), head(li_lr$rho,1000)))

theme_set(theme_bw())
par(mfrow=c(2,2))



ggplot(rho,aes(x=value,fill=Likelihood )) + geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#0C7EF5","#41F518", "#F5A318"),labels=c(bquote('\u2113'[i]), bquote('\u2113'[r]), 
                                                                      bquote('\u2113'[i]~+~'\u2113'[r])))+
  #labels=unname(TeX(c("$l_i$","$l_r$","l_i+l_r")))
  labs(x = TeX(r'($\\rho$)'), y = "density")+     
 # geom_density(aes(x=li$rho, y=..density..),size=2,colour=(rgb(163/255,35/255,142/255))) + 
  theme(text=element_text(size=12, family="LM Roman 10"),axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
  size = 1, linetype = "solid") ,panel.border = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+  
                        scale_y_continuous(expand = expansion(mult = c(0, .1)))+
                        scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+  
 #geom_segment(aes(x = mean(li$rho), y = 0, xend =  mean(li$rho), yend = 120), color="#0C7EF5")+
geom_point(aes(x=mean(li$rho), y=200), size=4, shape=22, color="black",fill="#0C7EF5",alpha=0.6)+
  geom_point(aes(x=mean(lr$rho), y=200), size=4, shape=23, color="black",fill="#41F518",alpha=0.6)+
geom_point(aes(x=mean(li_lr$rho), y=200), size=4, shape=24, color="black",fill="#F5A318",alpha=0.6)+
  geom_point(aes(x=50/9950, y=200), size=4, shape=25, color="black",fill="black",alpha=0.6)  




beta <- data.frame(Likelihood = c(rep("value_1",len), rep("value_2",len), rep("value_3",len) ) ,
                  value = c(head(li$beta,1000), head(lr$beta,1000), head(li_lr$beta,1000)))


ggplot(beta,aes(x=value,fill=Likelihood )) + geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#0C7EF5","#41F518", "#F5A318"),labels=c(bquote('\u2113'[i]), bquote('\u2113'[r]), 
                                                                      bquote('\u2113'[i]~+~'\u2113'[r])))+
  #scale_color_hue(labels = c(TeX(r'($\\ell_i$)'), "T888","mannaggia la madonna"))+
  labs(x = TeX(r'($\\beta$)'), y = "density")+     
  # geom_density(aes(x=li$rho, y=..density..),size=2,colour=(rgb(163/255,35/255,142/255))) + 
  theme(text=element_text(size=10, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                 size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  xlim(0.18,0.4)+
  #geom_segment(aes(x = mean(li$rho), y = 0, xend =  mean(li$rho), yend = 120), color="#0C7EF5")+
  geom_point(aes(x=mean(li$beta), y=200), size=4, shape=22, color="black",fill="#0C7EF5",alpha=0.6)+
  geom_point(aes(x=mean(lr$beta), y=200), size=4, shape=23, color="black",fill="#41F518",alpha=0.6)+
  geom_point(aes(x=mean(li_lr$beta), y=200), size=4, shape=24, color="black",fill="#F5A318",alpha=0.6)+
  geom_point(aes(x=0.25, y=200), size=4, shape=25, color="black",fill="black",alpha=0.6)  




mea <- data.frame(Likelihood = c(rep("value_1",len), rep("value_2",len), rep("value_3",len) ) ,
                  value = c(head(li$mean,1000), head(lr$mean,1000), head(li_lr$mean,1000)))
vara <- data.frame(Likelihood = c(rep("value_1",len), rep("value_2",len), rep("value_3",len) ) ,
                  value = c(head(li$var,1000), head(lr$var,1000), head(li_lr$var,1000)))




ggplot(mea,aes(x=value,fill=Likelihood )) + geom_histogram(alpha=0.6, position = 'identity',biwidth=0.1) +
  scale_fill_manual(values=c("#0C7EF5","#41F518", "#F5A318"),labels=c(bquote('\u2113'[i]), bquote('\u2113'[r]), 
                                                                      bquote('\u2113'[i]~+~'\u2113'[r])))+
  theme(text=element_text(size=14, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  #geom_segment(aes(x = mean(li$rho), y = 0, xend =  mean(li$rho), yend = 120), color="#0C7EF5")+
  geom_point(aes(x=mean(li$mean), y=800), size=4, shape=22, color="black",fill="#0C7EF5",alpha=0.6)+
  geom_point(aes(x=mean(lr$mean), y=800), size=4, shape=23, color="black",fill="#41F518",alpha=0.6)+
  geom_point(aes(x=mean(li_lr$mean), y=800), size=4, shape=24, color="black",fill="#F5A318",alpha=0.6)+
  geom_point(aes(x=9, y=800), size=4, shape=25, color="black",fill="black",alpha=0.6)  



ggplot(vara,aes(x=value,fill=Likelihood )) + geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#0C7EF5","#41F518", "#F5A318"),labels=c(bquote('\u2113'[i]), bquote('\u2113'[r]), 
                                                                      bquote('\u2113'[i]~+~'\u2113'[r])))+
  #scale_color_hue(labels = c(TeX(r'($\\ell_i$)'), "T888","mannaggia la madonna"))+
  labs(x = TeX(r'($\\variance$)'), y = "density")+     
  # geom_density(aes(x=li$rho, y=..density..),size=2,colour=(rgb(163/255,35/255,142/255))) + 
  theme(text=element_text(size=14, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  #geom_segment(aes(x = mean(li$rho), y = 0, xend =  mean(li$rho), yend = 120), color="#0C7EF5")+
  geom_point(aes(x=mean(li$var), y=150), size=4, shape=22, color="black",fill="#0C7EF5",alpha=0.6)+
  geom_point(aes(x=mean(lr$var), y=150), size=4, shape=23, color="black",fill="#41F518",alpha=0.6)+
  geom_point(aes(x=mean(li_lr$var), y=150), size=4, shape=24, color="black",fill="#F5A318",alpha=0.6)+
  geom_point(aes(x=6, y=150), size=4, shape=25, color="black",fill="black",alpha=0.6)  


  


#Rho_0
gg<- ggplot(li, aes(x=li$rho,y=..density..)) + geom_histogram(color="black", fill="white") +
                   labs(x = TeX(r'($\\rho$)'), y = "density")+     
                   geom_density(aes(x=li$rho, y=..density..),size=2,colour=(rgb(163/255,35/255,142/255))) + 
                   theme(axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                  size = 1, linetype = "solid", ),  panel.border = element_blank(),
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
                   scale_y_continuous(expand = expansion(mult = c(0, .1)))+
                   scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))
                 
gg+   geom_segment(aes(x = mean(li$rho), y = 0, xend =  mean(li$rho), yend = 120),colour = rgb(163/255,35/255,142/255))+
                   geom_point(aes(x=mean(li$rho), y=120,color=rgb(1,0,0)), size=3, shape=23)



#Histogram recovery mean vs std


rec_i = read.table("../Synthetic_data/mean_i.txt", col.names = c("mean_i","var_i"))
rec_r = read.table("../Synthetic_data/mean_r.txt", col.names = c("mean_r","var_r"))
rec_ir = read.table("../Synthetic_data/mean_ir.txt", col.names = c("mean_ir","var_ir"))



base<- ggplot(recdist, aes(x=meani,y=sqrt(var), color=cols))+
  geom_point()+
  theme(text=element_text(size=10, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.1,0.8))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5),limits=c(0,16))+
  scale_color_manual(values=c('#0C7EF5','#41F518','#F5A318'),labels=c(bquote('\u2113'[i]),bquote('\u2113'[r]), 
                                                                      bquote('\u2113'[i]~+~'\u2113'[r])), name=NULL)+
  labs(x = 'mean', y = "standard deviation")+
  geom_point(x=9,y=sqrt(6),color='black',size=2, fill='black')



dens1 <- ggplot(Density, aes(x=meani, fill=cols)) + 
  geom_density(alpha=0.6)+ 
  theme_void() + 
  #scale_fill_manual(values = unique(Density$cols))+
  scale_fill_manual(values=c('#0C7EF5','#41F518','#F5A318'))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5),limits=c(0,16))+
  geom_point(x=9,y=0.1,color='black',size=2)+
  theme(legend.position = "none")
dens2 <- ggplot(Density2, aes(x = sqrt(var),fill=cols)) + 
  geom_density(alpha=0.6)+ 
  scale_fill_manual(values=c('#0C7EF5','#41F518','#F5A318'))+
  theme_void() + 
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5),limits=c(1,4))+
  geom_point(x=sqrt(6),y=0.03,color='black',size=2)+
  theme(legend.position = "none")+
  coord_flip()

dens1 + plot_spacer() + base + dens2 + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(4, 4)
  ) 



