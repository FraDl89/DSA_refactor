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
li <- read.table('/home/fra/Projects/DSA_inference/Synthetic_data/li.txt',col.names=colnames)
lr <- read.table('/home/fra/Projects/DSA_inference/Synthetic_data/lr.txt',col.names=colnames)
li_lr <- read.table('/home/fra/Projects/DSA_inference/Synthetic_data/li_lr.txt',col.names=colnames)

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

ggsave("outputFile",plot=last_plot(), "pdf", units="in", width=5, height = 4)



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
  #scale_color_hue(labels = c(TeX(r'($\\ell_i$)'), "T888","mannaggia la madonna"))+
  labs(x = 'mean'), y = "density")+     
  # geom_density(aes(x=li$rho, y=..density..),size=2,colour=(rgb(163/255,35/255,142/255))) + 
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














colnames <- c("data_2_S","data_25_S","data_50_S","data_75_S","data_97_S")
data_S = read.table("/home/fra/Projects/DSA_inference/data_S.txt", col.names = colnames)

colnames <- c("data_2_I","data_25_I","data_50_I","data_75_I","data_97_I")
data_I = read.table("/home/fra/Projects/DSA_inference/data_I.txt", col.names = colnames)

colnames <- c("data_2_dS","data_25_dS","data_50_dS","data_75_dS","data_97_dS")
data_dS = read.table("/home/fra/Projects/DSA_inference/data_dS.txt", col.names = colnames)

colname<-c("dS")
true_dS=read.table("/home/fra/Projects/DSA_inference/true_ds", col.names = colname)

gridspacing=1000
time = seq(0,80-80/gridspacing,80/gridspacing,)
data_S$time = time


colnames <- c("li_lr_2_S","li_lr_25_S","li_lr_50_S","li_lr_75_S","li_lr_97_S")
li_lr_S = read.table("/home/fra/Projects/DSA_inference/li_lr_S.txt", col.names = colnames)

colnames <- c("li_lr_2_I","li_lr_25_I","li_lr_50_I","li_lr_75_I","li_lr_97_I")
li_lr_I = read.table("/home/fra/Projects/DSA_inference/li_lr_I.txt", col.names = colnames)

colnames <- c("li_lr_2_dS","li_lr_25_dS","li_lr_50_dS","li_lr_75_dS","li_lr_97_dS")
li_lr_dS = read.table("/home/fra/Projects/DSA_inference/li_lr_dS.txt", col.names = colnames)




colnames <- c("li_2_S","li_25_S","li_50_S","li_75_S","li_97_S")
li_S = read.table("/home/fra/Projects/DSA_inference/li_S.txt", col.names = colnames)

colnames <- c("li_2_I","li_25_I","li_50_I","li_75_I","li_97_I")
li_I = read.table("/home/fra/Projects/DSA_inference/li_I.txt", col.names = colnames)

colnames <- c("li_2_dS","li_25_dS","li_50_dS","li_75_dS","li_97_dS")
li_dS = read.table("/home/fra/Projects/DSA_inference/li_dS.txt", col.names = colnames)




colnames <- c("lr_2_S","lr_25_S","lr_50_S","lr_75_S","lr_97_S")
lr_S = read.table("/home/fra/Projects/DSA_inference/lr_S.txt", col.names = colnames)

colnames <- c("lr_2_I","lr_25_I","lr_50_I","lr_75_I","lr_97_I")
lr_I = read.table("/home/fra/Projects/DSA_inference/lr_I.txt", col.names = colnames)

colnames <- c("lr_2_dS","lr_25_dS","lr_50_dS","lr_75_dS","lr_97_dS")
lr_dS = read.table("/home/fra/Projects/DSA_inference/lr_dS.txt", col.names = colnames)



li_2_S = (li_S$li_2_S-li_S$li_2_S[5000])/(1-li_S$li_2_S[5000])
li_97_S = (li_S$li_97_S-li_S$li_97_S[5000])/(1-li_S$li_97_S[5000])

lr_2_S = (lr_S$lr_2_S-lr_S$lr_2_S[5000])/(1-lr_S$lr_2_S[5000])
lr_97_S = (lr_S$lr_97_S-lr_S$lr_97_S[5000])/(1-lr_S$lr_97_S[5000])

li_lr_2_S = (li_lr_S$li_lr_2_S-li_lr_S$li_lr_2_S[5000])/(1-li_lr_S$li_lr_2_S[5000])
li_lr_97_S = (li_lr_S$li_lr_97_S-li_lr_S$li_lr_97_S[5000])/(1-li_lr_S$li_lr_97_S[5000])





gg<- ggplot(data_S, aes(x=time,y=data_50_S)) + geom_line() +
      #geom_line(y=data_S$data_25_S,linetype = "dashed") +
      #geom_line(y=data_S$data_75_S,linetype = "dashed") + 
      geom_line(y=data_S$data_2_S,linetype = "dotted") +
      geom_line(y=data_S$data_97_S,linetype = "dotted")+
  labs(x = TeX(r'($time$)'), y = TeX(r"($S$)"))+     
  theme(text=element_text(size=10, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5), limits=c(0,50))+
 
  geom_ribbon(aes(ymin=li_2_S[seq(1, length(li_S$li_2_S),5)],ymax=li_97_S[seq(1, length(li_S$li_97_S), 5)]), fill="#0C7EF5",alpha=.5)+     
geom_ribbon(aes(ymin=lr_2_S[seq(1, length(li_S$li_2_S),5)],ymax=lr_97_S[seq(1, length(li_S$li_97_S), 5)]), fill="#41F518",alpha=.5)+     
geom_ribbon(aes(ymin=li_lr_2_S[seq(1, length(li_S$li_2_S),5)],ymax=li_lr_97_S[seq(1, length(li_S$li_97_S), 5)]), fill="#F5A318",alpha=.5)  
  
   
  




len2=999

li_2_dS   = li_dS$li_2_dS[seq(1, length(li_dS$li_2_dS),5)]
li_50_dS  = li_dS$li_50_dS[seq(1, length(li_dS$li_2_dS),5)]
li_97_dS  = li_dS$li_97_dS[seq(1, length(li_dS$li_2_dS),5)]

lr_2_dS   = lr_dS$lr_2_dS[seq(1, length(lr_dS$lr_2_dS),5)]
lr_97_dS  = lr_dS$lr_97_dS[seq(1, length(lr_dS$lr_2_dS),5)]

li_lr_2_dS  = li_lr_dS$li_lr_2_dS[seq(1, length(lr_dS$lr_2_dS),5)]
li_lr_97_dS = li_lr_dS$li_lr_97_dS[seq(1, length(lr_dS$lr_2_dS),5)]

datads = true_dS$dS[seq(1, length(lr_dS$lr_2_dS),5)]

time2=head(time,len2)


empirical <- data.frame( time= time2, li_2_dS=head(li_2_dS,len2),li_97_dS=head(li_97_dS,len2),
                         li_2_dS=head(li_2_dS,len2),li_97_dS=head(li_97_dS,len2),
                         lr_2_dS=head(lr_2_dS,len2),lr_97_dS=head(lr_97_dS,len2),
                         li_lr_2_dS=head(li_lr_2_dS,len2),li_lr_97_dS=head(li_lr_97_dS,len2),
                         d=head(datads,len2),
                         li_50_dS=head(li_50_dS,len2))





ggplot(empirical,aes(x=time)) +
  geom_ribbon(aes(ymin=li_2_dS,ymax=li_97_dS), fill="#0C7EF5",alpha=0.9)+
  geom_ribbon(aes(ymin=lr_2_dS,ymax=lr_97_dS, fill="#41F518"),alpha=0.9)+
  geom_ribbon(aes(ymin=li_lr_2_dS,ymax=li_lr_97_dS, fill="#F5A318"),alpha=0.9)+
  geom_line(aes(y=empirical$d, color='black'))+
  scale_linetype_manual(values=c("dashed"), labels=c(bquote('\u2113'[i])), aes(color='#0C7EF5'), name=NULL)+
  
  scale_fill_manual(values=c("#41F518", "#F5A318"),labels=c(bquote('\u2113'[r]), 
                                                                      bquote('\u2113'[i]~+~'\u2113'[r])), name=NULL)+
  scale_color_manual(values=c("black"), labels=c(bquote('True')), name=NULL)+
  labs(x = 'time', y = "empirical density")+     
  # geom_density(aes(x=li$rho, y=..density..),size=2,colour=(rgb(163/255,35/255,142/255))) + 
  theme(text=element_text(size=10, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  geom_line(aes(y=li_2_dS, linetype='dashed'), color='#0C7EF5')+
  geom_line(aes(y=li_97_dS), color='#0C7EF5', linetype='dashed')
 





rec_i = read.table("/home/fra/Projects/DSA_inference/mean_i.txt", col.names = c("mean_i","var_i"))
rec_r = read.table("/home/fra/Projects/DSA_inference/mean_r.txt", col.names = c("mean_r","var_r"))
rec_ir = read.table("/home/fra/Projects/DSA_inference/mean_ir.txt", col.names = c("mean_ir","var_ir"))



Density <- data.frame(meani=c(rec_i$mean_i,rec_r$mean_r,rec_ir$mean_ir), 
                      likelihood=c(rep("li",907),rep("lr",907),rep("lir",907)),
                      cols= c(rep("blue",907),rep("green",907),rep("orange",907)))

Density2 <- data.frame(var=c(rec_i$var_i,rec_r$var_r,rec_ir$var_ir), 
                       likelihood=c(rep("li",907),rep("lr",907),rep("lir",907)),
                       cols= c(rep("blue",907),rep("green",907),rep("orange",907)))

recdist <- data.frame(meani=rec_i$mean_i, vari=rec_i$var_i, meanr=rec_r$mean_r, varr=rec_r$var_r, 
                      meanir=rec_ir$mean_ir, varir=rec_ir$var_ir)

recdist <- data.frame(meani=c(rec_i$mean_i,rec_r$mean_r,rec_ir$mean_ir), 
                      var=c(rec_i$var_i,rec_r$var_r,rec_ir$var_ir),
                      cols= c(rep("#0C7EF5",907),rep("#41F518",907),rep("#F5A318",907)))




base<- ggplot(recdist, aes(x=meani,y=vari,color='#0C7EF5'))+geom_point(alpha=0.6)+
geom_point(x=recdist$meanir,y=recdist$varir, aes(color='#41F518'),alpha=0.6)+
  geom_point(x=recdist$meanr,y=recdist$varr, aes(color='#F5A318'),alpha=0.6)+
  theme(text=element_text(size=10, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.1,0.8))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5),limits=c(0,16))+
  scale_color_manual(values=c('#0C7EF5','#41F518','#F5A318'),labels=c(bquote('\u2113'[i]),bquote('\u2113'[r]), 
                                                            bquote('\u2113'[i]~+~'\u2113'[r])), name=NULL)+
  labs(x = 'mean', y = "variance")+
  geom_point(x=9,y=6,color='black',size=3, fill='black')
  


base<- ggplot(recdist, aes(x=meani,y=var, color=cols))+
      geom_point()+
  theme(text=element_text(size=10, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.1,0.8))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5),limits=c(0,16))+
  scale_color_manual(values=c('#0C7EF5','#41F518','#F5A318'),labels=c(bquote('\u2113'[i]),bquote('\u2113'[r]), 
                                                                      bquote('\u2113'[i]~+~'\u2113'[r])), name=NULL)+
  labs(x = 'mean', y = "variance")+
  geom_point(x=9,y=6,color='black',size=2, fill='black')



  dens1 <- ggplot(Density, aes(x=meani, fill=cols)) + 
  geom_density(alpha=0.6)+ 
  theme_void() + 
    #scale_fill_manual(values = unique(Density$cols))+
  scale_fill_manual(values=c('#0C7EF5','#41F518','#F5A318'))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5),limits=c(5.5,11.7))+
  geom_point(x=9,y=0.1,color='black',size=2)+
  theme(legend.position = "none")
dens2 <- ggplot(Density2, aes(x = var,fill=cols)) + 
  geom_density(alpha=0.6)+ 
  scale_fill_manual(values=c('#0C7EF5','#41F518','#F5A318'))+
  theme_void() + 
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5),limits=c(0,16))+
  geom_point(x=5.5,y=0.03,color='black',size=2)+
  theme(legend.position = "none")+
  coord_flip()

dens1 + plot_spacer() + base + dens2 + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(4, 4)
  ) 



x <- rnorm(300)
y <- rt(300, df = 2)
xy <- data.frame(x, y)

plot1 <- ggplot(xy, aes(x = x, y = y)) + 
  geom_point() 

dens1 <- ggplot(xy, aes(x = x)) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_void()
dens2 <- ggplot(xy, aes(x = y)) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_void() + 
  coord_flip()

