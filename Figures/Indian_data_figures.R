library(ggplot2)
library(latex2exp)

library(extrafont)
library(patchwork)

colnames<-c('day','cases')

dataset <- read.table("../India/Data/COVID19_timeseries_India_1227.csv", sep=',', header=TRUE)

real_cases <- read.table("../India/Data/infected_cases.txt",col.names=colnames)

real_rec <- read.table("../India/Data/recovered_recoveries.txt",col.names=colnames)
library(gghighlight)

options(stringsAsFactors = FALSE)
dataset$Date <- as.Date(dataset$Date_YMD, "%Y-%m-%d")

ggplot(dataset, aes(x=Date),color='black')  +
  geom_point(aes(y=Daily.Confirmed),size=2,colour=(rgb(163/255,35/255,142/255))) + 
  gghighlight(Date>="2021-02-15",Date<="2021-07-01",unhighlighted_params = list(color='black'))+
  theme_bw()+
  theme(text=element_text(size=12, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(x = 'time (days)', y = "daily new cases")     
  

#R_0 plots
r_0 <- read.table('../India/Data/R_0_distr.txt', col.names='r0')

ggplot(r_0, aes(x=r0))+
  geom_histogram(alpha=0.6,fill='#F5A318',color="#F5A318")+
  theme_bw()+
  theme(text=element_text(size=12, family="LM Roman 10"),axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                  size = 1, linetype = "solid") ,panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  labs(x = TeX(r'($R_0$)'), y='ćounnts')



result_x <- c(4.06914629e-04, 2.81346094e+00, 1.57537281e+00, 5.69340232e+00,
          1.89918160e+01)


x_lower <- 0

x_upper <- 25

shape_gam_R = (5.69340232e+00)**2/1.89918160e+01
scale_gam_R = 1.89918160e+01/(5.69340232e+00) 

shape_gam_I = 2.81346094e+00
scale_gam_I = 1.57537281e+00


max_height_I <- max(dgamma(x_lower:x_upper, shape = shape_gam_I, scale=scale_gam_I, log = FALSE))+0.03

#Pdfs

ggplot(data.frame(x = c(x_lower, x_upper)), aes(x = x)) + xlim(x_lower, x_upper) + 
  ylim(0, max_height_I) +
  theme_bw()+
  stat_function(fun = dgamma, args = list(shape = shape_gam_I, scale=scale_gam_I), size=2, aes(colour="orange")) + 
  stat_function(fun = dgamma, args = list(shape = shape_gam_R , scale=scale_gam_R),size=2, aes(colour="blue"))+
  labs(x = "\n time (days)", y = "pdf \n") + 
  theme(text=element_text(size=12, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  scale_colour_manual(breaks=c("orange","blue"), values=c("#0C7EF5","#F5A318"),labels=c("Infectiousness", "Recovery time"), name="Distributions")



library(gridExtra)

#Confidence intervals

day = head(day, T_f)
day = tail(day, -6)
colnames=c("t",'two','five','nine')

conf_int_i = read.table("../India/Data/infected_confidence_intervals.txt",col.names=colnames)
conf_int_r = read.table("../India/Data/recovered_confidence_intervals.txt",col.names=colnames)

colnames=c("t",'i')
cases = read.table("../India/Data/infected_cases.txt",col.names=colnames)
recoveries = read.table("../India/Data/recovered_recoveries.txt",col.names=colnames)


pp1 <- ggplot(conf_int_i, aes(x=t,y=five))  + geom_line(aes(color='conf'))+
  labs(x = 'days since 2021-02-15', y ='cases')+     
  theme_bw()+
  theme(text=element_text(size=12, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top")+  
  geom_ribbon(aes(ymin=two,ymax=nine, fill="#0C7EF5"),alpha=.5)+     
  geom_point(data=cases,x=cases$t, y=cases$i,size= 1, color=rgb(163/255,35/255,142/255)) +
  scale_fill_manual(values="#0C7EF5",labels="95 CI", name=NULL)+
  scale_colour_manual("", breaks=c("conf"), values=c('blue'), labels=c("median"))+
  theme(legend.spacing.y = unit(0.00, 'cm'))+
  guides(size= guide_legend(order = 1),
         colour = guide_legend(order = 2), 
         shape = guide_legend(order = 3))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5), limits=c(0,4.3e5))+
  scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5),limits = c(0,136))

pp2 <- ggplot(conf_int_r, aes(x=t,y=five))  + geom_line(aes(color='conf'))+
  labs(x = 'days since 2021-02-15', y ='removals')+     
  theme_bw()+
  theme(text=element_text(size=12, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid"),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top")+  
  geom_ribbon(aes(ymin=two,ymax=nine, fill="#F5A318"),alpha=.5)+     
  geom_point(data=recoveries,x=recoveries$t, y=recoveries$i,size= 1, shape=0, color=rgb(163/255,35/255,142/255)) +
  scale_fill_manual(values="#F5A318",labels="95 CI", name=NULL)+
  scale_colour_manual("", breaks=c("conf"), values=c('red'), labels=c("median"))+
  theme(legend.spacing.y = unit(0.00, 'cm'))+
  guides(size= guide_legend(order = 1),
         colour = guide_legend(order = 2), 
         shape = guide_legend(order = 3))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5), limits=c(0,4.3e5))+
  scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5),limits = c(0,136))

grid.arrange(pp1,pp2,ncol=2)

