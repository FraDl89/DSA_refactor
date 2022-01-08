library(ggplot2)
library(latex2exp)

library(extrafont)
library(patchwork)

colnames<-c('day','cases')
real_cases <- read.table("../foot_and_mouth/Data/processed_time_daily_cases.txt",col.names=colnames)
empirical <- read.table('../foot_and_mouth/Data/empirical_density_foot_and_mouth.txt',col.names=colnames)
fitempirical <- read.table('../foot_and_mouth/Data/fit_empirical_density_foot_and_mouth.txt',col.names=colnames)
library(gghighlight)


#Cases
ggplot(real_cases, aes(x=day),color='black')  +
  labs(x = 'time (days)', y = "daily new cases")+     
  geom_point(aes(y=cases),size=2,colour=(rgb(163/255,35/255,142/255))) + 
  gghighlight(day<84,unhighlighted_params = list(color='black'))+
  theme_bw()+
  theme(text=element_text(size=12, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))



#Cases fitted (not published)

gg<- ggplot(empirical, aes(x=day))  +
  labs(x = 'time', y = "empirical density")+     
  geom_point(aes(y=cases),size=2,colour=(rgb(163/255,35/255,142/255))) + 
  theme(axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                 size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  geom_line(data=fitempirical,aes(y=fitempirical$cases),size=2,colour="black")



colnames<-c('time','weib','gamma')

x_lower <- 0

x_upper <- 20

shape_gam = (4.97839683e+00)**2/1.08327439e+01
scale_gam = 1.08327439e+01/(4.97839683e+00) 

max_height <- max(dweibull(x_lower:x_upper, shape = 2.13623063e+00, scale=4.75098558e+00, log = FALSE))+0.1


#PDFs

ggplot(data.frame(x = c(x_lower, x_upper)), aes(x = x)) + xlim(x_lower, x_upper) + 
  ylim(0, max_height) +
  theme_bw()+
  stat_function(fun = dweibull, args = list(shape = 2.13623063e+00, scale=4.75098558e+00), size=2, aes(colour="orange")) + 
  stat_function(fun = dgamma, args = list(shape = shape_gam , scale=scale_gam),size=2, aes(colour="blue"))+
  labs(x = "\n time (days)", y = "pdf \n") + 
  theme(text=element_text(size=12, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.8,0.6))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  scale_colour_manual(breaks=c("orange","blue"), values=c("#0C7EF5","#F5A318"),labels=c("Infectiousness", "Removal time"), name="Distributions")




T_f = 80
day = real_cases$day-real_cases$day[1] #make time begin from day 0

day = head(day, T_f)
day = tail(day, -6)


cases = real_cases$cases
cases = head(cases, T_f)



mov_avg=apply(embed(cases, 7), 1, mean)

mov_avg_class<- data.frame(day,mov_avg)

differences=mov_avg-tail(cases, -6)
var_diff= var(differences)

colnames<-c('time','lower','central','upper')
conf_intervals <- read.table('../foot_and_mouth/Data/inflted_confidence_intervals.txt',col.names=colnames)


centralpoint= (conf_intervals$lower + conf_intervals$upper)/2
distance= (conf_intervals$upper - conf_intervals$lower)*sqrt(var_diff)


#Inflated confidence interval
ggplot(conf_intervals, aes(x=time,y=central)) + geom_line(aes(color='conf')) +
  labs(x = 'days since 2001-02-04', y ='cases')+     
  theme_bw()+
  theme(text=element_text(size=12, family="LM Roman 10"), axis.line = element_line(colour = rgb(163/255,35/255,142/255), 
                                                                                   size = 1, linetype = "solid", ),  panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="top")+  

  geom_ribbon(aes(ymin=pmax(centralpoint-0.5*distance,rep(0,1599)),ymax=centralpoint+0.5*distance, fill="#0C7EF5"),alpha=.5)+     
  geom_point(data=mov_avg_class,x=day, y=mov_avg,aes(size= "data"), color=rgb(163/255,35/255,142/255)) +
   scale_size_manual("", breaks=c("data"), values=c(1), labels=c("7-day moving average"))+
scale_fill_manual(values="#0C7EF5",labels="95 CI", name=NULL)+
  scale_colour_manual("", breaks=c("conf"), values=c('blue'), labels=c("median"))+
theme(legend.spacing.y = unit(0.00, 'cm'))+
  guides(size= guide_legend(order = 1),
        colour = guide_legend(order = 2), 
         shape = guide_legend(order = 3))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = scales::pretty_breaks(n = 5))+
  scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5),limits = c(0,80))


