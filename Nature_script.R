#librarys
###########################################################################################
library(lubridate)
library(dplyr)
library(dbplyr)
library(RPostgreSQL)
library(lamW)
library(minpack.lm)
library(nls.multstart)
library(boot)
library(tidyr)
library(purrr)
library(broom)
library(reshape2)

library(MASS)
library(ggplot2)
library(viridis)
library(grid)
library(latex2exp)
theme_set(theme_bw(base_size = 16))



get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

bootmean <- function(d, i) mean(d[i])
bootsd <- function(d, i) sd(d[i])


library(devtools)

###########################################################################################


# load data:
# exercises of marathon runners 
# user profiles of marathon runners

activities <- read.csv2("activities_used_180days.csv",stringsAsFactors = FALSE)
profiles <-  read.csv2("marathon_runners_prof.csv",stringsAsFactors = FALSE)

# group activities by runner_id and marathon_day, order by date and then by distance,

activities <- activities %>% group_by(h_user_id,marathon_day) %>% arrange(distance)

# summarize into training seasons, per runner and marathon

training <- summarise(activities,runs=n(),total_dist=sum(distance),total_time_min=sum(duration)/60000,mean_vel=mean(velocity),Mdist=last(distance),Mvel=last(velocity))

# exclude too long distance (>43000) and "new world records"
# lower limit for number of total runs in 6 months before marathon = 50

training <- training %>% filter(Mdist<43000,Mvel<42195/(121*60))
training <- training %>% filter(runs >= 50)

### 
### select RACE events ###
###

# consider only standard / typical race distances: 5k, 10k, HM=21.1k, M=42.195k
# with an error of +/- 3%

d_bins <- c(1000,0.97*5000,1.03*5000,0.97*10000,1.03*10000,0.97*21097.5,1.03*21097.5,0.97*42195,1.03*42195,100000)

d_labs <- c(0.,5000.,0.,10000.,0.,21097.5,0.,42195.,0.)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

activities$bin <- cut(activities$distance,breaks=d_bins,labels=d_labs) 
activities <- activities %>% mutate(vel_actual=1000*as.numeric.factor(bin)/duration)
activities_bins <- activities %>% group_by(h_user_id,marathon_day,bin)  %>% summarise(vel_max=max(vel_actual))

# remove races that are "new world record", and also 0 distance bin!

activities_bins <- activities_bins %>% filter((bin==5000 & vel_max < 6.60) | (bin==10000 & vel_max < 6.34) | (bin==21097.5 & vel_max < 6.03) | (bin==42195 & vel_max < 5.78))

races_by_vel <- semi_join(activities,activities_bins,by=c("h_user_id"="h_user_id","marathon_day"="marathon_day","vel_actual"="vel_max","bin"="bin")) %>% filter(bin!=0) %>%
  dplyr::select(h_user_id,marathon_day,bin,vel_actual,duration)

races_by_vel <- races_by_vel %>% group_by(h_user_id,marathon_day)

races_by_vel <- races_by_vel %>% arrange(h_user_id,marathon_day,desc(bin)) %>% 
  mutate(cur_max = ave(vel_actual, h_user_id, FUN=cummax))

races_by_vel <- races_by_vel %>% filter(vel_actual == cur_max) %>% dplyr::select(h_user_id,marathon_day,bin,vel_actual,duration) %>% mutate(time=duration/1000)

Real_M_times <- races_by_vel %>% filter(bin==42195) %>% 
summarise(Mdist=first(bin),Mtime=first(time)) %>% filter(Mtime>121*60+39)

# theoretical model fitting: select fixed tc = 6min, only training seasons with at least 2 races

time_fct = function(d, tc, vm, gl, gs) {(d < tc * vm) * (-d / (vm * gs * lambertWm1(-d * exp(-1 / gs) / (vm * tc *gs)))) + (d >= tc * vm) * (-d / (vm * gl * lambertWm1(-d * exp(-1 / gl) / (vm * tc * gl))))
}
Tc <- 360

# all races including marathon
races_by_vel <- races_by_vel %>% group_by(h_user_id,marathon_day)

# races WITHOUT marathon:

races_by_vel_woM <- races_by_vel %>% filter(bin != 42195)

races_by_vel <- races_by_vel %>% filter(n() > 1)
races_by_vel_woM <- races_by_vel_woM %>% filter(n() > 1)

races_by_vel <- races_by_vel %>% mutate(dist=as.numeric.factor(bin))
races_by_vel_woM <- races_by_vel_woM %>% mutate(dist=as.numeric.factor(bin))

model_fit <- races_by_vel_woM %>% nest() %>% 
mutate(fit = purrr::map(data, ~ nls_multstart(time ~ time_fct(d = dist, tc = 360, vm, gl, gs = 0.1),
                                             data = .,
                                             iter=50,
                                             start_lower = c(
                                               vm = 2, #* mean(.$velocity),
                                               gl = 0.05
                                             ),
                                             start_upper = c(
                                               vm = 7, #* mean(.$velocity),
                                               gl = 0.08
                                             ),
                                             upper = c(vm = 7, gl = 0.135),
                                             lower = c(vm = 2, gl = 0.039),
                                             supp_errors = 'N')))



model_fit_2 <- model_fit %>% mutate(n_short = purrr::map_dbl(data, ~ sum(.$time<=360)),
                                            n_long = purrr::map_dbl(data, ~ sum(.$time>360)),
                                            vm = purrr::map_dbl(fit, ~ broom::tidy(.,quick = TRUE)$estimate[1]),
                                            gl = purrr::map_dbl(fit, ~ broom::tidy(.,quick = TRUE)$estimate[2]))

model_fit <- model_fit[-c(11554), ]

model_fit_2 <- model_fit %>% mutate(n_short = purrr::map_dbl(data, ~ sum(.$time<=360)),
                                    n_long = purrr::map_dbl(data, ~ sum(.$time>360)),
                                    vm = purrr::map_dbl(fit, ~ broom::tidy(.,quick = TRUE)$estimate[1]),
                                    gl = purrr::map_dbl(fit, ~ broom::tidy(.,quick = TRUE)$estimate[2]))
  
model_fit_3 <- model_fit_2 %>% unnest(fit %>% purrr::map(broom::glance))

races_by_vel_pred <- model_fit %>% unnest(fit %>% purrr::map(broom::augment))
races_by_vel_pred <- races_by_vel_pred %>% group_by(h_user_id,marathon_day)
mean_race_errors_inperc <- races_by_vel_pred %>% summarise(mean_err_perc=100*sum(abs(.resid/time))/n())
model_fit_3 <- inner_join(model_fit_3,mean_race_errors_inperc,by=c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day"))

races_by_vel_pred <- model_fit %>% unnest(fit %>% purrr::map(broom::augment))
races_by_vel_pred <- races_by_vel_pred %>% group_by(h_user_id,marathon_day)
mean_race_errors_inperc <- races_by_vel_pred %>% summarise(mean_err_perc=100*sum(abs(.resid/time))/n())
model_fit_3 <- inner_join(model_fit_2,mean_race_errors_inperc,by=c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day"))

# we need at least two runs longer than time tc:

model_fit_3 <- model_fit_3 %>% filter(n_long>1)
model_fit_3_woM <- model_fit_3_woM %>% filter(n_long>1)

# remove results with endurance outside physiolgical range (fit converged to boundary values)
model_fit_real <- model_fit_3 %>% filter(gl>0.039,gl<0.135,vm<7)
model_fit_woM_real <- model_fit_3_woM %>% filter(gl>0.039,gl<0.135,vm<7)

# optional: select fits with 4 runs longer than time tc:
model_fit_woM_real <- model_fit_woM_real %>% filter(n_long>3)

sum(model_fit_real$n_long)  
sum(model_fit_woM_real$n_long)  

# number of users in model_fit_real:
length(unique(model_fit_real[["h_user_id"]]))

### save / load results

write.csv2(within(model_fit_real,rm(data,fit)), file="model_fit_real_woM_4races.csv", row.names=F, na="")
write.csv2(within(model_fit_woM_real,rm(data,fit)), file="model_fit_woM_real_v2.csv", row.names=F, na="")

model_fit_real <- read.csv2("model_fit_real_4races.csv")
model_fit_woM_real <- read.csv2("model_fit_real_woM_4races.csv")

###
### check results for consistency
###

# define function for converting vm to different tc>original tc:

vm_new = function(vm,tc,tcnew,gl) {vm*(1-gl*log(tcnew/tc))}

Mtime_selected <- races_by_vel %>% filter(distance>41500 & distance<42900 & 
                                          time>(2*60+29)*60 & time < (2*60+35)*60) %>% 
                  select(h_user_id,marathon_day,time)

Param_for_selected <- inner_join(Mtime_selected,model_fit_real,by=c("h_user_id", "marathon_day")) %>% select(h_user_id,marathon_day,time,vm,gl,mean_err_perc) %>% mutate(vm_JD = vm_new(vm,6*60,12.7*60,gl))

mean(exp(0.1/Param_for_selected$gl))

mean(Param_for_selected$vm_JD)
1000/mean(Param_for_selected$vm_JD)

mean(filter(Param_for_selected,exp(0.1/gl)>=6.7)$vm_JD)
1000/mean(filter(Param_for_selected,exp(0.1/gl)>=6.7)$vm_JD)


###
### ANALYSE ALL RESULTS (inclduing the ones without proper training as define above)
###

# histogram for mean error of race time prediction (n_long > 2)
# (note: for density mutiply count by 1/binwidth and 1/sum(count))

model_fit_real_2 <- model_fit_real %>% filter(n_long>2)

mean(filter(model_fit_real_2)$mean_err_perc)
median(filter(model_fit_real_2)$mean_err_perc)

# number of users in model_fit_real_2:
length(unique(model_fit_real_2[["h_user_id"]]))

# number of users in model_fit_real:
length(unique(model_fit_woM_real_2[["h_user_id"]]))

# histogram for error (only for more than two races since otherwise error is zero)

ggplot() + geom_histogram(data=filter(model_fit_real,n_long>2),aes(x=mean_err_perc,y=2*..count../sum(..count..)),binwidth = .5, color="white", fill="blue",alpha=.3)+ 
#geom_histogram(data=filter(model_fit_real,n_long>2),aes(x=mean_err_perc,y=2*..count../sum(..count..)),binwidth = .5, color="white", fill="red",alpha=.3)+
xlim(-1,12)+labs(y="probability density",x="mean error for race time [%]")+
theme(axis.text = element_text(size=18),axis.title = element_text(size=18))

# histogram for vm

ggplot() + 
geom_histogram(data=model_fit_real, aes(x=vm,y=5* ..count../sum(..count..)),binwidth = 0.2, color="white", fill="blue",alpha=0.3)+
#geom_histogram(data=model_fit_woM_real, aes(x=vm,y=5* ..count../sum(..count..)),binwidth = 0.2, color="white", fill="red",alpha=0.3)+
labs(y="probability density",x=expression('v'['m']*' [m/sec]'))+theme(axis.text = element_text(size=18),axis.title = element_text(size=18))+scale_x_continuous(breaks=seq(2,7,1))


# histogram for endurance 

ggplot() + 
geom_histogram(data=model_fit_real,aes(x=exp(0.1/gl),y=4* ..count../sum(..count..)),breaks=seq(1, 13, by = .25), color="white", fill="blue",alpha=0.3)+
#geom_histogram(data=model_fit_woM_real,aes(x=exp(0.1/gl),y=4* ..count../sum(..count..)),breaks=seq(1, 13, by = .25), color="white", fill="red",alpha=0.3)+
  labs(y="probability density",x=expression('endurance E'['l']))+theme(axis.text = element_text(size=18),axis.title = element_text(size=18))+scale_x_continuous(breaks=seq(1,13,2))+geom_line(data=model_fit_real, aes(x=exp(0.1/gl),y=1.98*exp(-.6*exp(0.1/gl))), colour="red") 


# correlation between vm and endurance?

# endurance E_l(v_m) function for fixed marathon time T_M:
El_vm <- function(TM,vm) {
  res <- (TM/Tc)^(0.1/(1-42195/(vm*TM)))
  return(res)
  }

# model parameters and real marathon time (time < 6h)
# use either model_fit_real (n_long>1) or model_fit_real_2 (n_long>2)

M_times_VS_vm_El <- Real_M_times %>% filter(Mdist == 42195, Mtime<360*60) %>% inner_join(.,model_fit_real,by = c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day")) %>% dplyr::select(h_user_id,marathon_day,Mtime,Mdist,vm,gl) %>% transmute(time_act=Mtime/60,time_com=time_fct(as.numeric.factor(Mdist),360,vm,gl,0.1)/60,diff_time=time_act-time_com,vm,El=exp(0.1/gl))

# model woM parameters and real marathon time (time < 6h, distance p/m 2%)

M_times_VS_vm_El <- Real_M_times %>% filter(Mdist == 42195, Mtime<360*60) %>% inner_join(.,model_fit_real_2,by = c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day")) %>% dplyr::select(h_user_id,marathon_day,Mtime,Mdist,vm,gl) %>% transmute(time_act=Mtime/60,time_com=time_fct(as.numeric.factor(Mdist),360,vm,gl,0.1)/60,diff_time=time_act-time_com,vm,El=exp(0.1/gl))

mean(M_times_VS_vm_El$diff_time)
mean(abs(M_times_VS_vm_El$diff_time))

# histogram for Mtime prediction error

M_times_VS_vm_El <- M_times_VS_vm_El %>% mutate(diff_rel=100*diff_time/time_act)

MT160 <- filter(M_times_VS_vm_El,time_act > 0, time_act <= 160)$diff_time
MT180 <- filter(M_times_VS_vm_El,time_act > 160, time_act <= 180)$diff_time
MT200 <- filter(M_times_VS_vm_El,time_act > 180, time_act <= 200)$diff_time
MT220 <- filter(M_times_VS_vm_El,time_act > 200, time_act <= 220)$diff_time
MT240 <- filter(M_times_VS_vm_El,time_act > 220, time_act <= 240)$diff_time

MT160r <- filter(M_times_VS_vm_El,time_act > 0, time_act <= 160)$diff_rel
MT180r <- filter(M_times_VS_vm_El,time_act > 160, time_act <= 180)$diff_rel
MT200r <- filter(M_times_VS_vm_El,time_act > 180, time_act <= 200)$diff_rel
MT220r <- filter(M_times_VS_vm_El,time_act > 200, time_act <= 220)$diff_rel
MT240r <- filter(M_times_VS_vm_El,time_act > 220, time_act <= 240)$diff_rel


bucket<-list(MT160=MT160,MT180=MT180,MT200=MT200,MT220=MT220,MT240=MT240)
bucket<-list(MT160r=MT160r,MT180r=MT180r,MT200r=MT200r,MT220r=MT220r,MT240r=MT240r)

ggplot() + 
  geom_histogram(data=melt(bucket), aes(x=value,fill=L1,y=0.2* ..count../sum(..count..)), breaks=seq(-60, 60, by = 5), color="white",alpha=1.,position = "stack")+
  labs(y="probability density",x=TeX("$\\Delta T_{m} \\, \\[min\\]$"))+theme(axis.text = element_text(size=18),axis.title = element_text(size=18))+scale_x_continuous(breaks=seq(-60,60,20)) + scale_fill_brewer(palette="Spectral",name=TeX("Actual $T_{m} \\, \\[min\\]$"),labels=c("..160","160..180","180..200","200..220","220..240"))

ggplot() + 
  geom_histogram(data=melt(bucket), aes(x=value,fill=L1,y=1* ..count../sum(..count..)), breaks=seq(-20, 20, by = 1), color="white",alpha=1.,position = "stack")+
  labs(y="probability density",x=TeX("$\\Delta T_{m} \\, \\[\\%\\]$"))+theme(axis.text = element_text(size=14),axis.title = element_text(size=16))+scale_x_continuous(breaks=seq(-20,20,5)) + scale_fill_brewer(palette="Spectral",name=TeX("Actual $T_{m} \\, \\[min\\]$"),labels=c("..160","160..180","180..200","200..220","220..240"))


#ggplot() + 
#  geom_histogram(data=filter(M_times_VS_vm_El,time_act < 240), aes(x=diff_time,y=0.2* ..count../sum#(..count..)), breaks=seq(-60, 60, by = 5), color="white", fill="blue",alpha=0.3)+
#  labs(y="probability density",x=expression(' [min]'))+theme(axis.text = element_text(size=18),axis.title = element_text(size=18))+scale_x_continuous(breaks=seq(-60,60,20))


bin_x <- 0:102
bins_x <- 2+0.05*bin_x
bin_y <- 0:122
bins_y <- 1+0.1*bin_y
M_times_VS_vm_El$bin_x <- cut(M_times_VS_vm_El$vm,breaks = bins_x,labels = FALSE)
M_times_VS_vm_El$bin_y <- cut(M_times_VS_vm_El$El,breaks = bins_y,labels = FALSE) 
M_times_VS_vm_El <- group_by(M_times_VS_vm_El,bin_x,bin_y)
M_times_VS_vm_El_binned <- M_times_VS_vm_El %>% summarise(n=n(),mean=mean(time_act),sd=sd(time_act))


# plot for mean act.time in each (vm,El)-bin

dat<-data.frame(x=2.025+0.05*M_times_VS_vm_El_binned$bin_x,y=1.05+0.1*M_times_VS_vm_El_binned$bin_y,z=M_times_VS_vm_El_binned$mean,n=M_times_VS_vm_El_binned$n)

ggplot(dat) + geom_point(aes(x, y, color = z),size = 1.,shape=15,alpha=1.) + 
  scale_color_distiller(direction = 1,palette = "Spectral", breaks = seq(130,360,by=40),name=TeX("$T_{marathon}\\,\\[min\\]$"))+ 
  scale_x_continuous(limits=c(2,7),breaks = seq(2,7,1))+
  scale_y_continuous(limits=c(2,14),breaks = seq(2,12,2))+
  ylab(TeX("Endurance $E_l$"))+
  xlab(TeX("velocity $v_m \\,\\[m/sec\\]"))+
  stat_function(fun=El_vm, n=1000, args = list(TM=120*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=140*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=160*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=180*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=200*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=220*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=240*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=270*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=300*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=330*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=360*60),
    colour="black",alpha=0.5,linetype="twodash",size=.5)+
  geom_text(x=6.6,y=13.8,label="2:00",angle=90)+
  geom_text(x=5.7,y=13.8,label="2:20",angle=90)+
  geom_text(x=5,y=13.8,label="2:40",angle=90)+
  geom_text(x=4.5,y=13.8,label="3:00",angle=90)+
  geom_text(x=4.1,y=13.8,label="3:20",angle=90)+
  geom_text(x=3.75,y=13.8,label="3:40",angle=90)+
  geom_text(x=3.4,y=13.8,label="4:00",angle=90)+
  geom_text(x=3.,y=13.8,label="4:30",angle=90)+
  geom_text(x=2.75,y=13.8,label="5:00",angle=90)+
  geom_text(x=2.5,y=13.8,label="5:30",angle=90)+
  geom_text(x=2.3,y=13.8,label="6:00",angle=90)
  
# plot for standard error of the mean in each (vm,El)-bin

dat<-data.frame(x=2.025+0.05*M_times_VS_vm_El_binned$bin_x,y=1.05+0.1*M_times_VS_vm_El_binned$bin_y,z=M_times_VS_vm_El_binned$sd/sqrt(M_times_VS_vm_El_binned$n))

ggplot(dat) + geom_point(aes(x, y, color = z),size = 1.,shape=15,alpha=1.) + 
  scale_color_distiller(direction = -1, palette = "RdBu",breaks = seq(0,60,by=10.),name="stand. error of the\nmean marathon\ntime [min]")+
  scale_x_continuous(limits=c(2,7),breaks = seq(2,7,1.))+
  scale_y_continuous(limits=c(2,13),breaks = seq(2,13,2))+
  ylab(TeX("Endurance $E_l$"))+
  xlab(TeX("velocity $v_m \\,\\[m/sec\\]"))

# plot for n in each (vm,El)-bin

dat<-data.frame(x=2.025+0.05*M_times_VS_vm_El_binned$bin_x,y=1.05+0.1*M_times_VS_vm_El_binned$bin_y,z=M_times_VS_vm_El_binned$n)

ggplot(dat) + geom_point(aes(x, y, color = z),size = 1.,shape=15,alpha=1.) + 
  scale_color_distiller(direction = 1, palette = "Oranges",breaks = seq(0,50,by=10),name=TeX("n"))+ 
  scale_x_continuous(limits=c(2,7),breaks = seq(2,7,1))+
  scale_y_continuous(limits=c(2,14),breaks = seq(2,12,2))+
  ylab(TeX("Endurance $E_l$"))+
  xlab(TeX("velocity $v_m \\,\\[m/sec\\]"))+
  stat_function(fun=El_vm, n=1000, args = list(TM=160*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=180*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=200*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=220*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=240*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=270*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=300*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=330*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  stat_function(fun=El_vm, n=1000, args = list(TM=360*60),
                colour="black",alpha=0.5,linetype="twodash",size=.5)+
  geom_text(x=5,y=13.8,label="2:40",angle=90)+
  geom_text(x=4.5,y=13.8,label="3:00",angle=90)+
  geom_text(x=4.1,y=13.8,label="3:20",angle=90)+
  geom_text(x=3.75,y=13.8,label="3:40",angle=90)+
  geom_text(x=3.4,y=13.8,label="4:00",angle=90)+
  geom_text(x=3.,y=13.8,label="4:30",angle=90)+
  geom_text(x=2.75,y=13.8,label="5:00",angle=90)+
  geom_text(x=2.5,y=13.8,label="5:30",angle=90)+
  geom_text(x=2.3,y=13.8,label="6:00",angle=90)


# plot for Mtime error (marathon time < 6h)

mean(M_times_VS_vm_El$diff_time)
mean(abs(M_times_VS_vm_El$diff_time))

bin_x <- 0:102
bins_x <- 2+0.05*bin_x
bin_y <- 0:122
bins_y <- 1+0.1*bin_y
M_times_VS_vm_El$bin_x <- cut(M_times_VS_vm_El$vm,breaks = bins_x,labels = FALSE)
M_times_VS_vm_El$bin_y <- cut(M_times_VS_vm_El$El,breaks = bins_y,labels = FALSE) 
M_times_VS_vm_El <- group_by(M_times_VS_vm_El,bin_x,bin_y)
# for absolute time difference:
M_times_VS_vm_El_binned <- M_times_VS_vm_El %>% summarise(n=n(),mean=mean(diff_time),sd=sd(diff_time))
# for relative time difference:
M_times_VS_vm_El_binned <- M_times_VS_vm_El %>% summarise(n=n(),mean=mean(diff_rel),sd=sd(diff_rel))
#bins_error <- c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
#bins_error <- c(-100,-10,-5,-3,-1,0,1,3,5,10,100)
bins_error <- c(-200,-20,-15,-10,-5,0,5,10,15,20,200)
M_times_VS_vm_El_binned$bin_error <- cut(M_times_VS_vm_El_binned$mean,breaks = bins_error,
#labels = c("< -4 min","-4 to -3 min","-3 to -2 min","-2 to -1 min","-1 to 0 min","0 to 1 min","1 to 2 min","2 to 3 min","3 to 4 min","> 4 min")) 
#labels = c("< -10 min","-10 to -5 min","-5 to -3 min","-3 to -1 min","-1 to 0 min","0 to 1 min","1 to 3 min","3 to 5 min","5 to 10 min","> 10 min")) 
labels = c("< -20 %","-20 to -15 %","-15 to -10 %","-10 to -5 %","-5 to 0 %","0 to 5 %","5 to 10 %","10 to 15 %","15 to 20 %","> 20 %")) 

dat<-data.frame(x=2.025+0.05*M_times_VS_vm_El_binned$bin_x,y=1.05+0.1*M_times_VS_vm_El_binned$bin_y,z=M_times_VS_vm_El_binned$bin_error,n=M_times_VS_vm_El_binned$n)
#dat$density <- get_density(dat$x, dat$y, n = 300)

ggplot(dat) + geom_point(aes(x, y, color = z),size=1,shape=15,alpha=1.) + 
  #scale_colour_brewer(direction=-1,palette = "RdYlBu",name=TeX("$\\Delta T_{marathon}$"))+ 
  scale_color_brewer(direction = 1,palette = "Spectral",name=TeX("$\\Delta T_{marathon}$"))+ 
  scale_size(range = c(1.5,7),breaks = seq(0,150,by=25))+
  scale_x_continuous(limits=c(2.0,7),breaks = seq(2.0,7,1.))+
  scale_y_continuous(limits=c(2,14),breaks = seq(2,12,2))+
  ylab(TeX("Endurance $E_l$"))+
  xlab(TeX("velocity $v_m \\,\\[m/sec\\]"))



# bootstrap curve for endurance vs. vm:

model_fit_El <- model_fit_real  %>% mutate(El = exp(0.1/gl)) 
dat<-data.frame(x=model_fit_El$vm,y=model_fit_El$El)
bin <- 0:32
bins <- 1.+0.2*bin 
dat$bin <- cut(dat$x,breaks = bins,labels = FALSE) 
dat_binned <- group_by(dat,bin) 

plot_dat <- dat_binned %>% filter(n()>7, !is.na(bin)) %>% summarise(n=n(),x_av=mean(x),y_av=mean(y),y_SEM=sd(boot(y, bootmean, R=1000, stype="i")$t),y_sd=sd(y),y_SEsd=sd(boot(y, bootsd, R=1000, stype="i")$t)) %>% filter(x_av>=0)

ggplot(plot_dat,aes(x=x_av))+
  geom_line(aes(y=y_av),color="red",size=.4)+
  geom_point(aes(y=y_av),color="red",size=3., shape=21, fill="white")+
  geom_ribbon(aes(ymin=y_av-y_sd,ymax=y_av+y_sd),fill="gray",alpha=0.5)+
  geom_ribbon(aes(ymin=y_av-y_SEM,ymax=y_av+y_SEM),fill="red",alpha=0.2)+
  geom_line(aes(y=y_av+y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av+y_sd-y_SEsd,ymax=y_av+y_sd+y_SEsd),fill="blue",alpha=0.2)+
  geom_line(aes(y=y_av-y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av-y_sd-y_SEsd,ymax=y_av-y_sd+y_SEsd),fill="blue",alpha=0.2)+
  scale_x_continuous(breaks=seq(0,7,1))+
  scale_y_continuous(breaks=seq(0,13,2))+
  ylab(TeX("endurance $E_l$"))+xlab(TeX("velocity $v_{m}\\, \\[m/sec\\]$"))+
  theme(axis.text = element_text(size=18),axis.title = element_text(size=18))



### 
### JOIN WITH TRAINING (when available)
###

# redefine training table without filters 
# add vm and sex to activities, and compute TRIMP

activities_2 <- inner_join(activities,model_fit_real,by = c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day")) %>% dplyr::select(h_user_id,marathon_day,distance,duration,velocity,hr_avg,vm)

activities_2 <- inner_join(activities_2,profiles,by = "h_user_id") %>% dplyr::select(h_user_id,marathon_day,distance,duration,velocity,hr_avg,vm,sex)

# duration is ms, hence devide by 60000 to get duration in min:
activities_2 <- activities_2 %>% mutate(trimp = if (first(sex)=="MALE") 
{duration/60000*0.64*velocity/vm*exp(1.92*velocity/vm)} else
{duration/60000*0.86*velocity/vm*exp(1.67*velocity/vm)}) %>% filter(!is.na(vm))

length(unique(activities_2[["h_user_id"]]))

training <- summarise(activities_2,runs=n(),total_dist=sum(distance),total_time_min=sum(duration)/60000,mean_vel=mean(velocity),vm=first(vm),mean_rel_vel=sum(duration/60000*velocity)/(total_time_min*vm),Mdist=last(distance),MHRav=last(hr_avg),Mvel=last(velocity),sex=first(sex),sum_trimp=sum(trimp),
#frac_fast=sum(trimp[velocity>0.98*vm])/sum_trimp
frac_fast=sum(distance[velocity>0.80*vm])/total_dist
#frac_fast=sum(duration[velocity>0.95*vm])/60000/total_time_min
)

# training <- training %>% filter(total_dist < 6000*1000)

# exclude too long distance (>43000) and "new world records"
# (lower limit for number of total runs in 6 months before marathon = 50)

# training <- training %>% filter(Mdist>41000,Mdist<43000,Mvel<42195/(121*60),vm<420/60)

# only model results from more than 1 long race, and different mean errors < 5%, at least 30 runs
# CHOOSE MALE or FEMALE or both:

test <- inner_join(model_fit_real,activities_2,by = c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day"))

test <- inner_join(test,training,by = c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day")) %>% filter(n_long>1,mean_err_perc<5,runs>30)
  
model_training <- inner_join(model_fit_real,training,by = c("h_user_id" = "h_user_id", "marathon_day" = "marathon_day")) %>% mutate(vm=vm.x,total_dist=total_dist/1000,El = exp(0.1/gl)) %>% filter(n_long>1,mean_err_perc<5,runs>30)#sex=="FEMALE")

length(unique(model_training[["h_user_id"]]))
sum(model_training$n_long)

# compare endurance and weighted mean training velocity / vm

dat<-data.frame(x=model_training$mean_rel_vel,y=model_training$El) %>% filter(!is.na(y),y<13) 

# compare endurance and total training trimp 

dat<-data.frame(x=model_training$sum_trimp,y=model_training$El) %>% filter(!is.na(y),y<13) 

# create distance bins, and compute mean endurance and SEM per bin

# for av.v / vm
bin <- 0:40
dist_bins <- 0.375+0.025*bin 

# for trimp
bin <- 0:25
dist_bins <- 0+2000*bin 

dat$bin <- cut(dat$x,breaks = dist_bins,labels = FALSE) 
dat_binned <- group_by(dat,bin) 

# SEM from bootstrap

El_vs_dist <- dat_binned %>% filter(n()>7, !is.na(bin)) %>% summarise(n=n(),x_av=mean(x),y_av=mean(y),y_SEM=sd(boot(y, bootmean, R=1000, stype="i")$t),y_sd=sd(y),y_SEsd=sd(boot(y, bootsd, R=1000, stype="i")$t)) #%>% filter(x_av>=.5, x_av<0.85)

# fit for av.v over vm only:
fit<-nlsLM(
  y_av ~ a*x_av*exp(b*x_av)+c,
  data = El_vs_dist,
  start = list( a=1.0, b=1.0,c=0.)
)
coef(summary(fit))
sigma(fit)

# for av.v over vm:
ggplot(El_vs_dist,aes(x=x_av))+
  geom_line(aes(y=y_av),color="red",size=.4)+
  geom_point(aes(y=y_av),color="red",size=3., shape=21, fill="white")+
  geom_ribbon(aes(ymin=y_av-y_sd,ymax=y_av+y_sd),fill="gray",alpha=0.5)+
  geom_ribbon(aes(ymin=y_av-y_SEM,ymax=y_av+y_SEM),fill="red",alpha=0.2)+
  geom_line(aes(y=y_av+y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av+y_sd-y_SEsd,ymax=y_av+y_sd+y_SEsd),fill="blue",alpha=0.2)+
  geom_line(aes(y=y_av-y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av-y_sd-y_SEsd,ymax=y_av-y_sd+y_SEsd),fill="blue",alpha=0.2)+
  scale_x_continuous(breaks=seq(0.5,1,0.1))+
  scale_y_continuous(breaks=seq(1,14,2))+
  ylab(TeX("Endurance $E_l$"))+xlab(TeX("Training intensity $\\bar{p}_{train} = \\frac{\\bar{v}_{train}}{v_m}$"))+
  theme(axis.text = element_text(size=18),axis.title = element_text(size=18))+
  stat_function(fun=function(x) 
    coef(fit)[1]*x*exp(coef(fit)[2]*x)+coef(fit)[3],
    colour="black",alpha=0.5,linetype="twodash",size=1.)+
  annotate("text",label=TeX("$E_l =  1.64 + 0.123\\,\\bar{p}_{train}\\,e^{5.02 \\,\\bar{p}_{train}}$"),x=.65,y=10,size=6)


# for trimp:
ggplot(El_vs_dist,aes(x=x_av/1000))+
  geom_line(aes(y=y_av),color="red",size=.2)+
  geom_point(aes(y=y_av),color="red",size=2., shape=21, fill="white")+
  geom_ribbon(aes(ymin=y_av-y_sd,ymax=y_av+y_sd),fill="gray",alpha=0.5)+
  geom_ribbon(aes(ymin=y_av-y_SEM,ymax=y_av+y_SEM),fill="red",alpha=0.2)+
  geom_line(aes(y=y_av+y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av+y_sd-y_SEsd,ymax=y_av+y_sd+y_SEsd),fill="blue",alpha=0.2)+
  geom_line(aes(y=y_av-y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av-y_sd-y_SEsd,ymax=y_av-y_sd+y_SEsd),fill="blue",alpha=0.2)+
  scale_x_continuous(breaks=seq(0,32,5))+
  scale_y_continuous(breaks=seq(1,14,2))+
  labs(y=TeX("Endurance $E_l$"),x=TeX("Training TRIMP \\[$\\times 10^3\\, \\]$"))+
  theme(axis.text = element_text(size=18),axis.title = element_text(size=18))

# compare vm and (avergane training speed / vm)

dat<-data.frame(x=model_training$mean_rel_vel,y=model_training$vm)

# compare vm and (distance @ fast speed > 0.8vm / total distance)

dat<-data.frame(x=model_training$frac_fast,y=model_training$vm)

# compare vm and total distance

dat<-data.frame(x=model_training$total_dist,y=model_training$vm)  

# create  bins, and compute mean vm and SEM per bin

# for mean_rel_vel / frac dist fast
bin <- 0:40
bins <- 0+0.025*bin 

# for total_dist
bin <- 0:25
bins <- 200+300*bin 

dat$bin <- cut(dat$x,breaks = bins,labels = FALSE) 
dat_binned <- group_by(dat,bin) 

# SEM from bootstrap

plot_dat <- dat_binned %>% filter(n()>7, !is.na(bin)) %>% summarise(n=n(),x_av=mean(x),y_av=mean(y),y_SEM=sd(boot(y, bootmean, R=1000, stype="i")$t),y_sd=sd(y),y_SEsd=sd(boot(y, bootsd, R=1000, stype="i")$t)) %>% filter(x_av>=.0)

# for av.v over vm:
ggplot(plot_dat,aes(x=x_av))+
  geom_line(aes(y=y_av),color="red",size=.2)+
  geom_point(aes(y=y_av),color="red",size=2., shape=21, fill="white")+
  geom_ribbon(aes(ymin=y_av-y_sd,ymax=y_av+y_sd),fill="gray",alpha=0.5)+
  geom_ribbon(aes(ymin=y_av-y_SEM,ymax=y_av+y_SEM),fill="red",alpha=0.2)+
  geom_line(aes(y=y_av+y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av+y_sd-y_SEsd,ymax=y_av+y_sd+y_SEsd),fill="blue",alpha=0.2)+
  geom_line(aes(y=y_av-y_sd),color="gray",size=.5)+
  geom_ribbon(aes(ymin=y_av-y_sd-y_SEsd,ymax=y_av-y_sd+y_SEsd),fill="blue",alpha=0.2)+
  scale_x_continuous(limits = c(400,3700),breaks=seq(0,4000,1000))+
  #scale_x_continuous(limits = c(0.0,1.), breaks=seq(0.3,1.,.1))+
  scale_y_continuous(limits = c(3.,6.5),breaks=seq(3.,6.5,.5))+
  labs(y=TeX("$v_m \\, \\[m/sec\\]$"),x=TeX("$d_{train}\\, \\[km\\]$"))+
  #labs(y=TeX("$v_m \\, \\[m/sec\\]$"),x=TeX("Training intensity$\\, \\bar{p}_{train}=\\bar{v}_{train}\\,/\\,v_m$"))+
  #labs(y=TeX("$v_m \\, \\[m/sec\\]$"),x=TeX("$(\\, d_{train} \\, at \\, \\bar{v}_{train} > 0.8 v_m \\,) / d_{train}$"))+
  theme(axis.text = element_text(size=18),axis.title = element_text(size=18))
  

############################  END  #######################


