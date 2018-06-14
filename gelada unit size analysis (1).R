## load required packages
require(plyr);require(dplyr);require(lme4); require(lmerTest); require(MuMIn); require(ggplot2);require(ggplot2);require(ggthemes);require(ggbeeswarm)
theme_set(theme_minimal(base_size=20))

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


## load in data
#load(url("https://github.com/GeladaResearchProject/Tinsley_Johnson_2018_Behav_Ecol/blob/master/seasonality_in_geladas_data.RData?raw=true"))

load("~/Documents/UMICH/Gelada Lab/5. Projects/3. Unit Size/Data/gelada_unit_size_data.RData")



## There are 12 R objects:
# repro_data: 

# repro_data_3mo:

# repro_data_6mo:

# repro_plot:

# repro_plot_summary:

# mortality_data:

# takeover_data:

# infanticide_data:

# feeding_data:

# gc_data:

# grooming_data:

# clustering_data:



#############################################################################
###                 Reproductive Performance                              ###
#############################################################################

## 1-month windows ----

mod1<-glmer(prod_mo~poly(AGE,2)+(max_males)+poly(max_females,2)+(1|anon)+(1|mth)+(1|year),data=repro_data,family="binomial", control=glmerControl(optimizer="bobyqa")) 
summary(mod1)

## 3-month windows ----

mod2<-glmer(cbind(prod_mos,total_mos-prod_mos)~poly(age,2)+(max_male)+poly(max_females,2)+(1|anon)+(1|cum_3mo),data=repro_data_3mo,family="binomial")
summary(mod2)


## 6-month windows ----

mod3<-glmer(cbind(prod_mos,total_mos-prod_mos)~poly(age,2)+max_male+poly(max_females,2)+(1|anon)+(1|cum_6mo),data=repro_data_6mo,family="binomial")
summary(mod3)


## PANEL 2a & 2b: REPRODUCTIVE PERFORMANCE (6 mo windows) ----

ggplot(repro_plot, aes(x=as.factor(max_females),y=resid_prod_mos, fill=unit.cat)) + 
  geom_beeswarm(size=2,alpha=0.75,color="black",shape=21,priority='ascending',cex=0.5, dodge.width=.1) +
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("residual productive months")+
  scale_y_continuous(breaks=c(-6,-3,0,3,6),labels=c(-6,"",0,"",6)) +
  xlab("")+
  stat_summary(fun.data=data_summary, shape = 23, size = 0.75, color = "black", fill = "black", position = position_nudge(x = 0.2))+
  geom_smooth(aes(x=as.integer(max_females), y=resid_prod_mos),method=lm,formula=y ~ poly(x,2),color="#5a5a5a", fill="#cccccc", alpha = 0.25)+
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))


ggplot(repro_plot_summary, aes(x=as.factor(unit.cat),y=prod_ratio, fill=unit.cat)) + 
  geom_hline(yintercept=0.5247238, color="#5a5a5a", linetype="dotted",size=0.8)+
  geom_errorbar(aes(ymin=prod_ratio-stdev, ymax=prod_ratio+stdev), width=.1, position=position_dodge(0.5)) +
  geom_point(shape=21,size=5, position=position_dodge(0.5), fill=c("#993404","#ec7014","#fec44f"), alpha=0.95) +
  ylab("productive months/\ntotal months observed")+
  scale_y_continuous(breaks=c(0.50,0.55,0.60),labels=c(" 0.50",""," 0.60")) +
  xlab("")+
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))



#############################################################################
###                       Female Mortality                                ###
#############################################################################

mod4<-glmer(cbind(n_deaths,(n_females-n_deaths)) ~ poly(max_fem,2)+(mean_age_yrs)+offset(n_fem_yrs) + (1|unit)+(1|year), family="binomial", data=mortality_data, control=glmerControl(optimizer="bobyqa")) 
summary(mod4)


## PANEL 2c & 2d: FEMALE MORTALITY ----

ggplot(mortality_data, aes(x=as.factor(max_fem),y=n_deaths/n_fem_yrs, fill=unit.cat)) + 
  geom_jitter(data = mortality_data[mortality_data$n_deaths>0,], aes(fill=unit.cat),shape=21,size=3, height=0, alpha=0.8, width=0.2) +
  stat_summary(data = mortality_data[mortality_data$n_deaths>0,], fun.data=data_summary, shape = 23, size = 0.75, color = "black", fill = "black", position = position_nudge(x = 0.25))+
  geom_smooth(data = mortality_data[mortality_data$n_deaths>0,], aes(x=as.integer(max_fem), y=n_deaths/n_fem_yrs),method=lm,formula=y ~ poly(x,2),color="#5a5a5a", fill="#cccccc")+
  geom_rug(data = mortality_data[mortality_data$n_deaths==0,], aes(x=max_fem), position = "jitter", sides="b", alpha=0.5) +
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("deaths per female-year")+
  scale_y_continuous(breaks=c(0,1.5,3,4.5,6),labels=c(0,"   ",3,"   ",6)) +
  xlab("number of females in unit")+
  theme(plot.title = element_text(hjust = 0,face = "bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))


mortality_data_summary<-ddply(mortality_data,. (unit.cat), summarize, ttl_fem = sum(n_females), ttl_fem_yrs=sum(n_fem_yrs), ttl_deaths = sum(n_deaths), death_rate = ttl_deaths/ttl_fem_yrs, stdev = sqrt(death_rate*(abs(1-death_rate))/ttl_fem_yrs), avg_age=mean(mean_age/365), stdev_age=sd(mean_age/365))


ggplot(mortality_data_summary, aes(x=as.factor(unit.cat),y=death_rate, fill=unit.cat)) + 
  geom_hline(yintercept=0.08275989, color="#5a5a5a", linetype="dotted",size=0.8)+
  geom_errorbar(aes(ymin=death_rate-stdev, ymax=death_rate+stdev), width=.1, position=position_dodge(0.5)) +
  geom_point(shape=21,size=5, position=position_dodge(0.5), fill=c("#993404","#ec7014","#fec44f"), alpha=0.95) +
  ylab("\ndeaths per female-year")+
  scale_y_continuous(breaks=c(0,0.03,0.06,0.09,0.12),labels=c(0,"",0.06,"",0.12)) +
  xlab("unit size category")+
  theme(plot.title = element_text(hjust = 0,face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))


#############################################################################
###                       Takeover Rates                                  ###
#############################################################################

mod5<-glmer(ttl_takeovers ~ scale(max_fem) + scale(max_males) + offset(ttl_yrs) + (1|unit), data=takeover_data, family="poisson")
summary(mod5)


ggplot(takeover_data, aes(x=as.factor(unit.cat),y=takeover_rate, fill=unit.cat)) + 
  geom_hline(yintercept=0.5498373, color="grey", linetype="dotted")+
  geom_jitter(aes(fill=unit.cat),shape=21,size=3, height=0, alpha=0.8, width=0.2) +
  stat_summary(fun.data=data_summary, shape = 23, size = 0.8, color = "black", fill = "black", position = position_nudge(x = 0.35))+
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("takeovers per unit-year\n")+
  scale_y_continuous(breaks=c(0,3,6,9,12),labels=c(0,"   ",6,"   ",12)) +
  xlab("unit size category")+
  theme(plot.title = element_text(hjust = 0,face = "bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))




#############################################################################
###                       Infanticide Rates                               ###
#############################################################################

mod6<-glmer(cbind(n_infanticides,(n_births-n_infanticides)) ~ poly(max_fem,2) + offset(n_fem_yrs) + (1|unit)+(1|year),family="binomial",data=infanticide_data)
summary(mod6)


ggplot(infanticide_data[infanticide_data$n_births>0,], aes(x=as.factor(unit.cat),y=infanticide_rate, fill=unit.cat)) + 
  geom_hline(yintercept=0.09566198, color="grey", linetype="dotted")+
  geom_jitter(aes(fill=unit.cat),shape=21, size=4, height=0, alpha=0.8, width=0.2) +
  scale_size_continuous(range = c(2, 9), name="Total infant deaths") +
  stat_summary(fun.data=data_summary, shape = 23, size = 0.8, color = "black", fill = "black", position = position_nudge(x = 0.35))+
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("infanticides per birth")+
  scale_y_continuous(breaks=c(0,0.3,0.6,0.9,1.2),labels=c(0,"   ",0.6,"   ",1.2)) +
  xlab("unit size category")+
  guides(size = guide_legend(override.aes = list(fill = "grey"))) +
  theme(plot.title = element_text(hjust = 0,face = "bold"), legend.text=element_text(size=11),legend.title=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "top", legend.background = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))


#############################################################################
###                       Feeding Rates                                   ###
#############################################################################

mod7<-lmer(ttl_f ~ scale(I(max_fem^2)) + scale(max_male) + scale(I(temp^2)) + scale(max_age) + repro + offset(ttl_focal) + (1|anon) + (1|month) + (1|year) + (1|unit), data=feeding_data)
summary(mod7)


mod7a<-lmer(ttl_f ~ unit.cat + scale(max_male) + scale(I(temp^2)) + scale(max_age) + repro + offset(ttl_focal) + (1|anon) + (1|month) + (1|year) + (1|unit), data=feeding_data)
summary(mod7a)

exp(4.0053)-1 #53.88829


ggplot(feeding_data, aes(x=as.factor(max_fem),y=ttl_f/focalhrs, fill=unit.cat)) + 
  geom_beeswarm(alpha=0.75,color="black",shape=21,size=2,priority='ascending',cex=0.35, dodge.width=.1) +
  stat_summary(fun.data=data_summary, shape = 23, size = 0.75, color = "black", fill = "black", position = position_nudge(x = 0.2))+
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("feeding rate\n(minutes per focal hour)")+
  scale_y_continuous(limits=c(0,60),breaks=c(0,15,30,45,60),labels=c(0,"   ",30,"   ",60)) +
  xlab("number of females in unit")+
  geom_smooth(aes(x=as.integer(max_fem), y=ttl_f/focalhrs),method=lm,formula=y ~ poly(x,2),color="#5a5a5a", fill="#cccccc")+
  theme(plot.title = element_text(hjust = 0,face = "bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))



#############################################################################
###                       Glucocorticoids                                 ###
#############################################################################

mod8<-lmer(log(GC) ~ poly(max_females,2) + scale(max_males) + scale(maxt_1) + scale(mint_1) + to_yn * repro + poly(AGE,2) * repro + poly(AGE,2) * repro + (1|anon)+(1|Month)+(1|UNIT)+(1|Year),data=gc_data)
summary(mod8)


ggplot(gc_data, aes(x=as.factor(unit.cat),y=GC, fill=unit.cat)) + 
  geom_hline(yintercept=25.10794, color="grey",linetype="dotted")+
  geom_beeswarm(alpha=0.75,color="black",shape=21,size=2,priority='ascending',cex=0.35, dodge.width=.1) +
  stat_summary(fun.data=data_summary, shape = 23, size = 0.9, color = "black", fill = "black", position = position_nudge(x = 0.3))+
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("glucocorticoid metabolites (ng/g)")+
  scale_y_continuous(breaks=c(0,15,30,45,60,75,90),labels=c(0,"   ",30,"   ",60,"   ",90)) +
  xlab("unit size category")+
  theme(plot.title = element_text(hjust = 0,face = "bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))



#############################################################################
###                       Grooming Time                                   ###
#############################################################################

mod9<-lmer(g_rate~poly(max_females,2)+scale(max_males)+scale(age)+(1|year)+(1|anon)+(1|id1unit),data=grooming_data)
summary(mod9)


ggplot(grooming_data,aes(x=as.factor(max_females),y=g_rate, fill=unit.cat))+
  geom_beeswarm(alpha=0.75,color="black",shape=21,size=3,priority='ascending',cex=0.75, dodge.width=.1) +
  stat_summary(fun.data=data_summary, shape = 23, size = 0.75, color = "black", fill = "black", position = position_nudge(x = 0.2))+
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("grooming rate\n(minutes per focal hour)")+
  scale_y_continuous(breaks=c(0,4,8,12,16),labels=c(0,"    ",8,"    ",16)) +
  xlab("")+
  geom_smooth(aes(x=as.integer(as.factor(max_females)), y=g_rate), method=lm, formula=y ~ x, color="#5a5a5a", fill="#cccccc")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))


#############################################################################
###                       Grooming Partners                               ###
#############################################################################

mod10<-lmer(g_partners~poly(max_females,2)+scale(max_males)+scale(age)+(1|year)+(1|anon)+(1|id1unit),data=grooming_data)
summary(mod10)

ggplot(grooming_data, aes(x=as.factor(max_females),y=g_partners, fill=unit.cat)) +
  geom_jitter(aes(fill=unit.cat),shape=21,size=3, height=0, alpha=0.8, width=0.2) +
  stat_summary(fun.data=data_summary, shape = 23, size = 0.75, color = "black", fill = "black", position = position_nudge(x = 0.2))+
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("female grooming partners\n")+
  scale_y_continuous(breaks=c(0,3,6),labels=c(0,"  3",6)) +
  xlab("number of females in unit")+
  scale_x_discrete(breaks = c(2:13))+
  geom_smooth(aes(x=as.integer(max_females-0.5), y=g_partners),method=lm,formula=y ~ poly(x,2),color="#5a5a5a", fill="#cccccc")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))  


#############################################################################
###                   Clustering Coefficient                              ###
#############################################################################

mod11<-lmer(mean_clust~poly(max_females,2)+(1|unit)+(1|yr), data=clustering_data)
summary(mod11)
 

ggplot(clustering_data, aes(x=as.factor(unit.cat),y=mean_clust, fill=unit.cat)) + 
  geom_hline(yintercept=0.3139578, color="grey", linetype="dotted")+
  geom_jitter(aes(fill=unit.cat),shape=21,size=3, height=0, alpha=0.8, width=0.2) +
  stat_summary(fun.data=data_summary, shape = 23, size = 0.75, color = "black", fill = "black", position = position_nudge(x = 0.25))+
  scale_fill_manual(aes(fill=unit.cat),values=c("#993404","#ec7014","#fec44f"),guide=FALSE) +
  ylab("network average\nclustering coefficient")+
  scale_y_continuous(breaks=c(0,0.25,0.50,0.75,1),labels=c(0,"",0.5,"",1)) +
  xlab("unit size category")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size=0.25), axis.ticks = element_line(size=0.25))

