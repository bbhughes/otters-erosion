###Data analysis for testing sea otter impacts on crabs, salt marsh and tidal creek bank erosion####
#By Brent B. Hughes et al. 2022#

#Load libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(glmmTMB)
library(sdmTMB)

####Before/After Comparison during the most recent otter expansion in Elkhorn Slough####
####GLMM tweedie analysis of sea otters in creeks, crab consumption, and erosion across 13 creeks from 2011-2017####
before_after_sea_otter <- read.csv("before_after_sea_otter.csv")
str(before_after_sea_otter)

#convert year to factor
before_after_sea_otter$Year <- as.factor(as.character(before_after_sea_otter$Year))

#Sea Otters
#subset data to run B/A tweedie analysis
baso_2012_17 <- subset(before_after_sea_otter, Year=="2012"|Year=="2017")
str(baso_2012_17)

#Run tweedie analysis
library(glmmTMB)
fit_otter <- glmmTMB(
  Otter_per_ha ~ Year + (1 | Creek),
  data = baso_2012_17,
  #dispformula = ~ Year,
  family = tweedie(),
  #verbose = TRUE
)
summary(fit_otter)
#dispersion parameter = 0.0844, df residuals = 21, p-value < 0.0005, significant difference between 2012 and 2017.
r_otter <- DHARMa::simulateResiduals(fit_otter)
plot(r_otter)

#figure of otters in creeks in 2012 and 2017
min_y<-min(0)
max_y<-max(0.60)
baso_plot<-ggplot(data=baso_2012_17, aes(x=Year,y=Otter_per_ha, fill=Year))+
  theme_bw(base_size=20)+theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_boxplot()+
  geom_point(alpha = 0.5, position = position_jitter(w=0.02, h=0.00), size=4) +
  annotate("text",x=0.75, y=0.55, label= "***", size=14)+
  ggtitle('A')+theme(plot.title=element_text(hjust=0))+
  labs(y = "Sea otters per ha of tidal creek\n", x = "")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
baso_plot

#Pachygrapsus consumed per day
fit_pachy <- glmmTMB(
  Pachy_consumed_per_day ~ Year + (1 | Creek),
  data = baso_2012_17,
  family = tweedie()
  # verbose = TRUE
)
summary(fit_pachy)
#dispersion parameter = 0.674, df residuals = 21, p-value < 0.0005, significant difference between 2012 and 2017.
r_pachy <- DHARMa::simulateResiduals(fit_pachy)
plot(r_pachy)

#Befor-After figure of Pachygrapsus crab consumed
min_y<-min(0)
max_y<-max(45)
baso_pachy_plot<-ggplot(data=before_after_sea_otter, aes(x=Year,y=Pachy_consumed_per_day, fill=Year))+
  theme_bw(base_size=20)+theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_boxplot()+
  geom_point(alpha = 0.5, position = position_jitter(w=0.02, h=0.00), size=4) +
  annotate("text",x=0.75, y=40, label= "***", size=14)+
  ggtitle('B')+theme(plot.title=element_text(hjust=0))+
  labs(y = "Crab consumed per ha per day\n", x = "")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())
baso_pachy_plot

#Before After comparison of annual erosion in 13 creek banks
#load data
before_after_erosion <- read.csv("before_after_erosion.csv")
str(before_after_erosion)

#convert year to factor
before_after_erosion$Year <- as.factor(as.character(before_after_erosion$Year))

#Run LMM
fit_erosion <- glmmTMB(
  erosion_m_yr ~ Year + (1 | Creek),
  data = before_after_erosion,
  family = gaussian()
  # verbose = TRUE
)
summary(fit_erosion)
#dispersion parameter for gaussian (sigma^2) = 0.0076, df residuals = 22, p-value < 0.0005, significant difference between 2011 and 2017.
r_pachy <- DHARMa::simulateResiduals(fit_pachy)
plot(r_pachy)

#figure of erosion in 2011 and 2017 at 13 creeks
min_y<-min(-0.3)
max_y<-max(0.60)
ba_erosion_plot<-ggplot(data=before_after_erosion, aes(x=Year,y=erosion_m_yr, fill=Year))+
  theme_bw(base_size=20)+theme(legend.position="none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_boxplot()+
  geom_point(alpha = 0.5, position = position_jitter(w=0.02, h=0.00), size=4) +
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.75, y=0.55, label= "***", size=14)+
  ggtitle('C')+theme(plot.title=element_text(hjust=0))+
  labs(y = "Bank erosion (m per yr)\n", x = "Year") +
  scale_x_discrete(labels = c("Before Otter\nExpansion", "After Otter\nExpansion"))
ba_erosion_plot

#Save Fig. Before After
pdf("ba_plot.pdf", width=6, height=13)
gridExtra::grid.arrange(baso_plot, baso_pachy_plot, ba_erosion_plot, ncol = 1, nrow = 3)
dev.off()


####TRANSECTS, No assumptions of when sea otters occupied tidal creeks prior to 2013.####
#Crab differences in marsh transects between 2014 and 2016####
crab_transect <- read.csv("otter_marsh_transects.csv")
str(crab_transect)
crab_transect$Year <- as.factor(as.integer(crab_transect$Year))
crab_transect$Otters_ha <- as.numeric(as.character(crab_transect$Otters_ha))
crab_transect$bg_biomass_2015<-as.numeric(as.character(crab_transect$bg_biomass_2015))
crab_transect$bank_retreat_2013_2015<-as.numeric(as.character(crab_transect$bank_retreat_2013_2015))
str(crab_transect)

#Relationship between sea otters and crabs in marsh for 2014
crab_marsh14 <- subset(crab_transect, Year=="2014")
str(crab_marsh14)

#GLM model using log link
crab_marsh14_glm1 <- glm(Density ~ usgs_oph_2015, data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))
crab_marsh14_glm2 <- glm(Density ~ log(usgs_oph_2015), data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))
AIC(crab_marsh14_glm1, crab_marsh14_glm2)
summary(crab_marsh14_glm1)

#P = 0.0024

min_x <- min(crab_marsh14$usgs_oph_2015)
max_x <- max(crab_marsh14$usgs_oph_2015)
min_y<-min(0)
max_y<-max(5)
crab_otter14_SI_plot<-ggplot(data=crab_marsh14, aes(x=usgs_oph_2015, y=Density))+
  theme_bw(base_size=20)+theme(legend.position="none")+
  geom_point(alpha = 0.5, position = position_jitter(w=0, h=0.02), size=4) +
  geom_smooth(method = "glm", formula = y~x, se = TRUE, colour='black',
              size=1.5, alpha = 0.3, method.args = list(family = Gamma(link = "log"))) +
  scale_x_continuous(limits=c(min_x,max_x)) +
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.05, y=4.5, label= "*", size=14)+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(y = "Crab density per trap\n", x = "\nOtters per ha")
crab_otter14_SI_plot

#Comparison of otter abundance and Pachygrapsus consumed per year
#GLM model using log tranformation
pachy_otter_glm <- glm(usgs_pachy_eaten_ha_yr ~ usgs_oph_2015, data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))
summary(pachy_otter_glm)

pachy_otter_glm2 <- glm(usgs_pachy_eaten_ha_yr ~ log(usgs_oph_2015), data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))

AIC(pachy_otter_glm, pachy_otter_glm2)

summary(pachy_otter_glm2)
#P=0.00346

#Plot the data
pachy_otter_plot<-ggplot(data=crab_marsh14, aes(x=usgs_oph_2015, y=(usgs_pachy_eaten_ha_yr/1000)))+
  theme_bw(base_size=20)+theme(legend.position="none")+
  geom_point(alpha = 0.5, position = position_jitter(w=0, h=0.02), size=4) +
  geom_smooth(method = "glm", formula = y~log(x), se = TRUE,colour='black', size=1.5, alpha = 0.3, method.args = list(family = Gamma(link = "log"), na.action=na.omit)) +
  scale_x_continuous(limits=c(min_x,max_x)) +
  annotate("text",x=0.005, y=5.5, label= "**", size=14)+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ggtitle('A')+theme(plot.title=element_text(hjust=0))+
  labs(y = "1000s shore crab consumed/ha/yr\n", x = "\nOtters per ha")
pachy_otter_plot

#Comparison of bank retreat and Pachygrapsus consumed per year,
#GLM model using log tranformation
# Sean TODO: log(usgs_pachy_eaten_ha_yr) or usgs_pachy_eaten_ha_yr?
pachy_retreat_glm <- glm(bank_retreat_2013_2015 ~ usgs_pachy_eaten_ha_yr, data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))
summary(pachy_retreat_glm)

pachy_retreat_glm2 <- glm(bank_retreat_2013_2015 ~ log(usgs_pachy_eaten_ha_yr), data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))

AIC(pachy_retreat_glm, pachy_retreat_glm2)
summary(pachy_retreat_glm)
#P= 0.00953

#Plot the data
min_y <- min(crab_marsh14$bank_retreat_2013_2015)
max_y <- max(90)
min_x<-min(0)
max_x<-max(4)
pachy_retreat_plot<-ggplot(data=crab_marsh14, aes(x=(usgs_pachy_eaten_ha_yr/1000), y=bank_retreat_2013_2015))+
  theme_bw(base_size=20)+theme(legend.position="none")+
  geom_point(alpha = 0.5, position = position_jitter(w=0.02, h=0.02), size=4) +
  geom_smooth(method = "glm", formula = y~x, se = TRUE,colour='black', size=1.5, alpha = 0.3, method.args = list(family = Gamma(link = "log"))) +
  scale_x_continuous(limits=c(min_x,max_x)) +
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.5, y=82, label= "*", size=14)+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ggtitle('C')+theme(plot.title=element_text(hjust=0))+
  labs(y = "Marsh edge retreat (cm/yr)\n", x = "\n1000s shore crab consumed/ha/yr")
pachy_retreat_plot

#Comparison of marsh biomass and Pachygrapsus consumed per year,
#GLM model using log tranformation
pachy_bg_glm <- glm(bg_biomass_2015~usgs_pachy_eaten_ha_yr, data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))
summary(pachy_bg_glm)
#P = 0.139

pachy_bg_glm2 <- glm(bg_biomass_2015~log(usgs_pachy_eaten_ha_yr), data=crab_marsh14, na.action=na.omit, family = Gamma(link = "log"))
AIC(pachy_bg_glm, pachy_bg_glm2)
summary(pachy_bg_glm2)
#P = 0.0322

#Plot the data
min_y <- min(crab_marsh14$bg_biomass_2015)
max_y <- max(30)
min_x<-min(0)
max_x<-max(4)
pachy_bg_plot<-ggplot(data=crab_marsh14, aes(x=(usgs_pachy_eaten_ha_yr/1000), y=bg_biomass_2015))+
  theme_bw(base_size=20)+theme(legend.position="none")+
  geom_point(alpha = 0.5, position = position_jitter(w=0, h=0.02), size=4) +
  geom_smooth(method = "glm", formula = y~log(x), se = TRUE,colour='black', size=1.5, alpha = 0.3, method.args = list(family = Gamma(link = "log"))) +
  scale_x_continuous(limits=c(min_x,max_x)) +
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.5, y=28, label= "*", size=14)+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ggtitle('B')+theme(plot.title=element_text(hjust=0))+
  labs(y = "Marsh belowground mass (g/core)\n", x = "\n1000s shore crab consumed/ha/yr")
pachy_bg_plot

#Save Spatial Analysis Figure
pdf("Spatial_Figure.pdf", width=7, height=15)
gridExtra::grid.arrange(pachy_otter_plot,pachy_bg_plot,pachy_retreat_plot, ncol = 1, nrow = 3)
dev.off()

####Experimental Data analyses####
#Load data for t-tests
Otter_marsh_data<-read.table(file="Otters_Marsh_2015_R.csv", header=TRUE, sep=",")
str(Otter_marsh_data)
Otter_marsh_data$AGMass_OE_kg_msq <- as.numeric(as.character(Otter_marsh_data$AGMass_OE_kg_msq))

# Otter_marsh_data$AG_Consum_d <- as.numeric(as.character(Otter_marsh_data$AG_Consum_d))
# Otter_marsh_data$BG_Consum_d <- as.numeric(as.character(Otter_marsh_data$BG_Consum_d))
Otter_marsh_data$AG_Consum_Crab <- as.numeric(as.character(Otter_marsh_data$AG_Consum_Crab))
Otter_marsh_data$AG_Consum_NoCrab <- as.numeric(as.character(Otter_marsh_data$AG_Consum_NoCrab))
Otter_marsh_data$BG_Consum_Crab <- as.numeric(as.character(Otter_marsh_data$BG_Consum_Crab))
Otter_marsh_data$BG_Consum_NoCrab <- as.numeric(as.character(Otter_marsh_data$BG_Consum_NoCrab))
Otter_marsh_data$biomass_consumed <- as.numeric(as.character(Otter_marsh_data$biomass_consumed))

####Feeding Assay####
#Aboveground feeding trial
var.test(Otter_marsh_data$AG_Consum_Crab, Otter_marsh_data$AG_Consum_NoCrab, data=two)
t.test(Otter_marsh_data$AG_Consum_Crab, Otter_marsh_data$AG_Consum_NoCrab, data=two, var.equal=T)
#t = -0.30593, df = 18, p-value = 0.7632, not significant

#Belowground feeding trial
var.test(Otter_marsh_data$BG_Consum_Crab, Otter_marsh_data$BG_Consum_NoCrab, data=two)
t.test(Otter_marsh_data$BG_Consum_Crab, Otter_marsh_data$BG_Consum_NoCrab, data=two, var.equal=T)
#t = 3.3884, df = 18, p-value = 0.003275, significant

#Plot the feeding assay data
#Load Figure data
Otter_marsh_data_assay_fig<-read.table(file="Otters_Marsh_2015_R_assay_fig_data.csv", header=TRUE, sep=",")
str(Otter_marsh_data_assay_fig)

#plot the data
assay_plot <- ggplot(Otter_marsh_data_assay_fig, aes(x = Marsh_type, y = biomass_consumed, fill= Treatment, group = paste(Treatment, Marsh_type))) +
  geom_violin(position=position_dodge(1))+
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1), color="black", size=2.5)+
  theme_bw(base_size=15)+theme(legend.position= c(0.15, 0.85))+
  annotate("text", x=0.75, y=0.35, label= "a", size=7) +
  annotate("text", x = 1.25, y=0.45, label = "a", size=7)+
  annotate("text", x=1.75, y=0.37, label= "x", size=7) +
  annotate("text", x = 2.15, y=0.5, label = "y", size=7)+
  xlab("") + ylab(expression("Marsh consumed (g per trial)")) +
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.text.x=element_text())
assay_plot
#Fig. S2

####Field Experimental Data analyses 2015####

#glmm for burrows with site as a random factor
library(glmmTMB)
burrow_glmm_nb2<- glmmTMB(Burrows_OE_msq ~ treat + (1|plot_id) + (1|Site),data = Otter_marsh_data, family=nbinom2())
summary(burrow_glmm_nb2)

#Both variance components collapsed to zero
# diagnose(burrow_glmm)

burrow_glmm<- glmmTMB(Burrows_OE_msq ~ treat + (1|Site),data = Otter_marsh_data, family=nbinom2())
summary(burrow_glmm)

burrow_glm<- glmmTMB(Burrows_OE_msq ~ treat,data = Otter_marsh_data, family=nbinom2())
summary(burrow_glm)

burrow_glm_nb1<- glmmTMB(Burrows_OE_msq ~ treat,data = Otter_marsh_data, family=nbinom1())
summary(burrow_glm_nb1)

burrow_glm_p <- glmmTMB(Burrows_OE_msq ~ treat,data = Otter_marsh_data, family=poisson())
summary(burrow_glm_p)

burrow_glm_qp <- glm(Burrows_OE_msq ~ treat,data = Otter_marsh_data, family=quasipoisson())
summary(burrow_glm_qp)

AIC(burrow_glm, burrow_glm_nb1, burrow_glm_p)
summary(burrow_glm)
#non-Significant Otter effect, P = 0.134
#decrease in crab burrows with sea otters present

ggplot(Otter_marsh_data, aes(x = treat, y = Burrows_OE_msq)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

tapply(Otter_marsh_data$Burrows_OE_msq,Otter_marsh_data$treat, FUN=mean )
burrow_dec<-22.24/13.04*100
burrow_dec #171% increase in burrows without otters


#glmm for crab density with site as a random factor
crab_glmm<- glmmTMB(crab_density_OE ~ treat + (1|plot_id)+(1|Site),data = Otter_marsh_data,family=nbinom2(), na.action=na.omit)
summary(crab_glmm) #significant otter effect, P=0.0307

#Random variance nearly 0
crab_glm<- glm(crab_density_OE ~ treat,data = Otter_marsh_data,family=nbinom2(), na.action=na.omit)
summary(crab_glm)
#P = 0.03669

#decrease in crab density with sea otters present
tapply(Otter_marsh_data$crab_density_OE,Otter_marsh_data$treat, FUN=mean )
crab_dec<-2.12/1.04*100
crab_dec #204% increase in burrows without otters

ggplot(Otter_marsh_data, aes(x = treat, y = crab_density_OE)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

#glmm for bulk density
bulkdens_glmm<- glmmTMB(BulkDens_OE_kg_msq ~ treat + (1|plot_id)+(1|Site),data = Otter_marsh_data, family=Gamma(link="log"), na.action=na.omit)
summary(bulkdens_glmm)

#But Variance collapses to 0, so:
bulkdens_glm <- glmmTMB(BulkDens_OE_kg_msq ~ treat,data = Otter_marsh_data, family=Gamma(link="log"), na.action=na.omit)
summary(bulkdens_glm)
#P 0.0911


tapply(Otter_marsh_data$BulkDens_OE_kg_msq,Otter_marsh_data$treat, FUN=mean )
bulk_inc<-((136.96-126.18)/136.96)*100 #Decrease by 8%
bulk_inc

ggplot(Otter_marsh_data, aes(x = treat, y = BulkDens_OE_kg_msq)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

#glmm for belowground biomass

ggplot(Otter_marsh_data, aes(x = treat, y = BGMass_OE_kg_msq)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

BG_mass_glmm<- glmmTMB(BGMass_OE_kg_msq ~ treat + (1|plot_id)+(1|Site),data = Otter_marsh_data, family=Gamma(link="log"), na.action=na.omit)
summary(BG_mass_glmm)

#Significant otter effect, P = 0.0353
#Otters had 0.46 kg more BG biomass per msq
#increase in bulk density with sea otters present
tapply(Otter_marsh_data$BGMass_OE_kg_msq,Otter_marsh_data$treat, FUN=mean )
BG_inc<-((3.709-3.25)/3.709)*100
BG_inc #12% Decrease in BG biomass

#glmm for aboveground biomass

ggplot(Otter_marsh_data, aes(x = treat, y = AGMass_OE_kg_msq)) +
  geom_point(position = position_jitter(width = 0.1, height = 0))

AG_mass_glmm<- glmmTMB(AGMass_OE_kg_msq ~ treat + (1|plot_id)+(1|Site),data = Otter_marsh_data,family=Gamma(link="log"), na.action=na.omit)
summary(AG_mass_glmm)

#Variance collapses to 0
AG_mass_glm<- glmmTMB(AGMass_OE_kg_msq ~ treat,data = Otter_marsh_data,family=Gamma(link="log"), na.action=na.omit)
summary(AG_mass_glm)

#Significant otter effect, P = 0.001179
#Otters had 0.37 kg more BG biomass per msq
#increase in bulk density with sea otters present
tapply(Otter_marsh_data$AGMass_OE_kg_msq,Otter_marsh_data$treat, FUN=mean, na.action=na.omit)
AG_inc<-(1.07-0.7)/1.07*100
AG_inc #Decreased by 35%

#Dataframe for figures
Otter_marsh_data_figs2<-read.table(file="Otters_Marsh_2015_R_figdata2.csv", header=TRUE, sep=",")
str(Otter_marsh_data_figs2)

#Burrows
min_y<-min(0)
max_y<-max(100)
p1<-ggplot(data=Otter_marsh_data, aes(x=treat, y=Burrows_OE_msq, fill=treat)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize = .7, fill="black")+
  theme_bw() +
  xlab("") + ylab(expression(Burrows%.%m^-2))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.75, y=80, label= "NS", size=8)+
  theme_bw(base_size = 14, base_family = "")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(),legend.position="none") +
  ggtitle('A')+theme(plot.title=element_text(hjust=0))+
  theme(axis.text.x=element_text())+
  scale_x_discrete(limit = c("no_otter", "otter"), labels = c("No Otter","Otter"))
p1

#crab density
min_y<-min(0)
max_y<-max(10)
p2<-ggplot(data=Otter_marsh_data, aes(x=treat, y=crab_density_OE, fill=treat)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize = .7, fill="black")+
  theme_bw() +
  xlab("") + ylab(expression(Crabs%.%m^-2))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.75, y=9.5, label= "*", size=14)+
  theme_bw(base_size = 14, base_family = "")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), legend.position="none") +
  ggtitle('A')+theme(plot.title=element_text(hjust=0))+theme(axis.text.x=element_blank())
p2

#Aboveground biomass
min_y<-min(0)
max_y<-max(3)
p3<-ggplot(data=Otter_marsh_data, aes(x=treat, y=AGMass_OE_kg_msq, fill=treat)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize = .7, fill="black")+
  theme_bw() +
  xlab("") + ylab(expression("Aboveground mass " (kg%.%m^-2)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.75, y=2.8, label= "**", size=14)+
  theme_bw(base_size = 14, base_family = "")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(), legend.position="none") +
  ggtitle('B')+theme(plot.title=element_text(hjust=0))+theme(axis.text.x=element_blank())
p3

#belowground biomass
min_y<-min(0)
max_y<-max(8)
p4<-ggplot(data=Otter_marsh_data, aes(x=treat, y=BGMass_OE_kg_msq, fill=treat)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize = .7, fill="black")+
  theme_bw() +
  xlab("") + ylab(expression("Belowground mass " (kg%.%m^-2)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.75, y=7.5, label= "*", size=14)+
  theme_bw(base_size = 14, base_family = "")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(),legend.position="none") +
  ggtitle('C')+theme(plot.title=element_text(hjust=0))+
  theme(axis.text.x=element_text())+
  scale_x_discrete(limit = c("no_otter", "otter"), labels = c("No Otter","Otter"))
p4

#Bulk density
min_y<-min(0)
max_y<-max(180)
p5<-ggplot(data=Otter_marsh_data, aes(x=treat, y=BulkDens_OE_kg_msq, fill=treat)) +
  geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize = .7, fill="black")+
  theme_bw() +
  xlab("") + ylab(expression("Bulk density " (kg%.%m^-2)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_y_continuous(limits=c(min_y,max_y)) +
  annotate("text",x=0.75, y=180, label= "P = 0.091", size=6)+
  theme_bw(base_size = 14, base_family = "")+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(),legend.position="none") +
  ggtitle('B')+theme(plot.title=element_text(hjust=0))+
  theme(axis.text.x=element_text())+
  scale_x_discrete(limit = c("no_otter", "otter"), labels = c("No Otter","Otter"))
p5

#Make Experimental figure
cowplot::plot_grid(p2, p3, p4, ncol = 1, align = "v")
ggsave("experiment_figure_aligned.pdf", width = 6, height = 13)

#End