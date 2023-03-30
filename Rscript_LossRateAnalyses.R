######## libraries ########  
library(MASS)
library(lme4)
library(car) 
library(ggplot2)
library(broom)
library(gridExtra)
library(dplyr)
library(segmented)
library(lmerTest)
library(MuMIn)
library("scales")  
library('geoR')


quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

bc.transform <- function(x,L){(x^L-1)/L}
bc.backtransform <- function(y,L){(L*y+1)^(1/L)}

MLLoss <- read.table('/Users/jennyjames/Desktop/LossRate_10iterations/LossRates/PfamMLLikeBayesPipeline_MLLoss_AnimalOnly.txt', sep = '\t', header = TRUE)
AnimalCounts <- read.table('/Users/jennyjames/Desktop/LossRate_10iterations/CountPfam/Count_pfam_animal.tsv', sep = '\t', header = TRUE)
PfamUIDs <- read.csv('/Users/jennyjames/Desktop/LossRates/LossRatesForManuscript/PfamUIDsTable_EnsemblAndNCBI.csv', header = TRUE, na.strings = 'None')

MLLoss_Stats <- inner_join(MLLoss, PfamUIDs, by = c('Pfam' = 'PfamUID'))
MLLoss_Stats <- merge(x = MLLoss_Stats, y = AnimalCounts, by = "Pfam", all.x = TRUE)

### this is only needed when looking at number of pfam instances- slightly restricts the dataset (6 
### pfam alignments unable to be downloaded)
MLLoss_Stats <- inner_join(x = MLLoss_Stats, y = AnimalCounts, by = "Pfam")

attach(MLLoss_Stats)

summary(MLLossRate)

###Do all of our transforms, make supplementary figure 1
par(mfrow = c(4,2))

Loss.box = boxcoxfit(MLLoss_Stats$MLLossRate, lambda2=TRUE)
Loss.box
lambda = Loss.box$lambda[1]
lambda2 = Loss.box$lambda[2]
if(lambda==0){Loss.norm = log(MLLoss_Stats$MLLossRate + lambda2)}
if(lambda!=0){Loss.norm = ((MLLoss_Stats$MLLossRate + lambda2) ^ lambda - 1) / lambda}

hist(MLLoss_Stats$MLLossRate, col="gray", xlab = "Loss rate/MY", main = '')
text(0.6, 6000, 'A)', cex = 1.5)

hist(Loss.norm, col="gray", xlab = "Box-Cox transformed loss rate/MY", main = '')
text(-0.5, 1000, 'B)', cex = 1.5)

model <- lm(MLLoss_Stats$MLLossRate ~ MLLoss_Stats$MeanISD_AnimalSpecific)
res <- resid(model)
qqnorm(res, main = '')
qqline(res)
text(-1.8, 0.58, 'C) Loss rate ~ ISD', cex = 1.5)
model.box <- lm(Loss.norm ~ MLLoss_Stats$MeanISD_AnimalSpecific)
res <- resid(model.box)
qqnorm(res, main = '')
qqline(res)
text(-1.2, 4.6, 'D) Transformed loss rate 
     ~ ISD', cex = 1.5)

Count.box = boxcoxfit(MLLoss_Stats$Animal_mean, lambda2=TRUE)
Count.box
countlambda = Count.box$lambda[1]
countlambda2 = Count.box$lambda[2]
if(countlambda==0){Count.norm = log(MLLoss_Stats$Animal_mean + countlambda2)}
if(countlambda!=0){Count.norm = ((MLLoss_Stats$Animal_mean + countlambda2) ^ countlambda - 1) / countlambda}
hist(MLLoss_Stats$Animal_mean, col = "grey", xlab = "Mean instances", main = '')
text(560, 6000, 'E)', cex = 1.5)

hist(Count.norm, col="gray", xlab = "Box-Cox transformed mean instances", main = '')
text(1.5, 690, 'F)', cex = 1.5)

model <- lm(MLLoss_Stats$Animal_mean ~ MLLoss_Stats$MeanISD_AnimalSpecific)
res <- resid(model)
qqnorm(res, main = '')
qqline(res)
text(-0.32, 505, 'G) Transformed mean instances 
     ~ ISD', cex = 1.5)
model.box <- lm(Count.norm ~ MLLoss_Stats$MeanISD_AnimalSpecific)
res <- resid(model.box)
qqnorm(res, main = '')
qqline(res)
text(-1.2, 0.8, 'H) Box-Cox transformed mean 
     instances ~ ISD', cex = 1.5)

par(mfrow = c(1,1))

###Supplementary figure 2
boxplot_p <- ggplot()+  
  geom_violin(data = MLLoss_Stats, aes(group=cut_interval(as.numeric(Count.norm), n = 20), x=Count.norm, y = Loss.norm), 
               fill= "azure2", colour = "darkslategray4", varwidth = TRUE, position = position_nudge(-0.04), outlier.size = NULL) +
  geom_jitter(data = MLLoss_Stats, aes(x=Count.norm, y = Loss.norm), color="black", size=0.4, alpha=0.3) +
  geom_smooth(method = "lm", aes(x=Count.norm, y = Loss.norm), color = 'midnightblue') +
  scale_y_continuous(breaks = c(0, -3, -6, -9), labels = signif(bc.backtransform(c(0, -3, -6, -9), lambda), digits = 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5), labels = signif(bc.backtransform(c(0, 0.5, 1.0, 1.5), countlambda), digits = 1))



plot <- boxplot_p + labs(x = "Mean pfam instances") +
  labs(y = "Loss rate") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))
plot  

summary(lm(Count.norm ~ Loss.norm))



### ISD analysis

ISDMod_BC <- lm(Loss.norm ~ MeanISD_AnimalSpecific)
summary(ISDMod_BC)
ISD_notransform <- lm(MLLossRate ~ MeanISD_AnimalSpecific)

### Include a break in ISD
ISDBreak <- segmented(ISDMod_BC)
print.segmented(ISDBreak)
summary(ISDBreak)
confint.segmented(ISDBreak)
slope(ISDBreak)
intercept(ISDBreak)
plot.segmented(ISDBreak)

anova(ISDMod_BC, ISDBreak)

MLLoss_Stats$high_ISD <- ifelse(MLLoss_Stats$MeanISD_AnimalSpecific-0.179>0, MLLoss_Stats$MeanISD_AnimalSpecific-0.179, 0)
MLLoss_Stats$low_ISD<-ifelse(0.179-MLLoss_Stats$MeanISD_AnimalSpecific>0, 0.179-MLLoss_Stats$MeanISD_AnimalSpecific, 0)
attach(MLLoss_Stats)

low <- lm(Loss.norm ~ low_ISD)
high <- lm(Loss.norm ~ high_ISD)
anova(ISDBreak, high)
anova(ISDBreak, high)$Pr
anova(ISDBreak, low)
anova(ISDBreak, low)$Pr

ISDBreakMod <- lm(Loss.norm ~ high_ISD + low_ISD)
lmp(ISDBreakMod)


### no we want to look at how variance changes with pfam age- but we need to re-date some pfams 
### based on their updated species composition dates, although we ignore any pfams that predate the
### animal plant split- 1404 MY. Age_Oldest_MY lists the James et al. dates- and Date are the newly
### calculated dates using the ML pipeline.

MLLoss_Stats$Updated_Age <- with(MLLoss_Stats, ifelse(MLLoss_Stats$Age_Oldest_MY >= 1404, MLLoss_Stats$Age_Oldest_MY,
                                 ifelse(MLLoss_Stats$Date < MLLoss_Stats$Age_Oldest_MY, MLLoss_Stats$Date, MLLoss_Stats$Date)))

summary(lm(MLLoss_Stats$ISDVariance_AnimalSpecific ~ MLLoss_Stats$Updated_Age))
lmp(lm(MLLoss_Stats$ISDVariance_AnimalSpecific ~ MLLoss_Stats$Updated_Age))
Ancient <- MLLoss_Stats[which(MLLoss_Stats$Updated_Age > 2000),]
summary(Ancient$MeanISD_AnimalSpecific)

### ISD standard deviation per age category, and loss rate
ISDVarAge <- aggregate(MLLoss_Stats$MeanISD_AnimalSpecific, by=list(MLLoss_Stats$Updated_Age), FUN = sd, na.rm = TRUE)
ISDCountAge <- MLLoss_Stats %>% count(Updated_Age)

ISDVarAgeDf <- data.frame(cbind(ISDVarAge, ISDCountAge))
colnames(ISDVarAgeDf) <- c("Age_Oldest_MY", "Var", "Age_Oldest_My", "n")

summary(lm(ISDVarAgeDf$Var ~ ISDVarAgeDf$Age_Oldest_My, weights = ISDVarAgeDf$n))
lmp(lm(ISDVarAgeDf$Var ~ ISDVarAgeDf$Age_Oldest_My, weights = ISDVarAgeDf$n))


### Figures 3A and 3B
boxplot_p <- ggplot()+  
  geom_violin(varwidth = TRUE, outlier.shape = NA, data = MLLoss_Stats, aes(group=cut_interval(as.numeric(MeanISD_AnimalSpecific), n = 18), x=as.numeric(MeanISD_AnimalSpecific), y = Loss.norm), 
               fill= "azure2", colour = "darkslategray4") + 
  scale_y_continuous(breaks = c(-9, -6, -3, 0), labels = c(8e-05, 0.002, 0.05, 1)) +
  geom_jitter(data = MLLoss_Stats, aes(x=as.numeric(MeanISD_AnimalSpecific), y = Loss.norm), size=0.4, alpha = 0.3) +
  geom_text(label = "A)", aes(x = 0.9, y = -0.7), size = 10) +
  geom_smooth(method = "loess", aes(x=as.numeric(MeanISD_AnimalSpecific), y = Loss.norm), color = 'midnightblue')

Loss_ISD <- boxplot_p + labs(x = "ISD") +
  labs(y = "Losses/Million years") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30)) + theme(axis.title.x = element_blank())


plot_p <- ggplot(data = ISDVarAgeDf)+  
  geom_point(aes(x = Age_Oldest_MY, y = Var, size = n), colour = "darkslategray4")+ scale_size(range = c(2,10))+ 
  #  geom_smooth(method = "lm", aes(x = Age_Oldest_MY, y = Sd), colour = "midnightblue")+
  geom_abline(intercept = 2.344e-01, slope = -4.100e-05, colour = "midnightblue", size = 1)+
  #+- 2 standard errors on slope
  geom_abline(intercept = 2.344e-01 + 5.095e-03+ 5.095e-03, slope = -4.100e-05+2.237e-06+2.237e-06, colour = "grey", size = 1, alpha = 0.7) +
  geom_abline(intercept = 2.344e-01- 5.095e-03- 5.095e-03, slope = -4.100e-05-2.237e-06-2.237e-06, colour = "grey", size = 1, alpha = 0.7)

Loss_ISDstdev <- plot_p + labs(x = "Pfam Age, MY") +
  labs(y = "ISD,\nstandard deviation") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30)) + theme(axis.text.x = element_text(angle = 90))
 


ISD_ageplot <- ggplot(data = MLLoss_Stats)+  
  geom_boxplot(aes(group=Updated_Age, x=Updated_Age, y = MeanISD_AnimalSpecific), 
              fill= "azure2", colour = "darkslategray4", varwidth = TRUE, width = 200, outlier.shape = NA) +
  geom_jitter(aes(x = Updated_Age , y = MeanISD_AnimalSpecific), size=0.4, alpha = 0.3) +
  geom_smooth(method = "loess", aes(x = Updated_Age , y = MeanISD_AnimalSpecific), color = 'midnightblue')

ISD_ageplot<- ISD_ageplot + labs(x = "Pfam Age, MY") +
  labs(y = "ISD") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30)) + theme(axis.text.x = element_text(angle = 90))


### Clustering analysis

ClusteringMod_BC <- lm(Loss.norm ~ MeanClustering_AnimalSpecific)
lmp(ClusteringMod_BC)
ClusteringBreak <- segmented(ClusteringMod_BC)
print.segmented(ClusteringBreak)
summary(ClusteringBreak)
confint.segmented(ClusteringBreak)
slope(ClusteringBreak)
intercept(ClusteringBreak)
plot.segmented(ClusteringBreak)

anova(ClusteringBreak, ClusteringMod_BC)
anova(ClusteringBreak, ClusteringMod_BC)$Pr

MLLoss_Stats$high_Clustering <- ifelse(MLLoss_Stats$MeanClustering_AnimalSpecific-0.81>0,  MLLoss_Stats$MeanClustering_AnimalSpecific, 0)
MLLoss_Stats$low_Clustering<-ifelse(0.81-MLLoss_Stats$MeanClustering_AnimalSpecific>0, MLLoss_Stats$MeanClustering_AnimalSpecific, 0)
attach(MLLoss_Stats)

low <- lm(Loss.norm ~ low_Clustering)
high <- lm(Loss.norm ~ high_Clustering)
anova(ClusteringBreak, high)
anova(ClusteringBreak, high)$Pr
anova(ClusteringBreak, low)
anova(ClusteringBreak, low)$Pr
summary(Ancient$MeanClustering_AnimalSpecific)

summary(lm(PfamUIDs$MeanClustering_AnimalSpecific ~ PfamUIDs$Age_Oldest_MY))


###Clustering variance per age category, and loss rate
summary(lm(MLLoss_Stats$ClusteringVariance_AnimalSpecific ~ MLLoss_Stats$Updated_Age))
lmp(lm(MLLoss_Stats$ClusteringVariance_AnimalSpecific ~ MLLoss_Stats$Updated_Age))

ClusteringVarAge <- aggregate(MLLoss_Stats$MeanClustering_AnimalSpecific, by=list(MLLoss_Stats$Updated_Age), FUN = sd, na.rm = TRUE)
ClusteringCountAge <- MLLoss_Stats %>% count(Updated_Age)

ClusteringVarAgeDf <- data.frame(cbind(ClusteringVarAge, ClusteringCountAge))
colnames(ClusteringVarAgeDf) <- c("Age_Oldest_MY", "Var", "Age_Oldest_spare", "n")

summary(lm(ClusteringVarAgeDf$Var ~ ClusteringVarAgeDf$Age_Oldest_MY, weights = ClusteringVarAgeDf$n))
lmp(lm(ClusteringVarAgeDf$Var ~ ClusteringVarAgeDf$Age_Oldest_MY, weights = ClusteringVarAgeDf$n))


### Figures 
boxplot_p <- ggplot(data = MLLoss_Stats)+  
  geom_violin(varwidth = TRUE, outlier.shape = NA, aes(group=cut_interval(as.numeric(MeanClustering_AnimalSpecific), n = 25), x=as.numeric(MeanClustering_AnimalSpecific), y = Loss.norm), 
               fill= "azure2", colour = "darkslategray4") + 
  xlim(c(0,2.5)) +
  scale_y_continuous(breaks = c(-9, -6, -3, 0), labels = c(8e-05, 0.002, 0.05, 1)) +
  geom_text(label = "B)", aes(x = 2.25, y = -0.7), size = 10) +
  geom_jitter(aes(x=as.numeric(MeanClustering_AnimalSpecific), y = Loss.norm), size=0.4, alpha = 0.3) +
  geom_smooth(method = "loess", aes(x=as.numeric(MeanClustering_AnimalSpecific), y = Loss.norm), color = 'midnightblue', span = 0.99)

Loss_clustering <- boxplot_p + labs(x = "Clustering") +
  labs(y = "Losses/Million years") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))+ theme(axis.title.x = element_blank())

plot_p <- ggplot(data = ClusteringVarAgeDf) +  
  geom_point(aes(x = Age_Oldest_MY, y = Var, size = n), colour = "darkslategray4") + scale_size(range = c(2,10)) +
  #  geom_smooth(method = "lm", aes(x = Age_Oldest_MY, y = Sd), colour = "midnightblue") +
  geom_abline(intercept = 3.406e-01, slope = -5.394e-05, colour = "midnightblue", size = 1)+
  #+- 2 standard errors on slope: 6.194e-06
  geom_abline(intercept = 3.406e-01+1.506e-02+1.506e-02, slope = -5.394e-05+6.612e-06+6.612e-06, colour = "grey", size = 1, alpha = 0.7) +
  geom_abline(intercept = 3.406e-01-1.506e-02-1.506e-02, slope = -5.394e-05-6.612e-06-6.612e-06, colour = "grey", size = 1, alpha = 0.7)

Loss_clusteringstdev <- plot_p + labs(x = "Pfam Age, MY") +
  labs(y = "Clustering,\n standard deviation") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30)) + theme(axis.text.x = element_text(angle = 90))

Clustering_ageplot <- ggplot(data = MLLoss_Stats)+  
  geom_boxplot(aes(group=Updated_Age, x=Updated_Age, y = MeanClustering_AnimalSpecific), 
               fill= "azure2", colour = "darkslategray4", varwidth = TRUE, width = 200, outlier.shape = NA) +
  geom_jitter(aes(x = Updated_Age , y = MeanClustering_AnimalSpecific), size=0.4, alpha = 0.3) +
  geom_smooth(method = "loess", aes(x = Updated_Age , y = MeanClustering_AnimalSpecific), color = 'midnightblue')

Clustering_ageplot<- Clustering_ageplot + labs(x = "Pfam Age, MY") +
  labs(y = "Clustering") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30)) + theme(axis.text.x = element_text(angle = 90))


### Age analysis

AgeLoss <- lm(Loss.norm ~ Updated_Age)
summary(AgeLoss)
lmp(AgeLoss)

### supplementary figure 3:  Age and loss rate 
boxplot_p <- ggplot()+  
  geom_boxplot(data = MLLoss_Stats, aes(group=Updated_Age, x=Updated_Age, y = Loss.norm), 
               fill= "azure2", colour = "darkslategray4", position = "dodge2", varwidth = TRUE, width = 200, outlier.shape = NA) +
  scale_y_continuous(breaks = c(-9, -6, -3, 0), labels = c(8e-05, 0.002, 0.05, 1)) +
  geom_jitter(data = MLLoss_Stats, aes(x = Updated_Age , y = Loss.norm), size=0.4, alpha = 0.3) +
  geom_smooth(method = "lm", aes(x = Updated_Age , y = Loss.norm), color = 'midnightblue')

boxplot_p<- ggplot(data = MLLoss_Stats)+  
  geom_violin(aes(group=Updated_Age, x=Updated_Age, y = Loss.norm), 
              fill= "azure2", colour = "darkslategray4", width = 200, outlier.shape = NA) +
  scale_y_continuous(breaks = c(-9, -6, -3, 0), labels = c(8e-05, 0.002, 0.05, 1)) +
  geom_jitter(aes(x = Updated_Age , y = Loss.norm), size=0.4, alpha = 0.3) +
  geom_smooth(method = "lm", aes(x = Updated_Age , y = Loss.norm), color = 'midnightblue')

boxplot_plot <- boxplot_p + labs(x = "Pfam Age, MY") +
  labs(y = "Loss Rate") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))
boxplot_plot 


### defining all our models here for clarity
ISDAgeModBreak <- lm(Loss.norm ~ high_ISD + low_ISD + Updated_Age)
ISDModBreak <- lm(Loss.norm ~ high_ISD + low_ISD)
ISDAgeModBreakRan <- lmer(Loss.norm ~ high_ISD + low_ISD + (1|Updated_Age))
ISDAgeModRan <- lmer(Loss.norm ~ MeanISD_AnimalSpecific + (1|Updated_Age))

ClusteringAgeModBreak <- lm(Loss.norm ~ high_Clustering + low_Clustering + Updated_Age)
ClusteringModBreak <- lm(Loss.norm ~ high_Clustering + low_Clustering)
ClusteringAgeModBreakRan <- lmer(Loss.norm ~ high_Clustering + low_Clustering + (1|Updated_Age))
ClusteringAgeModRan <- lmer(Loss.norm ~ MeanClustering_AnimalSpecific + (1|Updated_Age))

### comparing random and fixed effects models
r.squaredGLMM(ISDAgeModBreakRan)
r.squaredGLMM(ISDAgeModBreak)
r.squaredGLMM(ClusteringAgeModBreakRan)
r.squaredGLMM(ClusteringAgeModBreak)
AIC(ISDAgeModBreak, ISDAgeModBreakRan)
AIC(ClusteringAgeModBreak, ClusteringAgeModBreakRan)


### checking statistical support for breakpoints in random effects models

r.squaredGLMM(ISDAgeModBreakRan)
anova(ISDAgeModBreakRan, ISDAgeModRan)

anova(ISDAgeModBreak, ISDModBreak)
anova(ISDAgeModBreak, ISDModBreak)$Pr

anova(ClusteringAgeModBreakRan, ClusteringAgeModRan)

anova(ClusteringAgeModBreak, ClusteringModBreak)
anova(ClusteringAgeModBreak, ClusteringModBreak)$Pr


####plot layout ideas?
gl <- c(plot_3A_x, plot_3B_x, plot_3C_x, plot_3D_x)


grid.arrange(plot_3A_x, plot_3B_x, plot_3C_x, plot_3D_x, ncol = 2, 
             layout_matrix = rbind(c(1,1,2),
                                   c(3,3,4)), left = 'Loss rate / MY')



### Pfam instances analysis
### Mean number of Pfam instances results, calculated as a mean that has been transformed, Count.norm

summary(lm(Count.norm ~ Age_Oldest_MY))
summary(lm(Count.norm ~ MeanISD_AnimalSpecific))
summary(lm(Count.norm ~ MeanClustering_AnimalSpecific))

CopyISD <- lm(Count.norm ~ MeanISD_AnimalSpecific)
CopyISDBreak <- segmented(CopyISD)
print.segmented(CopyISDBreak)
summary(CopyISDBreak)
confint.segmented(CopyISDBreak)
slope(CopyISDBreak)

anova(CopyISD, CopyISDBreak)

CopyClustering <- lm(Count.norm ~ MeanClustering_AnimalSpecific)
CopyClusteringBreak <- segmented(CopyClustering)
print.segmented(CopyClusteringBreak)
summary(CopyClusteringBreak)
confint.segmented(CopyClusteringBreak)
slope(CopyClusteringBreak)

anova(CopyClustering, CopyClusteringBreak)
anova(CopyClustering, CopyClusteringBreak)$Pr

### and plotting results





Instances_ISD <- ggplot()+
  geom_violin(varwidth = TRUE, outlier.shape = NA, data = MLLoss_Stats, aes(group=cut_interval(as.numeric(MeanISD_AnimalSpecific), n = 18), x=as.numeric(MeanISD_AnimalSpecific), y = Count.norm),
               fill= "azure2", colour = "darkslategray4") +
  geom_smooth(method='loess', aes(MeanISD_AnimalSpecific, Count.norm), color = 'midnightblue') +
#  scale_y_continuous(breaks = c(-9, -6, -3, 0), labels = c(8e-05, 0.002, 0.05, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), labels = c(1, 2, 5, 60)) + 
  geom_text(label = "C)", aes(x = 0.9, y = 1.5), size = 10) +
  geom_jitter(data = MLLoss_Stats, aes(x=as.numeric(MeanISD_AnimalSpecific), y = Count.norm), size=0.4, alpha = 0.3) +
  labs(x = "ISD") +
  labs(y = "Mean instances") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))


Instances_clustering <- ggplot()+
  geom_violin(varwidth = TRUE, outlier.shape = NA, aes(group=cut_interval(as.numeric(MeanClustering_AnimalSpecific), n = 25), x=as.numeric(MeanClustering_AnimalSpecific), y = Count.norm), 
               fill= "azure2", colour = "darkslategray4") + 
  geom_jitter(data = MLLoss_Stats, aes(MeanClustering_AnimalSpecific, Count.norm), size=0.4, alpha = 0.3) +
  geom_smooth(method='loess', aes(MeanClustering_AnimalSpecific, Count.norm), color = 'midnightblue', span = 0.99) +
  #  scale_y_continuous(breaks = c(-9, -6, -3, 0), labels = c(8e-05, 0.002, 0.05, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), labels = c(1, 2, 5, 60)) + 
  geom_text(label = "D)", aes(x = 2.25, y = 1.5), size = 10) +
  xlim(c(0,2.5)) +
  labs(x = "Clustering") +
  labs(y = "Mean instances") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))



grid.arrange(Loss_ISD, Instances_ISD, Loss_clustering, Instances_clustering, ncol = 2, 
             layout_matrix = rbind(c(1,3),
                                   c(2,4)))


plots <- list(Loss_ISD, Loss_clustering, Instances_ISD, Instances_clustering)
grobs <- list()
widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

do.call("grid.arrange", c(grobs, ncol = 2))

