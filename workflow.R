#Basic analysis of amplicon data

setwd("/Users/Admin/Desktop/CNR/Chemostat Marcels/")

OTUList<-read.csv("zotutab.csv")
row.names(OTUList)<-OTUList[,1]
OTUList<-OTUList[,-1]


taxa<-read.csv("taxonomy.csv")
row.names(taxa)<-taxa[,1]
taxaall<-taxa[,-1]


taxaNC<-subset(taxaall, taxaall$c!="c:Chloroplast")
taxaNM<-subset(taxaNC, taxaNC$f!="f:Mitochondria")
taxa<-subset(taxaNM, taxaNM$p!="")
ttaxa<-as.data.frame(t(taxa))
OTU<-OTUList[row.names(OTUList)%in%row.names(taxa),]
OTUList<-as.data.frame(OTU)
tOTUList<-as.data.frame(t(OTU))
SumSample<-colSums(OTUList)
taxa<-taxa[row.names(taxa)%in%row.names(OTUList),]

library("GUniFrac")
rareOTU<-Rarefy(tOTUList, depth = min(SumSample))  #Rarefaction
rffOTU<-as.data.frame(rareOTU$otu.tab.rff)
totalSeq<-rowSums(rffOTU)
trffOTU<-t(rffOTU)
raOTU <- rffOTU  #[ ,colSums(rffOTU)!=0] 
traOTU<-t(raOTU)
paOTU<-raOTU
paOTU[paOTU>0]=1

variables<-read.csv("variables.csv")
row.names(variables)<-variables[,1]
variables<-variables[,-1]

#Alpha diversity

alpha<-rowSums(paOTU) #Number of OTUs
x<-variables$SOURCE

library("ggplot2")
alphaV<-as.data.frame(cbind(rowSums(paOTU), variables))
colnames(alphaV)<-c("sum", "SOURCE", "REPLICATE", "TIME", "AB", "PAAB", "TT", "post")
pR<-ggplot(alphaV, aes(x=factor(TIME), y=sum, fill = factor(SOURCE)))  + geom_dotplot(binaxis = "y", stackdir = "center")
pRR<-pR+labs(x = "time", y="richness", fill="source")+theme(text = element_text(size=20), legend.position = "bottom")+ scale_fill_manual(name = NULL, labels = c("Orta", "Varese", "Effluent"),values=c("snow1", "grey43", "grey1" ))
pRR
library(lmerTest)
library(car)
Source<-variables$SOURCE
time<-variables$TIME
glmR<-glm(alphaV$sum~Source+factor(time), family=poisson)
summary(glmR)
Anova(glmR)
library(emmeans)
marginal = emmeans(glmR, ~ Source)
marginal = emmeans(glmR, ~ Source+factor(time))
marginal = emmeans(glmR, ~ factor(time))

pairs(marginal)


model1<-aov(log(rowSums(paOTU)) ~ x)
posthoc <- TukeyHSD(x=model1) #which ones are significantly different from each other
posthoc
out<-capture.output(posthoc)
write.csv(out,"alpha_HSD.csv")

y<-variables$TIME
plot(as.factor(y),alpha, main="number of OTUs")
model1<-aov(rowSums(paOTU) ~ as.factor(y))
posthoc <- TukeyHSD(x=model1) #which ones are significantly different from each other
posthoc
out<-capture.output(posthoc)
write.csv(out,"alpha_x.csv")

alpha<-rowSums(paOTU) #Number of OTUs
x<-as.factor(variables$TT)
plot(x,alpha, main="number of OTUs")
model1<-aov(log(rowSums(paOTU)) ~ x)
posthoc <- TukeyHSD(x=model1) #which ones are significantly different from each other
posthoc
out<-capture.output(posthoc)
write.csv(out,"alpha_HSD.csv")





library(vegan)


#Beta diversity 
library("betapart")

library("RColorBrewer")
library(vegan)
#ADONIS in vegan: y~a where y is a dist matrix (dissimilarity) 
#Strata: If the experimental design has nestedness, then use strata to test hypotheses. For instance, imagine we are testing the whether a plant community is influenced by nitrate amendments, and we have two replicate plots at each of two levels of nitrate (0, 10 ppm). We have replicated the experiment in three fields with (perhaps) different average productivity. In this design, we would need to specify strata = field so that randomizations occur only within each field and not across all fields
betabray<-vegdist(raOTU,method="bray")
plot(hclust(betabray, method="complete"),hang=-1, main='Bray-Curtis Bacteria', sub='', xlab='', cex=1) #plot cluster analysis of betapair


adonis1<-adonis(betabray~variables$SOURCE*variables$TIME, permutations=9999, strata=variables$VESSEL)
out<-capture.output(adonis1)
write(out, "10_adonisall.txt", append=TRUE)
adonis1

taxa<-read.csv("taxonomy.csv")
row.names(taxa)<-taxa[,1]
taxa<-taxa[,-1]
taxan0<-taxa[row.names(taxa)%in%row.names(traOTU),] 
taxa<-taxan0
library("reshape2")
otu<-traOTU

vari<-read.csv("variables.csv")



#plot combined data

library(ggplot2)

sumTaxa1sample<-do.call("rbind", as.list(by(otu[,], taxa$p, colSums)))
sumTaxa2sample<-do.call("rbind", as.list(by(otu[,], taxa$c, colSums)))
sumTaxa3sample<-do.call("rbind", as.list(by(otu[,], taxa$o, colSums)))
sumTaxa4sample<-do.call("rbind", as.list(by(otu[,], taxa$f, colSums)))
sumTaxa5sample<-do.call("rbind", as.list(by(otu[,], taxa$g, colSums)))

alldataT1<-rowSums(sumTaxa1sample)
alldataT2<-rowSums(sumTaxa2sample)
alldataT3<-rowSums(sumTaxa3sample)
alldataT4<-rowSums(sumTaxa4sample)
alldataT5<-as.data.frame(rowSums(sumTaxa5sample))

write.csv(alldataT5, "sum_taxa.csv")

#library("dplyr")
tsumtaxa<-as.data.frame(t(sumTaxa5sample))

abundance<-(colSums(tsumtaxa))

rareST2<-subset(sumTaxa5sample, abundance<10000)
abundST2<-subset(sumTaxa5sample, abundance>=10000)
sumrareST2<-as.data.frame(colSums(rareST2))
sumrareST2<-t(sumrareST2)
row.names(sumrareST2)<-"other"
colnames(sumrareST2)<-colnames(abundST2)
allA2<-as.data.frame(rbind(abundST2, sumrareST2))
abundST2<-abundST2[-1,]
rel<-prop.table(abundST2)

heatmap(as.matrix(abundST2,col.regions=heat.colors))

library("RColorBrewer")
library("apcluster")

library("randomForest")
set.seed(2014)
ntree<-100000
group<-variables$TIME
sample.size.per.treatment=table(group)
#predictor<-randomForest(y=as.factor(group), x=raOTU, ntree=ntree, importance=TRUE)
predictor
group<-variables$SOURCE
#predictor<-randomForest(y=as.factor(group), x=raOTU, ntree=ntree, importance=TRUE)
predictor
group<-variables$post
#predictor<-randomForest(y=as.factor(group), x=raOTU, ntree=ntree, importance=TRUE)
predictor
#write.csv(predictor$importance, "meandecrease.csv")
imp<-importance(predictor)
varImpPlot(predictor)






t4<-as.data.frame(t(sumTaxa4sample))
t4<-t4[,-1]



library("RColorBrewer")
NMDS<-metaMDS(raOTU, distance = "bray", k=2)
NMDS

col8=brewer.pal(12, "Paired")
colNew<-c('#4363d8', '#3cb44b', '#e6194b' )# '#911eb4''#f58231') #'#ffe119', '#f58231', , '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
#colNew<-c("snow1", "grey43", "grey1")
species.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)


data.scores = as.data.frame(scores(NMDS))
data.scores$Site=row.names(data.scores)
data.scores$Sample = variables$SOURCE
data.scores$Time = as.factor(variables$TIME)

head(data.scores)

xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Time, colour = Sample))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12,  colour ="black"), 
        legend.position = "right", axis.title.y = element_text(size = 20), 
        axis.title.x = element_text( size = 20, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "grey", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Sample", y = "NMDS2", shape = "Time")  + 
  scale_colour_manual(values = c("grey77", "grey40", "grey1")) 

xx

###UNIFRAC
library(ape)
require(ade4)
tree1<-read.nexus("tree_new")
tree<-as.phylo(tree1$tree_1)
tipsL<-rownames(traOTU)
pruned.tree<-drop.tip(tree,tree$tip.label[-match(tipsL, tree$tip.label)])
GUni<-GUniFrac(raOTU, pruned.tree, alpha=c(0,0.5,1))
unifracs<-GUni$unifracs
dv <- unifracs[, , "d_VAW"]
adonis(as.dist(dv) ~ variables$SOURCE*variables$TIME)
plot(hclust(as.dist(dv), method="average"),hang=-1, main='unifrac', sub='', xlab='', cex=1) #plot cluster analysis of betapair

#Specialists
traOTU<-t(raOTU)
sumhabit <-as.data.frame(subset(traOTU,rowSums(traOTU)>6))
sumOTU<-rowSums(sumhabit)
relhabit<-sumhabit/sumOTU
squ<-relhabit*relhabit
B=1/(rowSums(squ)) #formula according to Pandit et al 2009 as in Szekely & Langenheder 2014
Bn=B
rBn<-round(Bn)  #round B to full numbers
rB=as.data.frame(rBn)

summary(B)

Blow<-as.matrix(t(subset(sumhabit, rB<=4)))
Bhigh<-as.matrix(t(subset(sumhabit, rB>10)))

par(mfrow=c(1,5))
boxplot(Bhigh~as.factor(variables$TIME))

summary(lm(Bhigh~variables$TIME))

par(mfrow=c(2,1))

pasum<-sumhabit
pasum[pasum>0]=1
paBlow<-as.matrix(t(subset(pasum, rB<=4)))
paBhigh<-as.matrix(t(subset(pasum, rB>4)))
boxplot(rowSums(paBlow)~as.factor(variables$TIME))
boxplot(rowSums(paBhigh)~as.factor(variables$TIME))

summary(lm((rowSums(paBlow)/rowSums(paOTU))~variables$TIME))
summary(lm((rowSums(paBhigh)/rowSums(paOTU))~variables$TIME))
glmB<-glm((rowSums(paBhigh)/rowSums(paOTU))~variables$TIME+variables$SOURCE)
Anova(glmB)
summary(lm((rowSums(paBhigh)/rowSums(paOTU))~variables$TIME))

par(mfrow=c(2,2))
dotplot((rowSums(paBlow)/rowSums(paOTU))~variables$TIME)
dotplot((rowSums(paBhigh)/rowSums(paOTU))~variables$TIME)
paBOTU<-as.data.frame(cbind((rowSums(paBlow)/rowSums(paOTU)), (rowSums(paBhigh)/rowSums(paOTU)), variables$TIME))
colnames(paBOTU)<-c("low","high","time")
pH<-ggplot(paBOTU, aes(x=factor(time), y=high, fill = factor(variables$SOURCE)))  + geom_dotplot(binaxis = "y", stackdir = "center")
pHH<-pH+scale_fill_manual(name = NULL, labels = c("Orta", "Varese", "Effluent"),values=c("snow1", "grey43", "grey1" ))+labs(x = "time", fill="source")+theme(text = element_text(size=20),legend.position = "bottom")

grid.arrange(pRR,pHH, nrow
             =2)

grid.arrange(arrangeGrob(pRR + theme(legend.position="none"),
                               pHH ,
                               nrow=1))
#                   mylegend, nrow=2,heights=c(10, 1))
pHH
pRR

paBhigh<-as.matrix(t(subset(pasum, rB>12)))
pH<-ggplot(paBOTU, aes(x=factor(time), y=high, fill = factor(variables$SOURCE)))  + geom_dotplot(binaxis = "y", stackdir = "center", alpha=0.5)
pHH<-pH+scale_fill_manual(name = NULL, labels = c("Orta", "Varese", "Effluent"),values=c("#edae49", "violetred4", "turquoise4" ))+labs(x = "time", fill="source")+theme(text = element_text(size=20),legend.position = "bottom")

m2<-ggplot(paBOTU, aes(x=factor(time), y=high, fill = factor(variables$SOURCE))) + 
  geom_dotplot(binaxis = "y", stackdir = "center", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values=c("#edae49", "violetred4", "turquoise4"))+labs(x = "",y="", fill="source")+theme(text = element_text(size=10))+
  theme(legend.position = "none")

distN<-read.csv("dist_overtime.csv")
distoly<-read.csv("dist_0112430.csv")
tdist<-as.matrix(t(distoly))
count<-c(0,11,24,30)
plot(tdist~as.factor(count))
summary(lm(tdist~count))
count<-as.data.frame(count)
dic<-as.matrix(tdist)
rownames(dic)<-count$count
library("reshape2")

mdic<-melt(dic)
di<-ggplot(mdic, aes(x=factor(Var1), y=value))+ 
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=Var2)) +
  labs(y = "distance between samples from different sources",x="time")+
  theme(legend.position = "none", text = element_text(size=20))
plot(distoly)
plot(distN$dist0,distN$dist30)
plot(distN$dist0,distN$dist24)
plot(distN$dist24,distN$dist30)
library("gridExtra")
grid.arrange(xx,di, ncol=2)
grid.arrange(pR,pH, ncol=2)
lmdist<-lm(mdic$value~factor(mdic$Var1))
summary(lmdist)
Anova(lmdist)
marginal = emmeans(lmdist, ~ Var1)

pairs(marginal)

cor.test(distN$dist11,distN$dist30)
cor.test(distN$dist24,distN$dist30)

distoly<-read.csv("dist_within.csv")
tdist<-as.matrix(t(distoly))
count<-c(0,11,24,30)
plot(tdist~as.factor(count))
summary(lm(tdist~count))
count<-as.data.frame(count)
dic<-as.matrix(tdist)
rownames(dic)<-count$count
library("reshape2")

mdic<-melt(dic)
diR<-ggplot(mdic, aes(x=factor(Var1), y=value))+ 
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=Var2)) +
  labs(y = "distance between replicates",x="time")+
  theme(legend.position = "none", text = element_text(size=20))

grid.arrange(di,diR, ncol=2)


lmdist<-lm(mdic$value~factor(mdic$Var1))
summary(lmdist)
Anova(lmdist)
marginal = emmeans(lmdist, ~ Var1)

pairs(marginal)



###Cell numbers

cells<-read.csv("cells.csv")
row.names(cells)<-cells[,1]
cells<-cells[,-1]
##Test for stability
library(trend)

#stability all vessels from day 5-11

tsfc<-as.ts(cells[,5:8])
mk.test(log(tsfc+1))
tsfc<-as.ts(cells[1:3,5:8])
mk.test(log(tsfc+1))
tsfc<-as.ts(cells[4:6,5:8])
mk.test(log(tsfc+1))
tsfc<-as.ts(cells[7:9,5:8])
mk.test(log(tsfc+1))


###relative abundace of specific zOTUs from random forest
decr<- c("Zotu198", "Zotu365", "Zotu16","Zotu341","Zotu510","Zotu717")
decrOTU<-raOTU[,colnames(raOTU)%in%decr]
decrOTU$sum<-rowSums(decrOTU)
decrOTU$total<-rowSums(raOTU)
decrOTU$perc<-100*(decrOTU$sum/decrOTU$total)

incr<- c("Zotu1", "Zotu44", "Zotu262","zOTU566","Zotu180","Zotu390","Zotu192","Zotu32","Zotu307")
incrOTU<-raOTU[,colnames(raOTU)%in%incr]
incrOri<-as.matrix(incrOTU)
incrOTU$sum<-rowSums(incrOTU)
incrOTU$total<-rowSums(raOTU)
incrOTU$perc<-100*(incrOTU$sum/incrOTU$total)


plot(factor(variables$TIME),incrOTU$perc, main="incresing OTUs")
incrO<-cbind(incrOri,variables)
mincrori<-melt(incrO)
otOTU<-subset(mincrori, mincrori$variable%in%c("Zotu1", "Zotu44", "Zotu262","zOTU566", "Zotu180","Zotu390","Zotu192","Zotu32","Zotu307"))
otOTU<-subset(mincrori, mincrori$variable%in%c("Zotu1", "Zotu44", "Zotu262","zOTU566", "Zotu180","Zotu390","Zotu192","Zotu32","Zotu307"))

ZOTU1<-subset(mincrori, mincrori$variable=="Zotu1")

library("reshape2")
incrOrel<-100*(incrOTU/rowSums(raOTU))
varSou<-list(variables$TT)
agginc<-aggregate(incrOrel, by=varSou, mean)  
rownames(agginc)<-agginc$Group.1
agginc<-agginc[,-1]
meanAgginc<-rowMeans(agginc)
minAgginc<-rowMins(as.matrix(agginc))
maxAgginc<-rowMaxs(as.matrix(agginc))
range<-cbind(meanAgginc,minAgginc,maxAgginc)
colnames(range)<-c("mean","min","max")
abure<-cbind(variables$TIME, variables$SOURCE, incrOrel[,1:8])
colnames(abure)[1]<-"time"
colnames(abure)[2]<-"treatment"

mabure<-melt(abure, id=c("time","treatment"))
plot(as.factor(mabure$time),mabure$value)
p<-ggplot(mabure, aes(x=factor(time), y=as.numeric(value), color = factor(treatment)))  + 
  geom_boxplot(binaxis = "y", stackdir = "center")+ 
  scale_y_continuous("percentage reads")
W<-p + scale_color_manual(values=c("#edae49", "violetred4", "turquoise4" ))+labs(x = "days", color="treatment")+theme(text = element_text(size=20)) 
W + geom_boxplot(outlier.colour = NULL)

  
#Graph for families: Incr and decr
FAin<- c("f:Sphingomonadaceae")
         #, "f:Verrucomicrobiaceae", "f:Cytophagaceae", "f:Rhodospirillaceae") 
FAdi<- c("f:Pseudomonadaceae", "f:Comamonadaceae")
famdi<-t4[,colnames(t4)%in%FAin]
famdi<-as.matrix(famdi)
vafam<-cbind(famdi, variables)
mfam<-melt(vafam)
mfam<-subset(mfam, mfam$variable%in%c("f:Sphingomonadaceae", "f:Verrucomicrobiaceae", "f:Cytophagaceae", "f:Pseudomonadaceae", "f:Comamonadaceae", "f:Rhodospirillaceae"))


t41000<-t4[,colSums(t4)>1000]
vat4<-cbind.data.frame(famdi,variables$TIME,variables$SOURCE)
colnames(vat4)<-c("value","TIME","SOURCE")

###plots
FAin<-"f:Verrucomicrobiaceae"
famdi<-t4[,colnames(t4)%in%FAin]
vat4<-cbind.data.frame(famdi,variables$TIME,variables$SOURCE)
colnames(vat4)<-c("value","TIME","SOURCE")

m5<-ggplot(vat4, aes(x=TIME, y=value, color=SOURCE)) +
  geom_point(aes(size=1), alpha=0.5) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+
  scale_color_manual(values=c("#edae49", "violetred4", "turquoise4"))+labs(x = "days",y="", fill="source")+theme(text = element_text(size=10))+ggtitle("Verrucomicrobiacea")+
  theme(legend.position = "none")

FAin<-"f:Sphingomonadaceae"
famdi<-t4[,colnames(t4)%in%FAin]
vat4<-cbind.data.frame(famdi,variables$TIME,variables$SOURCE)
colnames(vat4)<-c("value","TIME","SOURCE")

m4<-ggplot(vat4, aes(x=TIME, y=value, color=SOURCE)) +
  geom_point(aes(size=1), alpha=0.5) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+
  scale_color_manual(values=c("#edae49", "violetred4", "turquoise4"))+labs(x = "",y="", fill="source")+theme(text = element_text(size=10))+ggtitle("Sphingomonadaceae")+
  theme(legend.position = "none")

FAin<-"f:Cytophagaceae"
famdi<-t4[,colnames(t4)%in%FAin]
vat4<-cbind.data.frame(famdi,variables$TIME,variables$SOURCE)
colnames(vat4)<-c("value","TIME","SOURCE")

m2<-ggplot(vat4, aes(x=TIME, y=value, color=SOURCE)) +
  geom_point(aes(size=1), alpha=0.5) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+
  scale_color_manual(values=c("#edae49", "violetred4", "turquoise4"))+labs(x = "",y="", fill="source")+theme(text = element_text(size=10))+ggtitle("Cytophagaceae")+
  theme(legend.position = "none")

FAin<-"f:Pseudomonadaceae"
famdi<-t4[,colnames(t4)%in%FAin]
vat4<-cbind.data.frame(famdi,variables$TIME,variables$SOURCE)
colnames(vat4)<-c("value","TIME","SOURCE")

m3<-ggplot(vat4, aes(x=TIME, y=value, color=SOURCE)) +
  geom_point(aes(size=1), alpha=0.5) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+
  scale_color_manual(values=c("#edae49", "violetred4", "turquoise4"))+labs(x = "",y="reads", fill="source")+theme(text = element_text(size=10))+ggtitle("Pseudomonadaceae")+
  theme(legend.position = "none")

FAin<-"f:Comamonadaceae"
famdi<-t4[,colnames(t4)%in%FAin]
vat4<-cbind.data.frame(famdi,variables$TIME,variables$SOURCE)
colnames(vat4)<-c("value","TIME","SOURCE")

m1<-ggplot(vat4, aes(x=TIME, y=value, color=SOURCE)) +
  geom_point(aes(size=1), alpha=0.5) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+
  scale_color_manual(values=c("#edae49", "violetred4", "turquoise4"))+labs(x = "",y="", fill="source")+theme(text = element_text(size=10))+ggtitle("Comamonadaceae")+
  theme(legend.position = "none")

FAin<-"f:Rhodospirillaceae"
famdi<-t4[,colnames(t4)%in%FAin]
vat4<-cbind.data.frame(famdi,variables$TIME,variables$SOURCE)
colnames(vat4)<-c("value","TIME","SOURCE")

m6<-ggplot(vat4, aes(x=TIME, y=value, color=SOURCE)) +
  geom_point(aes(size=1), alpha=0.5) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  theme_classic()+
  scale_color_manual(values=c("#edae49", "violetred4", "turquoise4"))+labs(x = "days",y="", fill="source")+theme(text = element_text(size=10))+ggtitle("Rhodospirillaceae")+
  theme(legend.position = "none")
grid.arrange(m1,m2,m3,m4,m5,m6,ncol=2,nrow=3)

