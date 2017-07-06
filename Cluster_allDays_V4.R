rm(list=ls())
setwd("/Users/angels/Google Drive/Sozzani's lab/Franks' project/clustering_miguel")
getwd()

library (pheatmap)
library (NMF)
library(dendextend)
library(d3heatmap)
library(ape)
library(gplots)
library(ggplot2)
library(cluster)
library(fpc)
library(reshape)
library(matrixStats)

# get data
Day00= read.table(file="d00_FDR01_RPKM.txt",header=TRUE, sep="\t")
Day00$d00 <- (Day00$Day00_BR1 + Day00$Day00_BR4 + Day00$Day00_BR5 + Day00$Day00_BR6  + Day00$Day00_BR7)/5
Day00=Day00[,c(-2:-6)]

Day02= read.table(file="d02_FDR01_RPKM.txt",header=TRUE, sep="\t")
Day02$d02 <- (Day02$Day02d_BR4 + Day02$Day02d_BR5 + Day02$Day02d_BR6 + Day02$Day02d_BR7)/4
Day02=Day02[,c(-2:-5)]

Day04= read.table(file="d04_FDR01_RPKM.txt",header=TRUE, sep="\t")
Day04$d04 <- (Day04$Day04_BR2 +  Day04$Day04_BR3 +  Day04$Day04_BR4)/3
Day04=Day04[,c(-2:-4)]

Day06= read.table(file="d06_FDR01_RPKM.txt",header=TRUE, sep="\t")
Day06$d06 <- (Day06$Day06_BR1 +Day06$Day06_BR2 + Day06$Day06_BR3 + Day06$Day06_BR4)/4
Day06=Day06[,c(-2:-5)]

Day08= read.table(file="d08_FDR01_RPKM.txt",header=TRUE, sep="\t")
Day08$d08 <- (Day08$Day08_BR1+ Day08$Day08_BR2+Day08$Day08_BR3+Day08$Day08_BR4)/4
Day08=Day08[,c(-2:-5)]

#merge all data
all = merge(Day00,Day02,by=c(1), all=TRUE)
all1 =merge(all,Day04,by=c(1), all=TRUE)
all2=merge(all1,Day06,by=c(1), all=TRUE)
all3=merge(all2,Day08,by=c(1), all=TRUE)

all3[is.na(all3)] <- 0
all3 = unique(all3)


# read RPKM values all genes
Day00_all= read.table(file="d00_ALL_RPKM.txt",header=TRUE, sep="\t")
Day00_all$d00 <- (Day00_all$Day00_BR1 + Day00_all$Day00_BR4 + Day00_all$Day00_BR5 + Day00_all$Day00_BR6  + Day00_all$Day00_BR7)/5
Day00_all=Day00_all[,c(-2:-6)]

Day02_all= read.table(file="d02_ALL_RPKM.txt",header=TRUE, sep="\t")
Day02_all$d02 <- (Day02_all$Day02d_BR4 + Day02_all$Day02d_BR5 + Day02_all$Day02d_BR6 + Day02_all$Day02d_BR7)/4
Day02_all=Day02_all[,c(-2:-5)]

Day04_all= read.table(file="d04_ALL_RPKM.txt",header=TRUE, sep="\t")
Day04_all$d04 <- (Day04_all$Day04_BR2 +  Day04_all$Day04_BR3 +  Day04_all$Day04_BR4)/3
Day04_all=Day04_all[,c(-2:-4)]

Day06_all= read.table(file="d06_ALL_RPKM.txt",header=TRUE, sep="\t")
Day06_all$d06 <- (Day06_all$Day06_BR1 +Day06_all$Day06_BR2 + Day06_all$Day06_BR3 + Day06_all$Day06_BR4)/4
Day06_all=Day06_all[,c(-2:-5)]

Day08_all= read.table(file="d08_ALL_RPKM.txt",header=TRUE, sep="\t")
Day08_all$d08 <- (Day08_all$Day08_BR1+ Day08_all$Day08_BR2+Day08_all$Day08_BR3+Day08_all$Day08_BR4)/4
Day08_all=Day08_all[,c(-2:-5)]


#merge all data
RPKM_all = merge(Day00_all,Day02_all,by=c(1), all=TRUE)
RPKM_all = unique(RPKM_all)
RPKM_all2 =merge(RPKM_all,Day04_all,by=c(1), all=TRUE)
RPKM_all2 = unique(RPKM_all2)
RPKM_all3=merge(RPKM_all2,Day06_all,by=c(1), all=TRUE)
RPKM_all3 = unique(RPKM_all3)
RPKM_all4=merge(RPKM_all3,Day08_all,by=c(1), all=TRUE)
RPKM_all4 = unique(RPKM_all4)
RPKM_all4[is.na(RPKM_all4)] <- 0

all3 = all3[,1, drop=FALSE]
all4 =merge(all3,RPKM_all4,by=c(1), all=FALSE)


###############################################################
# k-means cluster
###############################################################


###remove dup rows ###
all4 = unique(all4)
########remove a cluter that only have two genes: Migut.E01377, Migut.M01044, it seems they are outlier.
#all4 = all4[-1*which(all4[,1]=="Migut.E01377" | all4[,1]=="Migut.M01044"),]
#all4 = all4[-1*which(all4[,1]=="Migut.M01349" | all4[,1]=="Migut.J00046"),]

###scale data#######
all4_num = t(all4[,2:6])
all4_scaled = matrix(1,dim(all4_num)[1],dim(all4_num)[2])

for (i in 1:15097){
  m = mean(all4_num[,i])
  s = sd(all4_num[,i])
  all4_scaled[,i] = (all4_num[,i]-m)/(s+0.000001)
}
all4_scaled = t(all4_scaled)  
all4_scaled = as.data.frame(all4_scaled)

# all4_scaled = as.data.frame(scale(all4[,2:6]))
# 
all4_scaled = cbind(all4[,1],all4_scaled)
names(all4_scaled) = names(all4)


####select best k for clustering #############
#library(fpc)
#library(cluster)
#asw <- numeric(20)
#for (k in 2:20)
#  asw[[k]] <- pam(all4_scaled[2:6], k) $ silinfo $ avg.width
#k.best <- which.max(asw)

#plot(asw[1:20], xlab="K",ylab="Average Silhouette",type="b",pch=16)

###
cl <- kmeans(all4_scaled[,2:6], 30, iter.max = 25, nstart = 100)
cl.cluster<-as.factor(cl$cluster)

clustT1WIN2<- data.frame (all4_scaled[,2:6], cl.cluster) 
clustT1WIN2$.row=rownames(all4_scaled)


molten2 <- melt(clustT1WIN2, id = c(".row", "cl.cluster") )
names(molten2) = c("gp","cl.cluster","variable","value")

pcp_cl = ggplot(molten2,aes(variable, value,group= gp)) + geom_line() + stat_summary(fun.y = mean,geom = "line",colour= "red",size=1.5,aes(group=cl.cluster)) + facet_wrap(~cl.cluster) 

pcp_cl


pcp_cl+ theme_minimal(base_size = 11, base_family = "")

####select genes in a specific cluster########
###e.g. cluster 2
cl1_gene = all4_scaled[cl.cluster==1,]
cl2_gene = all4_scaled[cl.cluster==2,]
cl3_gene = all4_scaled[cl.cluster==3,]
cl4_gene = all4_scaled[cl.cluster==4,]
cl5_gene = all4_scaled[cl.cluster==5,]
cl6_gene = all4_scaled[cl.cluster==6,]
cl7_gene = all4_scaled[cl.cluster==7,]
cl8_gene = all4_scaled[cl.cluster==8,]
cl9_gene = all4_scaled[cl.cluster==9,]
cl10_gene = all4_scaled[cl.cluster==10,]
cl11_gene = all4_scaled[cl.cluster==11,]
cl12_gene = all4_scaled[cl.cluster==12,]
cl13_gene = all4_scaled[cl.cluster==13,]
cl14_gene = all4_scaled[cl.cluster==14,]
cl15_gene = all4_scaled[cl.cluster==15,]
cl16_gene = all4_scaled[cl.cluster==16,]
cl17_gene = all4_scaled[cl.cluster==17,]
cl18_gene = all4_scaled[cl.cluster==18,]
cl19_gene = all4_scaled[cl.cluster==19,]
cl20_gene = all4_scaled[cl.cluster==20,]
cl21_gene = all4_scaled[cl.cluster==21,]
cl22_gene = all4_scaled[cl.cluster==22,]
cl23_gene = all4_scaled[cl.cluster==23,]
cl24_gene = all4_scaled[cl.cluster==24,]
cl25_gene = all4_scaled[cl.cluster==25,]
cl26_gene = all4_scaled[cl.cluster==26,]
cl27_gene = all4_scaled[cl.cluster==27,]
cl28_gene = all4_scaled[cl.cluster==28,]
cl29_gene = all4_scaled[cl.cluster==29,]
cl30_gene = all4_scaled[cl.cluster==30,]

write.table(cl1_gene,"cl1")
write.table(cl2_gene,"cl2")
write.table(cl3_gene,"cl3")
write.table(cl4_gene,"cl4")
write.table(cl5_gene,"cl5")
write.table(cl6_gene,"cl6")
write.table(cl7_gene,"cl7")
write.table(cl8_gene,"cl8")
write.table(cl9_gene,"cl9")
write.table(cl10_gene,"cl10")
write.table(cl11_gene,"cl11")
write.table(cl12_gene,"cl12")
write.table(cl13_gene,"cl13")
write.table(cl14_gene,"cl14")
write.table(cl15_gene,"cl15")
write.table(cl16_gene,"cl16")
write.table(cl17_gene,"cl17")
write.table(cl18_gene,"cl18")
write.table(cl19_gene,"cl19")
write.table(cl20_gene,"cl20")
write.table(cl21_gene,"cl21")
write.table(cl22_gene,"cl22")
write.table(cl23_gene,"cl23")
write.table(cl24_gene,"cl24")
write.table(cl25_gene,"cl25")
write.table(cl26_gene,"cl26")
write.table(cl27_gene,"cl27")
write.table(cl28_gene,"cl28")
write.table(cl29_gene,"cl29")
write.table(cl30_gene,"cl30")

