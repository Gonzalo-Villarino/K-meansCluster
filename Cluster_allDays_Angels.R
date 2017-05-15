rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_Mimulus/HeatMaps/FDR001_up_dn/")
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
all2 =merge(all,Day04,by=c(1), all=TRUE)
all3=merge(all2,Day06,by=c(1), all=TRUE)
all4=merge(all3,Day08,by=c(1), all=TRUE)

all4[is.na(all4)] <- 0


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

# with 15 cluster one is empty 
cl <- kmeans(all4_scaled[,2:6], 14, iter.max = 25, nstart = 100)
cl.cluster<-as.factor(cl$cluster)

clustT1WIN2<- data.frame (all4_scaled[,2:6], cl.cluster) 
clustT1WIN2$.row=rownames(all4_scaled)


molten2 <- melt(clustT1WIN2, id = c(".row", "cl.cluster") )
names(molten2) = c("gp","cl.cluster","variable","value")

pcp_cl = ggplot(molten2,aes(variable, value,group= gp)) + geom_line() + stat_summary(fun.y = mean,geom = "line",colour= "red",size=1.5,aes(group=cl.cluster)) + facet_wrap(~cl.cluster) 

pcp_cl


pcp_cl+ theme_minimal(base_size = 11, base_family = "")


pcp_cl + theme_minimal()axis.text.x = element_text(angle = 90, hjust = 1, size=13,color="darkred"))
pcp_cl + scale_color_manual(values=c("Red"))

####select genes in a specific cluster########

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


 test=cl4_gene[which(cl4_gene[,1]=="Migut.A00116"),]


#which(cl1_gene$gene_short_name == "Migut.A00116")
#which(cl2_gene$gene_short_name == "Migut.A00116")
#which(cl3_gene$gene_short_name == "Migut.A00116")
test=which(cl4_gene$gene_short_name == "Migut.A00116")
#which(cl5_gene$gene_short_name == "Migut.A00116")
#which(cl6_gene$gene_short_name == "Migut.A00116")
#which(cl7_gene$gene_short_name == "Migut.A00116")
#which(cl8_gene$gene_short_name == "Migut.A00116")
#which(cl9_gene$gene_short_name == "Migut.A00116")
#which(cl10_gene$gene_short_name == "Migut.A00116")
#which(cl11_gene$gene_short_name == "Migut.A00116")
#which(cl12_gene$gene_short_name == "Migut.A00116")
#which(cl13_gene$gene_short_name == "Migut.A00116")
#which(cl14_gene$gene_short_name == "Migut.A00116")
