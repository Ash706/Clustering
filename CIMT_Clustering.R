set.seed(17)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(corrplot)
library(limma)
library(tidyr)
library(gtools)
Metabolomics_data <- read.csv("Z:/AQAI - Ashfaq Ali/CIMT_trial/Metabolomics_data.csv")
Metabolomics_data_2<-Metabolomics_data[!(is.na(Metabolomics_data$grp)),]

##### Extract data for clustering 
## Select baseline data, metformin arm
CIMT_for_clustering<-filter(Metabolomics_data_2, grp==2) %>% filter(.,visit=="v1")
rownames(CIMT_for_clustering)<-CIMT_for_clustering$subjid
number<-c(1:205)
disc<-sample(number, 68)
CIMT_for_clustering_discovery<- CIMT_for_clustering[ number %in% disc,]
CIMT_for_clustering_discovery<-CIMT_for_clustering_discovery[,51:148]
Clust_clinical_data_CIMT<-cbind(CIMT_for_clustering[rownames(X.cluster.mean),], X.cluster.mean)
write.csv(CIMT_for_clustering_discovery,"Z:/AQAI - Ashfaq Ali/CIMT_trial/Data_for_clustering_CIMT.csv" , row.names = TRUE)


CIMT_visit2_clustered<-filter(Metabolomics_data_2, grp==2) %>% filter(visit=="v7") 
rownames(CIMT_visit2_clustered)<-CIMT_visit2_clustered$subjid
CIMT_visit2_clustered<-CIMT_visit2_clustered[rownames(X.cluster.mean),]
colnames(CIMT_visit2_clustered)<-paste(colnames(CIMT_visit2_clustered), "V7", sep="")

CIMT_v1v2_clusters<-(data.frame(CIMT_visit2_clustered,Clust_clinical_data_CIMT))
CIMT_v1v2_clusters$detlatHbA1c<-CIMT_v1v2_clusters$hba1cV7-CIMT_v1v2_clusters$hba1c
CIMT_v1v2_clusters$
CIMT_v1v2_clusters$percent_response<-((CIMT_v1v2_clusters$detlatHbA1c)/(CIMT_v1v2_clusters$hba1c))*100
CIMT_v1v2_clusters$response_cat<-ifelse(CIMT_v1v2_clusters$percent_response >= (-10), 0, 1 )

################# Test association with CIMT_v1v2_clusters
#remove NAs containing rows in y varibles (detlatHbA1c)
data<-CIMT_v1v2_clusters[!is.na(CIMT_v1v2_clusters$detlatHbA1c),]
Datam<-t(data[,grep("Cluster",names(data))])
design1<-model.matrix(~detlatHbA1c+hba1c+pre.ins.1,data=data) #data without nas
fit1<-lmFit(Datam,design1) #data is metabolotes
fit1 <- eBayes(fit1)
deltahba1c_sig_clust<-topTable(fit1,coef=2,number=20, p.value = 1)  #significant clusters

data1<-CIMT_v1v2_clusters[!is.na(CIMT_v1v2_clusters$detlatHbA1c),]
Datam1<-t(data[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~detlatHbA1c+hba1c,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
deltahba1c_sig_met<-topTable(fit2,coef=2,number=20, p.value = 1)  #significant clusters

data1<-CIMT_v1v2_clusters[!is.na(CIMT_v1v2_clusters$hba1c),]
Datam1<-t(data1[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~hba1c,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
blhba1c_sig_met_pl<-topTable(fit2,coef=2,number=20, p.value = 1)


data1<-CIMT_v1v2_clusters[!is.na(CIMT_v1v2_clusters$hba1cV7),]
Datam1<-t(data1[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~hba1c,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
fuphba1c_sig_met<-topTable(fit2,coef=2,number=20, p.value = 1)

##### Analyses of the whole dataset
CIMT_whole<-Metabolomics_data_2
CIMT_v1<-filter(Metabolomics_data_2, visit=="v1") 
CIMT_v7<-filter(Metabolomics_data_2,visit=="v7")
colnames(CIMT_v7)<-paste(names(CIMT_v7), "v7", sep = "")
CIMT_both<-merge(CIMT_v7, CIMT_v1, by.x = "subjidv7",by.y = "subjid")
CIMT_both$deltahba1c<-CIMT_both$hba1cv7-CIMT_both$hba1c
CIMT_both$percent_response<-((CIMT_both$deltahba1c)/(CIMT_both$hba1c))*100
CIMT_both$response_cat<-ifelse(CIMT_both$percent_response >= (-10), 0, 1 )
CIMT_both$non_response_cat<-ifelse(CIMT_both$percent_response >= (-10), 1, 0 )
#Analyses of the whole metabolomics data

data1<-CIMT_v7[!is.na(CIMT_v7$hba1cv7),]
data1<-filter(data1, grp=="2")
Datam1<-t(data1[,paste(names(CIMT_for_clustering_discovery),"v7",sep="")])
design2<-model.matrix(~log(hba1cv7),data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
fuphba1c_sig_met<-topTable(fit2,coef=2,number=50, p.value = 1)

######## Re analyses
###Test for associations with delta hba1c with baseline metabolites overall

data1<-CIMT_both[!is.na(CIMT_both$deltahba1c),]
#data1<-filter(data1, grp=="1")
Datam1<-t(data1[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~deltahba1c+hba1c+pre.ins.1,data=data1) #data without nas+pre.ins.1
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
deltahba1c_sig_met_grp12<-topTable(fit2,coef=2,adjust.method="fdr",number=nrow(fit2), p.value = 1)

data1<-CIMT_both[!is.na(CIMT_both$deltahba1c),]
data1<-filter(data1, grp=="2")
Datam1<-t(data1[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~deltahba1c+hba1c,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
deltahba1c_sig_met_grp2<-topTable(fit2,coef=2,adjust.method="fdr",number=nrow(fit2), p.value = 1)

data1<-CIMT_both[!is.na(CIMT_both$deltahba1c),]
data1<-filter(data1, grp=="2")
Datam1<-t(data1[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~deltahba1c+hba1c+pre.ins.1,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
deltahba1c_sig_met_grp2_pre<-topTable(fit2,coef=2,adjust.method="fdr",number=nrow(fit2), p.value = 1)

#### Test 
data1<-CIMT_both[!is.na(CIMT_both$response_cat),]
data1<-filter(data1, grp=="1")
Datam1<-t(data1[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~response_cat+hba1c,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
percent_sig_met_grp1<-topTable(fit2,coef=2,number=50, p.value = 1)

data1<-CIMT_both[!is.na(CIMT_both$hba1c),]
#data1<-filter(data1, grp=="1")
Datam1<-t(data1[,names(CIMT_for_clustering_discovery)])
design2<-model.matrix(~hba1c,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
blhba1c_all<-topTable(fit2,coef=2,number=50, p.value = 1)


#Reanalyses of the clustering data
data<-CIMT_v1v2_clusters[!is.na(CIMT_v1v2_clusters$hba1c),]
#data<-filter(data, grp=="1")
Datam<-t(data[,grep("Cluster", names(data))])
design1<-model.matrix(~log(hba1c),data=data) #data without nas+pre.ins.1
fit1<-lmFit(Datam,design1) #data is metabolotes
fit1 <- eBayes(fit1)
blhba1c_sig_clust<-topTable(fit1,coef=2,number=20, p.value = 1)  #significant clusters


data1<-CIMT_both[!is.na(CIMT_both$hba1cv7),]
#data1<-filter(data1, grp=="1")
Datam1<-t(data1[,paste(names(CIMT_for_clustering_discovery),"v7",sep="")])
design2<-model.matrix(~log(hba1cv7),data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
v7hba1c_all<-topTable(fit2,coef=2,number=50, p.value = 1)

data1<-CIMT_both[!is.na(CIMT_both$grp),]
#data1<-filter(data1, grp=="1")
Datam1<-t(data1[,paste(names(CIMT_for_clustering_discovery),"v7",sep="")])
design2<-model.matrix(~grp,data=data1) #data without nas
fit2<-lmFit(Datam1,design2) #data is metabolotes
fit2 <- eBayes(fit2)
grphba1c_all<-topTable(fit2,coef=2,number=50, p.value = 1)

#Reanalyses of the clustering data
data<-CIMT_v1v2_clusters[!is.na(CIMT_v1v2_clusters$hba1c),]
#data<-filter(data, grp=="1")
Datam<-t(data[,grep("Cluster", names(data))])
design1<-model.matrix(~log(hba1c),data=data) #data without nas+pre.ins.1
fit1<-lmFit(Datam,design1) #data is metabolotes
fit1 <- eBayes(fit1)
blhba1c_sig_clust<-topTable(fit1,coef=2,number=20, p.value = 1)  #significant clusters



number<-!number[disc]

sample_n(CIMT_for_clustering, )
plot_1<-ggplot(data=(Metabolomics_data_2), aes(x=visit,y=hba1c, group=subjid))+geom_point()+geom_line()+facet_grid(.~grp)
plot_1

plot_1<-ggplot(data=(Metabolomics_data_2), aes(x=visit,y=hba1c,group=subjid))+geom_point()+geom_line()+facet_grid(pre.ins ~grp )
plot_1
