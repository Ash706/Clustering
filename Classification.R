
library(glmnet)
library(Hmisc)
library(corrplot)
library(limma)
library(dplyr)
library(tidyr)
library(gtools)
library(caret)
library(gbm)
set.seed(1)
load("Z:/AQAI - Ashfaq Ali/MEDIM/Medim_new_data/Clustering_work/Clustering_workspace.RData")
All_data<-Medim_new
All_data<-subset(Medim_new, select=-c(MEDIM_IDnr,Batch,Sample,Compound.Name, MEDIM_IDnr1))

inTrainingSet<-createDataPartition(All_data$Country, p= .5, list=FALSE)
Country_train<-All_data[inTrainingSet,]
Country_test<- All_data[-inTrainingSet,]
predictros<-subset(Country_train, select=-c(Country,MEDIM_IDnr,Batch,Sample,Compound.Name, MEDIM_IDnr1, Country))

###Model tuning using train

ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE, summaryFunction = twoClassSummary)

grid <- expand.grid(interaction.depth = seq(1,7, by = 2),
                    n.trees = seq(100, 1000, by = 50),
                    shrinkage = c(0.01,0.1))

gbmTune<- train(x = Country_train[,names(predictros)],
               y = factor(Country_train$Country, labels = c("Iraq", "Sweden")),
                method = "gbm",
               metric = "ROC",
               tuneGrid = grid,
               verbose = FALSE,
               trControl = ctrl
              )
