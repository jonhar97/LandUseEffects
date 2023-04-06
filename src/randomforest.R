library(dplyr)
library(randomForest)
library(caret)
library(e1071)
require(ggplot2)
require(reshape2)
require(ggpubr)
require(readr)
library(doParallel)
require(ranger)
require(tibble)
library("cowplot")
set.seed(1234)

# Step 1: Import the data
path <- "C:/Users/joah/OneDrive - Skogforsk/Documents/Projekt/Old/Reservvatten/LandUse/"

# load OTU and taxonomy tables
load(file = paste(path,"otu_table_psno.rds",sep=""))
load(file = paste(path,"tax_table_psno.rds",sep=""))
metadata <- read_tsv(file = paste(path,"metadata.tsv",sep="")) # load sample metadata
# add columns
metadata <- metadata %>%
  mutate(SourceSink=rep("source",nrow(metadata))) %>%
  mutate(Env= if_else(.$Type == "Point", 'Inflow', 'Lake'))
# rename otus
rn <- rownames(otu_full)
nn <- paste("asv_",seq(1,nrow(otu_full)),sep="")# new name
df <-cbind(rn,nn)
colnames(df) <- c("Original","Recoded")
write.table(file = paste(path,"key_recoded_asvs_MST.txt",sep=""),x = df,quote=F,sep="\t")
rownames(otu_full) <- nn
rownames(tax_full) <- nn
######## keep top abundant asvs ############
top_limit <- 100
otus <- otu_full[order(rowSums(otu_full),decreasing = TRUE)[1:top_limit],]
tax2 <- tax_full[order(rowSums(otu_full),decreasing = TRUE)[1:top_limit],]
otus <- as.matrix(t(otus))

# extract only those samples in common between the two tables
common.sample.ids <- intersect(metadata$SampleID, rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[metadata$SampleID %in% common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                    'between the metadata file and data table')
    stop(message)
}
# loop across sink samples (i.e. lake samples)
############## step 1 optimize hyperparameters #################
samplenames <- metadata %>%
  filter(Env == "Lake") %>%
  select(SampleID) 
write.table(x = metadata,file = paste(path,"metadata_MST.txt",sep = ""),quote = FALSE,sep = "\t") # Write edited metadata to file

data.MST <- data.frame(otus,Env = metadata$Env) # MST
data.MST <- data.MST %>%
  mutate(across(Env,as.factor))
write.table(x = data.MST,file = paste(path,"data_MST.txt",sep=""),sep="\t",quote=F,row.names = TRUE)
data.CST <- metadata %>%
  select(c(SampleID,Turbidity,ColorVal,pH,Alkalinity,CODMn,Conductivity,TotalHardness,Cloride,Sulfate,Ammonium,         
         Ammonium_nitrogen,Nitrate,Nitrate_nitrogen,Sodium,Potassium,Calcium,Iron,Magnesium,Manganese,
         Aluminium,Copper,SourceSink,Env)) %>%
  na.omit()
data.CST[data.CST == 0] <- 0.001
data.CST <- data.CST %>%
  mutate_at(vars(Turbidity,ColorVal,pH,Alkalinity,CODMn,Conductivity,TotalHardness,Sulfate,Ammonium,         
                 Ammonium_nitrogen,Nitrate,Nitrate_nitrogen,Calcium,Iron,Magnesium,Manganese,
                 Aluminium,Copper), log) %>%
  mutate_at(vars(Turbidity,ColorVal,pH,Alkalinity,CODMn,Conductivity,TotalHardness,Cloride,Sulfate,Ammonium,         
                 Ammonium_nitrogen,Nitrate,Nitrate_nitrogen,Sodium,Potassium,Calcium,Iron,Magnesium,Manganese,
                 Aluminium,Copper),scale) %>%
  column_to_rownames('SampleID')  %>%
  select(-c('SourceSink')) %>%
  mutate(across(Env,as.factor)) %>%
  as.data.frame(.)
write.table(x = data.CST,file = paste(path,"data_CST.txt",sep=""),sep="\t",quote=F,row.names = TRUE)

sensitivity.df <- expand.grid(
  num.trees = c(250, 500, 1000, 1500),
  mtry = c(2,3,4,5),
  min.node.size = c(1, 5, 10, 20)
)
sensitivity.df$path.data.MST = paste(path,"data_MST.txt",sep="")
sensitivity.df$path.data.CST = as.character(paste(path,"data_CST.txt",sep=""))

kableExtra::kable(sensitivity.df)
# set up parallel backend
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#fitting each rf model with different hyperparameters
prediction.error <- foreach(
  num.trees = sensitivity.df$num.trees,
  mtry = sensitivity.df$mtry,
  min.node.size = sensitivity.df$min.node.size,
  path.data.MST = sensitivity.df$path.data.MST,
  path.data.CST = sensitivity.df$path.data.CST,
  .verbose = TRUE,
  .combine = 'rbind', 
  .packages = "randomForest" #"ranger"
) %dopar% {

  data.MST <- read.table(file = path.data.MST,sep="\t",header=T,row.names = 1)
  data.CST <- read.table(file = path.data.CST,sep="\t",header=T,row.names = 1)
  data.MST$Env <- as.factor(data.MST$Env)
  data.CST$Env <- as.factor(data.CST$Env)
  #fit the random forest model
  m.i <- randomForest(
    formula = Env ~ .,
    ntree = num.trees,
    mtry = mtry,
    nodesize = min.node.size,
    data = data.MST
  )  
  m.ii <- randomForest(
    formula = Env ~ .,
    ntree = num.trees,
    mtry = mtry,
    nodesize = min.node.size,
    data = data.CST
  )

  #returning prediction error as percentage
  return(c(m.i$err.rate[num.trees,1] * 100,m.ii$err.rate[num.trees,1] * 100)) # 
  
}
# stop cluster
parallel::stopCluster(cl = my.cluster) 
#adding the prediction error column
sensitivity.df <- cbind(sensitivity.df,prediction.error)
colnames(sensitivity.df)[(ncol(sensitivity.df)-1):ncol(sensitivity.df)] <- c("OOB_MST","OOB_CST")
cat(sprintf('Total number of optima for MST: %d\n',length(which(sensitivity.df$OOB_MST==min(sensitivity.df$OOB_MST)))))
cat(sprintf('Total number of optima for CST: %d\n',length(which(sensitivity.df$OOB_CST==min(sensitivity.df$OOB_CST)))))
cat(sprintf('Lowest MST prediction error: %f\n',min(sensitivity.df$OOB_MST)))
cat(sprintf('Lowest CST prediction error: %f\n',min(sensitivity.df$OOB_CST)))
cat(sprintf('Optimal setting MST\n'))
print(sensitivity.df[which(sensitivity.df$OOB_MST==min(sensitivity.df$OOB_MST)),c(1,2,3,6,7)])
cat(sprintf('Optimal setting CST\n'))
print(sensitivity.df[which(sensitivity.df$OOB_CST==min(sensitivity.df$OOB_CST)),c(1,2,3,6,7)])
cat(sprintf('Checking for joint optima\n'))
if(length(intersect(which(sensitivity.df$OOB_MST==min(sensitivity.df$OOB_MST)),which(sensitivity.df$OOB_CST==min(sensitivity.df$OOB_CST))))>0){
  print(sensitivity.df[intersect(which(sensitivity.df$OOB_MST==min(sensitivity.df$OOB_MST)),which(sensitivity.df$OOB_CST==min(sensitivity.df$OOB_CST))),c(1,2,3,6,7)])
}
############## Step 2: Use optimized hyperparameter set and loop across samples ############
# MST
sensitivity.df.new <- expand.grid(
  sample = as.data.frame(samplenames)
)
index <- which(sensitivity.df$OOB_MST==min(sensitivity.df$OOB_MST))
sensitivity.df.new$mtry <- sensitivity.df$mtry[index[1]]
sensitivity.df.new$num.trees <- sensitivity.df$num.trees[index[1]]
sensitivity.df.new$min.node.size <- sensitivity.df$min.node.size[index[1]]
sensitivity.df.new$path.data.MST = paste(path,"data_MST.txt",sep="")
# set up parallel backend
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#fitting each rf model with different hyperparameters
res.mst <- foreach(
  sample = sensitivity.df.new$SampleID,
  num.trees = sensitivity.df.new$num.trees,
  min.node.size = sensitivity.df.new$min.node.size,
  mtry = sensitivity.df.new$mtry,
  path.data.MST = sensitivity.df$path.data.MST,
  .verbose = TRUE,
  .combine = 'rbind', 
  .packages = "randomForest"
) %dopar% {
  # read data
  data.MST <- read.table(path.data.MST,sep='\t', header=T,row.names=1)

  # extract the source environments and source/sink indices
  train.ix <- which(rownames(data.MST)!=sample)
  test.ix <- which(rownames(data.MST)==sample)
  Y <- as.factor(data.MST$Env)
  X <- data.MST[,-ncol(data.MST)]
  #fit model
  rf <- randomForest(x = X[train.ix,], # predictors (i.e. OTUs)
        y = Y[train.ix], # A response vector of sources
        ntree=num.trees, # Number of trees to grow
        keep.forest=TRUE,
        importance=FALSE,
        nodesize = min.node.size, # Minimum size of terminal nodes
        mtry = mtry # Number of variables randomly sampled as candidates at each split
  )
  
  # predict sink sample source
  # Use the model on test data
  pred <- predict(object = rf,X[c(1,test.ix),],type="prob")
  
  # Save results
  #returning prediction estimates
  return(pred[2,]) # ,ntree,mtry,nodesize,notu) c(rownames(pred)[2],
  
}
# stop cluster
parallel::stopCluster(cl = my.cluster) 
sensitivity.df.new$Inflow <- res.mst[,1]
sensitivity.df.new$Lake <- res.mst[,2]
##########################################
# Chem phys data
# first remove non-overlapping samples
rn <- as.data.frame(rownames(data.CST)) 
colnames(rn) <- "temp"
print(samplenames)
samplenames <- samplenames$SampleID[-which(samplenames$SampleID %in% rn$temp == FALSE)]
intersect(as.data.frame(samplenames$SampleID),as.data.frame(rn$temp))
sensitivity.df.new2 <- expand.grid(
  sample = as.data.frame(samplenames)
)
index <- which(sensitivity.df$OOB_CST==min(sensitivity.df$OOB_CST))
sensitivity.df.new2$mtry <- sensitivity.df$mtry[index[1]]
sensitivity.df.new2$num.trees <- sensitivity.df$num.trees[index[1]]
sensitivity.df.new2$min.node.size <- sensitivity.df$min.node.size[index[1]]
sensitivity.df.new2$path.data.CST = paste(path,"data_CST.txt",sep="")
# set up parallel backend
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#fitting each rf model with different hyperparameters
res.cst <- foreach(
  sample = sensitivity.df.new2$samplenames,
  num.trees = sensitivity.df.new2$num.trees,
  min.node.size = sensitivity.df.new2$min.node.size,
  mtry = sensitivity.df.new2$mtry,
  path.data.CST = sensitivity.df.new2$path.data.CST,
  .verbose = TRUE,
  .combine = 'rbind', 
  .packages = "randomForest"
) %dopar% {
  # read data
  data.CST <- read.table(path.data.CST,sep='\t', header=T,row.names=1)
  
  # extract the source environments and source/sink indices
  train.ix <- which(rownames(data.CST)!=sample)
  test.ix <- which(rownames(data.CST)==sample)
  Y <- as.factor(data.CST$Env)
  X <- data.CST[,-ncol(data.CST)]
  #fit model
  rf <- randomForest(x = X[train.ix,], # predictors (i.e. physiochem variables)
                     y = Y[train.ix], # A response vector of sources
                     ntree=num.trees, # Number of trees to grow
                     keep.forest=TRUE,
                     importance=FALSE,
                     nodesize = min.node.size, # Minimum size of terminal nodes
                     mtry = mtry # Number of variables randomly sampled as candidates at each split
  )
  
  # predict sink sample source
  # Use the model on test data
  pred <- predict(object = rf,X[c(1,test.ix),],type="prob")
  
  # Save results
  #returning prediction estimates
  return(pred[2,]) # ,ntree,mtry,nodesize,notu) c(rownames(pred)[2],
  
}
# stop cluster
parallel::stopCluster(cl = my.cluster) 
sensitivity.df.new2$Inflow <- res.cst[,1]
sensitivity.df.new2$Lake <- res.cst[,2]



###################################################
# Find variable importance

#fit the optimal random forest model
m.1 <- randomForest(
  formula = Env ~ .,
  ntree = sensitivity.df.new$num.trees[1],
  mtry = sensitivity.df.new$mtry[1],
  nodesize = sensitivity.df.new$min.node.size[1],
  importance=TRUE,
  data = data.MST
)  
m.2 <- randomForest(
  formula = Env ~ .,
  ntree = sensitivity.df.new2$num.trees[1],
  mtry = sensitivity.df.new2$mtry[1],
  nodesize = sensitivity.df.new2$min.node.size[1],
  importance=TRUE,
  data = data.CST
)
varImpPlot(m.1,main = "Imporance of ASVs for predicting Lake/Inflow environments",n.var = 15)
varImpPlot(m.2,main = "Imporance of Chem variables for predicting Lake/Inflow environments",n.var = 15)
png("variable_importance_top15_ASVs.png")
varImpPlot(m.1,main = "Imporance of ASVs for predicting Lake/Inflow environments",n.var = 15)
dev.off()

png("variable_importance_top15_Chems.png")
varImpPlot(m.2,main = "Imporance of Chem variables for predicting Lake/Inflow environments",n.var = 15)
dev.off()
asvs <- as.data.frame(m.1$importance)
asvs <- cbind(asvs,m.1$importanceSD[,'MeanDecreaseAccuracy'])
colnames(asvs)[5] <- "AccuracySD"
tax.asvs <- tax_full[rownames(asvs),]
imp.asvs <- cbind(asvs,tax.asvs)
imp.asvs <- imp.asvs %>%
  rownames_to_column(.)
imp.asvs <- imp.asvs[order(imp.asvs$MeanDecreaseAccuracy,decreasing = TRUE)[1:15],]
imp.asvs$rowname <- factor(imp.asvs$rowname,                                    # Factor levels in decreasing order
                  levels = imp.asvs$rowname[order(imp.asvs$MeanDecreaseAccuracy, decreasing = FALSE)])
p<-ggplot(data=imp.asvs, aes(x=rowname,y=MeanDecreaseAccuracy,fill=Genus)) +
  geom_bar(stat="identity", color="black") + 
  coord_flip() + 
  theme(legend.position = "top",legend.text = element_text(size = 10, color = "black"),axis.text.x = element_text(size=10,color="black"),axis.title.y=element_blank(),axis.text.y = element_blank()) +
  geom_errorbar(aes(ymin=MeanDecreaseAccuracy-AccuracySD, max=MeanDecreaseAccuracy+AccuracySD), width=.2, position=position_dodge(width=0.7)) +
  ylab("Variable importance"); p # + 
#  scale_x_discrete(limits = imp.asvs$Position); p # Horizontal bar plot
chem <- as.data.frame(m.2$importance)
chem <- cbind(chem,m.2$importanceSD[,'MeanDecreaseAccuracy'])
colnames(chem)[5] <- "AccuracySD"
chem <- chem %>%
  rownames_to_column(.)
chem <- chem[order(chem$MeanDecreaseAccuracy,decreasing = TRUE)[1:15],]
chem$rowname <- factor(chem$rowname,                                    # Factor levels in decreasing order
                           levels = chem$rowname[order(chem$MeanDecreaseAccuracy, decreasing = FALSE)])
p2<-ggplot(data=chem, aes(x=rowname,y=MeanDecreaseAccuracy)) +
  geom_bar(stat="identity", color="black") + 
  coord_flip() + 
  theme(axis.title.y=element_blank(),axis.text.x = element_text(size=10,color="black"),axis.text.y = element_text(color="black")) +
  geom_errorbar(aes(ymin=MeanDecreaseAccuracy-AccuracySD, max=MeanDecreaseAccuracy+AccuracySD), width=.2, position=position_dodge(width=0.7)) +
  ylab("Variable importance"); p2 # + 
ggsave(filename = "variable_importance_top20_ASVs_genus.pdf",plot = p)
ggsave(filename = "variable_importance_top20_chem.pdf",plot = p2)
ggsave(filename = "variable_importance_top20_ASVs_genus.png",plot = p)
ggsave(filename = "variable_importance_top20_chem.png",plot = p2)
ggsave(filename = "variable_importance_top20_ASVs_genus.svg",plot = p)
ggsave(filename = "variable_importance_top20_ASVs_genus.tiff",plot = p,width = 8,height = 4)
# Summary statistics
# Chemistry
data.CST %>%
  group_by(Env) %>%
  summarize(
    count = n(),
    mean.pH = mean(pH, na.rm = TRUE),
    sd.pH = sd(pH, na.rm = TRUE),
    mean.ALK = mean(Alkalinity, na.rm = TRUE),
    sd.ALK = sd(Alkalinity, na.rm = TRUE),
    mean.AL = mean(Aluminium, na.rm = TRUE),
    sd.AL = sd(Aluminium, na.rm = TRUE)
  )
# Microbiology
# 1
data.MST %>%
  group_by(Env) %>%
  summarize(
    count = n(),
    mean = mean(asv_31, na.rm = TRUE)/sum(asv_31),
    sd = sd(asv_31/sum(asv_31), na.rm = TRUE)
  )
# 2
data.MST %>%
  group_by(Env) %>%
  summarize(
    count = n(),
    mean = mean(asv_4, na.rm = TRUE)/sum(asv_4),
    sd = sd(asv_4/sum(asv_4), na.rm = TRUE)
  )
# 4
data.MST %>%
  group_by(Env) %>%
  summarize(
    count = n(),
    mean = mean(asv_39, na.rm = TRUE)/sum(asv_39),
    sd = sd(asv_39/sum(asv_39), na.rm = TRUE)
  )
########### barplot ################
# add lake as factor
meta <- metadata %>%
column_to_rownames("SampleID")
sensitivity.df.new$LakeSystem <- meta[as.character(sensitivity.df.new$SampleID),"Lake"]  # which(colnames(meta) == "Lake")
m1.long <- melt(sensitivity.df.new,
                id.vars = c("SampleID","mtry","num.trees", "min.node.size","path.data.MST","LakeSystem"))
m1.long <- m1.long %>%
  rename(Environment = variable,Proportion = value)
# Stacked barplot with multiple groups
bp1 <- ggplot(data=m1.long, aes(x=SampleID, y=Proportion, fill=Environment)) +
  geom_bar(stat="identity",color = "black") +
  facet_grid(cols = vars(LakeSystem),scales = "free") +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 12),axis.title.x = element_blank()) +
  scale_fill_manual(values=c( "#56B4E9", "#009E73")); bp1
ggsave(filename = "barplot_ASV_lakes_prediction.pdf",bp1)
  #facet_wrap(Lake ~ .,scale = "free) +
sensitivity.df.new2$LakeSystem <- meta[as.character(sensitivity.df.new2$samplenames),"Lake"]  # which(colnames(meta) == "Lake")
m2.long <- melt(sensitivity.df.new2,
                id.vars = c("samplenames","mtry","num.trees", "min.node.size","path.data.CST","LakeSystem"))
m2.long <- m2.long %>%
  rename(Environment = variable,Proportion = value)
# Stacked barplot with multiple groups
bp2 <- ggplot(data=m2.long, aes(x=samplenames, y=Proportion, fill=Environment)) +
  geom_bar(stat="identity",color = "black") +
  facet_grid(cols = vars(LakeSystem),scales = "free") +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 12),axis.title.x = element_blank()) +
  scale_fill_manual(values=c( "#56B4E9", "#009E73")); bp2
ggsave(filename = "barplot_ASV_lakes_prediction.pdf",bp2)
# joint barplot
sensitivity.df.new3 <- sensitivity.df.new2 %>% 
  rename(SampleID = samplenames,Physiochem = Inflow) %>%
  select(-c(Lake,path.data.CST,min.node.size,num.trees,mtry)) #%>%
  #mutate(DataType = "Physiochem")
sensitivity.df.new4 <- sensitivity.df.new %>% 
  select(-c(Lake,path.data.MST,min.node.size,num.trees,mtry )) %>%
  rename(V3V4_16SrRNA = Inflow)
  #mutate(DataType = "16rRNA_V3V4")
data_ST <- merge(sensitivity.df.new4,sensitivity.df.new3,by=c("SampleID","LakeSystem"))
data_ST.long <- data_ST %>% 
  melt(.) %>%
  rename(DataType = variable, Proportion = value,Sample = SampleID)
bpd <- ggplot(data = data_ST.long,aes(x = Sample,y = Proportion,fill = DataType)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 12,color = "black"),legend.text = element_text(size = 12,color = "black"),strip.text = element_text(size = 12),legend.position = "top") +
  #theme_minimal() + 
  scale_fill_manual(values=c('#999999','#E69F00')) + # Use custom colors
  facet_wrap(vars(LakeSystem), scales = "free"); bpd 
ggsave(filename = "joint_barplot_ST_RF.pdf",plot = bpd)
ggsave(filename = "joint_barplot_ST_RF.tiff",plot = bpd,scale = 2.5)
# Heatmap

bph <- ggplot(data_ST.long, aes(Sample, DataType)) +
  geom_tile(aes(fill = Proportion), colour = "white") +
  scale_fill_gradient(low = "darkblue", high = "red") + 
  theme(axis.text.x = element_blank(),strip.text = element_text(size = 10),axis.text.y = element_text(size = 8,color = "black"),legend.text = element_text(size = 10,color = "black")) + 
  facet_wrap(nrow = 5,ncol = 1,vars(LakeSystem), scales = "free"); bph
ggsave(filename = "joint_heatmap_ST_RF.pdf",plot = bph)
##### Joint figure ######

plot_grid(p, p2, bpd,  ncol=2)
