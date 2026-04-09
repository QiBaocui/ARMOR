source("model.R")
library(dplyr)
library(data.table)
library(ggplot2)
library(optparse)
library(caret)
library(h2o)

option_list <- list(
  make_option(c("-o", "--outdir"), type = "character", default = NULL),
  make_option(c("-g", "--group"), type = "character", default = NULL),
  make_option(c("-G", "--gene"), type = "character", default = NULL),
  make_option(c("-s", "--snp"), type = "character", default = NULL),
  make_option(c("-k", "--kmer"), type = "character", default = NULL),
  make_option(c("-t", "--trainNum"), type = "numeric", default = 0.6),
  make_option(c("-S", "--split_type"), type = "character", default = "proportional"),
  make_option(c("-p", "--ip"), type = "numeric", default = 54321)
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

# group <- read.table("group_cefoxitin_matrix.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
# gene_dat <- read.table("gene_matrix.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
# snp_dat <- read.table("SNP_matrix.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
# kmer_dat <- read.table("kmer_matrix.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
# trainRadio=0.6
# outdir <- "/gpfs01/home/qibc/project/h2o/test/cefoxitin"
# IPaddress <- 54321

group <- read.table(opt$group,sep = "\t",header = T,stringsAsFactors = F,check.names = F)
gene_dat <- fread(opt$gene,sep = "\t",header = T) %>% as.data.frame()
snp_dat <- fread(opt$snp,sep = "\t",header = T) %>% as.data.frame()
kmer_dat <- fread(opt$kmer,sep = "\t",header = T) %>% as.data.frame()
trainNum <- opt$trainNum
outdir <- opt$outdir
IPaddress <- opt$ip
split_type <- opt$split_type
if(!dir.exists(outdir)){
  dir.create(outdir)
}

group_split <- splitSample(group = group,split_type = split_type,ratio=trainNum,num=trainNum)
train_group <- group_split$train_group
test_group <- group_split$test_group

gene_train <- left_join(train_group,gene_dat,by="ID")
gene_test <-left_join(test_group,gene_dat,by="ID")
snp_train <- left_join(train_group,snp_dat,by="ID")
snp_test <-left_join(test_group,snp_dat,by="ID")
kmer_train <- left_join(train_group,kmer_dat,by="ID")
kmer_test <-left_join(test_group,kmer_dat,by="ID")
train <- left_join(gene_train,snp_train[,-2],by="ID") %>% left_join(.,kmer_train[,-2],by="ID")
test <- left_join(gene_test,snp_test[,-2],by="ID") %>% left_join(.,kmer_test[,-2],by="ID")

rownames(gene_train) <- gene_train[,1]
gene_train <- gene_train[,-1]
rownames(gene_test) <- gene_test[,1]
gene_test <- gene_test[,-1]
rownames(snp_train) <- snp_train[,1]
snp_train <- snp_train[,-1]
rownames(snp_test) <- snp_test[,1]
snp_test <- snp_test[,-1]
rownames(kmer_train) <- kmer_train[,1]
kmer_train <- kmer_train[,-1]
rownames(kmer_test) <- kmer_test[,1]
kmer_test <- kmer_test[,-1]
rownames(train) <- train[,1]
train <- train[,-1]
rownames(test) <- test[,1]
test <- test[,-1]

dat <- list(train_group=train_group,test_group=test_group,train=train,test=test,gene_train=gene_train,gene_test=gene_test,snp_train=snp_train,snp_test=snp_test,kmer_train=kmer_train,kmer_test=kmer_test)
save(dat,file = paste0(outdir,"/train_data.Rdata"))

h2o.init(ip = "localhost", port = IPaddress)
modelFun(train.R=gene_train,test.R=gene_test,outdir=outdir,feature="Gene",feature_type_list=c("GBM","GLM","RF","DL"))
modelFun(train.R=snp_train,test.R=snp_test,outdir=outdir,feature="SNP",feature_type_list=c("GBM","GLM","RF","DL"))
modelFun(train.R=kmer_train,test.R=kmer_test,outdir=outdir,feature="Kmer",feature_type_list=c("GBM","GLM","RF","DL"))
modelFun(train.R=train,test.R=test,outdir=outdir,feature="All",feature_type_list=c("GBM","GLM","RF","DL"))
#h2o.shutdown()
h2o.shutdown(prompt = F)

#loaded_model <- h2o.loadModel(model_path)

