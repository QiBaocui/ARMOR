
splitSample <- function(group,split_type="proportional",ratio=0.6,num=100){
  colnames(group) <- c("ID","Group")
  if(split_type == "proportional"){
    library(caret)
    train_index <- createDataPartition(group[,2], p=ratio, list=F)
    train_group <- group[train_index,]
    test_group <- group[-train_index,]
  }else if(split_type == "sample_size"){
    per_group_num <- floor(num / 2)  # 两组各取该数量，保证总数≤num且均等
    group_levels <- unique(group$Group)
    if(length(group_levels) != 2){
      stop("split_type='sample_size' 仅支持2个分组的样本拆分！当前分组数：", length(group_levels))
    }
    train_list <- list()
    for(i in 1:length(group_levels)){
      # 获取当前分组的所有样本
      sub_group <- group[group$Group == group_levels[i], ]
      # 校验：当前分组的样本量是否足够抽取per_group_num个
      if(nrow(sub_group) < per_group_num){
        stop("分组 ", group_levels[i], " 样本量不足（仅", nrow(sub_group), "个），无法抽取", per_group_num, "个训练样本！")
      }
      # 随机抽取per_group_num个样本（无放回抽样）
      sub_train_index <- sample(1:nrow(sub_group), size = per_group_num, replace = FALSE)
      train_list[[i]] <- sub_group[sub_train_index, ]
    }
    
    # 4. 合并两个分组的训练样本，形成最终训练集
    train_group <- do.call(rbind, train_list)
    
    # 5. 提取测试集：不在训练集中的所有样本（通过ID匹配去重）
    test_group <- group[!group$ID %in% train_group$ID, ]
    
    # 可选：重置行名（避免行名混乱）
    rownames(train_group) <- NULL
    rownames(test_group) <- NULL
  }else{
    print("Sample Split Fail")
    break
  }
  df <- list(train_group=train_group,test_group=test_group)
  return(df)
}


get_outputperformance <- function(perf){
  outputperformance <- perf  %>% .@metrics %>% .$thresholds_and_metric_scores
  outputperformance.R <- as.data.frame(outputperformance)
  youdenIndex <- outputperformance.R$tpr - outputperformance.R$fpr #计算cutoff，选择Sensitivity+Specificity最大的值,youden指数=Sensitivity+Specificity-1=TPR−FPR最大值
  cutoff <- max(youdenIndex)
  site <- which(youdenIndex==cutoff)
  cutoffpermance <- outputperformance.R[site,]
  return(cutoffpermance)
}
get_predict <- function(model,modelName,test.R,test.h2o,outdir,feature){
  pred_matri <- h2o.predict(object = model, newdata = test.h2o) %>% as.data.frame()
  pred_matri <- cbind(Group=test.R[,1],pred_matri)
  pred_matri <- cbind(sample_code=rownames(test.R),pred_matri)
  test_perf <- h2o.performance(model = model,newdata = test.h2o)
  cutoffpermance <- get_outputperformance(test_perf)
  cutoffpermance$auc <- h2o.auc(test_perf)
  cutoffpermance$model <- modelName
  rocData <- test_perf %>% .@metrics %>% .$thresholds_and_metric_scores %>% .[c('tpr','fpr')] %>% add_row(tpr=0,fpr=0,.before=T)
  rocData <- rocData[order(rocData$fpr),]
  rocData$label <- paste0(modelName,"_AUC: ",signif(h2o.auc(test_perf),digits = 3))
  rocCurve <- ggplot(rocData,aes(x=fpr,y=tpr,color=label))+
    geom_line(size=1.2)+
    theme_bw()+
    xlab("False positive rate")+ylab("True positive rate")+
    labs(color='') +
    theme(legend.position = c(.8,.2))+
    scale_color_manual(values=c("#4e639e"))+
    theme(axis.text.x = element_text(size = 10,  face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          axis.title.y = element_text(size = 15,face = "bold"),
          legend.title = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), #去掉背景线和网格线
          plot.title = element_text(hjust = 0.5, vjust=1), ##更改标题位置
          plot.margin =  unit(c(1,1,1,1), "cm"))+
    geom_abline(intercept=0,slope=1,colour="grey")+
    ggtitle('ROC Curve')
  
  if(modelName=="stacked"){
    var <- h2o.permutation_importance(model, newdata = test.h2o) %>% as.data.frame() #置换重要性通过打乱每个特征的数据并观察模型性能下降的程度来评估该特征的重要性，其评估结果通常比模型自带的变量重要性更具解释性。模型无关，结果更稳健，能更好地反映特征对模型预测的真实贡献，但计算成本较高
  }else{
    var <- h2o.varimp(model) %>% as.data.frame()
  }
  options(bitmapType = 'cairo')
  #ggsave(paste0(outdir,"/",feature,"_",modelName,"_ROC.png"),rocCurve,width = 800,height = 800,units = "px")
  png(paste0(outdir, "/", feature, "_", modelName, "_ROC.png"), width = 800, height = 800)
  print(rocCurve)  # 必须显式打印图形
  dev.off()
  ggsave(paste0(outdir,"/",feature,"_",modelName,"_ROC.pdf"),rocCurve,width = 8,height = 8)
  write.table(pred_matri,paste0(outdir,"/",feature,"_",modelName,"_pred.tsv"),sep="\t",quote=F,row.names = F,col.names = T,na = "")
  write.table(cutoffpermance,paste0(outdir,"/",feature,"_",modelName,"_cutoffpermance.tsv"),sep="\t",quote=F,row.names = F,col.names = T,na = "")
  write.table(rocData,paste0(outdir,"/",feature,"_",modelName,"_rocData.tsv"),sep="\t",quote=F,row.names = F,col.names = T,na = "")
  write.table(var,paste0(outdir,"/",feature,"_",modelName,"_varimp.tsv"),sep="\t",quote=F,row.names = F,col.names = T,na = "")
}


modelFun <- function(train.R,test.R,outdir,feature,feature_type_list,nfolds=10){
  #h2o.init()
  train.h2o <- as.h2o(train.R)
  train.h2o$Group <- h2o.asfactor(train.h2o$Group)
  test.h2o <- as.h2o(test.R)
  test.h2o$Group <- h2o.asfactor(test.h2o$Group)
  
  y <- "Group"
  x <- setdiff(names(train.R),y)
  model_list <- list()
  if("GBM" %in% feature_type_list){
    model_list[["GBM"]] <- h2o.gbm(x = x,y = y,training_frame = train.h2o,nfolds = nfolds,keep_cross_validation_predictions = TRUE,fold_assignment = "Stratified",keep_cross_validation_fold_assignment = TRUE,seed = 1)
    get_predict(model=model_list[["GBM"]],modelName="GBM",test.R=test.R,test.h2o=test.h2o,outdir=outdir,feature=feature)
  }
  if("GLM" %in% feature_type_list){
    model_list[["GLM"]] <- h2o.glm(x = x,y = y,training_frame = train.h2o,nfolds = nfolds,keep_cross_validation_predictions = TRUE,fold_assignment = "Stratified",keep_cross_validation_fold_assignment = TRUE,seed = 1)
    get_predict(model=model_list[["GLM"]],modelName="GLM",test.R=test.R,test.h2o=test.h2o,outdir=outdir,feature=feature)
  }
  if("RF" %in% feature_type_list){
    model_list[["RF"]] <- h2o.randomForest(x = x,y = y,training_frame = train.h2o,nfolds = nfolds,keep_cross_validation_predictions = TRUE,fold_assignment = "Stratified",seed = 1)
    get_predict(model=model_list[["RF"]],modelName="RF",test.R=test.R,test.h2o=test.h2o,outdir=outdir,feature=feature)
  }
  if("DL" %in% feature_type_list){
    model_list[["DL"]] <- h2o.deeplearning(x = x,y = y,training_frame = train.h2o,nfolds = nfolds,keep_cross_validation_predictions = TRUE,fold_assignment = "Stratified",seed = 1)
    get_predict(model=model_list[["DL"]],modelName="DL",test.R=test.R,test.h2o=test.h2o,outdir=outdir,feature=feature)
  }
  if("xGboost" %in% feature_type_list){
    model_list[["xGboost"]] <- h2o.xgboost(x = x,y = y,training_frame = train.h2o,nfolds = nfolds,keep_cross_validation_predictions = TRUE,fold_assignment = "Stratified",seed = 1)
    get_predict(model=model_list[["xGboost"]],modelName="xGboost",test.R=test.R,test.h2o=test.h2o,outdir=outdir,feature=feature)
  }
  if("Bayes" %in% feature_type_list){
    model_list[["Bayes"]] <- h2o.naiveBayes(x = x,y = y,training_frame = train.h2o,nfolds = nfolds,keep_cross_validation_predictions = TRUE,fold_assignment = "Stratified",seed = 1)
    get_predict(model=model_list[["Bayes"]],modelName="Bayes",test.R=test.R,test.h2o=test.h2o,outdir=outdir,feature=feature)
  }
  stacked <- h2o.stackedEnsemble(x = x,y = y,training_frame = train.h2o,base_models = model_list)
  get_predict(model=stacked,modelName="stacked",test.R=test.R,test.h2o=test.h2o,outdir=outdir,feature=feature)
  h2o.saveModel(object = stacked, path = outdir, filename = paste0(feature,"_stacked_model"), force = TRUE)
  for (i in feature_type_list) {
    h2o.saveModel(object = model_list[[i]], path = outdir, filename = paste0(feature,"_",i,"_model"), force = TRUE)
  }
}
