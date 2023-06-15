###
######      Prognosis Signature Select     ######      
##  Supporter: Yongwei_Zhu Wei Zhang
##  E-mail: zhuyongwei@csu.edu.cn  2204170416@csu.edu.cn
##  Vision: V1.0(2022-11-01)


### 策略： UniCox + Lasso 
## 准备两个必需文件
# gene_list：基因集
# inputSet：数据集的表达谱矩阵(行：基因ID；列：样本ID) +预后信息；
################第一列为样本名称为ID， 第二列为OS.time（day）， 第三列为OS,后面的列为基因名。不允许空值，OS为1/0。


# # 参数设置示例
# #unicox
# unicox_pcutoff <- 0.05 #单因素回归的筛选阈值
# #LASSO
# iter.times <- 100 #设置LASSO迭代次数，一般为1000
# nfolds = 10 # 10-fold交叉验证选取最优lambda
# alpha = 1 # alpha = 1 意味着 lasso
# family = "cox" # 依赖cox模型

SigSelectUL <- function(gene_list,
                        inputSet,
                        unicox_pcutoff, #单因素回归的筛选阈值0.05
                        iter.times, #设置LASSO迭代次数，一般为1000
                        nfolds,  # 10-fold交叉验证选取最优lambda, 
                        alpha,   # alpha = 1 意味着 lasso
                        family, # 依赖cox模型
                        pred.time){
  ##############加载相关包##########
  if(!require(multtest))  BiocManager::install("clusterProfiler")#ID转化
  if(!require(multtest))install.packages("survival")#生存分析
  if(!require(multtest))install.packages("glmnet")#LASSO回归
  if(!require(multtest))install.packages("pbapply")#批量处理LASSO回归
  
  if(!require(pbapply)){ 
    install.packages("pbapply")
  } else {library(pbapply)}  #批量处理LASSO回归
  
  
  if(!require(tidyverse)){ 
    install.packages("tidyverse")
  } else {library(tidyverse)} 
  
  if(!require(dplyr)){ 
    install.packages("dplyr")
  } else {library(dplyr)} 
  
  if(!require(clusterProfiler)){ 
    BiocManager::install("clusterProfiler")
  } else {library(clusterProfiler)} ##ID转化
  
  if(!require(survival)){ 
    install.packages("survival")
  } else {library(survival)} #生存分析
  
  if(!require(glmnet)){ 
    install.packages("glmnet")
  } else {library(glmnet)} #LASSO回归
  
  if(!require(org.Hs.eg.db)){ 
    BiocManager::install("org.Hs.eg.db")
  } else {library(org.Hs.eg.db)} #计算AUC
  
  
  
  print("Starting the data preprocess")
  ###############数据预处理#######
  inputSet <- inputSet %>% as.data.frame()
  inputSet$OS.time <- as.numeric(inputSet$OS.time)#时间数值化
  inputSet <- inputSet[inputSet$OS.time>0,]#剔除时间为0
  
  print("Rejecting a null value")
  
  #剔除空值的基因
  table(is.na(inputSet))
  inputSet <- t(inputSet) %>% na.omit()
  inputSet <- t(inputSet) %>% as.data.frame()
  table(is.na(inputSet))
  
  print("Correcting gene set")
  #将基因集矫正 策略为转化从ALIAS为ENTREZID， 然后再转化为SYMBOL,避免别名
  gene.df <- bitr(gene_list, fromType = "ALIAS", 
                  toType = c("ENTREZID"), 
                  OrgDb = org.Hs.eg.db)
  gene.df1 <- bitr(gene.df$ENTREZID, fromType = "ENTREZID", 
                   toType = c("SYMBOL"), 
                   OrgDb = org.Hs.eg.db)
  gene_list <- gene.df1$SYMBOL
  write.table(gene_list, "ID_transformed_genelist.txt",row.names = F, quote = F)
  
  #将genelist和表达矩阵的基因名称格式统一
  gene_list <- gsub("-","_",gene_list)
  colnames(inputSet)[4:ncol(inputSet)] <- gsub("-","_",colnames(inputSet)[4:ncol(inputSet)])
  
  print("Gets the intersection of genelist and expression profile")
  #获取genelist和表达谱的交集
  comsa1 <- intersect(colnames(inputSet)[4:ncol(inputSet)],gene_list)
  write.table(comsa1,"intersection_genelist_exprSet_gene.txt", row.names = F, quote = F)
  
  print("Processing the  input representation matrix")
  #对输入的表达矩阵进行处理
  inputSet <- inputSet[,c("ID","OS.time","OS",comsa1)]
  
  inputSet[,c(1:2)] <- apply(inputSet[,c(1:2)],2,as.factor)
  inputSet[,c(2:ncol(inputSet))] <- apply(inputSet[,c(2:ncol(inputSet))],2,as.numeric)
  inputSet <- as.data.frame(inputSet)
  rownames(inputSet) <- inputSet$ID
  
  print("Data preprocessing completed")
  #自定义显示进程函数
  display.progress = function (index, totalN, breakN=20) {
    if ( index %% ceiling(totalN/breakN)  ==0  ) {
      cat(paste(round(index*100/totalN), "% ", sep=""))
    }
  }  
  
  
  ###############单变量cox#######
  print("Stating the univariable cox regression")
  
  unicox <- data.frame()
  for(i in 1:ncol(inputSet[,4:ncol(inputSet)])){
    
    display.progress(index = i, totalN = ncol(inputSet[,4:ncol(inputSet)]))
    gene <- colnames(inputSet[,4:ncol(inputSet)])[i]
    tmp <- data.frame(expr = as.numeric(inputSet[,4:ncol(inputSet)][,i]),
                      futime = inputSet$OS.time,
                      fustat = inputSet$OS,
                      stringsAsFactors = F)
    cox <- coxph(Surv(futime, fustat) ~ expr, data = tmp)
    coxSummary <- summary(cox)
    unicox <- rbind.data.frame(unicox,
                               data.frame(gene = gene,
                                          HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                          z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                          pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                          lower = as.numeric(coxSummary$conf.int[,3][1]),
                                          upper = as.numeric(coxSummary$conf.int[,4][1]),
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
  }
  
  write.csv(unicox,"unicox_results.csv")
  
  print("Finished the univariable cox regression")
  
  
  #进行变量筛选
  selgene <- unicox[which(unicox$pvalue < unicox_pcutoff), "gene"]
  write.table(selgene,paste("unicox_selected_cutoff_",unicox_pcutoff,"_genes.txt",sep = ""))
  
  
  print("Starting the lasso regression 1000 times")
  # 运行1000次 lasso penalty
  
  candidate.gene <- selgene
  surv.obj <- Surv(inputSet$OS.time, inputSet$OS)
  lasso_fea_list <- list()
  
  list.of.seed <- 1:iter.times
  lasso_fea_list <- pblapply(list.of.seed, function(x){ # 大概运行2天
    set.seed(list.of.seed[x])
    cvfit = cv.glmnet(x = as.matrix(inputSet[,candidate.gene]), 
                      y = surv.obj, 
                      nfolds = nfolds, # 10-fold交叉验证选取最优lambda
                      alpha = alpha, # alpha = 1 意味着 lasso
                      family = family, # 依赖cox模型
                      maxit = 1000) 
    
    # 取出最优lambda
    fea <- rownames(coef(cvfit, s = 'lambda.min'))[coef(cvfit, s = 'lambda.min')[,1]!= 0]
    if(is.element("(Intercept)", fea)) {
      lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
    } else {
      lasso_fea <- sort(fea)
    }
    return(lasso_fea)
  })
  
  save(lasso_fea_list,file = "lasso_fea_list.rda") # 保存该结果
  
  print("Finish running 1000 lasso regressions")
  
  print("Output the set of genes for each run")
  
  # 输出每次运行的基因集合
  lasso_res <- NULL
  for(i in 1:iter.times) {
    lasso_res <- rbind.data.frame(lasso_res,
                                  data.frame(iteration = i,
                                             n.gene = length(lasso_fea_list[[i]]),
                                             genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
                                             stringsAsFactors = F),
                                  stringsAsFactors = F)
  }
  
  print("Strategy 1: Determine the model according to the occurrence times of the obtained gene model")
  
  uniquelist <- unique(lasso_res$genelist) #注意这里就已经按照频数从高到低排序了
  uniquelab <- LETTERS[1:length(uniquelist)]
  lasso_res$uniquelab <- NA
  for (i in 1:length(uniquelist)) {
    lasso_res[which(lasso_res$genelist == uniquelist[i]),"uniquelab"] <- uniquelab[i]
  }
  lasso_res$label <- paste(lasso_res$n.gene,"genes",lasso_res$uniquelab,sep = "_") # 最终模型标签
  
  #将最终得到的模型标签输出
  write.csv(lasso_res,"model_label.csv")
  
  genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # 根据基因出现的频次排序
  # 如果觉得出现次数较少的基因是不鲁棒的，也可以仅选择top基因或者出现次数大于5%的，例如运行一千次如果出现次数小于五十
  #那么被认为是不重要的
  freq.cutoff <- 0.05*iter.times
  genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50=1000*0.05的基因，认为是多次lasso的共识基因
  write.table(genes, "top50gene.txt", row.names = F, quote = F) #输出到文件
  
}
