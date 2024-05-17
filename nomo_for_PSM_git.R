setwd("location")
library(readxl)
library(reportROC)
library(pROC)
library(plyr)
library(gmodels)
library(generics)
library(glmnet)
library(MASS)
library(tableone)
library(tidyverse)
library(table1)
library(rms)
library(DescTools)
library(ResourceSelection)
library(boot)

data1 <- data.frame(read_xlsx("file.xlsx", sheet=1))

factor_vars <- c("bGG4", "PIRADS", "PIRADS3", "fEPE", "cT_summary", "pSM", "m_psm", "robot")
data1[factor_vars] <- lapply(data1[factor_vars], as.factor)

#Shapiro-Wilk test
variable_names <- c("age", "tPSA_ng.ml", "PSAD", "posi.core.percent", "prostate_volume", "MR_size", "CCL", "UDL", "depth", "annual_surg_vol")
shapiro_results <- lapply(variable_names, function(var) shapiro.test(data1[[var]]))
names(shapiro_results) <- variable_names
shapiro_results

#demographic characteristics
calculate_statistics <- function(data) {
  list(
    summary = summary(data),
    sd = if (is.numeric(data)) sd(data) else NA,
    count = if (is.factor(data) ) table(data) else NA,
    proportion = if (is.factor(data)) prop.table(table(data)) * 100 else NA
  )
}

variables <- c("age", "tPSA_ng.ml", "PSAD", "bGG4", "posi.core.percent", "prostate_volume", 
               "PIRADS", "MR_size", "CCL", "fEPE", "UDL", "depth", "cT_summary", 
               "annual_surg_vol", "pSM", "m_psm", "robot")

statistics <- lapply(variables, function(var) calculate_statistics(data1[[var]]))
names(statistics) <- variables
statistics

#######################################################
#univariate analysis for oPSM
fit_and_show <- function(formula, data) {
  lm <- glm(formula, family = "binomial", data = data)
  ShowRegTable(lm, exp = TRUE, digits = 3, pDigits = 3, printToggle = TRUE, quote = FALSE, ciFun = confint)
}

formulas <- list(
  pSM ~ age,
  pSM ~ tPSA_ng.ml,
  pSM ~ PSAD,
  pSM ~ bGG4,
  pSM ~ posi.core.percent,
  pSM ~ prostate_volume,
  pSM ~ PIRADS3,
  pSM ~ MR_size,
  pSM ~ fEPE,
  pSM ~ CCL,
  pSM ~ UDL,
  pSM ~ depth,
  pSM ~ annual_surg_vol,
  pSM ~ cT_summary,
  pSM ~ robot
)

results <- lapply(formulas, fit_and_show, data = data1)

#lasso
y<-as.matrix(data1[,"pSM"])
data1[factor_vars] <- lapply(data1[factor_vars], as.numeric)
x<-as.matrix(data1[,c("age","tPSA_ng.ml","PSAD","bGG4","posi.core.percent","prostate_volume", "PIRADS3","MR_size","CCL","fEPE",
                      "UDL","depth","cT_summary", "annual_surg_vol","robot")])

#select variables
set.seed(123)
cvfit <- cv.glmnet(x, y, alpha = 1, family = "binomial")
lambda_1se <- cvfit$lambda.1se
#coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact = F) #min
lasso_model_psm <- glmnet(x, y, alpha = 1, lambda = lambda_1se, family = "binomial")
coef_lasso_psm<- coef(lasso_model_psm)
selected_variables <- rownames(coef_lasso_psm)[which(coef_lasso_psm != 0)]
selected_variables <- selected_variables[selected_variables != "(Intercept)"]
print(selected_variables)

#build the model
formula_psm <- as.formula(paste("pSM ~", paste(selected_variables, collapse = " + ")))
final_model_psm <- glm(formula_psm, data = data1, family = "binomial")
table1<-ShowRegTable(final_model_psm, exp = TRUE, digits= 3, pDigits = 3, printToggle=TRUE, quote=FALSE,ciFun=confint)

#roc
data1$prepsm <-predict(final_model_psm, newdata = data1)
reportROC(gold=data1$pSM,data1$prepsm,important="se", plot=FALSE)
roc1 <- roc(data1$pSM, data1$prepsm)
plot(roc1, col="darkred", identity.col="grey",identity.lty=2, identity.lwd=1)
legend("bottomright", legend="AUC for PSM: 0.717 (0.648–0.785)   ",
       col="darkred",
       lty=1)

#nomo
data.psm = data1 %>% select(pSM=pSM, Maximum_tumor_size = MR_size, Frank_EPE = fEPE, Annual_surgery_volume = annual_surg_vol)
dd<-datadist(data.psm)
options(datadist='dd')
lmpsmnomo<-lrm(pSM~Maximum_tumor_size+Frank_EPE+Annual_surgery_volume, data.psm,x=T,y=T)
print(lmpsmnomo,digits=3)
coef(lmpsmnomo)#b值
exp(coef(lmpsmnomo))#OR
#summary(lmpsmnomo, Maximum_tumor_size=c(0,1), Frank_EPE=c(0,1),Annual_surgery_volume=c(0,1))
nompsm <- nomogram(lmpsmnomo, fun=plogis,
                   lp=F, funlabel="Risk of PSM")
plot(nompsm)

#bootstrap
f<-"pSM~Maximum_tumor_size+Frank_EPE+Annual_surgery_volume"
c<-nomogram_Bootstrap_validation(f=f,dataset=data.psm,Iterations=1000,seed=123)
cat(c$print_report)

#calibration
calpsm<-calibrate(lmpsmnomo,  method = "boot", B = 1000)
plot(calpsm,lwd=2,lty=1,  
     xlab="Predicted risk for PSM",  
     ylab="Observed risk",  
     xlim = c(0,1),ylim= c(0,1),
     col=c(rgb(192,98,83,maxColorValue=255)))
#Hosmer-Lemeshow test
f.glm <- glm(pSM~.,data=data.psm,family = binomial(link = "logit"))
P1 <- predict(f.glm,type = 'response')
val.prob(P1,data1$pSM)

#######################################################
#univariate analysis for mPSM
formulas2 <- list(
  m_psm ~ age,
  m_psm ~ tPSA_ng.ml,
  m_psm ~ PSAD,
  m_psm ~ bGG4,
  m_psm ~ posi.core.percent,
  m_psm ~ prostate_volume,
  m_psm ~ PIRADS3,
  m_psm ~ MR_size,
  m_psm ~ fEPE,
  m_psm ~ CCL,
  m_psm ~ UDL,
  m_psm ~ depth,
  m_psm ~ annual_surg_vol,
  m_psm ~ cT_summary,
  m_psm ~ robot
)

results <- lapply(formulas2, fit_and_show, data = data1)

#lasso
y2<-as.matrix(data1[,"m_psm"])
data1[factor_vars] <- lapply(data1[factor_vars], as.numeric)
x<-as.matrix(data1[,c("age","tPSA_ng.ml","PSAD","bGG4","posi.core.percent","prostate_volume", "PIRADS3","MR_size","CCL","fEPE",
                      "UDL","depth","cT_summary", "annual_surg_vol","robot")])

#select variables
set.seed(1234)
cvfit2 <- cv.glmnet(x, y2, alpha = 1, family = "binomial")
lasso_model_mpsm <- glmnet(x, y2, alpha = 1, lambda = cvfit2$lambda.1se, family = "binomial")
coef_lasso_mpsm<- coef(lasso_model_mpsm)
coef(cvfit2$glmnet.fit,s=cvfit2$lambda.1se,exact = F) 
selected_variables <- rownames(coef_lasso_mpsm)[which(coef_lasso_mpsm != 0)]
selected_variables <- selected_variables[selected_variables != "(Intercept)"]
print(selected_variables)

#build the model
formula_mpsm <- as.formula(paste("m_psm ~", paste(selected_variables, collapse = " + ")))
final_model_mpsm <- glm(formula_mpsm, data = data1, family = "binomial")
table1<-ShowRegTable(final_model_mpsm, exp = TRUE, digits= 3, pDigits = 3, printToggle=TRUE, quote=FALSE,ciFun=confint)

#roc
data1$prempsm <-predict(final_model_mpsm, newdata = data1)
reportROC(gold=data1$m_psm,data1$prempsm,important="se", plot=FALSE)
roc2 <- roc(data1$m_psm, data1$prempsm)
plot(roc2, col="darkred", identity.col="grey",identity.lty=2, identity.lwd=1)
legend("bottomright", legend="AUC for mPSM: 0.790 (0.719–0.861)   ",
       col="darkred",
       lty=1)

#nomo
data.mpsm = data1 %>% select(mpsm=m_psm, Positive_core_percentage=posi.core.percent, Maximum_tumor_size = MR_size, Apex_depth = depth, Annual_surgery_volume = annual_surg_vol)
dd2<-datadist(data.mpsm)
options(datadist='dd2')
lmmpsmnomo<-lrm(mpsm~Positive_core_percentage+Maximum_tumor_size+Apex_depth+Annual_surgery_volume, data.mpsm,x=T,y=T)
lmmpsmnomo
coef(lmmpsmnomo)#b值
exp(coef(lmmpsmnomo))#OR
summary(lmmpsmnomo, Positive_core_percentage=c(0,1), Maximum_tumor_size=c(0,1), Apex_depth=c(0,1),Annual_surgery_volume=c(0,1))
nommpsm <- nomogram(lmmpsmnomo, fun=plogis,
                    lp=F, funlabel="Risk of multifocal PSM")
plot(nommpsm)

#bootstrap
f<-"mpsm~Positive_core_percentage+Maximum_tumor_size+Apex_depth+Annual_surgery_volume"
c<-nomogram_Bootstrap_validation(f=f,dataset=data.mpsm,Iterations=1000,seed=NULL)
cat(c$print_report)

#calibration
calmpsm<-calibrate(lmmpsmnomo,  method = "boot", B = 1000)
plot(calmpsm,lwd=2,lty=1,  
     xlab="Predicted risk for PSM",  
     ylab="Observed risk",  
     xlim = c(0,1),ylim= c(0,1),
     col=c(rgb(192,98,83,maxColorValue=255)))
#Hosmer-Lemeshow test
f.glm2 <- glm(mpsm~.,data=data.mpsm,family = binomial(link = "logit"))
P2 <- predict(f.glm2,type = 'response')  ##获得预测概率值
val.prob(P2,data.mpsm$mpsm)


##########################################################
#Bootstrap_validation
library(car)
library(survival)
library(pROC)
library(rms)
library(tcltk)
#下方为内部验证相关预设函数,请直接运行加载
{
  nomogram_cross_validation<-function(f,dataset,N_fold=10,Iterations=200,seed=NULL){
    
    #设置随机数种子
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    #提取公式中的因变量和自变量
    f<-as.formula(f)
    y_list<-as.character(f[[2]])
    x_list<-as.character(unlist(strsplit(as.character(f[3]),split=" \\+ ")))
    
    #检查列名对不对
    if(!all(c(y_list,x_list) %in% colnames(dataset))){
      output_result<-list()
      print("方程中有变量不在数据列名中！请检查后再运行！")
      return(output_result)
    }
    
    #筛选出式子中包含的列，并筛选数据完整的行
    clean_dataset<-dataset[,(c(y_list,x_list))]
    clean_dataset<-clean_dataset[(complete.cases(clean_dataset)),]
    #生成空的输出对象
    output_result<-list()
    output_result$f<-f
    output_result$y_list<-y_list
    output_result$x_list<-x_list
    output_result$N_fold<-N_fold
    output_result$Iterations<-Iterations
    output_result$clean_dataset<-clean_dataset
    
    output_result$AUCs<-data.frame(matrix(NA,output_result$Iterations,3),stringsAsFactors = FALSE)
    colnames(output_result$AUCs)<-c("Iterations","AUC","C_index")
    output_result$AUCs$Iterations<-1:output_result$Iterations
    
    #计算需要随机分组的次数
    round<-ceiling(Iterations/N_fold)
    
    #生成初始化的迭代计数器
    inter_counter=0
    
    #进度条
    pb <- tkProgressBar("进度","已完成 %", 0, 100)
    
    #开始计算
    if(output_result$N_fold>nrow(output_result$clean_dataset)){
      print("交叉验证折数大于数据中变量完整的行数！")
    }else{
      for(index_round in 1:round){
        #随机分组
        rand_seq<-rep(1:output_result$N_fold,ceiling(nrow(output_result$clean_dataset)/output_result$N_fold))[1:nrow(output_result$clean_dataset)]
        rand_seq<-sample(rand_seq, nrow(output_result$clean_dataset), replace = FALSE)
        #计算该round内各折的结果
        for(index_fold in 1:output_result$N_fold){
          #检验迭代次数是否足够
          if(inter_counter==output_result$Iterations){
            break
          }else{
            
            #迭代次数不够继续计算
            #划分训练集和验证集
            temp_training_dataset<-output_result$clean_dataset[(rand_seq!=index_fold),]
            temp_validation_dataset<-output_result$clean_dataset[(rand_seq==index_fold),]
            ddist <<- datadist(temp_training_dataset)
            options(datadist='ddist')
            #建模
            if(length(table(temp_validation_dataset[,output_result$y_list]))==1){
              #验证集的结局事情全为0或全为1时无法计算AUC，默认为NA
            }else{
              f_reg <- tryCatch(lrm(output_result$f, data=temp_training_dataset, x=TRUE, y=TRUE,maxit=1000),error=function(e){return(list(fail=TRUE))} )
              if(f_reg$fail | 'try-error' %in% class(f_reg)){           # 判断当前循环的try语句中的表达式是否运行正确
                # 若该回归分析无法进行则该轮循环不赋值，默认为NA
              }else{
                #验证集中计算AUC
                pred_f_validation<-predict(f_reg,temp_validation_dataset)
                modelroc <- roc(temp_validation_dataset[,output_result$y_list],pred_f_validation,quiet=TRUE)
                output_result$AUCs[(inter_counter+1),"AUC"]<-modelroc$auc[1]
                #验证集中计算C-index
                temp_pred_list<-data.frame(grop=temp_validation_dataset[,output_result$y_list],pred=pred_f_validation)
                #所有样本两两配对
                valid_num=0
                concordance_num=0
                for(index_temp_pred_list_1 in 1:nrow(temp_pred_list)){
                  for(index_temp_pred_list_0 in 1:nrow(temp_pred_list)){
                    if((temp_pred_list$grop[index_temp_pred_list_1]==1)&(temp_pred_list$grop[index_temp_pred_list_0]==0)){
                      #当匹配的对子中第一个为group=1（实验组）;第二个为group=0（对照组）时为有效对子
                      valid_num<-valid_num+1
                      if(temp_pred_list$pred[index_temp_pred_list_1]>temp_pred_list$pred[index_temp_pred_list_0]){
                        concordance_num<-concordance_num+1
                      }
                    }
                  }
                }
                #计算C-index
                output_result$AUCs[(inter_counter+1),"C_index"]<-concordance_num/valid_num
              }
            }
            #进度条
            info <- sprintf("已完成 %d%%",round(inter_counter*100/output_result$Iterations))
            setTkProgressBar(pb, inter_counter*100/output_result$Iterations, sprintf("进度 (%s)", info), info)
            inter_counter<-inter_counter+1
          }
        }
      }
    }
    output_result$mean_AUC<-mean(output_result$AUCs$AUC,na.rm = TRUE)
    output_result$ci_left_AUC<-t.test(output_result$AUCs$AUC,na.rm = TRUE)$conf.int[1]
    output_result$ci_right_AUC<-t.test(output_result$AUCs$AUC,na.rm = TRUE)$conf.int[2]
    output_result$mean_C_index<-mean(output_result$AUCs$C_index,na.rm = TRUE)
    output_result$ci_left_C_index<-t.test(output_result$AUCs$C_index,na.rm = TRUE)$conf.int[1]
    output_result$ci_right_C_index<-t.test(output_result$AUCs$C_index,na.rm = TRUE)$conf.int[2]
    output_result$print_report<-paste0("本次验证方法: ",output_result$N_fold," 折交叉验证\n"
                                       ,"迭代次数: ",output_result$Iterations," 次\n"
                                       ,"其中迭代失败次数: ",sum(is.na(output_result$AUCs$AUC))," 次\n"
                                       ,"计算所得的平均AUC: ",output_result$mean_AUC,"\n"
                                       ,"计算所得的AUC的95%置信区间为: ",output_result$ci_left_AUC,", ",output_result$ci_right_AUC,"\n"
                                       ,"计算所得的平均C_index: ",output_result$mean_C_index,"\n"
                                       ,"计算所得的C_index的95%置信区间为: ",output_result$ci_left_C_index,", ",output_result$ci_right_C_index,"\n")
    return(output_result)
    
  }
  nomogram_Jackknife_validation<-function(f,dataset){
    
    
    #提取公式中的因变量和自变量
    f<-as.formula(f)
    y_list<-as.character(f[[2]])
    x_list<-as.character(unlist(strsplit(as.character(f[3]),split=" \\+ ")))
    #检查列名对不对
    if(!all(c(y_list,x_list) %in% colnames(dataset))){
      output_result<-list()
      print("方程中有变量不在数据列名中！请检查后再运行！")
      return(output_result)
    }
    #筛选出式子中包含的列，并筛选数据完整的行
    clean_dataset<-dataset[,(c(y_list,x_list))]
    clean_dataset<-clean_dataset[(complete.cases(clean_dataset)),]
    
    N_fold=nrow(clean_dataset)
    Iterations=nrow(clean_dataset)
    
    #生成空的输出对象
    output_result<-list()
    output_result$f<-f
    output_result$y_list<-y_list
    output_result$x_list<-x_list
    output_result$N_fold<-N_fold
    output_result$Iterations<-Iterations
    output_result$clean_dataset<-clean_dataset
    
    output_result$predicted_value<-data.frame(matrix(NA,output_result$Iterations,3),stringsAsFactors = FALSE)
    colnames(output_result$predicted_value)<-c("Iterations","actual_group","predicted_value")
    output_result$predicted_value$Iterations<-1:output_result$Iterations
    
    #计算需要随机分组的次数
    round<-ceiling(Iterations/N_fold)
    
    #进度条
    pb <- tkProgressBar("进度","已完成 %", 0, 100)
    
    #生成初始化的迭代计数器
    inter_counter=0
    
    #开始计算
    if(output_result$N_fold>nrow(output_result$clean_dataset)){
      print("交叉验证折数大于数据中变量完整的行数！")
    }else{
      for(index_round in 1:round){
        #无需随机分组
        rand_seq<-1:nrow(output_result$clean_dataset)
        #计算该round内各折的结果
        for(index_fold in 1:output_result$N_fold){
          #检验迭代次数是否足够
          if(inter_counter==output_result$Iterations){
            break
          }else{
            #迭代次数不够继续计算
            #划分训练集和验证集
            temp_training_dataset<-output_result$clean_dataset[(rand_seq!=index_fold),]
            temp_validation_dataset<-output_result$clean_dataset[(rand_seq==index_fold),]
            ddist <<- datadist(temp_training_dataset)
            options(datadist='ddist')
            #建模
            f_reg <- tryCatch(lrm(output_result$f, data=temp_training_dataset, x=TRUE, y=TRUE,maxit=1000),error=function(e){return(list(fail=TRUE))} )
            if(f_reg$fail | 'try-error' %in% class(f_reg)){           # 判断当前循环的try语句中的表达式是否运行正确
              # 若该回归分析无法进行则该轮循环不赋值，默认为NA
            }else{
              #验证集中计算预测值
              pred_f_validation<-predict(f_reg,temp_validation_dataset)
              output_result$predicted_value[(inter_counter+1),"predicted_value"]<-pred_f_validation
              output_result$predicted_value[(inter_counter+1),"actual_group"]<-temp_validation_dataset[,output_result$y_list]
            }
            #进度条
            info <- sprintf("已完成 %d%%",round(inter_counter*100/output_result$Iterations))
            setTkProgressBar(pb, inter_counter*100/output_result$Iterations, sprintf("进度 (%s)", info), info)
            inter_counter<-inter_counter+1
          }
        }
      }
    }
    #计算AUC
    modelroc <- roc(output_result$predicted_value$actual_group,output_result$predicted_value$predicted_value,quiet=TRUE)
    output_result$AUC_value<-modelroc$auc[1]
    #计算C-index
    temp_pred_list<-data.frame(grop=output_result$predicted_value$actual_group,pred=output_result$predicted_value$predicted_value)
    #所有样本两两配对
    valid_num=0
    concordance_num=0
    for(index_temp_pred_list_1 in 1:nrow(temp_pred_list)){
      for(index_temp_pred_list_0 in 1:nrow(temp_pred_list)){
        if((temp_pred_list$grop[index_temp_pred_list_1]==1)&(temp_pred_list$grop[index_temp_pred_list_0]==0)){
          #当匹配的对子中第一个为group=1（实验组）;第二个为group=0（对照组）时为有效对子
          valid_num<-valid_num+1
          if(temp_pred_list$pred[index_temp_pred_list_1]>temp_pred_list$pred[index_temp_pred_list_0]){
            concordance_num<-concordance_num+1
          }
        }
      }
    }
    output_result$C_index_value<-concordance_num/valid_num
    output_result$print_report<-paste0("本次验证方法: Jackknife验证(留一法交叉验证，leave-one-out cross validation)\n"
                                       ,"迭代次数: ",output_result$Iterations," 次\n"
                                       ,"计算所得的AUC: ",output_result$AUC_value,"\n"
                                       ,"计算所得的C_index: ",output_result$C_index_value,"\n")
    return(output_result)
  }
  nomogram_Bootstrap_validation<-function(f,dataset,Iterations=200,seed=NULL){
    
    #设置随机数种子
    if(!is.null(seed)){
      set.seed(seed)
    }
    
    #提取公式中的因变量和自变量
    f<-as.formula(f)
    y_list<-as.character(f[[2]])
    x_list<-as.character(unlist(strsplit(as.character(f[3]),split=" \\+ ")))
    
    #检查列名对不对
    if(!all(c(y_list,x_list) %in% colnames(dataset))){
      output_result<-list()
      print("方程中有变量不在数据列名中！请检查后再运行！")
      return(output_result)
    }
    
    #筛选出式子中包含的列，并筛选数据完整的行
    clean_dataset<-dataset[,(c(y_list,x_list))]
    clean_dataset<-clean_dataset[(complete.cases(clean_dataset)),]
    #生成空的输出对象
    output_result<-list()
    output_result$f<-f
    output_result$y_list<-y_list
    output_result$x_list<-x_list
    
    output_result$Iterations<-Iterations
    output_result$clean_dataset<-clean_dataset
    
    output_result$AUCs<-data.frame(matrix(NA,output_result$Iterations,3),stringsAsFactors = FALSE)
    colnames(output_result$AUCs)<-c("Iterations","AUC","C_index")
    output_result$AUCs$Iterations<-1:output_result$Iterations
    
    #进度条
    pb <- tkProgressBar("进度","已完成 %", 0, 100)
    
    #生成初始化的迭代计数器
    inter_counter=0
    
    #开始计算
    
    for(index_iterations in 1:output_result$Iterations){
      #随机分组
      rand_seq<-1:nrow(output_result$clean_dataset)
      rand_seq<-sample(rand_seq, nrow(output_result$clean_dataset), replace = TRUE)
      
      if(inter_counter==output_result$Iterations){
        break
      }else{
        #迭代次数不够继续计算
        #划分训练集和验证集
        temp_training_dataset<-output_result$clean_dataset
        temp_validation_dataset<-output_result$clean_dataset[c(rand_seq),]
        ddist <<- datadist(temp_training_dataset)
        options(datadist='ddist')
        #建模
        if(length(table(temp_validation_dataset[,output_result$y_list]))==1){
          #验证集的结局事情全为0或全为1时无法计算AUC，默认为NA
        }else{
          f_reg <- tryCatch(lrm(output_result$f, data=temp_training_dataset, x=TRUE, y=TRUE,maxit=1000),error=function(e){return(list(fail=TRUE))} )
          if(f_reg$fail | 'try-error' %in% class(f_reg)){           # 判断当前循环的try语句中的表达式是否运行正确
            # 若该回归分析无法进行则该轮循环不赋值，默认为NA
          }else{
            #验证集中计算AUC
            pred_f_validation<-predict(f_reg,temp_validation_dataset)
            modelroc <- roc(temp_validation_dataset[,output_result$y_list],pred_f_validation,quiet=TRUE)
            output_result$AUCs[(inter_counter+1),"AUC"]<-modelroc$auc[1]
            #验证集中计算C-index
            temp_pred_list<-data.frame(grop=temp_validation_dataset[,output_result$y_list],pred=pred_f_validation)
            #所有样本两两配对
            valid_num=0
            concordance_num=0
            for(index_temp_pred_list_1 in 1:nrow(temp_pred_list)){
              for(index_temp_pred_list_0 in 1:nrow(temp_pred_list)){
                if((temp_pred_list$grop[index_temp_pred_list_1]==1)&(temp_pred_list$grop[index_temp_pred_list_0]==0)){
                  #当匹配的对子中第一个为group=1（实验组）;第二个为group=0（对照组）时为有效对子
                  valid_num<-valid_num+1
                  if(temp_pred_list$pred[index_temp_pred_list_1]>temp_pred_list$pred[index_temp_pred_list_0]){
                    concordance_num<-concordance_num+1
                  }
                }
              }
            }
            #计算C-index
            output_result$AUCs[(inter_counter+1),"C_index"]<-concordance_num/valid_num
          }
        }
        #进度条
        info <- sprintf("已完成 %d%%",round(inter_counter*100/output_result$Iterations))
        setTkProgressBar(pb, inter_counter*100/output_result$Iterations, sprintf("进度 (%s)", info), info)
        inter_counter<-inter_counter+1
      }
      
    }
    
    output_result$mean_AUC<-mean(output_result$AUCs$AUC,na.rm = TRUE)
    output_result$ci_left_AUC<-t.test(output_result$AUCs$AUC,na.rm = TRUE)$conf.int[1]
    output_result$ci_right_AUC<-t.test(output_result$AUCs$AUC,na.rm = TRUE)$conf.int[2]
    output_result$mean_C_index<-mean(output_result$AUCs$C_index,na.rm = TRUE)
    output_result$ci_left_C_index<-t.test(output_result$AUCs$C_index,na.rm = TRUE)$conf.int[1]
    output_result$ci_right_C_index<-t.test(output_result$AUCs$C_index,na.rm = TRUE)$conf.int[2]
    output_result$print_report<-paste0("本次验证方法: Bootstrap验证\n"
                                       ,"迭代次数: ",output_result$Iterations," 次\n"
                                       ,"其中迭代失败次数: ",sum(is.na(output_result$AUCs$AUC))," 次\n"
                                       ,"计算所得的平均AUC: ",output_result$mean_AUC,"\n"
                                       ,"计算所得的AUC的95%置信区间为: ",output_result$ci_left_AUC,", ",output_result$ci_right_AUC,"\n"
                                       ,"计算所得的平均C_index: ",output_result$mean_C_index,"\n"
                                       ,"计算所得的C_index的95%置信区间为: ",output_result$ci_left_C_index,", ",output_result$ci_right_C_index,"\n")
    
    return(output_result)
    
  }
}
