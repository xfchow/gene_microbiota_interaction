#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 09:06:32 2023

@author: zhou
"""

############Machine learning classification models############
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix,accuracy_score,f1_score,roc_auc_score,recall_score,precision_score,roc_curve
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import KFold
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import cross_val_score
from matplotlib import pyplot
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import cross_validate
from sklearn.metrics import roc_curve, auc
from sklearn import metrics 
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from collections import Counter
from sklearn.model_selection import StratifiedKFold

LUSC_gene = pd.read_csv(r'LUSC_fpkm.csv',index_col=0)
LUSC_microbe = pd.read_csv(r'LUSC_microbe.csv',index_col=0)
label = pd.read_csv(r'label.csv',index_col = 0)

LUSC_all = pd.merge(LUSC_gene, LUSC_microbe, left_index = True, right_index = True)
x = LUSC_all
y = label

scores_1=[]
scores_2=[]
scores_3=[]
scores_4=[]
scores_5=[]
scores_6=[]
scores_7=[]
scores_8=[]

tprs=[]
aucs=[]

mean_fpr_BV=np.linspace(0,1,100)
fig,ax=plt.subplots(1,1,figsize=[8,8],dpi=300)
skf=StratifiedKFold(n_splits=5, shuffle=True, random_state=44)

for i,(train_index,test_index) in  enumerate(skf.split(x,y)):
    print("i",i)
    x_train,x_test=x.iloc[train_index,:],x.iloc[test_index,:]
    y_train,y_test=y.iloc[train_index],y.iloc[test_index]
   
   
    
    model=RandomForestClassifier(random_state=18)
    #model=AdaBoostClassifier(random_state=18)
    #model = LogisticRegression()
    #model = SVC(probability=True)
    #model = GaussianNB()
    #parameters = {'max_depth':np.arange(1,10,1)}
    #grid_search = GridSearchCV(RandomForestClassifier(), parameters,  verbose=0, scoring='roc_auc', cv=15)
    #grid =grid_search.fit(x_train, y_train)
    #best_parameters= grid_search.best_estimator_.get_params()
    #model = RandomForestClassifier(max_depth=best_parameters['max_depth'])
    
    model.fit(x_train, y_train)   
    y_proba=model.predict_proba(x_test)[:,1]
    y_pred = model.predict(x_test)
    
    auc = metrics.roc_auc_score(y_test, y_proba) 
    acc = accuracy_score(y_test, y_pred)
    
    precision, recall, _thresholds = metrics.precision_recall_curve(y_test, y_proba)
    pr_auc = metrics.auc(recall, precision)
    mcc = matthews_corrcoef(y_test, y_pred)
    f1=f1_score(y_test, y_pred, average='macro')  
    
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    conf=confusion_matrix(y_test, y_pred)
    
    total=tn+fp+fn+tp
    recall = float(tp)/float(tp+fn)
    sps = float(tn)/float((tn+fp))
    precision =float(tp)/float(tp+fp)
   
    print("auc",auc)

    scores_1.append(auc)
    scores_2.append(acc)
    scores_3.append(pr_auc)
    scores_4.append(precision)
    scores_5.append(recall)
    scores_6.append(sps)
    scores_7.append(mcc)
    scores_8.append(f1)
    
    #i=i+1
     
    viz=plot_roc_curve(model, x_test, y_test,
                       name='ROC fold {}'.format(i+1),
                       alpha=0.3,lw=1,ax=ax)
    interp_tpr=np.interp(mean_fpr_BV,viz.fpr,viz.tpr)
    interp_tpr[0]=0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
    
ax.plot([0, 1],[0, 1],linestyle='-',lw=2,color='r',label='Chance', alpha=.8 )    


mean_tpr_BV=np.mean(tprs,axis=0)
mean_tpr_BV[-1]=1.0
mean_auc_BV=metrics.auc(mean_fpr_BV,mean_tpr_BV)
std_auc=np.std(aucs)
ax.plot(mean_fpr_BV,mean_tpr_BV,color='b',
        label=r'Mean ROC (AUC= %0.2f $\pm$ %0.2f)' % ( mean_auc_BV,std_auc),
        lw=2,alpha=.8)

std_tpr=np.std(tprs,axis=0)
tprs_upper=np.minimum(std_tpr + mean_tpr_BV,1)
tprs_lower=np.maximum(mean_tpr_BV - std_tpr,0)


ax.set(xlim=[-0.05,1.05],ylim=[-0.05,1.05])

plt.tick_params(labelsize=20)
font2={'size':23}

plt.xlabel('False Positive Rate',font2)

plt.ylabel('True Positive Rate',font2)
ax.legend(fontsize=10,loc='lower right')

plt.savefig('ROC-RF曲线.pdf')
plt.show()
                                                  
print('micro-5-CV-AUC:',np.mean(scores_1))
print('micro-5-CV-accuracy:',np.mean(scores_2))
print('micro-5-CV-pr_auc:',np.mean(scores_3))
print('micro-5-CV-precision:',np.mean(scores_4))
print('micro-5-CV-recall:',np.mean(scores_5))
print('micro-5-CV-sps:',np.mean(scores_6))
print('micro-5-CV-mcc:',np.mean(scores_7))
print('micro-5-CV-f1:',np.mean(scores_8))

#############Calculating feature importance#############
from sklearn import tree
from sklearn import datasets, model_selection

y = pd.read_csv('label.csv',index_col = 0)
X = pd.read_csv('LUSC_microbe.csv', index_col = 0)
#Dividing the training set and test set
x_train,x_test, y_train, y_text = model_selection.train_test_split(X, y, test_size=0.3, random_state=35)

classification_tree = tree.DecisionTreeClassifier()
classification_tree.fit(x_train, y_train)

#Export the importance of each feature
print(classification_tree.feature_importances_)
importance = classification_tree.feature_importances_

############# Survival analysis##############
#DPI
data <- read.csv('survival_OS.csv',row.names = 1)
#data <- read.csv('survival_PFI.csv',row.names = 1)
#data <- read.csv('survival_DSS.csv',row.names = 1)
#data <- read.csv('survival_DFI.csv',row.names = 1)

res.cox <- coxph(Surv(DFI.time, DFI) ~ MEG3 + MIR3153 + RP11.568J23.6 + C9orf24 + C1orf111
                 + C20orf85 + MORN5 + MYBPC1 + RPS6P22 + CATSPERD + AC004158.3 + RP11.22B23.2
                 + RP11.177G23.1 + RPL15P18 + MKRN9P + HGFAC + RP4.756H11.4 + HSD3BP4 + RP11.684B2.3
                 + RP11.556K13.1 + TEKT1 + LINC00675 + MLIP + THSD7A + RP11.981G7.2 + RPSAP56 + RP11.350G13.1
                 + g__Shewanella + g__Pseudogulbenkiania + f__Comamonadaceae + g__Shigella + f__Chromobacteriaceae
                 + f__Shewanellaceae + g__Acidovorax,data =  data)

summary(res.cox)
ph_hypo_multi <- cox.zph(res.cox)
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
formula_for_multivariate <- as.formula(paste0('Surv(DFI.time, DFI)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep ='',collapse = '+')))
multi_variate_cox_2 <- coxph(formula_for_multivariate, data = data)
survival_cancer = read.csv('survival_DFI1.csv')
riskscore <- function(survival_cancer, candidate_genes_for_cox, cox_report){
  library('dplyr')
  risk_score_table <- survival_cancer[,candidate_genes_for_cox]
  for(each_sig_gene in colnames(risk_score_table)){
    risk_score_table[,each_sig_gene] <- risk_score_table[,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table, 'total_risk_score'=exp(rowSums(risk_score_table)))%>%
    cbind(survival_cancer[c('X','DFI.time','DFI')])
  risk_score_table <- risk_score_table[,c('X','DFI.time','DFI',candidate_genes_for_cox,'total_risk_score')]
  risk_score_table
}

candidate_genes_for_cox2 <- c(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(survival_cancer, candidate_genes_for_cox2,multi_variate_cox_2)

multi_ROC <- function(time_vector, risk_score_table){
  library('survivalROC')
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$DFI.time,
                           status = risk_score_table$DFI,
                           marker = risk_score_table$total_risk_score,
                           predict.time = single_time,method='KM')
    data.frame('True_positive'=for_ROC$TP,'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values,'Time_point'=rep(single_time,length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC,length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector,single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC(time_vector = 365*5,risk_score_table = risk_score_table_multi_cox2)
AUC_max <- max(for_multi_ROC$AUC)
AUC_max_time <- for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
AUC_max_time <- AUC_max_time[!duplicated(AUC_max_time)]
AUC_max_time <- AUC_max_time[length(AUC_max_time)]
for_multi_ROC$Time_point <- as.factor(for_multi_ROC$Time_point)
#optimal_time_ROC_df <- for_multi_ROC[which(for_multi_ROC$Time_point == AUC_max_time),]
optimal_time_ROC_df <- for_multi_ROC[which(for_multi_ROC$Time_point == 1825),]
cut.off <- optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive-optimal_time_ROC_df$False_positive)]
high_low <- (risk_score_table_multi_cox2$total_risk_score > cut.off)
high_low[high_low == TRUE] <- 'high'
high_low[high_low == FALSE] <- 'low'
risk_score_table_multi_cox2 <- cbind(risk_score_table_multi_cox2, high_low)
risk_score_table_multi_cox2$DFI[which(risk_score_table_multi_cox2$DFI.time > 1825)] <- 0#AUC_max_time
risk_score_table_multi_cox2$DFI.time[which(risk_score_table_multi_cox2$DFI.time > 1825)] <- AUC_max_time#AUC_max_time
fit_km <- survfit(Surv(DFI.time, DFI) ~high_low, data=risk_score_table_multi_cox2)

pdf('KM-curve-OS.pdf')
#pdf('KM-curve-PFI.pdf')
#pdf('KM-curve-DSS.pdf')
#pdf('KM-curve-DFI.pdf')
ggsurvplot(fit_km,data=risk_score_table_multi_cox2,
           font.legend=16, font.x=16, font.y=16,font.tickslab=14,
           pval = TRUE, conf.int = F,
           legend.title='total risk score',
           legend.labs=c(paste0('High risk'),paste0('Low risk')), 
           risk.table = T,
           #risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           #linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(),
           pval.method=TRUE,# Change ggplot2 theme
           palette = c('red', 'blue'),xlim = c(0, 2000))
dev.off()

