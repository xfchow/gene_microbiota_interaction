############# Survival analysis##############
#DFI
#data <- read.csv('survival_OS.csv',row.names = 1)
#data <- read.csv('survival_PFI.csv',row.names = 1)
#data <- read.csv('survival_DSS.csv',row.names = 1)
data <- read.csv('survival_DFI.csv',row.names = 1)

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