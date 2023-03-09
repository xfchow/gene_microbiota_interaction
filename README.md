R and python framework for identifying associations between lung squamous carcinoma and lung adenocarcinoma microbiome composition and host gene expression, 
and predicting recurrence metastasis and survival in LUSC patients. These scripts were used in the study by Xiangfeng Zhou et al. 
"Intratumoral microbiota-host interactions shape the variability of LUAD and LUSC in recurrence and metastasis" ( Under review)

Description of scripts 

-diversity_analysis.R
The script was used to analyze the diversity of tissue microorganisms in LUSC versus LUAD patients. Among them, alpha diversity was measured by Simpson's index, 
Richness's index, Chao's index and Shannon's index. The beta diversity of the tissue microorganisms of LUSC versus LUAD patients was analyzed by PCoA.

-Procrustes_analysis.R
By performing this script, we implemented Procrustes analysis and Mantel test for host gene expression and tissue microbita abundance.

-collection_microbiota_at_all_levels.py
This script is used to collate microbial abundance data at the phylum, class, order, family, and genus levels for LUSC and LUAD patients.

-Spearman_correlation_analysis.py
We used the script to calculate the correlation between genes and microbiota.

-RM_classification.py
In this script, we predict the recurrence and metastasis of LUSC patients by three machine learning classification models, Adaboost, RandomForest and NaiveBayes, 
and evaluate the generalization ability of the classification models by a five-fold cross-validation.

-Survival_analysis.R
This script was used to analyze the postoperative survival of LUSC patients.The survival analysis regression model was constructed by the 'coxph' function of the R package, 
with biomarker as the model-independent variable and survival time as the model-dependent variable. 
The survival curve was performed by using the Kaplan–Meier method and the log-rank test was used to compare the difference in survival probability with the R package “survival”.
