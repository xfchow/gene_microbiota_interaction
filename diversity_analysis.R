###############Diversity analysis###############
#######alpha diversity analysis#######
BiocManager::install("rlang")
install.packages('tidyverse')
library('rlang')
library('tidyverse')
install.packages("vegan")
install.packages("ggpubr")
install.packages('RColorBrewer')
install.packages('dplyr')
install.packages('ape')
library(ape)
library(vegan)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
input = 'LUAD_microbe.csv'
#input = 'LUSC_microbe.csv'
df = read.csv(input,row.names = 1,check.names=F)#header=T,row.names = 
1,check.names=F
dt = t(df)
#compute alpha diversity index
Shannon = diversity(df,index = 'shannon',MARGIN=2,base=exp(1))
Simpson = diversity(df,index = 'simpson',MARGIN=2,base=exp(1))
Richness = specnumber(df,MARGIN=2)
index = as.data.frame(cbind(Shannon,Simpson,Richness))
tdf = t(df)
tdf = ceiling(as.data.frame(t(df)))
obs_chao_ace = t(estimateR(tdf))
obs_chao_ace = obs_chao_ace[rownames(index),]
index$Chao = obs_chao_ace[,2]
index$Ace = obs_chao_ace[,4]
index$Sobs = obs_chao_ace[,1]
index$Pielou = Shannon / log(Richness,2)
index$Goods_converage = 1 - colSums(df == 1) / colSums(df)
write.table(index,'diversity_index_LUAD.csv',row.names = T,sep = ',')
#write.table(index,'diversity_index_LUSC.csv',row.names = T,sep = ',')
#plot diversity index plot
df = read.csv('diversity_index_LUAD.csv',row.names = 1,check.names = F)
#df = read.csv('diversity_index_LUSC.csv',row.names = 1,check.names = F)
a<-ggplot(df,aes(x=group,y=Shannon,fill=group)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("Non￾recurrence","Recurrence")),map_signif_level = TRUE,test = t.test,y_position = 
                c(4.3,2),tip_length = c(0.13,0.13))+theme_bw()+ 
  scale_fill_manual(values = c("#DE6757","#5B9BD5"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = 
          element_blank())+
  labs(y="Shannon's Diversity Index")+theme(axis.title.x 
                                            =element_text(size=12,face = "bold"),axis.title.y=element_text(size=12,face = 
                                                                                                             "bold"),axis.text = element_text(size = 12,face = "bold"),legend.text = 
                                              element_text(size = 11))
b<-ggplot(df,aes(x=group,y=Simpson,fill=group)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("Non￾recurrence","Recurrence")),map_signif_level = TRUE,test = t.test,y_position = 
                c(1.05,5),tip_length = c(0.09,0.09))+theme_bw()+ 
  scale_fill_manual(values = 
                      c("#DE6757","#5B9BD5"))+theme(panel.grid.major = 
                                                      element_blank(),panel.grid.minor = element_blank())+
  labs(x="Group", y="Simpson's Diversity Index")+
  theme(axis.title.x =element_text(size=12,face = "bold"), 
        axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size 
                                                                                  = 12,face = "bold"),legend.text = element_text(size = 11))
c<-ggplot(df,aes(x=group,y=Richness,fill=group)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("Non￾recurrence","Recurrence")),map_signif_level = TRUE,test = t.test,y_position = 
                c(1100,0),tip_length = c(0.05,0.05))+
  theme_bw()+ 
  scale_fill_manual(values = 
                      c("#DE6757","#5B9BD5"))+theme(panel.grid.major = 
                                                      element_blank(),panel.grid.minor = element_blank())+
  labs(x="Group", y="Richness")+
  theme(axis.title.x =element_text(size=12,face = "bold"), 
        axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size 
                                                                                  = 12,face = "bold"),legend.text = element_text(size = 11))
d<-ggplot(df,aes(x=group,y=Chao,fill=group)) +
  geom_boxplot()+
  geom_signif(comparisons = list(c("Non￾recurrence","Recurrence")),map_signif_level = TRUE,test = t.test,y_position = 
                c(1230,0),tip_length = c(0.04,0.04))+theme_bw()+ 
  scale_fill_manual(values = c("#DE6757","#5B9BD5"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = 
          element_blank())+
  labs(x="Group", y="Chao")+
  theme(axis.title.x =element_text(size=12,face = "bold"), 
        axis.title.y=element_text(size=12,face = "bold"),axis.text = element_text(size 
                                                                                  = 12,face = "bold"),legend.text = element_text(size = 11))
pdf("LUAD_diersity.pdf",width = 12,height = 11)
#pdf("LUSC_diersity.pdf",width = 12,height = 11)
ggarrange(a,b,c,d, ncol = 2, nrow = 2)
dev.off()
######## beta diversity analysis#######
dune = read.csv('LUSC_microbe.csv',row.names = 1,check.names = F)
#dune = read.csv('LUAD_microbe.csv',row.names = 1,check.names = F)
dune = t(dune)
dune.env = read.csv('list.csv',row.names=1)
dune.dist <- vegdist(dune)#构造距离矩阵。
res=pcoa(dune.dist,correction="cailliez")
biplot(res)
dune.ano <- with(dune.env, anosim(dune.dist, label))
summary(dune.ano)
df.plot<-data.frame(res$vectors)
head(df.plot)
pcoa <- dudi.pco(dune.dist, scan=FALSE,nf=10)
pcoa_eig <- (pcoa$eig)[1:10]/sum(pcoa$eig)
pcoa_exp = pcoa$eig/sum(pcoa$eig)
x_label = round(100*pcoa_exp[1],2)
y_label = round(100*pcoa_exp[2],2)
p = ggplot(data=df.plot,aes(x=x,y=y,
                            color=type))+
  # geom_point(size=5)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  # geom_vline(xintercept = 0,lty="dashed")+
  # geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))+
  stat_ellipse(data=df.plot,
               # geom = "polygon",
               # aes(fill=type),
               alpha=0.8)+
  geom_point(aes(color = type),size=6)+
  scale_color_manual(values=c("#DE6757","#5B9BD5"))
p
ggsave(‘LUSC_PCoA_plot.pdf',p)
#ggsave(‘LUAD_PCoA_plot.pdf',p)