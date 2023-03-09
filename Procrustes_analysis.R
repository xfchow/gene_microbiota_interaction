############Procrustes analysis#################
geneset<-read.table("LUSC_fpkm.txt",header = T,sep='\t',row.names=1)
microbe<-read.csv("LUSC_microbe_relative.csv",header=T,row.names=1)
geneset.dist<-vegdist(geneset,method = “euclidean”)
geneset.dist
microbe <- decostand(microbe, method = 'hellinger')
microbe.dist <- vegdist(microbe, method = "bray")
mds.gene<-monoMDS(geneset.dist)
mds.microbe<-monoMDS(microbe.dist)
pro.g.s<-procrustes(mds.gene,mds.microbe)
prot <- protest(X = mds.gene, Y = mds.microbe, permutations = how(nperm = 999))
Y <- cbind(data.frame(pro.g.s$Yrot), data.frame(pro.g.s$X))
X <- data.frame(pro.g.s$rotation)
p <- ggplot(Y) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + 
                                                                   MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#B2182B", size = 1) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend 
                   = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#56B4E9", size = 1) +
  geom_point(aes(X1, X2), fill = "#B2182B", size = 4, shape = 23) +
  geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 4, shape = 21) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between Host gene expression and microbiome abundance 
in LUSC") + 
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\n M2 = 0.946, p-value = 0.004',x = 
             -2.2, y = 2.1, size = 4.2,hjust = 0) +
  theme(plot.title = element_text(size=14,colour = "black",hjust = 0,face = "bold"))
p
ggsave('LUSC_Procrustes.pdf',p)
#Mantel test
set.seed(520)
mantel(geneset.dist, microbe.dist, method="spearman", permutations=999)

