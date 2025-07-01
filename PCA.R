merge_one_batch<-read.table('/Users/kdkt095/Downloads/merge_one_batch_fivedata.txt',sep = '\t')
install.packages("factoextra")
library(factoextra)

d<-as.matrix(merge_one_batch)
dim(d)
10281*0.5 ###
mads=apply(d,1,mad) ### measured by median absolute deviation
d=d[rev(order(mads))[1:5141],] ### top 50% variable genes

data_t <- t(as.matrix(d))
variableL <- ncol(data_t)
pca <- prcomp(data_t[,1:variableL], scale=T)
fviz_pca_ind(pca, repel=T)
data_t<-as.data.frame(data_t)
data_t$conditions<-as.character(sapply(colnames(merge_one_batch),function(x){
  strsplit(x,'_')[[1]][1]
}))

p<-(fviz_pca_ind(pca, col.ind=data_t$conditions, mean.point=F, addEllipses = T, legend.title="Disease status",
                 geom.ind = "point")+scale_shape_manual(values=seq(0,17))+theme_classic())

p





library(ggplot2)
library(factoextra)

# 1. 先提取得分矩阵
pca_scores <- as.data.frame(pca$x) # 得到PC1-PCn

pca_scores$conditions <- data_t$conditions # 加上分组信息


pca_scores_scaled <- as.data.frame(scale(pca_scores[, paste0("PC", 1:20)])) # 标准化每列
pca_scores_scaled$conditions <- pca_scores$conditions

# 2. 设置你需要的主成分数，比如前20维
pcs <- paste0("PC", 1:10)

# 3. 循环两两组合
PCA_distance<-as.data.frame(matrix(data=NA,
                                   nrow =190,ncol=2 ))
colnames(PCA_distance)<-c('PC','Distance')
plots <- list()
counter <- 1

i=1
j=4

pi <- pcs[i]
pj <- pcs[j]
p <- ggplot(pca_scores, aes_string(x=pi, y=pj, color='conditions')) +
  geom_point() +
  stat_ellipse() +
  theme_classic() +
  labs(title=paste("PCA:", pi, "vs", pj))
ggsave(filename=paste0("./PCA_plot/PCA_", i,'_',j, ".pdf"), 
       plot=p, width=5, height=4)

for (i in 1:(length(pcs)-1)) {
  for (j in (i+1):length(pcs)) {
    pi <- pcs[i]
    pj <- pcs[j]
    p <- ggplot(pca_scores, aes_string(x=pi, y=pj, color='conditions')) +
      geom_point() +
      stat_ellipse() +
      theme_classic() +
      labs(title=paste("PCA:", pi, "vs", pj))
    plots[[counter]] <- p
    names(plots)[[counter]]<-paste0(i,'_',j)
    PCA_distance[counter,1]<-paste0('PC_',i,'_',j)
    
    aggregate_coords <- aggregate(. ~ conditions, 
                                  data = pca_scores_scaled[, c(paste0("PC", i),
                                                               paste0("PC", j), 
                                                               "conditions")], mean)
    # 取两个簇的中心
    center1 <- as.numeric(aggregate_coords[1, c(paste0("PC",i),
                                                paste0("PC",j))])
    center2 <- as.numeric(aggregate_coords[2, c(paste0("PC",i),
                                                paste0("PC",j))])
    # 计算欧氏距离
    dist <- sqrt(sum((center1 - center2)^2))
    
    PCA_distance[counter,2]<-dist
    counter <- counter + 1
  }
}

# 4. 按需展示或保存，例如一次展示4个
library(gridExtra)
grid.arrange(grobs=plots[1:4], ncol=2)

# 也可以自动批量保存
for(idx in seq_along(plots)) {
  ggsave(filename=paste0("./PCA_plot/PCA_", names(plots)[[idx]], ".png"), plot=plots[[idx]], width=5, height=4)
}



print(dist)



loadings <- pca$rotation  # 行是基因，列是PC

N <- 500  # 可以改为1即只找最大
# 绝对值排序，取前N个名字
top_genes_PC1 <- names(sort(abs(loadings[,pc]), decreasing=TRUE))[1:N]
print(top_genes_PC1)


N <- 100  # 最大的就是1
DEG_up<-read.table('/Users/kdkt095/Yilu_Work/endometrosis/result/microarray_DEG_up.txt')
DEG_down<-read.table('/Users/kdkt095/Yilu_Work/endometrosis/result/microarray_DEG_down.txt')

DEG<-unique(c(DEG_down$gene,DEG_up$gene))
PC_DEG<-as.data.frame(matrix(data=NA,
                             nrow = 10,
                             ncol=3))

colnames(PC_DEG)<-c('PC','N','Group')
PC_DEG$Group<-'Non-Separated PC'
PC_DEG$Group[c(1,2,4,5)]<-'Separated PC'
i=1
for (pc in c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')) {
  top_genes_PC1 <- names(sort(abs(loadings[,pc]), decreasing=TRUE))[1:N]
  print(paste0(pc,': ',length(intersect(top_genes_PC1,DEG))/length(top_genes_PC1 )*100))
  PC_DEG[i,1]<-pc
  PC_DEG[i,2]<-length(intersect(top_genes_PC1,DEG))/length(top_genes_PC1 )*100
  i=i+1
}

colnames(PC_DEG)

my_comparisons_train <- list(c('Separated PC','Non-Separated PC'))
ggviolin(PC_DEG, x = "Group", y = 'N',
         fill = "Group",
         add = "boxplot", add.params = list(fill = "white"),palette = c("#D22C6C", "#6688AB"))+
  stat_compare_means(comparisons = my_comparisons_train,size=5)+
  ylab("Number of Overlapping Genes\n") +
  xlab('')+
  #scale_y_continuous(breaks=seq(0, 1, 0.1))+
  #coord_cartesian(ylim = c(-2500,3000))+
  #geom_hline(yintercept = median(gsva_MYC_matrix[Purity_Test$Group=='MYCN_AMP',]$ImmuneScore), 
  #          size=2,linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        
        axis.text.x = element_text(  size=16,color = 'black',hjust = 0.5,vjust=0),
        axis.text.y = element_text( size=16,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=20,vjust=0),
        plot.margin = unit(c(1,1,1,1),"cm"))
ggsave('PC_DEG.pdf',width = 6,height = 6)



