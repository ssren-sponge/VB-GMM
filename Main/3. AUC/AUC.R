

###AUC
require(Hmisc)
library(pheatmap)

plotDir <- "./outputs"
outputDir <- "./outputs"

Bay_Label <- read.csv(paste0(outputDir,"/labels_KSneuro_k14.csv"),row.names = NULL)
Bay_Label <- as.matrix(Bay_Label +1)
load(paste0(outputDir,'/GW1722_typeMeans.rd'))

colIdx<-c("chartreuse4","slateblue3","tomato3","goldenrod3","thistle1","sienna1","skyblue1","turquoise3","violet","violetred","red","bisque3","aquamarine3","darkgrey","darkolivegreen3","darkseagreen3","navajowhite3","rosybrown3","dimgrey","darkorchid3","deepskyblue2","skyblue3","paleturquoise1","saddlebrown","tan3");#defining good colours for plotting
my_colors <- c( colorRampPalette(c("turquoise3", "white"))(7)[-7], colorRampPalette(c("white", "violetred"))(7)[-1])


##VI AUC
allAUCs<-sapply(sort(unique(Bay_Label)), function(tempKi){
  tempY<-as.vector(as.numeric(Bay_Label==tempKi))
  tempAUCs<-apply(typeMeans, 1, function(TempTMeans){
    return((somers2(TempTMeans, tempY) [["Dxy"]]+1)/2)
  })
  return(tempAUCs)
})

plotData <- t(allAUCs[rev(rownames(allAUCs)),])
rownames(plotData) <- paste0('c', 1:nrow(plotData))
tempk <- 14 #

pdf(file=paste0(plotDir,"/Figure8b.pdf"), useDingbats=F);{
  pheatmap(plotData, display_numbers = TRUE, number_color = "black", fontsize_number = 12,
           fontsize_row = 15,
           fontsize_col = 15,
           color = my_colors,
           cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 315,
           main = paste0("AUC Heatmap of VB (K = ",tempk,")"),
           width = 6, height = 4)
}
dev.off()
  



##GMM AUC
GMM_label <- read.csv(paste0(outputDir,'/gw1722_k14_GMM.csv'))
GMM_label <- as.matrix(GMM_label)

allAUCs<-sapply(sort(unique(GMM_label)), function(tempKi){
  tempY<-as.vector(as.numeric(GMM_label==tempKi))
  tempAUCs<-apply(typeMeans, 1, function(TempTMeans){
    return((somers2(TempTMeans, tempY) [["Dxy"]]+1)/2)
  })
  return(tempAUCs)
})

plotData <- t(allAUCs[rev(rownames(allAUCs)),])
rownames(plotData) <- paste0('c', 1:nrow(plotData))

pdf(file=paste0(plotDir,"/Figure8a.pdf"), useDingbats=F);{
  pheatmap(plotData, display_numbers = TRUE, number_color = "black", fontsize_number = 12,
           fontsize_row = 15,
           fontsize_col = 15,
           color = my_colors,
           cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 315,
           main = paste0("AUC Heatmap of GMM (K = ",tempk,")"),
           width = 6, height = 4)
}
dev.off()
