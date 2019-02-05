library(gplots)
plot_heatmap <- function (mat,outdir,outname){
  
  postscript(paste(outdir,"/",outname,".ps",sep=""), paper="letter");
  #pdf(paste(outdir,"/",outname,".pdf",sep=""))
  rowptsize <- "";
  if (dim(mat)[1] > 50) {
    rowptsize <- 0.1
  } else {
    rowptsize <- 0.4
  }
  
  heatmap.2(mat, distfun = function(x) dist(x,method = 'euclidean'), 
            hclustfun = function(x) hclust(x,method = 'average'), 
            dendrogram="both", 
            breaks=c(seq(-1,1,0.05)), 
            srtCol=45, 
            ##colsep=1:ncol(rpkm_table_for_heatmap_matrix_unique[c(1:400),]), 
            ##rowsep=1:nrow(rpkm_table_for_heatmap_matrix_unique[c(1:400),]), 
            ##sepcolor="black", sepwidth=c(0,0),
            col=colorRampPalette(c("blue",rep("white",8),"red"))(40), 
            symkey=F, 
            symbreaks=F,
            cexRow=rowptsize, 
            cexCol=0.6, 
            trace="none", 
            ###cellnote=round(as.matrix(ColonCa_cor_matrix),2), notecol="black", notecex=0.4,
            main="Heatmap of Pearson Correlation Values", ylab="",  key.xlab="Pearson Coefficient");
  dev.off();
}