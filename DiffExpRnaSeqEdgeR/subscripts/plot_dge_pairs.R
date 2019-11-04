plot_pairs <- function (cpms,outdir,outname){
  
  log2_cpms <- log2(cpms);
  
  pdf(paste(outdir,"/Pairwise_plots_",outname,".pdf",sep=""), width=22, height=22);
  #postscript(paste("QC/",outname,".ps",sep=""), paper="letter");
  pairs(log2_cpms,main=paste("Pairs: log2 CPMs  ",outname,sep=""), pch=16, cex=0.35, col="black" );
  dev.off();
}
