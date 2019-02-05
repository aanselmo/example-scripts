print_table <- function(df,outdir,outname) {
  
  write.table(df, file=paste(outdir,"/QC/",outname,".txt",sep=""), row.names=TRUE, quote=FALSE, sep="\t")
  
}