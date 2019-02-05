assemble_table <- function(files,samplenames,gene_id_colname,colname) {

i = 0
for (file in files) {
  mydata <- read.delim(file=file,header=TRUE,skip=0)
  if (i == 0) {
    df<- mydata[,c(gene_id_colname,colname)]
    
    i = i+1
  }
  else {df <- cbind(df,mydata[,colname])}
}

row.names(df) <- df[,gene_id_colname]
df[,gene_id_colname] <- NULL
colnames(df) <- samplenames
return(df)
}
