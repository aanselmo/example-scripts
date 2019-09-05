#!/usr/bin/Rscript

script.dir <- system("find ~/ -name \"DiffExpRnaSeqEdgeR\" 2>/dev/null", intern=TRUE)[1]
print(script.dir)
# Preferably, put DiffExpRnaSeqEdgeR directory in a home directory!!!

date <- Sys.Date()
print(date)

##### SUBS/FUNCTIONS ########
source(paste(script.dir,"/subscripts/get_args.R",sep=""))
source(paste(script.dir,"/subscripts/print_df_table.R",sep=""))
source(paste(script.dir,"/subscripts/plot_mat_heatmap.R",sep=""))
source(paste(script.dir,"/subscripts/plot_dge_pairs.R",sep=""))
source(paste(script.dir,"/subscripts/assemble_table.R",sep=""))

#source("subscripts/get_args.R")
#source("subscripts/print_df_table.R")
#source("subscripts/plot_mat_heatmap.R")
#source("subscripts/plot_dge_pairs.R")
#source("subscripts/assemble_table.R")
library(methods);
library(edgeR);

#ARGS

#--inputdir
#--outputdir
#--groupfile
#--compfile 
#--plot (0=NO or 1=YES)

args <- get_args()

#args <- read.table("temp.arg.txt", header=FALSE, quote="", comment.char="#", stringsAsFactors=FALSE)
#colnames(args) <- c("inputdir","outputdir","groupfile","comparisonfile")


#DIRECTORIES
datadir <- args$inputdir
outputdir <- args$outputdir
pairsbool <- args$plot
#pairsbool <- 1

dir.create(outputdir, showWarnings=FALSE, recursive=FALSE)

#GROUPINGS FILE
groupfile <- args$groupfile
compfile <- args$compfile
print(groupfile)
print(compfile)

groupings <- read.table(groupfile, header=FALSE, quote="", comment.char="#", stringsAsFactors=FALSE)
colnames(groupings) <- c("samp","grp")

comparisons <- read.table(compfile, header=FALSE, quote="", comment.char="#", stringsAsFactors=FALSE)
colnames(comparisons) <- c("group1","group2","logfold_cutoff","fdr_cutoff")

#outname <- paste(comparisons$group2,"_vs_",comparisons$group1,sep="")
outname <- "All_Samples"
fdr_cutoff <- comparisons$fdr_cutoff
log_cutoff <- comparisons$logfold_cutoff


##########_READ_IN_DATA_TABLES_#######################################
countfiles <- list.files(datadir, pattern="*.results*", full.names = TRUE)
filenames <- list.files(datadir, pattern="*.results*", full.names = FALSE)
sample_names <- unlist(lapply(strsplit(filenames,"\\."), function(x) x[1]))
print(sample_names)

######## Generate 1)CountsTable, 2) TPM Table, 3) FPKM Table where rownames are gene ids
expected_counts_df <- assemble_table(countfiles,sample_names,"gene_id", "expected_count")
TPM_df <- assemble_table(countfiles,sample_names,"gene_id", "TPM")
FPKM_df <- assemble_table(countfiles,sample_names,"gene_id","FPKM")


### Choose the samples to test from the Groupings txt file 
### and round to the nearest integer
chosen_count_data <-round(expected_counts_df[,c(groupings$samp)],digits=0)
chosen_tpm_data <- TPM_df[,c(groupings$samp)]
chosen_fpkm_data <- FPKM_df[,c(groupings$samp)]

#save.image(paste("Analysis_edgeR_",date,".RData",sep=""))

## WRITE COUNTS TABLE
counts_to_write <- cbind( row.names(as.data.frame(chosen_count_data)), as.data.frame(chosen_count_data)  )
colnames(counts_to_write)[1] <- "ID"
write.table(counts_to_write,file=paste(outputdir,"/",outname,"_complete_COUNTS_table.txt",sep=""),
            sep="\t",row.names=FALSE, quote=FALSE);

## WRITE TPM TABLE
tpm_to_write <- cbind( row.names(as.data.frame(chosen_tpm_data)), as.data.frame(chosen_tpm_data)  )
colnames(tpm_to_write)[1] <- "ID"
write.table(tpm_to_write,file=paste(outputdir,"/",outname,"_complete_TPM_table.txt",sep=""),
            sep="\t",row.names=FALSE, quote=FALSE);

## WRITE FPKM TABLE
fpkm_to_write <- cbind( row.names(as.data.frame(chosen_fpkm_data)), as.data.frame(chosen_fpkm_data)  )
colnames(fpkm_to_write)[1] <- "ID"
write.table(fpkm_to_write,file=paste(outputdir,"/",outname,"_complete_FPKM_table.txt",sep=""),
            sep="\t",row.names=FALSE, quote=FALSE);


############ DGEList
dge <- DGEList(chosen_count_data, group=groupings$grp)
dge$samples$lib.size <- colSums(dge$counts)

cpm.dge <- cpm(dge)


### WRITE CPM TABLE
cpm_to_write <- cbind( row.names(as.data.frame(cpm.dge)), as.data.frame(cpm.dge)  )
colnames(cpm_to_write)[1] <- "ID"
write.table(cpm_to_write,file=paste(outputdir,"/",outname,"_complete_CPM_table.txt",sep=""),
            sep="\t",row.names=FALSE, quote=FALSE);


### DISCARD GENES with LOW READ COUNTS
keep_genes_idx <- rowSums(chosen_fpkm_data>0) > 0 &
  rowSums(dge$counts) >= (5*ncol(cpm.dge)) &
  rowSums(cpm.dge>0) > 0 

chosen_fpkm_data_kept <- chosen_fpkm_data[keep_genes_idx,]
cpm.dge_kept <- cpm.dge[keep_genes_idx,]

dge <- dge[keep_genes_idx,]
dge$samples$lib.size <- colSums(dge$counts)
#dge <- calcNormFactors(dge);
cpm.dge <- cpm(dge)

save.image(paste("Analysis_edgeR_",date,".RData",sep=""))

######## SAMPLE CORRELATION (Pearson), PAIRWISE PLOTS
#qcdir <- paste(outputdir,"/QC",sep="")
#dir.create(qcdir, showWarnings=FALSE, recursive=FALSE)

correlation_matrix <- cor(cpm.dge)
correlation_df <- cbind(row.names(as.data.frame(correlation_matrix)), as.data.frame(correlation_matrix))
colnames(correlation_df)[1] <- "ID"

write.table(correlation_df,file=paste(outputdir,"/",outname,"_pearson_cor_table.txt",sep=""),
            sep="\t",row.names=FALSE, quote=FALSE);

plot_heatmap(correlation_matrix,outputdir,paste(outname,"_corr_matrix",sep=""))
if (pairsbool == 1) { plot_pairs(cpm.dge,outputdir,outname) }
    
#### MDS plot ##############
pdf(file=paste(outputdir,"/",outname,"_MDSplot.pdf",sep=""), height=6, width=6)
plotMDS(dge, main=paste("MDS Plot__ all genes ", outname, sep=""),cex=0.5)
dev.off()

    

#####
#####
for (i in 1:dim(comparisons)[1]) {
#i <- 1
#### Exact Tests
pairs_to_test <- c(comparisons$group1[i],comparisons$group2[i])

samples2test_idx <- groupings$grp %in% pairs_to_test
samples2test <- groupings$samp[samples2test_idx]


d2test <- dge[,samples2test]
d2_keep_genes_idx <- rowSums(chosen_fpkm_data_kept >0) > 0 &
  rowSums(d2test$counts) >= (5*length(samples2test)) &
  rowSums(cpm.dge_kept >0) > 0 
d2test <- d2test[d2_keep_genes_idx,]
d2test$samples$lib.size <- colSums(d2test$counts)
d2test <- calcNormFactors(d2test);

###### DISPERSION ######################

save.image(paste("Analysis_edgeR_",date,".RData",sep=""))
et <- exactTest(d2test,pair=pairs_to_test,dispersion=0.10)

save.image(paste("Analysis_edgeR_",date,".RData",sep=""))

#### EXACT TEST ########################
if (length(samples2test) > 5) {
  d2test <- estimateCommonDisp(d2test, verbose = TRUE)
  d2test <- estimateTagwiseDisp(d2test, verbose = TRUE) 
  et <- exactTest(d2test,pair=pairs_to_test)
} else {
   et <- exactTest(d2test,pair=pairs_to_test,dispersion=0.10)
}
etp <- topTags(et, n=nrow(et), adjust.method="BH", sort.by="none");

test_order <- sort(pairs_to_test,decreasing=FALSE)

	if (test_order[1] != comparisons$group1[i]){
		etp$logFC <- -1*etp$logFC
	}

compname <- paste(comparisons$group2[i],"_vs_",comparisons$group1[i], sep="")

#############_PRINT  MA_plots_###################################

pdf(paste(outputdir,"/MA_plot_", compname, "_fdr_",fdr_cutoff[i],".pdf", sep=""));
plot(etp$table$logCPM,etp$table$logFC, 
       xlab="Log2 CPM avg", 
       ylab=paste("Log2 FC: ",comparisons$group2[i],"/",comparisons$group1[i], sep=""), 
       main=paste("MA_plot: ", compname, sep=""), 
       xlim=c(-3, 20), 
       ylim=c(-10, 10),
       pch=20, 
       cex=ifelse(etp$table$FDR < fdr_cutoff[i], 0.5, 0.3 ),
       col=ifelse(etp$table$FDR < fdr_cutoff[i] & etp$table$logFC > 1 ,"red", 
                  ifelse(etp$table$FDR < fdr_cutoff[i] & etp$table$logFC < (-1), "green",
                         ifelse(etp$table$FDR > fdr_cutoff[i],"black","black"))))
dev.off()

############ SELECT UP and DOWN Genes : Log2FC >1 or <-1 ###########

logFC_fdr_table <- etp$table[order(etp$table$logFC,decreasing=TRUE),]
ID <- row.names(logFC_fdr_table)
logFC_fdr_table <- cbind(ID,logFC_fdr_table)


UP_genes_table <- logFC_fdr_table[logFC_fdr_table$FDR <= fdr_cutoff[i] & logFC_fdr_table$logFC >= log_cutoff[i],]
DOWN_genes_table <- logFC_fdr_table[logFC_fdr_table$FDR <= fdr_cutoff[i] & logFC_fdr_table$logFC <= (-1)*log_cutoff[i],]

#### WRITE TABLES ########
write.table(logFC_fdr_table,
            file=paste(outputdir,"/",compname, "_fdr_", fdr_cutoff[i], "_log2cutoff_",log_cutoff[i], "_expressed_genes_log2FC_table.txt", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE);

write.table(UP_genes_table[order(UP_genes_table$logFC,decreasing=TRUE),],
            file=paste(outputdir,"/",compname, "_fdr_", fdr_cutoff[i], "_log2cutoff_",log_cutoff[i], "_UP_genes_table.txt", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE)
              
  
write.table(DOWN_genes_table[order(DOWN_genes_table$logFC,decreasing=FALSE),],
              file=paste(outputdir,"/",compname, "_fdr_", fdr_cutoff[i], "_log2cutoff_",log_cutoff[i], "_DOWN_genes_table.txt", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE)

}

save.image(paste("Analysis_edgeR_",date,".RData",sep=""))


