library("optparse")

get_args <- function(){
option_list = list(
  make_option(c("-i", "--inputdir"), type="character", default=NULL, 
              help="data directory of count files", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default="OUTPUT", 
              help="output directory [default= %default]", metavar="character"),
  make_option(c("-g", "--groupfile"), type="character", default=NULL, 
              help="tab-delimited sample groupings .txt file; col1=samplename, col2=groupname", metavar="character"),
  make_option(c("-c", "--compfile"), type="character", default=NULL, 
              help="tab-delimited textfile of sample groups (col1=group1,col2=group2) to compare at \nspecified fdr (col3) and logFC (col4) cutoffs", metavar="character"),
  make_option(c("-p", "--plot"), type="integer", default=0, 
              help="plot 0 for no, 1 for yes [default=%default]", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$inputdir)){
  print_help(opt_parser)
  stop("An input directory containing count files MUST be specified.", call.=FALSE)
} 
if (is.null(opt$outputdir)) {
  print_help(opt_parser)
  stop("An output directory containing count files MUST be specified. Otherwise defaults to OUTPUT.", call.=FALSE)
}

if (is.null(opt$groupfile)) {
  print_help(opt_parser)
  stop("A tab-delimited file specifying sample groupings (col1 = samplename; col2=groupname) MUST be specified.", call.=FALSE)
}

if (is.null(opt$compfile)) {
  print_help(opt_parser)
  stop("A tab-delimited textfile of sample groups (cols 1&2) to compare at chosen fdr (col3) and logFC (col 4) cutoffs MUST be specified.", call.=FALSE)
}

print(paste("Processing ",opt$inputdir," as INPUT directory of count files",sep=""))
print(paste("Processing ",opt$outputdir," as OUTPUT directory",sep=""))

return(opt)

}
