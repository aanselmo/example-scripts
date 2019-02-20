library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("App"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "bandwidth",
                  label = "Bandwidth",
                  min = 1,
                  max = 50,
                  value = 1),
      
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "topgen",
                  label = "TopGenes",
                  min = 1,
                  max = 31,
                  value = 30)
      
    ),
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)

server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
 
  output$distPlot <- renderPlot({
    
    ### FUNCTIONS
    plot_pca <- function (mat,groupnames,mat_name){
      library(FactoMineR)
      tRPKM <- as.data.frame(t(mat))
      type <- groupnames
      tRPKM <- cbind(tRPKM,type)
      
      type_col <- dim(mat)[1]+1
      type_colors <- c("blue","red") 
      
      #postscript(paste("Plots/PCA_",mat_name,"_1v2.ps",sep=""),paper="letter",onefile=FALSE)
      tRPKM.fmr.1v2.label.pca <- PCA(tRPKM, quali.sup=type_col, axes=c(1,2))
      plot(tRPKM.fmr.1v2.label.pca, choix="ind", habillage=type_col,
           col.hab=type_colors,
           cex=2.0,label="none",
           title = paste("PCA: RNAseq ",mat_name," dim 1v2",sep=), axes=c(1,2))
      #legend("bottomright", legend=unique(groupnames),
      #lty=c(1,2), lwd=c(1,1),
      #text.col=c("blue","red"),cex=1.5,
      #col=c("blue","red"), border="black") # gives the legend lines the correct color and width
      #dev.off()
      
      
      #postscript(paste("Plots/PCA_",mat_name,"_1v3.ps",sep=""),paper="letter",onefile=FALSE)
      tRPKM.fmr.1v3.label.pca <- PCA(tRPKM, quali.sup=type_col, axes=c(1,3))
      pca_1v3 <- tRPKM.fmr.1v3.label.pca
      plot(tRPKM.fmr.1v3.label.pca, choix="ind", habillage=type_col,
           col.hab=type_colors,
           cex=2.0,label="none",
           title = paste("PCA: RNAseq ",mat_name," dim 1v3",sep=), axes=c(1,3))
      #legend("bottomright", legend=unique(groupnames),
      #lty=c(1,2), lwd=c(1,1),
      #  text.col=c("blue","red"),cex=1.5,
      # col=c("blue","red"), border="black") # gives the legend lines the correct color and width
      #dev.off()
      
      return(tRPKM.fmr.1v2.label.pca)
    }
    
    
    
    ##########_READ_IN_DATA_TABLES_#######################################
    library(FactoMineR)
    
    pca_1v2 <- ""
    pca_1v3 <- ""
    outname <- "Disease"
    
    inputdir <- "INPUT"
    groupfile <- "groupings.txt"
    groupings <- read.table(paste(inputdir,"/",groupfile,sep=""), header=FALSE, quote="", comment.char="", stringsAsFactors=FALSE)
    colnames(groupings) <- c("samp","grp")
    
    #Disease_TPM_table = input_table
    input_table <- read.delim("DiseasevCTRL_Salmon_TPM_table.txt", stringsAsFactors=FALSE)
    SYMBOL <- toupper(unlist(lapply(strsplit(input_table$ID,"_"), function(x) x[2])))
    input_table <- cbind(SYMBOL,input_table)
    tpm_mat <- expr_data_mat <- as.matrix(input_table[,c(groupings$samp)])
    row.names(tpm_mat) <- input_table$ID
    tpm_table <- input_table
    
    sig2diff <- read.table("sig2_and_diff.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
    
    ############################
    
    t_mat <- tpm_mat[rowSums(tpm_mat > 1) >= 1, ]
    t_table <- tpm_table[rowSums(tpm_mat > 1) >= 1,]
    t_mat_sig2diff <- t_mat[t_table$SYMBOL %in% sig2diff$V1,]
    pca1v2 <- plot_pca(t_mat_sig2diff,groupings$grp,"Disease_sig2diff")
    
    
    
    ### PLOT
    ####
    topgenes <- row.names(pca1v2$var$contrib[order(pca1v2$var$contrib[,"Dim.1"],decreasing=TRUE),])[1:input$topgen]
    topcontrib <- pca1v2$var$contrib[order(pca1v2$var$contrib[,"Dim.1"],decreasing=TRUE),"Dim.1"][1:input$topgen]
    
    top_mat <- as.matrix(topcontrib)
    colnames(top_mat) <- "contrib"
    
    t_mat_top <- t_mat_sig2diff[row.names(top_mat),]
    avg_Disease <- rowSums(t_mat_top[,1:99])/99
    avg_ctrl <- rowSums(t_mat_top[,100:117])/18
    
    t_mat_top <- t_mat_top/avg_ctrl
    t_mult_top <- top_mat[,"contrib"] * t_mat_top /100

    Sample_Scores <- colSums(t_mult_top)/max(colSums(t_mult_top))*100
    
    
    #########
    
    
    #pdf(file="Sample_scores_test_top.pdf")
    genexp <- Sample_Scores[1:99]
    genexp2 <- Sample_Scores[100:117]
    genexpall <- Sample_Scores[1:117]
    xhigh <- max(genexpall)
    bwidth = input$bandwidth
    xlow <- 0
    #yhigh <- 0.1
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
    plot(density(genexp2,bw=bwidth), lwd = 2, col = "black", main = "Distribution of Disease Signature Scores", xlab = "", ylab = "", #ylim = c(0, yhigh), 
         xlim = c(xlow, round(xhigh+1)), axes = TRUE)
    
    #axis(1, seq(xlow, xhigh, by = .1))
    #axis(2, labels = FALSE, lwd.ticks = 0)
    rug(jitter(genexp),col="red")
    rug(jitter(genexp2))
    mtext("Disease Score", side = 1, line = 2.5, cex = 1.5, font = 2)
    mtext("Density", side = 2, line = 4, cex = 1.5, font = 2, las = 0)
    lines(density(genexp,bw=bwidth), col="red", lwd= 2)
    legend("topright",c("Disease","CTRL"),col=c("red","black"),lty=1,lwd=2)
    
    
  })
  
}

shinyApp(ui, server)
