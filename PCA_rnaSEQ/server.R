#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

shinyServer(function(input, output) {
    
    ###############################################
    # Read Data
    Count_Data <- as.data.frame(read.table("../Data/Counts_Matrix.txt", header = T, sep = ","))
    row.names(Count_Data) <- Count_Data[,1]
    Count_Data <- Count_Data[,-1]
    
    ###############################################
    # Read Metadata
    
    MetaData <- read.csv("../Data/sample_Groups.txt", header = T)
    row.names(MetaData) <- MetaData[,1]
    MetaData <- MetaData[,-1] ; MetaData <- MetaData[order(row.names(MetaData)),]
    
    MetaData$Activation_State <- as.factor(MetaData$Activation_State)
    MetaData$Treatment <- as.factor(MetaData$Treatment)
    MetaData$Time_hrs <- as.factor(MetaData$Time_hrs)
    
    ###############################################
    # Rearranging the sample order in Metadata file according to the Count_Data matrix
    
    if(!all(rownames(MetaData) == colnames(Count_Data))) {
        if(all(rownames(MetaData) %in% colnames(Count_Data))){
            Count_Data <- Count_Data[,row.names(MetaData)]
            # Check
            all(rownames(MetaData) == colnames(Count_Data))
        }      
        
    }
    
    ###############################################
    # Normalize : The DESeq2 has functions to normalize, hence we use DESeq2
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = Count_Data,
                                  colData = MetaData,
                                  design = ~ Treatment)
    # Minimal pre-filtering to keep only rows that have at least 10 reads total.
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    # Doing rlog normalization offered in DESeq2
    rld <- rlog(dds, blind = FALSE)
    
    ###############################################
    library("pheatmap")
    
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix( sampleDists )
    rownames(sampleDistMatrix) <- paste( rld$Treatment, rld$Time_hrs, sep = " - " )
    colnames(sampleDistMatrix) <- rld$Activation_State
    
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists, main = "Heatmap of sample-to-sample distances using the variance stabilizing transformed values.")
    
    ###############################################
    
    Color_by <- reactive({
        col <- input$Col_grps
    })
    
    ###############################################
    output$Score_Plot <- renderPlot({
        
    })
    
    
    ###############################################
    
    output$Heatmap <- renderPlot({
        
    })
    
    ###############################################
    
    output$Scree_Plot <- renderPlot({
        
    })
    
    ###############################################
    
})
