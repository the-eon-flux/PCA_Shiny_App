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
    # Dimension reduction using PCA
    df_pca <- prcomp(assay(rld),center = TRUE,scale = TRUE)
    
    # $rotation matrix for sample-wise plot
    df_out_r <- as.data.frame(df_pca$rotation)
    
    ###############################################
    
    Color_by <- reactive({
        col <- as.character(input$Col_grps)
    })

    ###############################################
    
    Plot_By <- reactive({
        xPC <- as.character(input$xPC)
        yPC <- as.character(input$yPC)
        PCs <- c(xPC, yPC)
    })
    
    ############################################### Plot 1
    output$Score_Plot <- renderPlot({
        
        # Calculate % variance for each component
        percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
        percentage <- paste( colnames(df_out_r), "(", paste( as.character(percentage), "%", ")", sep="") )
        
        # Col Var Input
        Col <- Color_by()
        PCs <- Plot_By()
        
        #df_out_r$feature <- rld$Treatment
        X <- PCs[1]
        xI <- as.numeric(unlist(strsplit(X,"PC"))[2])
        
        Y <- PCs[2]
        yI <- as.numeric(unlist(strsplit(Y,"PC"))[2])
        
        plot_df <- data.frame(X = df_out_r[,PCs[1]], Y = df_out_r[,PCs[2]] )
        plot_df$Legend <- colData(rld)[,Col]
        
        p<-ggplot(plot_df,aes(x=X,y=Y,color=Legend ))
        p+geom_point(size=3, pch=15) +xlab(percentage[xI]) + ylab(percentage[yI])
    })
    
    
    ############################################### Plot 2
    
    output$Heatmap <- renderPlot({
        library("pheatmap")
        
        sampleDists <- dist(t(assay(rld)))
        sampleDistMatrix <- as.matrix( sampleDists )
        rownames(sampleDistMatrix) <- paste( rld$Treatment, rld$Time_hrs, sep = " - " )
        colnames(sampleDistMatrix) <- rld$Activation_State
        
        #pheatmap : not rendering without any error
        'pheatmap(sampleDistMatrix,
                 clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists, main = "Heatmap of sample-to-sample distances using the variance stabilizing transformed values.")
        '
        
        
        heatmap(sampleDistMatrix, scale = "none")
    })
    
    ############################################### Plot 3
    
    output$Scree_Plot <- renderPlot({
        
        # Scree plot
        
        Var <- ((round(df_pca$sdev / sum(df_pca$sdev),2)))*100
        
        plot(1:length(percentage), Var, xlab="Principal Components", ylab="% Variance Explained", pch=16, type="o", col="red", main = "Scree Plot", ylim=c(1,100))
        lines(1:length(percentage), cumsum(Var), type = "o", pch=16, col="dodgerblue")
        text(x=rep(10, 2), y=c(mean(Var), mean(cumsum(Var))), pos=4, labels = c("Variance", "Cumulative Variance"), col = c("red", "dodgerblue"))
        
    })
    
})
