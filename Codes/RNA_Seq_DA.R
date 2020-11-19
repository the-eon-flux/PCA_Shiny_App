'
Input : Reads count matrix. Each Row is gene & Column is a Sample
'

# Read Data
Count_Data <- as.data.frame(read.table("./Data/Counts_Matrix.txt", header = T, sep = ","))
row.names(Count_Data) <- Count_Data[,1]
Count_Data <- Count_Data[,-1]

MetaData <- read.csv("./Data/sample_Groups.txt", header = T)
row.names(MetaData) <- MetaData[,1]
MetaData <- MetaData[,-1] ; MetaData <- MetaData[order(row.names(MetaData)),]

MetaData$Activation_State <- as.factor(MetaData$Activation_State)
MetaData$Treatment <- as.factor(MetaData$Treatment)
MetaData$Time_hrs <- as.factor(MetaData$Time_hrs)

# Rearranging the sample order in Metadata file according to the Count_Data matrix

if(!all(rownames(MetaData) == colnames(Count_Data))) {
      if(all(rownames(MetaData) %in% colnames(Count_Data))){
            Count_Data <- Count_Data[,row.names(MetaData)]
            # Check
            all(rownames(MetaData) == colnames(Count_Data))
      }      
      
}


# Normalise : The DESeq2 model internally corrects for library size
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = Count_Data,
                              colData = MetaData,
                              design = ~ Treatment)
# Count data : assay(dds)
# Metadata info : colData(dds)

# minimal pre-filtering to keep only rows that have at least 10 reads total.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Normalisation by DESeq2::rlog() tranformation 
'PCA works best for data that generally has the same range of variance at different ranges of the mean values. When the expected amount of variance is approximately the same across different mean values, the data is said to be homoskedastic. For RNA-seq counts, however, the expected variance grows with the mean. For example, if one performs PCA directly on a matrix of counts or normalized counts (e.g. correcting for differences in sequencing depth), the resulting plot typically depends mostly on the genes with highest counts because they show the largest absolute differences between samples. '

library(vsn)

# Some pre-normalization diagnostic plots to show mean v/s variance
#meanSdPlot(assay(dds), ranks = FALSE)

' We can see that the variance shoots up with mean. There are a few high count genes also. Taking Log can help with that'

#Count_Data.log <- log2(assay(dds) + 1)
#meanSdPlot(Count_Data.log, ranks = FALSE) 

' Genes with small count amplifies differences when the values are close to 0. The low count genes with low signal-to-noise ratio will overly contribute to sample-sample distances and PCA plots. Hence we do some normalisations offered in DESeq2
'
# Plot post Rlog() normalization
rld <- rlog(dds, blind = FALSE)
#meanSdPlot(assay(rld), ranks = FALSE)

'blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).'

# Assess & visualize overall similarity between samples
sampleDists <- dist(t(assay(rld)))
library("pheatmap")
library(ggplot2)

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Treatment, rld$Time_hrs, sep = " - " )
colnames(sampleDistMatrix) <- rld$Activation_State

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, main = "Heatmap of sample-to-sample distances using the variance stabilizing transformed values.")

## Dimension reduction using PCA

# Score plot Direct
plotPCA(rld, intgroup = c("Treatment","Time_hrs"))

# Manual
df_pca <- prcomp(assay(rld),center = TRUE,scale = TRUE)

# $x for row-wise
gg <- ggplot(data.frame(df_pca$x), aes(x=PC1, y=PC2)) + geom_point()
gg

# $rotation for sample-wise
df_out_r <- as.data.frame(df_pca$rotation)

# Calculate Variance of each PCs
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out_r), "(", paste( as.character(percentage), "%", ")", sep="") )


df_out_r$feature <- paste( rld$Treatment, rld$Time_hrs, sep = " - " )
p<-ggplot(df_out_r,aes(x=PC1,y=PC2,color=feature ))
p+geom_point(size=3, pch=15) +xlab(percentage[1]) + ylab(percentage[2])

' Our PCA plot
pdf("./1.pdf")
plotPCA(rld, intgroup = c("Treatment","Time_hrs"))
p+geom_point(size=3, pch=18) +xlab(percentage[1]) + ylab(percentage[2])
dev.off()
'

# Scree plot

Var <- ((round(df_pca$sdev / sum(df_pca$sdev),2)))*100

plot(1:length(percentage), Var, xlab="Principal Components", ylab="% Variance Explained", pch=16, type="o", col="red", main = "Scree Plot", ylim=c(1,100))
lines(1:length(percentage), cumsum(Var), type = "o", pch=16, col="dodgerblue")
text(x=rep(10, 2), y=c(mean(Var), mean(cumsum(Var))), pos=4, labels = c("Variance", "Cumulative Variance"), col = c("red", "dodgerblue"))














