library(dplyr)
library(Seurat)
packageVersion(pkg = "Seurat")
devtools::install_github("satijalab/seurat", ref = "v4.3.0")
library(patchwork)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

install.packages("Matrix")
install.packages("irlba")

packageVersion("Seurat")

## Load the dataset
hc1.data <- Read10X(data.dir = "./GSM6634568/")
hc2.data <- Read10X(data.dir = "./GSM6634570/")



## Initialize the Seurat object with the raw (non-normalized data).
hc1 <- CreateSeuratObject(counts = hc1.data, project = "hc1", min.cells = 3, min.features = 200)
hc2 <- CreateSeuratObject(counts = hc2.data, project = "hc2", min.cells = 3, min.features = 200)


###QC and selecting cells for further analysis
## The [[ operator can add columns to object metadata. This is a great place to stash QC stats
hc1[["percent.mt"]] <- PercentageFeatureSet(hc1, pattern = "^mt-")
hc2[["percent.mt"]] <- PercentageFeatureSet(hc2, pattern = "^mt-")

## Visualize QC metrics as a violin plot
options(repr.plot.width = 8, repr.plot.height = 5) #size change
VlnPlot(hc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1)
VlnPlot(hc2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(hc1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


##select cells
hc1 <- subset(hc1, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10)
hc2 <- subset(hc2, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt <10)


###integration 
###normalization
##Identification of highly variable features (feature selection)
cell <- c(hc1,hc2)
for (i in 1:length(cell)) {
  cell[[i]] <- NormalizeData(cell[[i]], verbose = F)
  cell[[i]] <- FindVariableFeatures(cell[[i]], selection.method = "vst", 
                                    nfeature = 2000, verbose = F)
}

cell.anchors <- FindIntegrationAnchors(object.list = cell, dims = 1:30)
cell.int <- IntegrateData(anchorset = cell.anchors, dims = 1:30) 

### Scaling and clustering
##Scaling the data
cell.int <- ScaleData(cell.int, verbose = T) 
#Perform linear dimensional reduction
cell.int <- RunPCA(cell.int, npcs = 30, verbose = T)

#Elbow Plot
ElbowPlot(cell.int)

##Cluster the cells
cell.int1<- FindNeighbors(cell.int, dims = 1:15)
cell.int1 <- FindClusters(cell.int1, resolution = 0.5)

#UMAP
cell.int1 <- RunUMAP(cell.int1, dims = 1:15)
DimPlot(cell.int1, reduction = "umap", label = TRUE)

# Visualization
p1 <- DimPlot(cell.int1, reduction = "umap", group.by="orig.ident")
p2 <- DimPlot(cell.int1, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


DimPlot(cell.int1, reduction = "umap", split.by = "orig.ident")


# Run SingleR to infer cell types of pbmc dataset using reference data
results <- SingleR(test = as.SingleCellExperiment(cell.int1), ref = ref, labels = ref$label.main, assay="integrated")

# Add inferred cell type labels to pbmc object
cell.int1$singlr_labels <- results$labels

# Visualize cell types in a UMAP plot with labels
DimPlot(cell.int1, reduction = 'umap', label = TRUE)
DimPlot(cell.int1, reduction = "umap", split.by = "orig.ident")

# Visualize expression of selected genes split by stimulation status
FeaturePlot(cell.int1, features = c("RORC","TBX21","GATA3","CXCR3","CCR6","CCR4","CCR7","CD247"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red"))
VlnPlot(cell.int1, features = c("TBX21","GATA3","RORC","CXCR3","CCR4","CCR6","CCR7","CD247"),  split.by = "orig.ident",stack = T,flip = T, assay="RNA" )

VlnPlot(cell.int1, features = c("RORC","TBX21","GATA3","CXCR3","CCR6","CCR4","CD247"),  split.by = "orig.ident", assay="RNA" )
VlnPlot(cell.int1, features = c("RORC","TBX21","GATA3","CXCR3","CCR6","CCR4","CD247"),  split.by = "orig.ident",stack = T,flip = T, assay="RNA" )

#Dotplot
DotPlot(cell.int1, features = c("RORC","TBX21","GATA3","CXCR3","CCR6","CCR4","CD247"),  split.by = "orig.ident",assay="RNA",
        cols = c("lightgray", "orange", "red", "darkred", "black"))



###ばらして表示

VlnPlot(cell.int1, features = c("RORC"),  split.by = "orig.ident", assay="RNA" )
VlnPlot(cell.int1, features = c("TBX21"),  split.by = "orig.ident", assay="RNA" )
VlnPlot(cell.int1, features = c("GATA3"),  split.by = "orig.ident", assay="RNA" )
VlnPlot(cell.int1, features = c("CXCR3"),  split.by = "orig.ident", assay="RNA" )
VlnPlot(cell.int1, features = c("CCR6"),  split.by = "orig.ident", assay="RNA" )
VlnPlot(cell.int1, features = c("CCR4"),  split.by = "orig.ident", assay="RNA" )
VlnPlot(cell.int1, features = c("CD247"),  split.by = "orig.ident", pt.size = 0.1, assay="RNA" )

VlnPlot(cell.int, features = c("CD247"),  split.by = "orig.ident", pt.size = 0.1, assay="RNA" )


FeaturePlot(cell.int1, features = c("RORC"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey95", "red"))
FeaturePlot(cell.int1, features = c("TBX21"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey95", "red"))
FeaturePlot(cell.int1, features = c("GATA3"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey95", "red"))
FeaturePlot(cell.int1, features = c("CXCR3"), split.by = "orig.ident", max.cutoff = 3, cols = c("grey95", "red"))
FeaturePlot(cell.int1, features = c("CCR6"), split.by = "orig.ident", max.cutoff = 3, cols = c("gray95", "red"))
FeaturePlot(cell.int1, features = c("CCR4"), split.by = "orig.ident", max.cutoff = 3, cols = c("gray95", "red"))
FeaturePlot(cell.int1, features = c("CCR7"), split.by = "orig.ident", max.cutoff = 3, cols = c("gray95", "red"))
FeaturePlot(cell.int1, features = c("CD247"), split.by = "orig.ident", max.cutoff = 3, cols = c("gray95", "red"))


FeaturePlot(cell.int1, features = c("CD247"), split.by = "orig.ident", max.cutoff = 3, cols = c("gray99","blue4"))
