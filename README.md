# ST_snRNA-seq
Technical manuscript

## Step1: 
**a)** Data pre-processing.
**b)** Data normalization.
**c)** Dimentionality reduction and clustering.

```r

library(Seurat)
#####Importing the spatial data processed by spaceranger (spaceranger mkfastq/ spaceranger count)
## the slices were not integrated by spaceranger (spaceranger agr). 

setwd("Directory")
Directory="..."
CR_TN1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "CR_TN1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)
setwd("Directory")
Directory="..."
CR_TN2<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "CR_TN2",
                              filter.matrix = TRUE,
                              to.upper = FALSE)
setwd("Directory")
Directory="..."
HE_TN1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "HE_TN1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)
setwd("Directory")
Directory="..."
HE_TN2<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "HE_TN2",
                              filter.matrix = TRUE,
                              to.upper = FALSE)
setwd("Directory")
Directory="..."
CR_ON1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "CR_ON1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)
setwd("Directory")
Directory="..."
CR_ON2<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "CR_ON2",
                              filter.matrix = TRUE,
                              to.upper = FALSE)
setwd("Directory")
Directory="..."
HE_ON1<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "HE_ON1",
                             filter.matrix = TRUE,
                             to.upper = FALSE)
setwd("Directory")
Directory="..."
HE_ON2<-Load10X_Spatial(Directory,filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "HE_ON2",
                              filter.matrix = TRUE,
                              to.upper = FALSE)



##adding the percentage of MT genes to metadata
CR_ON2[["percent.mt"]] <- PercentageFeatureSet(CR_ON2, pattern = "^MT-")
CR_TN2[["percent.mt"]] <- PercentageFeatureSet(CR_TN2, pattern = "^MT-")
CR_ON1[["percent.mt"]] <- PercentageFeatureSet(CR_ON1, pattern = "^MT-")
CR_TN1[["percent.mt"]] <- PercentageFeatureSet(CR_TN1, pattern = "^MT-")
HE_ON2[["percent.mt"]] <- PercentageFeatureSet(HE_ON2, pattern = "^MT-")
HE_TN2[["percent.mt"]] <- PercentageFeatureSet(HE_TN2, pattern = "^MT-")
HE_ON1[["percent.mt"]] <- PercentageFeatureSet(HE_ON1, pattern = "^MT-")
HE_TN1[["percent.mt"]] <- PercentageFeatureSet(HE_TN1, pattern = "^MT-")



##Filtering the metadata basd on various criteria
CR_TN2 <- subset(CR_TN2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
CR_TN1 <- subset(CR_TN1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
CR_ON2 <- subset(CR_ON2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
CR_ON1 <- subset(CR_ON1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
HE_ON2 <- subset(HE_ON2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_ON1 <- subset(HE_ON1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_TN2 <- subset(HE_TN2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_TN1 <- subset(HE_TN1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)

##Data normalization
CR_ON2 <- NormalizeData(CR_ON2)
CR_ON1 <- NormalizeData(CR_ON1)
CR_TN2 <- NormalizeData(CR_TN2)
CR_TN1 <- NormalizeData(CR_TN1)
HE_ON2 <- NormalizeData(HE_ON2)
HE_ON1 <- NormalizeData(HE_ON1)
HE_TN2 <- NormalizeData(HE_TN2)
HE_TN1 <- NormalizeData(HE_TN1)


##Finidng variable features
CR_ON2 <- FindVariableFeatures(CR_ON2, selection.method = "vst")
CR_ON1 <- FindVariableFeatures(CR_ON1, selection.method = "vst")
CR_TN2 <- FindVariableFeatures(CR_TN2, selection.method = "vst")
CR_TN1 <- FindVariableFeatures(CR_TN1, selection.method = "vst")
HE_ON2 <- FindVariableFeatures(HE_ON2, selection.method = "vst")
HE_ON1 <- FindVariableFeatures(HE_ON1, selection.method = "vst")
HE_TN2 <- FindVariableFeatures(HE_TN2, selection.method = "vst")
HE_TN1 <- FindVariableFeatures(HE_TN1, selection.method = "vst")


##Scaling the data and regressing out the MT genes
HE_TN1<-ScaleData(HE_TN1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
HE_TN2<-ScaleData(HE_TN2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_ON1 <- ScaleData(CR_ON1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_ON2<-ScaleData(CR_ON2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
HE_ON1<-ScaleData(HE_ON1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
HE_ON2<-ScaleData(HE_ON2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_TN1 <- ScaleData(CR_TN1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_TN2 <-ScaleData(CR_TN2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")


## Preparing the required information for clustering and PCA/UMAP plots
HE_TN1 <- RunPCA(HE_TN1, verbose = FALSE)
HE_TN1 <- RunUMAP(HE_TN1, dims = 1:30, verbose = FALSE)
HE_TN1 <- FindNeighbors(HE_TN1, dims = 1:30, verbose = FALSE)
HE_TN1 <- FindClusters(HE_TN1, verbose = FALSE,resolution = 0.3) 

HE_TN2 <- RunPCA(HE_TN2, verbose = FALSE)
HE_TN2 <- RunUMAP(HE_TN2, dims = 1:30, verbose = FALSE)
HE_TN2 <- FindNeighbors(HE_TN2, dims = 1:30, verbose = FALSE)
HE_TN2 <- FindClusters(HE_TN2, verbose = FALSE, resolution = 0.3)

CR_ON1 <- RunPCA(CR_ON1, verbose = FALSE)
CR_ON1 <- RunUMAP(CR_ON1, dims = 1:30, verbose = FALSE)
CR_ON1 <- FindNeighbors(CR_ON1, dims = 1:30, verbose = FALSE)
CR_ON1 <- FindClusters(CR_ON1, verbose = FALSE, resolution = 0.3)

CR_ON2 <- RunPCA(CR_ON2, verbose = FALSE)
CR_ON2 <- RunUMAP(CR_ON2, dims = 1:30, verbose = FALSE)
CR_ON2 <- FindNeighbors(CR_ON2, dims = 1:30, verbose = FALSE)
CR_ON2 <- FindClusters(CR_ON2, verbose = FALSE, resolution = 0.3) 

HE_ON1 <- RunPCA(HE_ON1, verbose = FALSE)
HE_ON1 <- RunUMAP(HE_ON1, dims = 1:30, verbose = FALSE)
HE_ON1 <- FindNeighbors(HE_ON1, dims = 1:30, verbose = FALSE)
HE_ON1 <- FindClusters(HE_ON1, verbose = FALSE,resolution = 0.3) 

HE_ON2 <- RunPCA(HE_ON2, verbose = FALSE)
HE_ON2 <- RunUMAP(HE_ON2, dims = 1:30, verbose = FALSE)
HE_ON2 <- FindNeighbors(HE_ON2, dims = 1:30, verbose = FALSE)
HE_ON2 <- FindClusters(HE_ON2, verbose = FALSE, resolution = 0.3)

CR_TN1 <- RunPCA(CR_TN1, verbose = FALSE)
CR_TN1 <- RunUMAP(CR_TN1, dims = 1:30, verbose = FALSE)
CR_TN1 <- FindNeighbors(CR_TN1, dims = 1:30, verbose = FALSE)
CR_TN1 <- FindClusters(CR_TN1, verbose = FALSE, resolution = 0.3)

CR_TN2 <- RunPCA(CR_TN2, verbose = FALSE)
CR_TN2 <- RunUMAP(CR_TN2, dims = 1:30, verbose = FALSE)
CR_TN2 <- FindNeighbors(CR_TN2, dims = 1:30, verbose = FALSE)
CR_TN2 <- FindClusters(CR_TN2, verbose = FALSE, resolution = 0.3) 

##clusters are shown prior to consequtive data integration
DimPlot(HE_ON1, reduction = "umap")
DimPlot(HE_TN1, reduction = "umap")
DimPlot(CR_ON1, reduction = "umap")
DimPlot(CR_TN1, reduction = "umap")
DimPlot(HE_ON2, reduction = "umap")
DimPlot(HE_TN2, reduction = "umap")
DimPlot(CR_ON2, reduction = "umap")
DimPlot(CR_TN2, reduction = "umap")


###Adding extra information to metadata 
HE_ON1$orig.ident <- "HE_ON1"
CR_ON1$orig.ident <- "CR_ON1"
HE_TN1$orig.ident <- "HE_TN1"
CR_TN1$orig.ident <- "CR_TN1"

HE_ON2$orig.ident <- "HE_ON2"
CR_ON2$orig.ident <- "CR_ON2"
HE_TN2$orig.ident <- "HE_TN2"
CR_TN2$orig.ident <- "CR_TN2"


HE_ON1$staining <- "HE"
CR_ON1$staining <- "CR"
HE_TN1$staining <- "HE"
CR_TN1$staining <- "CR"

HE_ON2$staining <- "HE"
CR_ON2$staining <- "CR"
HE_TN2$staining <- "HE"
CR_TN2$staining <- "CR"

###Integrating the consecutive slices
SP1<-list(HE_ON1, CR_ON1)
SP2<-list(HE_TN1, CR_TN1)
CT1<-list(HE_ON2, CR_ON2)
CT2<-list(HE_TN2, CR_TN2)


##DATA integration
#Finding anchors between consequtive slices
SP1.anchors <- FindIntegrationAnchors(object.list = SP1, dims = 1:30)
SP2.anchors <- FindIntegrationAnchors(object.list = SP2, dims = 1:30)
CT1.anchors <- FindIntegrationAnchors(object.list = CT1, dims = 1:30)
CT1.anchors <- FindIntegrationAnchors(object.list = CT2, dims = 1:30)

#Integratiing consequtive slices based on identified anchors
SP1.integrated <- IntegrateData(anchorset = SP1.anchors, dims = 1:30)
SP2.integrated <- IntegrateData(anchorset = SP2.anchors, dims = 1:30)
CT1.integrated <- IntegrateData(anchorset = CT1.anchors, dims = 1:30)
CT2.integrated <- IntegrateData(anchorset = CT2.anchors, dims = 1:30)


DefaultAssay(SP1.integrated) <- "integrated"
DefaultAssay(SP2.integrated) <- "integrated"
DefaultAssay(CT1.integrated) <- "integrated"
DefaultAssay(CT2.integrated) <- "integrated"


## Preparing the required information for clustering and PCA/UMAP plots
SP1.integrated <- ScaleData(SP1.integrated, verbose = FALSE)
SP1.integrated <- RunPCA(SP1.integrated, verbose = FALSE)
SP1.integrated <- RunUMAP(SP1.integrated, dims = 1:30, verbose = FALSE)
SP1.integrated <- FindNeighbors(SP1.integrated, dims = 1:30, verbose = FALSE)
SP1.integrated <- FindClusters(SP1.integrated, verbose = FALSE, resolution = 0.3)

SP2.integrated <- ScaleData(SP2.integrated, verbose = FALSE)
SP2.integrated <- RunPCA(SP2.integrated, verbose = FALSE)
SP2.integrated <- RunUMAP(SP2.integrated, dims = 1:30, verbose = FALSE)
SP2.integrated <- FindNeighbors(SP2.integrated, dims = 1:30, verbose = FALSE)
SP2.integrated <- FindClusters(SP2.integrated, verbose = FALSE, resolution = 0.3)

CT1.integrated <- ScaleData(CT1.integrated, verbose = FALSE)
CT1.integrated <- RunPCA(CT1.integrated, verbose = FALSE)
CT1.integrated <- RunUMAP(CT1.integrated, dims = 1:30, verbose = FALSE)
CT1.integrated <- FindNeighbors(CT1.integrated, dims = 1:30, verbose = FALSE)
CT1.integrated <- FindClusters(CT1.integrated, verbose = FALSE, resolution = 0.3)

CT2.integrated <- ScaleData(CT2.integrated, verbose = FALSE)
CT2.integrated <- RunPCA(CT2.integrated, verbose = FALSE)
CT2.integrated <- RunUMAP(CT2.integrated, dims = 1:30, verbose = FALSE)
CT2.integrated <- FindNeighbors(CT2.integrated, dims = 1:30, verbose = FALSE)
CT2.integrated <- FindClusters(CT2.integrated, verbose = FALSE, resolution = 0.3)

##cluster visualization
DimPlot(SP1.integrated, reduction = "umap")
DimPlot(SP2.integrated, reduction = "umap")
DimPlot(CT1.integrated, reduction = "umap")
DimPlot(CT2.integrated, reduction = "umap")


##Saving the integrated data for further analysis

#setwd()  Changing the directory to were you want to save the integrated objects
saveRDS(SP1.integrated, "SP1.integrated.rds")
saveRDS(SP1.integrated, "SP1.integrated.rds")
saveRDS(CT1.integrated, "CT1.integrated.rds")
saveRDS(CT2.integrated, "CT2.integrated.rds")


```
