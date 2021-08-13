# ST_snRNA-seq
Technical manuscript

Step1: 
a) Data pre-processing
b) data normalization
c) dimentionality reduction and clustering

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
CR_control_1[["percent.mt"]] <- PercentageFeatureSet(CR_control_1, pattern = "^MT-")
CR_control_2[["percent.mt"]] <- PercentageFeatureSet(CR_control_2, pattern = "^MT-")
CR_sample_1[["percent.mt"]] <- PercentageFeatureSet(CR_sample_1, pattern = "^MT-")
CR_sample_2[["percent.mt"]] <- PercentageFeatureSet(CR_sample_2, pattern = "^MT-")
HE_control_1[["percent.mt"]] <- PercentageFeatureSet(HE_control_1, pattern = "^MT-")
HE_sample_1[["percent.mt"]] <- PercentageFeatureSet(HE_sample_1, pattern = "^MT-")
HE_control_2[["percent.mt"]] <- PercentageFeatureSet(HE_control_2, pattern = "^MT-")
HE_sample_2[["percent.mt"]] <- PercentageFeatureSet(HE_sample_2, pattern = "^MT-")



##Filtering the metadata basd on various criteria
CR_control_2 <- subset(CR_control_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
CR_sample_2 <- subset(CR_sample_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
CR_control_1 <- subset(CR_control_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
CR_sample_1 <- subset(CR_sample_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15 )
HE_control_1 <- subset(HE_control_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_1 <- subset(HE_sample_1, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_control_2 <- subset(HE_control_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
HE_sample_2 <- subset(HE_sample_2, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)

##Data normalization
CR_control_1 <- NormalizeData(CR_control_1)
CR_sample_1 <- NormalizeData(CR_sample_1)
CR_control_2 <- NormalizeData(CR_control_2)
CR_sample_2 <- NormalizeData(CR_sample_2)
HE_control_1 <- NormalizeData(HE_control_1)
HE_sample_1 <- NormalizeData(HE_sample_1)
HE_control_2 <- NormalizeData(HE_control_2)
HE_sample_2 <- NormalizeData(HE_sample_2)


##Finidng variable features
CR_control_1 <- FindVariableFeatures(CR_control_1, selection.method = "vst")
CR_sample_1 <- FindVariableFeatures(CR_sample_1, selection.method = "vst")
CR_control_2 <- FindVariableFeatures(CR_control_2, selection.method = "vst")
CR_sample_2 <- FindVariableFeatures(CR_sample_2, selection.method = "vst")
HE_control_1 <- FindVariableFeatures(HE_control_1, selection.method = "vst")
HE_sample_1 <- FindVariableFeatures(HE_sample_1, selection.method = "vst")
HE_control_2 <- FindVariableFeatures(HE_control_2, selection.method = "vst")
HE_sample_2 <- FindVariableFeatures(HE_sample_2, selection.method = "vst")


##Scaling the data and regressing out the MT genes
HE_sample_2<-ScaleData(HE_sample_2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
HE_control_2<-ScaleData(HE_control_2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_sample_1 <- ScaleData(CR_sample_1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_control_1<-ScaleData(CR_control_1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
HE_sample_1<-ScaleData(HE_sample_1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
HE_control_1<-ScaleData(HE_control_1, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_sample_2 <- ScaleData(CR_sample_2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")
CR_control_2 <-ScaleData(CR_control_2, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")


## Preparing the required information for clustering and PCA/UMAP plots
HE_sample_2 <- RunPCA(HE_sample_2, verbose = FALSE)
HE_sample_2 <- RunUMAP(HE_sample_2, dims = 1:30, verbose = FALSE)
HE_sample_2 <- FindNeighbors(HE_sample_2, dims = 1:30, verbose = FALSE)
HE_sample_2 <- FindClusters(HE_sample_2, verbose = FALSE,resolution = 0.3) 

HE_control_2 <- RunPCA(HE_control_2, verbose = FALSE)
HE_control_2 <- RunUMAP(HE_control_2, dims = 1:30, verbose = FALSE)
HE_control_2 <- FindNeighbors(HE_control_2, dims = 1:30, verbose = FALSE)
HE_control_2 <- FindClusters(HE_control_2, verbose = FALSE, resolution = 0.3)

CR_sample_1 <- RunPCA(CR_sample_1, verbose = FALSE)
CR_sample_1 <- RunUMAP(CR_sample_1, dims = 1:30, verbose = FALSE)
CR_sample_1 <- FindNeighbors(CR_sample_1, dims = 1:30, verbose = FALSE)
CR_sample_1 <- FindClusters(CR_sample_1, verbose = FALSE, resolution = 0.3)

CR_control_1 <- RunPCA(CR_control_1, verbose = FALSE)
CR_control_1 <- RunUMAP(CR_control_1, dims = 1:30, verbose = FALSE)
CR_control_1 <- FindNeighbors(CR_control_1, dims = 1:30, verbose = FALSE)
CR_control_1 <- FindClusters(CR_control_1, verbose = FALSE, resolution = 0.3) 

HE_sample_1 <- RunPCA(HE_sample_1, verbose = FALSE)
HE_sample_1 <- RunUMAP(HE_sample_1, dims = 1:30, verbose = FALSE)
HE_sample_1 <- FindNeighbors(HE_sample_1, dims = 1:30, verbose = FALSE)
HE_sample_1 <- FindClusters(HE_sample_1, verbose = FALSE,resolution = 0.3) 

HE_control_1 <- RunPCA(HE_control_1, verbose = FALSE)
HE_control_1 <- RunUMAP(HE_control_1, dims = 1:30, verbose = FALSE)
HE_control_1 <- FindNeighbors(HE_control_1, dims = 1:30, verbose = FALSE)
HE_control_1 <- FindClusters(HE_control_1, verbose = FALSE, resolution = 0.3)

CR_sample_2 <- RunPCA(CR_sample_2, verbose = FALSE)
CR_sample_2 <- RunUMAP(CR_sample_2, dims = 1:30, verbose = FALSE)
CR_sample_2 <- FindNeighbors(CR_sample_2, dims = 1:30, verbose = FALSE)
CR_sample_2 <- FindClusters(CR_sample_2, verbose = FALSE, resolution = 0.3)

CR_control_2 <- RunPCA(CR_control_2, verbose = FALSE)
CR_control_2 <- RunUMAP(CR_control_2, dims = 1:30, verbose = FALSE)
CR_control_2 <- FindNeighbors(CR_control_2, dims = 1:30, verbose = FALSE)
CR_control_2 <- FindClusters(CR_control_2, verbose = FALSE, resolution = 0.3) 

##clusterS are shown prior to consequtive data integration
DimPlot(HE_sample_1, reduction = "umap")
DimPlot(HE_sample_2, reduction = "umap")
DimPlot(CR_sample_1, reduction = "umap")
DimPlot(CR_sample_2, reduction = "umap")
DimPlot(HE_control_1, reduction = "umap")
DimPlot(HE_control_2, reduction = "umap")
DimPlot(CR_control_1, reduction = "umap")
DimPlot(CR_control_2, reduction = "umap")


###Adding extra information to metadata 
HE_sample_1$orig.ident <- "HE_ON1"
CR_sample_1$orig.ident <- "CR_ON1"
HE_sample_2$orig.ident <- "HE_TN1"
CR_sample_2$orig.ident <- "CR_TN1"

HE_control_1$orig.ident <- "HE_ON2"
CR_control_1$orig.ident <- "CR_ON2"
HE_control_2$orig.ident <- "HE_TN2"
CR_control_2$orig.ident <- "CR_TN2"


HE_sample_1$staining <- "HE"
CR_sample_1$staining <- "CR"
HE_sample_2$staining <- "HE"
CR_sample_2$staining <- "CR"

HE_control_1$staining <- "HE"
CR_control_1$staining <- "CR"
HE_control_2$staining <- "HE"
CR_control_2$staining <- "CR"

###Integrating the consecutive slices
SP1<-list(HE_sample_1, CR_sample_1)
SP2<-list(HE_sample_2, CR_sample_2)
CT1<-list(HE_control_1, CR_control_1)
CT2<-list(HE_control_2, CR_control_2)


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

##clusterS are shown after consequtive data integration
DimPlot(SP1.integrated, reduction = "umap")
DimPlot(SP2.integrated, reduction = "umap")
DimPlot(CT1.integrated, reduction = "umap")
DimPlot(CT2.integrated, reduction = "umap")


##Saving the integrated data for further analysis

#setwd()  Changing the directory to were you want to save the integrated objects
saveRDS(SP1.integrated, "SP1.integrated.rds")
saveRDS(SP1.integrated, "SP1.integrated.rds")
saveRDS(SP1.integrated, "SP1.integrated.rds")
saveRDS(SP1.integrated, "SP1.integrated.rds")


```
