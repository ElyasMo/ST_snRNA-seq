# ST_snRNA-seq
Technical manuscript

## Step1: 
**a)** Data pre-processing.
**b)** Data normalization.
**c)** Dimentionality reduction and clustering.
**d)** Consecutive slices data integration

#### **a)** Data pre-processing.
```r

library(Seurat)
 
#Importing the required paths for the seurat objects which includes .h5 and spatial directory
path1 <- "//P1_ON1_A/"
path2 <- "//P1_ON2_A/"
path3 <- "//P2_ON1_B/"
path4 <- "//P2_ON2_B/"
path5 <- "//P3_TN1_A/"
path6 <- "//P3_TN2_A/"
path7 <- "//P4_TN1_B/"
path8 <- "//P4_TN2_B/"
  
#Making a list of the paths
paths <- c(path1, path2, path3, path4, path5, path6, path7, path8)  

#Creating the required lists to be added to metadata
slices <- as.character(c("P1_ON1_A", "P1_ON2_A", "P2_ON1_B", "P2_ON2_B", "P3_TN1_A", "P3_TN2_A", "P4_TN1_B", "P4_TN2_B"))
staining <- as.character(c("HE","CR", "HE","CR", "HE","CR", "HE","CR"))
Region <- as.character(c("Orbitofrontal", "Orbitofrontal", "Orbitofrontal", "Orbitofrontal", "Temporal", "Temporal", "Temporal", "Temporal"))
Order <- as.character(c("First_Slice", "Second_Slice", "First_Slice", "Second_Slice", "First_Slice", "Second_Slice", "First_Slice", "Second_Slice"))

#Reading in the Seurat objects
Seurat_Objects <-mapply(function(X,Y){
  Load10X_Spatial(X,filename = "filtered_feature_bc_matrix.h5",
                                                      assay = "Spatial",
                                                      slice = Y,
                                                      filter.matrix = TRUE,
                                                      to.upper = FALSE)},X=paths,Y=slices, USE.NAMES = FALSE) 

#Adding extra information to metadata 
nums=1
num <- c(1:8)
for (nums in num){
  Seurat_Objects[[nums]][["orig.ident"]] <- slices[nums]
  Seurat_Objects[[nums]][["staining"]] <- staining[nums]
  Seurat_Objects[[nums]][["Region"]] <- Region[nums]
  Seurat_Objects[[nums]][["Order"]] <- Order[nums] 
}


##adding the percentage of MT genes to metadata
Seurat_Objects <- lapply(Seurat_Objects, function(X){
  X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-");X
  })

##Filtering the metadata based on various criteria
Seurat_Objects <- lapply(Seurat_Objects, function(X){
  X <- subset(X, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15)
  })
  
```
#### **b)** Data normalization.
```r
##Data normalization
Seurat_Objects <- lapply(Seurat_Objects, function(X){
  X <- NormalizeData(X)
  })



##Finidng variable features
nums=1
num <- c(1:8)
for (nums in num){
  Seurat_Objects[[nums]] <- FindVariableFeatures(Seurat_Objects[[nums]], selection.method = "vst")
}


##Scaling the data and regressing out the MT genes, 
Seurat_Objects <- lapply(Seurat_Objects, function(X){
  ScaleData(X, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt",do.scale = FALSE, do.center = FALSE)
  })

 
```
#### **c)** Dimentionality reduction and clustering.
```r

## Preparing the required information for clustering and PCA/UMAP plots

Seurat_Objects_single <- lapply(Seurat_Objects, function(X){
  X <- RunPCA(X, verbose = FALSE)
  X <- RunUMAP(X, dims = 1:30, verbose = FALSE)
  X <- FindNeighbors(X, dims = 1:30, verbose = FALSE)
  X <- FindClusters(X, verbose = FALSE,resolution = 0.3) 
})


##clusters are shown prior to consequtive data integration
for (obj in Seurat_Objects_single){
  print(SpatialDimPlot(obj,label = TRUE, label.size = 3))
}

```
#### **d)** Consecutive slices data integration
```r
###Integrating the consecutive slices
SP1 <- list(Seurat_Objects[[1]], Seurat_Objects[[2]])
CT1 <- list(Seurat_Objects[[3]], Seurat_Objects[[4]])
SP2 <- list(Seurat_Objects[[5]], Seurat_Objects[[6]])
CT2 <- list(Seurat_Objects[[7]], Seurat_Objects[[8]])


##DATA integration
#Finding anchors between consequtive slices. Scale is TRUE because as we mentioned
# for the purpose of data integration, we leave the data unscaled in single slices.
SP1.anchors <- FindIntegrationAnchors(object.list = SP1, dims = 1:30, 
                                      scale = TRUE, normalization.method = "LogNormalize")
SP2.anchors <- FindIntegrationAnchors(object.list = SP2, dims = 1:30, 
                                      scale = TRUE, normalization.method = "LogNormalize")
CT1.anchors <- FindIntegrationAnchors(object.list = CT1, dims = 1:30, 
                                      scale = TRUE, normalization.method = "LogNormalize")
CT2.anchors <- FindIntegrationAnchors(object.list = CT2, dims = 1:30, 
                                      scale = TRUE, normalization.method = "LogNormalize")


#Integratiing consequtive slices based on identified anchors
integrated_obj <- list(
SP1.integrated <- IntegrateData(anchorset = SP1.anchors, dims = 1:30,normalization.method = "LogNormalize"),
SP2.integrated <- IntegrateData(anchorset = SP2.anchors, dims = 1:30,normalization.method = "LogNormalize"),
CT1.integrated <- IntegrateData(anchorset = CT1.anchors, dims = 1:30,normalization.method = "LogNormalize"),
CT2.integrated <- IntegrateData(anchorset = CT2.anchors, dims = 1:30,normalization.method = "LogNormalize"))

lapply(integrated_obj, function(X){
  DefaultAssay(X) <- "integrated"
})

## Preparing the required information for clustering and PCA/UMAP plots

integrated_obj <- lapply(integrated_obj, function(X){
  X <- ScaleData(X, verbose = FALSE,do.scale = FALSE, do.center = FALSE)
  X <- RunPCA(X, verbose = FALSE)
  X <- RunUMAP(X, dims = 1:30, verbose = FALSE)
  X <- FindNeighbors(X, dims = 1:30, verbose = FALSE)
  X <- FindClusters(X, verbose = FALSE,resolution = 0.3) 
})


##clusterS are shown after consequtive data integration
for (obj in integrated_obj){
  print(SpatialDimPlot(obj,label = TRUE, label.size = 3))
}


##Saving the integrated data for further analysis

#setwd()  Changing the directory to were you want to save the integrated objects
setwd("//Processed_data/")
saveRDS(integrated_obj[[1]], "SP1.integrated.rds")
saveRDS(integrated_obj[[2]], "SP2.integrated.rds")
saveRDS(integrated_obj[[3]], "CT1.integrated.rds")
saveRDS(integrated_obj[[4]], "CT2.integrated.rds")

```
# Step2
**a)** snRNA-seq data pre-processing.
**b)** snRNA-seq label transfering

#### **a)** snRNA-seq data pre-processing.
```r

library(Seurat)
##Changing the directory to where the integrated data is saved
setwd("//Processed_data/")


SP1.integrated <- readRDS("SP1.integrated.rds")
SP2.integrated <- readRDS("SP2.integrated.rds")
CT1.integrated <- readRDS("CT1.integrated.rds")
CT2.integrated <- readRDS("CT2.integrated.rds")


##Loading the labeled single nuclei datasets from prefrontal cortex of human non-AD and AD postmorem brain samples
AD_Prefrontal <- readRDS("//rds/AD01104_Prefrontal cortex_SP.rds")
H_Prefrontal <- readRDS("//rds/AD01102_Prefrontal cortex_CT.rds")

#normalizing the single nuclei datasets and providing the PCA/UMAP infromation
library(dplyr)
AD_Prefrontal <- NormalizeData(AD_Prefrontal, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)
H_Prefrontal <- NormalizeData(H_Prefrontal, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)

# DefaultAssay(SP1.integrated) <- "integrated"
# DefaultAssay(SP2.integrated) <- "integrated"
# DefaultAssay(CT1.integrated) <- "integrated"
# DefaultAssay(CT2.integrated) <- "integrated"
```
#### **b)** snRNA-seq label transfering
```r
#The same workflow can be applied to single slices
#Finding the anchors between single nuclei dataset and integrated spatial transcriptomics data.
SP1.anchors <- FindTransferAnchors(reference = AD_Prefrontal, query = SP1.integrated, normalization.method = "LogNormalize")
SP2.anchors <- FindTransferAnchors(reference = AD_Prefrontal, query = SP2.integrated, normalization.method = "LogNormalize")
CT1.anchors <- FindTransferAnchors(reference = H_Prefrontal, query = CT1.integrated, normalization.method = "LogNormalize")
CT2.anchors <- FindTransferAnchors(reference = H_Prefrontal, query = CT2.integrated, normalization.method = "LogNormalize")



##Transfering the identified anchors and consequently cell labels from single nuclei data and spatial transcriptomics data
SP1.predictions.assay1 <- TransferData(anchorset = SP1.anchors, refdata = AD_Prefrontal$predicted.id, 
                                   prediction.assay = TRUE, weight.reduction = SP1.integrated[["pca"]], dims = 1:30)
SP1.predictions <- TransferData(anchorset = SP1.anchors, refdata = AD_Prefrontal$predicted.id, 
                            dims = 1:50, weight.reduction = SP1.integrated[["pca"]])

SP2.predictions.assay1 <- TransferData(anchorset = SP2.anchors, refdata = AD_Prefrontal$predicted.id, 
                                       prediction.assay = TRUE, weight.reduction = SP2.integrated[["pca"]], dims = 1:30)
SP2.predictions <- TransferData(anchorset = SP2.anchors, refdata = AD_Prefrontal$predicted.id, 
                                dims = 1:50, weight.reduction = SP2.integrated[["pca"]])

CT1.predictions.assay1 <- TransferData(anchorset = CT1.anchors, refdata = H_Prefrontal$predicted.id, 
                                       prediction.assay = TRUE, weight.reduction = CT1.integrated[["pca"]], dims = 1:30)
CT1.predictions <- TransferData(anchorset = CT1.anchors, refdata = H_Prefrontal$predicted.id, 
                                dims = 1:50, weight.reduction = CT1.integrated[["pca"]])

CT2.predictions.assay1 <- TransferData(anchorset = CT2.anchors, refdata = H_Prefrontal$predicted.id, 
                                       prediction.assay = TRUE, weight.reduction = CT2.integrated[["pca"]], dims = 1:30)
CT2.predictions <- TransferData(anchorset = CT2.anchors, refdata = H_Prefrontal$predicted.id, 
                                dims = 1:50, weight.reduction = CT2.integrated[["pca"]])

##Ading the information regarding the label transfering to medatada
SP1.integrated[["prediction.labels"]] <- SP1.predictions.assay1
SP1.integrated[["predictions"]] <- SP1.predictions$predicted.id
SP2.integrated[["prediction.labels"]] <- SP2.predictions.assay1
SP2.integrated[["predictions"]] <- SP2.predictions$predicted.id

CT1.integrated[["prediction.labels"]] <- CT1.predictions.assay1
CT1.integrated[["predictions"]] <- CT1.predictions$predicted.id
CT2.integrated[["prediction.labels"]] <- CT2.predictions.assay1
CT2.integrated[["predictions"]] <- CT2.predictions$predicted.id

##Clusters befor label transfering
integrated_obj <- list(SP1.integrated,SP2.integrated,CT1.integrated,CT2.integrated)

for (objs in integrated_obj){
  print(SpatialDimPlot(objs))
}

for (objs in integrated_obj){
  print(DimPlot(objs, reduction = "umap"))
}

##Changing the clusters from automatic clustering to label transfering
num <- c(1:4)
for (nums in num){
  integrated_obj[[nums]] <- SetIdent(integrated_obj[[nums]], value = integrated_obj[[nums]]$predictions)
}

for (objs in integrated_obj){
  print(SpatialDimPlot(objs))
}

for (objs in integrated_obj){
  print(DimPlot(objs, reduction = "umap"))
}


##Ploting the transfered labels
for (objs in integrated_obj){
  print(SpatialFeaturePlot(objs, 
                           features = c("Astrocytes", "Oligodendrocytes", 
                                        "Excitatory neurons"),
                           pt.size.factor = 1.6, ncol = 2, crop = TRUE))
}

setwd("/mnt/data1/elyas/10X_Visuim/data/rds/")

saveRDS(SP1.integrated, "SP1.labels.rds")
saveRDS(SP2.integrated, "SP2.labels.rds")
saveRDS(CT1.integrated, "CT1.labels.rds")
saveRDS(CT2.integrated, "CT2.labels.rds")
```
