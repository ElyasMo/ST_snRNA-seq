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

library(Seurat)
 

#Importing the required paths for the seurat objects which includes .h5 and spatial directory
path1 <- "D:/Poland/PHD/spatial/Processed_data//P1_ON1_A/"
path2 <- "D:/Poland/PHD/spatial/Processed_data//P1_ON2_A/"
path3 <- "D:/Poland/PHD/spatial/Processed_data//P2_ON1_B/"
path4 <- "D:/Poland/PHD/spatial/Processed_data//P2_ON2_B/"
path5 <- "D:/Poland/PHD/spatial/Processed_data//P3_TN1_A/"
path6 <- "D:/Poland/PHD/spatial/Processed_data//P3_TN2_A/"
path7 <- "D:/Poland/PHD/spatial/Processed_data//P4_TN1_B/"
path8 <- "D:/Poland/PHD/spatial/Processed_data//P4_TN2_B/"
  
#Making a list of the paths
paths <- c(path1,path2,path3,path4,path5,path6,path7,path8)  

#Creating the required lists to be added to metadata
slices <- as.character(c("P1_ON1_A","P1_ON2_A", "P2_ON1_B", "P2_ON2_B", "P3_TN1_A", "P3_TN2_A", "P4_TN1_B", "P4_TN2_B"))
staining <- as.character(c("HE","CR", "HE","CR", "HE","CR", "HE","CR"))
Region <- as.character(c("Orbitofrontal","Orbitofrontal","Orbitofrontal","Orbitofrontal","Temporal","Temporal","Temporal","Temporal"))
Order <- as.character(c("First_Slice", "Second_Slice","First_Slice", "Second_Slice","First_Slice", "Second_Slice","First_Slice", "Second_Slice"))

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
  X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-");X})

##Filtering the metadata based on various criteria
Seurat_Objects <- lapply(Seurat_Objects, function(X){
  subset(X, subset = nFeature_Spatial > 200 & nFeature_Spatial < 7000 & percent.mt < 15);X})
```
#### **b)** Data normalization.
```r

##Data normalization
Seurat_Objects <- lapply(Seurat_Objects, function(X){
  NormalizeData(X);X})



##Finidng variable features
nums=1
num <- c(1:8)
for (nums in num){
  Seurat_Objects[[nums]] <-  FindVariableFeatures(Seurat_Objects[[nums]], selection.method = "vst")
}


##Scaling the data and regressing out the MT genes, 
# If you want to use the information from single slices for consecutive slices data integration,
# it is better to avoid scaling the data at this point and do this during the data integration
# To do so, connsider do.scale=FALSE and do.center=FALSE in the ScaleData function
Seurat_Objects <- lapply(Seurat_Objects, function(X){
  ScaleData(X, assay = "Spatial", verbose = FALSE, vars.to.regress = "percent.mt")})
  
```
#### **c)** Dimentionality reduction and clustering.
```r

## Preparing the required information for clustering and PCA/UMAP plots

Seurat_Objects <- lapply(Seurat_Objects, function(X){
  X <- RunPCA(X, verbose = FALSE)
  X <- RunUMAP(X, dims = 1:30, verbose = FALSE)
  X <- FindNeighbors(X, dims = 1:30, verbose = FALSE)
  X <- FindClusters(X, verbose = FALSE,resolution = 0.3) 
  
})


##clusters are shown prior to consequtive data integration
for (obj in Seurat_Objects){
  print(SpatialDimPlot(obj,label = TRUE, label.size = 3))
}

```
#### **d)** Consecutive slices data integration
```r

###Integrating the consecutive slices
SP1<-list(Seurat_Objects[[1]], Seurat_Objects[[2]])
CT1<-list(Seurat_Objects[[3]], Seurat_Objects[[4]])
SP2<-list(Seurat_Objects[[5]], Seurat_Objects[[6]])
CT2<-list(Seurat_Objects[[7]], Seurat_Objects[[8]])


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
SP1.integrated <- IntegrateData(anchorset = SP1.anchors, dims = 1:30),
SP2.integrated <- IntegrateData(anchorset = SP2.anchors, dims = 1:30),
CT1.integrated <- IntegrateData(anchorset = CT1.anchors, dims = 1:30),
CT2.integrated <- IntegrateData(anchorset = CT2.anchors, dims = 1:30))

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
saveRDS(integrated_obj[[1]], "SP1.integrated.rds")
saveRDS(integrated_obj[[2]], "CT1.integrated.rds")
saveRDS(integrated_obj[[3]], "SP2.integrated.rds")
saveRDS(integrated_obj[[4]], "CT2.integrated.rds")
```
# Step2
**a)** snRNA-seq data pre-processing.
**b)** snRNA-seq label transfering

#### **a)** snRNA-seq data pre-processing.
```r
library(Seurat)

##Changing the directory to where the integrated data is saved
SP1.integrated <- readRDS("SP1.integrated.rds")
SP2.integrated <- readRDS("SP2.integrated.rds")
CT1.integrated <- readRDS("CT1.integrated.rds")
CT2.integrated <- readRDS("CT2.integrated.rds")

##Loading the labeled single nuclei datasets from prefrontal cortex of human non-AD and AD postmorem brain samples
AD_Prefrontal <- readRDS(".../AD01104_Prefrontal cortex_SP.rds")
H_Prefrontal <- readRDS(".../AD01102_Prefrontal cortex_CT.rds")

#normalizing the single nuclei datasets and providing the PCA/UMAP infromation
library(dplyr)
AD_Prefrontal <- NormalizeData(AD_Prefrontal, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)
H_Prefrontal <- NormalizeData(H_Prefrontal, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)

DefaultAssay(SP1.integrated) <- "integrated"
DefaultAssay(SP2.integrated) <- "integrated"
DefaultAssay(CT1.integrated) <- "integrated"
DefaultAssay(CT2.integrated) <- "integrated"
```

#### **b)** snRNA-seq label transfering
```r
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


##Changing the clusters from automatic clustering to label transfering
num <- c(1:4)

for (nums in num){
  integrated_obj[[nums]] <- SetIdent(integrated_obj[[nums]], value = integrated_obj[[nums]]$predictions)
}


##Ploting the transfered labels
for (objs in integrated_obj){
  print(SpatialFeaturePlot(objs, 
                           features = c("Astrocytes", "Oligodendrocytes", 
                                                        "Excitatory neurons"),
                           pt.size.factor = 1.6, ncol = 2, crop = TRUE))
}

```

# Step3
**a)** Integrating the same cell types from two brain regions
**b)** comparing the transcriptomic profile of the same cell type in two different brain regions


#### **a)** Integrating the same cell types from two brain regions
```r
library(Seurat)
SP1 <-readRDS("D:/Poland/PHD/spatial/Second_set/Single_Nucli/SP1.integrated.rds")
SP2<-readRDS("D:/Poland/PHD/spatial/Second_set/Single_Nucli/SP2.integrated.rds")
CT1 <- readRDS("D:/Poland/PHD/spatial/Second_set/Single_Nucli/CT1.integrated.rds")
CT2<-readRDS("D:/Poland/PHD/spatial/Second_set/Single_Nucli/CT2.integrated.rds")

SP1@meta.data$orig.ident<- sub("CR_SP1","CR_ON1",SP1@meta.data$orig.ident)
SP1@meta.data$orig.ident<- sub("HE_SP1","HE_ON1",SP1@meta.data$orig.ident)
SP2@meta.data$orig.ident<- sub("CR_SP2","CR_TN1",SP2@meta.data$orig.ident)
SP2@meta.data$orig.ident<- sub("HE_SP2","HE_TN1",SP2@meta.data$orig.ident)

CT1@meta.data$orig.ident<- sub("CR_CT1","CR_ON2",CT1@meta.data$orig.ident)
CT1@meta.data$orig.ident<- sub("HE_CT1","HE_ON2",CT1@meta.data$orig.ident)
CT2@meta.data$orig.ident<- sub("CR_CT2","CR_TN2",CT2@meta.data$orig.ident)
CT2@meta.data$orig.ident<- sub("HE_CT2","HE_TN2",CT2@meta.data$orig.ident)


DefaultAssay(SP1)<- "integrated"
DefaultAssay(SP2)<- "integrated"
DefaultAssay(CT1)<- "integrated"
DefaultAssay(CT2)<- "integrated"

SP1$region <- "Prefrontal"
SP2$region <- "Temporal"
CT1$region <- "Prefrontal"
CT2$region <- "Temporal"

SP1<-SplitObject(SP1, split.by = "predictions")
SP2<-SplitObject(SP2, split.by = "predictions")
CT1<-SplitObject(CT1, split.by = "predictions.label")
CT2<-SplitObject(CT2, split.by = "predictions.label")


SP.Astro.list<-list(SP1$Astrocytes, SP2$Astrocytes)
SP.Neuron.list<-list(SP1$`Excitatory neurons`, SP2$`Excitatory neurons`)
SP.Oligo.list<-list(SP1$Oligodendrocytes, SP2$Oligodendrocytes)

CT.Astro.list<-list(CT1$Astrocytes, CT2$Astrocytes)
CT.Neuron.list<-list(CT1$`Excitatory neurons`, CT2$`Excitatory neurons`)
CT.Oligo.list<-list(CT1$Oligodendrocytes, CT2$Oligodendrocytes)



##DATA integration
SP.Astro.anchors <- FindIntegrationAnchors(object.list = SP.Astro.list, dims = 1:30, scale = FALSE, anchor.features = 3000)
SP.Neuron.anchors <- FindIntegrationAnchors(object.list = SP.Neuron.list, dims = 1:30, scale = FALSE, anchor.features = 3000)
SP.Oligo.anchors <- FindIntegrationAnchors(object.list = SP.Oligo.list, dims = 1:30, scale = FALSE, anchor.features = 3000)

CT.Astro.anchors <- FindIntegrationAnchors(object.list = CT.Astro.list, dims = 1:30, scale = FALSE, anchor.features = 3000)
CT.Neuron.anchors <- FindIntegrationAnchors(object.list = CT.Neuron.list, dims = 1:30, scale = FALSE, anchor.features = 3000)
CT.Oligo.anchors <- FindIntegrationAnchors(object.list = CT.Oligo.list, dims = 1:30, scale = FALSE, anchor.features = 3000)


SP.Astro.integrated <- IntegrateData(anchorset = SP.Astro.anchors, dims = 1:30)
SP.Neuron.integrated <- IntegrateData(anchorset = SP.Neuron.anchors, dims = 1:30)
SP.Oligo.integrated <- IntegrateData(anchorset = SP.Oligo.anchors, dims = 1:30)

CT.Astro.integrated <- IntegrateData(anchorset = CT.Astro.anchors, dims = 1:30)
CT.Neuron.integrated <- IntegrateData(anchorset = CT.Neuron.anchors, dims = 1:30)
CT.Oligo.integrated <- IntegrateData(anchorset = CT.Oligo.anchors, dims = 1:30)




SP.Astro.integrated <- ScaleData(SP.Astro.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
SP.Astro.integrated <- RunPCA(SP.Astro.integrated, verbose = FALSE)
SP.Astro.integrated <- RunUMAP(SP.Astro.integrated, dims = 1:30, verbose = FALSE)

SP.Neuron.integrated <- ScaleData(SP.Neuron.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
SP.Neuron.integrated <- RunPCA(SP.Neuron.integrated, verbose = FALSE)
SP.Neuron.integrated <- RunUMAP(SP.Neuron.integrated, dims = 1:30, verbose = FALSE)

SP.Oligo.integrated <- ScaleData(SP.Oligo.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
SP.Oligo.integrated <- RunPCA(SP.Oligo.integrated, verbose = FALSE)
SP.Oligo.integrated <- RunUMAP(SP.Oligo.integrated, dims = 1:30, verbose = FALSE)

CT.Astro.integrated <- ScaleData(CT.Astro.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
CT.Astro.integrated <- RunPCA(CT.Astro.integrated, verbose = FALSE)
CT.Astro.integrated <- RunUMAP(CT.Astro.integrated, dims = 1:30, verbose = FALSE)

CT.Neuron.integrated <- ScaleData(CT.Neuron.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
CT.Neuron.integrated <- RunPCA(CT.Neuron.integrated, verbose = FALSE)
CT.Neuron.integrated <- RunUMAP(CT.Neuron.integrated, dims = 1:30, verbose = FALSE)

CT.Oligo.integrated <- ScaleData(CT.Oligo.integrated, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
CT.Oligo.integrated <- RunPCA(CT.Oligo.integrated, verbose = FALSE)
CT.Oligo.integrated <- RunUMAP(CT.Oligo.integrated, dims = 1:30, verbose = FALSE)
```

#### **b)** comparing the transcriptomic profile of the same cell type in two different brain regions
```r
library(ggplot2)

##Looking at the UMAP of each cell type in each slice seperately
DimPlot(SP.Astro.integrated, group.by = "region", split.by = "orig.ident")+ ggtitle("Astrocytes_AD")
DimPlot(SP.Neuron.integrated, group.by = "region", split.by = "orig.ident")+ ggtitle("Excitatory Neurons_AD")
DimPlot(SP.Oligo.integrated, group.by = "region", split.by = "orig.ident")+ ggtitle("Oligodendrocytes_AD")

DimPlot(CT.Astro.integrated, group.by = "region", split.by = "orig.ident")+ ggtitle("Astrocytes")
DimPlot(CT.Neuron.integrated, group.by = "region", split.by = "orig.ident")+ ggtitle("Excitatory Neurons")
DimPlot(CT.Oligo.integrated, group.by = "region", split.by = "orig.ident")+ggtitle("Oligodendrocytes")


##Comparing cell types in consecutive slices
p1 <- DimPlot(SP.Astro.integrated, group.by = "orig.ident", split.by = "region")+ ggtitle("Astrocytes_AD")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p2 <-DimPlot(SP.Neuron.integrated, group.by = "orig.ident", split.by = "region")+ ggtitle("Excitatory Neurons_AD")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p3 <-DimPlot(SP.Oligo.integrated, group.by = "orig.ident", split.by = "region")+ ggtitle("Oligodendrocytes_AD")+theme(plot.title = element_text(size = 15, face = "bold"))

p4 <-DimPlot(CT.Astro.integrated, group.by = "orig.ident", split.by = "region")+ ggtitle("Astrocytes")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p5 <-DimPlot(CT.Neuron.integrated, group.by = "orig.ident", split.by = "region")+ ggtitle("Excitatory Neurons")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p6 <-DimPlot(CT.Oligo.integrated, group.by = "orig.ident", split.by = "region")+ggtitle("Oligodendrocytes")+theme(plot.title = element_text(size = 15, face = "bold"))

library(ggrepel)
library(ggpubr)

ggarrange(p1, p2,p3,p4, p5,p6,
                    ncol = 3, nrow = 2)


##Comparing cell types between brain regions
p1 <- DimPlot(SP.Astro.integrated, group.by = "region")+ ggtitle("Astrocytes1")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p2 <-DimPlot(SP.Neuron.integrated, group.by = "region")+ ggtitle("Excitatory Neurons1")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p3 <-DimPlot(SP.Oligo.integrated, group.by = "region")+ ggtitle("Oligodendrocytes1")+theme(plot.title = element_text(size = 15, face = "bold"))

p4 <-DimPlot(CT.Astro.integrated, group.by = "region")+ ggtitle("Astrocytes2")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p5 <-DimPlot(CT.Neuron.integrated, group.by = "region")+ ggtitle("Excitatory Neurons2")+theme(plot.title = element_text(size = 15, face = "bold"),legend.position = "none")
p6 <-DimPlot(CT.Oligo.integrated, group.by = "region")+ggtitle("Oligodendrocytes2")+theme(plot.title = element_text(size = 15, face = "bold"))

library(ggrepel)
library(ggpubr)

ggarrange(p1, p2,p3,p4, p5,p6,
                    ncol = 3, nrow = 2)

```
