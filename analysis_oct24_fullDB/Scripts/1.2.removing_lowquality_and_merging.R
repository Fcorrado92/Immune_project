library(qs)
library(Seurat)
integrated1<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/integrated1.qs")
integrated1[['ident']]<-Idents(integrated1)

FeaturePlot(integrated1, c("lognCount","percent.mt"), label=T,order=T)

integrated1<-RenameIdents(integrated1,
              c('7'="lowq",'16'="lowq",'36'="lowq",'22'="lowq",
                '21'="lowq",'32'="lowq",'33'="lowq",
                '37'="lowq",'6'="lowq",'10'="lowq"))
integrated1[['ident']]<-Idents(integrated1)
DimPlot(integrated1, label=T)
keep<-setdiff(unique(Idents(integrated1)), "lowq")
integrated1_filt<-subset(integrated1, idents = keep)

integrated2<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/integrated2.qs")

FeaturePlot(integrated2, c("lognCount","percent.mt"), label=T,order=T)
integrated2<-RenameIdents(integrated2,
                          c('22'="lowq",'3'="lowq",'35'="lowq",'34'="lowq",
                            '31'="lowq",'7'="lowq",'24'="lowq",
                            '20'="lowq",'15'="lowq"))
integrated2[['ident']]<-Idents(integrated2)
DimPlot(integrated2, label=T)
keep<-setdiff(unique(Idents(integrated2)), "lowq")
integrated2_filt<-subset(integrated2, idents = keep)


qsave(integrated1_filt,"/mnt/disks/cellranger/Accelerator_qs_files_oct24/integrated1_filt.qs")
qsave(integrated2_filt,"/mnt/disks/cellranger/Accelerator_qs_files_oct24/integrated2_filt.qs")


# ---------Accelerator_qs_files_oct24 this has been moved,now I move the two qs files to analysis_oct24_fullDB----------------------------------------------------------------
library(qs)
integrated1_filt<-qread("/mnt/disks/disk/full_dataset_qs/integrated1_filt.qs")
integrated2_filt<-qread("/mnt/disks/disk/full_dataset_qs/integrated2_filt.qs")
integrated_filt<-merge(integrated1_filt, integrated2_filt)

integrated_filt <- FindVariableFeatures(object=integrated_filt, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
integrated_filt <- ScaleData(object=integrated_filt, verbose = FALSE) 
print("Running PCA...")
integrated_filt <- RunPCA(object=integrated_filt, npcs = 30, verbose = FALSE)

##Integrate with Harmony
integrated_filt <- RunHarmony(integrated_filt, group.by.vars = "Pool")

##Cluster and embed corrected dataset
print("Clustering corrected dataset...")
integrated_filt <- RunUMAP(integrated_filt, reduction = "harmony", dims=1:30)
integrated_filt <- FindNeighbors(integrated_filt, reduction = "harmony", dims=1:30)
integrated_filt <- FindClusters(integrated_filt, resolution = 1)

qsave(integrated_filt, "/mnt/disks/disk/full_dataset_qs/integrated_filt_final_jan25.qs")

