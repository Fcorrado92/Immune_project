#Import dependencies
library(Seurat)
library(scran)
library(harmony)
library(scuttle)
library(scater)
library(bluster)
library(tidyverse)
library(qs)
library(ggplot2)
library(readxl)
library(purrr)
library(kableExtra)

# find TCR files
files <- list.files("~/Immune_project/TCR")
# rename to have same order fopr TCR and GEX
#new_files <- ifelse(grepl("^TCR", files), paste0("G_", files), files)
# rename
#file.rename(files, new_files)

##Set input directory for all samples
input_dir <- "~/Immune_project/COUNTS/SoupX_results/"
##Create sample list of GEX_5
samples<-list.files(input_dir)

new_samples<-samples[-c(39:56)]
  
df<-as.data.frame(new_samples, files) 


read_data <- function(sample) {
  rds_file_path <- paste0(input_dir, sample, "/SoupX/data_i_soup.rds")
  data <- read_rds(rds_file_path)
}
samplelist <- lapply(samples, read_data)
names(samplelist) <- samples

##For each sample, concatenate Scrublet, scDblFinder & SCDS results and add to metadata
##Do this before filtering souporcell singlets
dbl_cat <- function(dir,samplelist){
  res <- list()
  for(i in 1:length(samplelist)){
    scrublet <- read_csv(paste0(dir,samplelist[i],"/Scrublet/",samplelist[i],"_Scrublet_results.csv"), col_types = cols())
    scrublet$Cell <- gsub("\"","",scrublet$Cell)
    cDC <- read_csv(paste0(dir,samplelist[i],"/scDblFinder/",samplelist[i],"_computeDoubletDensity_results.csv"), col_types = cols())
    fDC <- read_csv(paste0(dir,samplelist[i],"/scDblFinder/",samplelist[i],"_findDoubletClusters_results.csv"), col_types = cols())
    dbl <- read_csv(paste0(dir,samplelist[i],"/scDblFinder/",samplelist[i],"_scDblFinder_results.csv"), col_types = cols())
    sc <- read_csv(paste0(dir,samplelist[i],"/SCDS/",samplelist[i],"_SCDS_results.csv"), col_types = cols())
    all_dbl <- merge(scrublet,cDC,"Cell")
    all_dbl <- merge(all_dbl,fDC,"Cell")
    all_dbl <- merge(all_dbl,dbl,"Cell")
    all_dbl <- merge(all_dbl,sc,"Cell")
    rownames(all_dbl) <- all_dbl$Cell
    res[[i]] <- all_dbl
  }
  return(res)
}
dbl <- dbl_cat(input_dir,samples)

##For each sample, concatenate VDJ data
##Do this before filtering souporcell singlets
vdj_folder<-"~/Immune_project/TCR/"
samples_VDJ<-  list.files(vdj_folder)
setwd("~/Immune_project/TCR/")
#addTCRdata
getTCR <- function(samples) {
  #This function will return a dataframe with unique cell barcodes, all TRA/TRB chain clonotype information and MAIT/iNKT evidence in a single row.
  ##Load filtered_contig_annotations.csv files for all samples
  df_l <- lapply(samples, function(x) {
    read.csv(paste0(x, "/", "/filtered_contig_annotations.csv"))
  })
  ## Load clonotypes.csv files for all samples
  cl_l <- lapply(samples, function(x) {
    read.csv(paste0(x, "/", "/clonotypes.csv"))
  })
  ##Parse filtered_contig_annotations.csv files for each sample
  for (i in seq_along(df_l)){
    df <- df_l[[i]]
    df <- df[(df$chain %in% c("TRA","TRB")) & (df$productive %in% c("true","TRUE","True")),]
    ##For every chain per cell barcode, keep the row with the most UMIs (remove multiple rows per chain)
    ###Sorting by UMIs rather than reads.
    df <- df %>% group_by(barcode,chain) %>% arrange(desc(umis)) %>% slice_head(n=1)
    ##For every cell barcode, keep the two chains with the most UMIs (remove multiple chains per barcode)
    ###Using slice_head instead of top_n to handle ties in the number of UMIs per chain.
    df <- df %>% group_by(barcode) %>% arrange(desc(umis)) %>% slice_head(n=2)
    ##Add columns with TRA/TRB chain gene-level clonotype
    df <- df %>%
      mutate(TRAct = ifelse(chain == "TRA", paste(interaction(v_gene,  j_gene, c_gene)), NA)) %>%
      mutate(TRBct = ifelse(chain == "TRB", paste(interaction(v_gene,  j_gene, c_gene)), NA))
    ##Create dataframe with unique cell barcodes in rows and TRA/TRB chain gene-, amino acid- and nucleotide-level clonotype information
    out <- parseTCR(df)
    ##Add MAIT/iNKT evidence per cell
    cl <- cl_l[[i]]
    out$MAIT <- cl$mait_evidence[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    out$iNKT <- cl$inkt_evidence[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    out$clonotype_freq <- cl$proportion[match(as.character(out$clonotypeID), as.character(cl$clonotype_id))]
    df_l[[i]] <- out
  }
  names(df_l) <- samples
  return(df_l)
}
parseTCR <- function(df){
  #This function will output a dataframe of unique cell barcodes and their TRA & TRB chain
  #clonotype information in a single row
  ##Create vector of unique cell barcodes in dataframe
  unique_cb <- unique(df$barcode)
  ##Create output dataframe where each cell barcode is matched to heavy and light chain clonotype information in a single row
  out <- data.frame(matrix(nrow=length(unique_cb),ncol=10))
  colnames(out) <- c("barcode","TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene","TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene","clonotypeID")
  out$barcode <- unique_cb
  ##For every unique cell barcode:
  for (i in seq_along(unique_cb)){
    barcode.i <- out$barcode[i]
    location.i <- which(df$barcode == barcode.i)
    ##Write in its clonotype ID
    ##Write in its clonotype ID
    if(length(location.i)>1){
      out[i,"clonotypeID"] <- df[location.i[1],"raw_clonotype_id"]
    } else {
      out[i,"clonotypeID"] <- df[location.i,"raw_clonotype_id"]
    }
    ##If there are two rows corresponding to it in the input data
    if (length(location.i) == 2){
      ##and the first row is not a TRA chain row
      if (is.na(df[location.i[1],c("TRAct")])) {
        ##then, write in the TRA chain clonotype info from the second row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i[2], c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##and write in the TRB chain clonotype info from the first row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i[1], c("TRBct","cdr3","cdr3_nt","v_gene")]
      } else {
        ##otherwise, write in the TRA chain clonotype info from the first row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i[1], c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##and write in the TRB chain clonotype info from the second row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i[2], c("TRBct","cdr3","cdr3_nt","v_gene")]
      }
      ##otherwise, if there is only one row for this cell barcode
    } else if (length(location.i) == 1) {
      chain.i <- df$chain[location.i]
      ##and it corresponds to a TRA chain
      if (chain.i == "TRA"){
        ##write in the TRA chain clonotype information from that row
        out[i,c("TRA_gene_ct","TRA_cdr3_aa","TRA_cdr3_nt","TRA_V_gene")]<-df[location.i, c("TRAct","cdr3","cdr3_nt","v_gene")]
        ##otherwise
      } else {
        ##write in the TRB chain clonotype information from that row
        out[i,c("TRB_gene_ct","TRB_cdr3_aa","TRB_cdr3_nt","TRB_V_gene")]<-df[location.i, c("TRBct","cdr3","cdr3_nt","v_gene")]
      }
    }
  }
  return(out)
}
vdj<-getTCR(samples_VDJ)



for(i in 1:38){
  rownames(dbl[[i]])<-dbl[[i]]$Cell
  samplelist[[i]] <- Seurat::AddMetaData(samplelist[[i]],dbl[[i]])
  rownames(vdj[[i]])<-vdj[[i]]$barcode
  samplelist[[i]] <- Seurat::AddMetaData(samplelist[[i]],vdj[[i]])
  
}

for(i in 39:56){
  rownames(dbl[[i]])<-dbl[[i]]$Cell
  samplelist[[i]] <- Seurat::AddMetaData(samplelist[[i]],dbl[[i]])
}
  
  
for(i in 57:131){
  rownames(dbl[[i]])<-dbl[[i]]$Cell
  samplelist[[i]] <- Seurat::AddMetaData(samplelist[[i]],dbl[[i]])
  rownames(vdj[[i-18]])<-vdj[[i-18]]$barcode
  samplelist[[i]] <- Seurat::AddMetaData(samplelist[[i]],vdj[[i-18]])
  
}


#how many clonotypes x samples##
  
  df <- data.frame(Sample=character(), NumRows=integer(), NumNA=integer(),Perc=integer(), stringsAsFactors=FALSE)
  #extract colnames of vdj
  colnames<-colnames(samplelist[[37]]@meta.data)[25:37]
  #then add them to those samples without vdj

  # Loop sui campioni senza dati VDJ per aggiungere le colonne mancanti
  for (i in 39:56) {
    missing_cols <- setdiff(colnames, colnames(samplelist[[i]]@meta.data))
    
    # Aggiungi le colonne mancanti con valori NA
    for (col in missing_cols) {
      samplelist[[i]]@meta.data[[col]] <- NA
    }
  }
  
  # Correggi il ciclo for per scorrere tutti gli elementi della lista
for (i in 1:length(samplelist)) {
  sample_name <- names(samplelist)[i]  # Prende il nome del campione dalla lista
  num_rows <- nrow(samplelist[[i]]@meta.data)  # Conta il numero di righe nel dataset corrente
  num_NA <- sum(is.na(samplelist[[i]]$clonotypeID))  # Conta il numero di NA in clonotypeID
  Perc_NA <- num_NA/num_rows  # Conta il numero di NA in clonotypeID
  
  # Aggiunge le informazioni al dataframe
  df <- rbind(df, data.frame(Sample=sample_name, NumRows=num_rows, NumNA=num_NA, Perc=Perc_NA))
}


demux_file<-read_csv("~/Immune_project/analysis_oct24_fullDB/Final_demux_acc_sheet_11.10.24.csv")
meta_file<-read_excel("~/Immune_project/analysis_oct24_fullDB/Meta_Immune_proj.xlsx")
colnames(demux_file)<-c("Pair" , "Prcnt_concordance", "Pool",  "Genotype" , "assignment", "Sample_ID")
colnames(meta_file)<-c("Sample_ID", "Pool" ,  "Disease" ,  "Timepoint", "Therapy" ,  "Tissue" ,   "Cohort" ,   "ID"    ,    "BOR" ,      "CAR" )
demux_file<-demux_file[-c(7,8)]


#remove pools which do not contain accelerator samples
samplelist2<-samplelist[-c(39,43,48,50,51,52,53,54)]
samples2<-names(samplelist2)
samples2

demux_dir<-"/mnt/disks/cellranger/Accelerator_souporcell_v8/"
demux_file$assignment<-as.character(demux_file$assignment)

#assign souporcell meta, filter singlets, match soupcluster and pool, assign Sample_ID
for (i in 1:48) {
  # Load soup clusters for each pool, take only singlets
  singlets <- read_tsv(file = paste0(demux_dir, samples2[i], "/clusters.tsv"))
  # Clean it
  singlets <- as.data.frame(singlets) %>%
    select(1:3)
  rownames(singlets) <- singlets$barcode
  # Add metadata and filter singlets
  samplelist2[[i]] <- Seurat::AddMetaData(samplelist2[[i]], singlets)
  samplelist2[[i]] <- subset(samplelist2[[i]], subset = status %in% "singlet")
  # Set Library metadata based on sample prefix
  if (startsWith(samples2[[i]], "Accelerator")) {
    samplelist2[[i]]@meta.data$Pool <- strsplit(samples2[[i]], "_")[[1]][2]
  } else if (startsWith(samples2[[i]], "B")) {
    samplelist2[[i]]@meta.data$Pool <- strsplit(samples2[[i]], "_")[[1]][1]}
  
  samplelist2[[i]]$original_barcodes<-rownames(samplelist2[[i]]@meta.data)
  samplelist2[[i]]@meta.data<-left_join(samplelist2[[i]]@meta.data,demux_file[c("Pool","assignment","Genotype","Sample_ID")], by=c("Pool", "assignment"))
  samplelist2[[i]]@meta.data<-left_join(samplelist2[[i]]@meta.data,meta_file, by=c("Pool", "Sample_ID"))
  rownames(samplelist2[[i]]@meta.data)<-samplelist2[[i]]$original_barcodes
}

##check results
samplelist3 <- samplelist2[1:48]
# Initialize an empty dataframe
df <- data.frame(Pool = character(), Sample_ID = character(), cluster=character(), Count = integer(), Disease=character())
# Loop through each Seurat object in samplelist3
for (i in 1:length(samplelist3)) {
  print(samplelist3[i])
  # Get the Sample_ID table for the current Seurat object
  
  sample_counts <- as.data.frame(table(samplelist3[[i]]$Sample_ID, samplelist3[[i]]$Disease))%>%filter(Freq>0)
  
  # Add a Pool column for the current Seurat object
  sample_counts$Pool <- unique(samplelist3[[i]]$Pool) # You can modify how Pool is named
  
  # Reorder columns to match desired format (Pool, Sample_ID, Count)
  sample_counts <- sample_counts[, c("Pool", "Var1", "Freq","Var2")]
  colnames(sample_counts) <- c("Pool", "Sample_ID", "Count","Disease")
  
  # Append to the df
  df <- rbind(df, sample_counts)
}
# View the resulting dataframe
df


#now assign a Pool label to the unpooled samples and then assign meta based on Pool
meta_file$Pool[meta_file$Pool=="CP1003"]<- "CP1003PB"
meta_file$Pool[meta_file$Pool=="MXMERZ002A_1_R"]<- "MXMERZ002A_01_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_2_R"]<- "MXMERZ002A_02_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_3_R"]<- "MXMERZ002A_03_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_4_R"]<- "MXMERZ002A_04_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_5_R"]<- "MXMERZ002A_05_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_6_R"]<- "MXMERZ002A_06_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_7_R"]<- "MXMERZ002A_07_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_8_R"]<- "MXMERZ002A_08_R" 
meta_file$Pool[meta_file$Pool=="MXMERZ002A_9_R"]<- "MXMERZ002A_09_R" 

for (i in 49:length(samples2)) {
  if (grepl("^X|^CP|^pM", samples2[[i]])) {  # If starts with X, C, or pM
    samplelist2[[i]]@meta.data$Pool <- strsplit(samples2[[i]], "_")[[1]][1]
  } 
  else if (startsWith(samples2[[i]], "PT")) {    # If starts with PT, take first two parts
    samplelist2[[i]]@meta.data$Pool <- paste(strsplit(samples2[[i]], "_")[[1]][1:2], collapse = "_")
  } 
  else if (startsWith(samples2[[i]], "MX") || startsWith(samples2[[i]], "G")) {  # If starts with MX or G, take the full name
    samplelist2[[i]]@meta.data$Pool <- samples2[[i]]
  } 
    } 

for (i in 49:length(samples2)) {
  samplelist2[[i]]$original_barcodes<-rownames(samplelist2[[i]]@meta.data)
  samplelist2[[i]]@meta.data<-left_join(samplelist2[[i]]@meta.data,meta_file, by=c("Pool"))
  rownames(samplelist2[[i]]@meta.data)<-samplelist2[[i]]$original_barcodes
}


#now check colnames discrepancies between 48 and 48 samples
colnames1<-colnames(samplelist2[[48]]@meta.data)
colnames2<-colnames(samplelist2[[49]]@meta.data)
diff<-setdiff(colnames1, colnames2)
print(diff)

#for the first 48 samples remove "status","assignment", "Genotype"
for (i in 48) {
  samplelist2[[i]]@meta.data <- samplelist2[[i]]@meta.data %>%
    select(-status, -assignment, -Genotype)
}

#now check colnames discrepancies between 48 and 48 samples
colnames1<-colnames(samplelist2[[48]]@meta.data)
colnames2<-colnames(samplelist2[[49]]@meta.data)
diff<-setdiff(colnames1, colnames2)
print(diff)

for (i in 2:length(samples2)) {
  samplelist2[[i]] <- Seurat::RenameCells(samplelist2[[i]], add.cell.id=i)
}
#save
qsave(samplelist2, "/mnt/disks/cellranger/Accelerator_qs_files_oct24/samplelist2.qs")
#reload
samplelist2<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/samplelist2.qs")
#Merge filtered Seurat objects
merged <- purrr::reduce(samplelist2,merge)
qsave(merged,"/mnt/disks/cellranger/Accelerator_qs_files_oct24/merged.qs" )


merged<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/merged.qs")
#remove JJ samples and remove cells assigned to unmatched soup clusters(B15, B16, B17)
metadata<-merged@meta.data
t<-as.data.frame(metadata%>%group_by(Pool, Disease, Cohort)%>%summarise(n=n_distinct(Sample_ID)))
keep<-setdiff(unique(merged$Cohort), "JJ")
merged_filt<-subset(merged, subset= Cohort%in%keep)


#quality checks
merged_filt[["percent.mt"]] <- Seurat::PercentageFeatureSet(merged_filt, pattern = "^MT-")
filtered <- subset(merged_filt, subset = nFeature_RNA > 200 & percent.mt < 10)
rm(merged_filt)

# Total number of cells
total_cells <- rownames(filtered@meta.data)

# Randomly sample half of the cells
set.seed(123)  # Set seed for reproducibility
split_cells <- sample(total_cells, size = length(total_cells) / 2)

# Subset the Seurat object into two parts
sobj_part1 <- subset(sobj, cells = split_cells)
sobj_part2 <- subset(sobj, cells = setdiff(total_cells, split_cells))

qsave(sobj_part1,"/mnt/disks/cellranger/Accelerator_qs_files_oct24/sobj1.qs")
qsave(sobj_part2,"/mnt/disks/cellranger/Accelerator_qs_files_oct24/sobj2.qs")

sobj1<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/sobj1.qs")
sobj1<-JoinLayers(sobj1)
##Normalize cells of the merged object with scran & create Seurat Object
print("Normalizing merged object with scran...")
sce <- as.SingleCellExperiment(sobj1)
rm(sobj1)
gc()
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- logNormCounts(sce)
norm <- as.Seurat(sce)
qsave(norm, "/mnt/disks/cellranger/Accelerator_qs_files_oct24/norm1.qs")









sobj2<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/sobj2.qs")
sobj2<-JoinLayers(sobj2)
##Normalize cells of the merged object with scran & create Seurat Object
print("Normalizing merged object with scran...")
sce <- as.SingleCellExperiment(sobj2)
rm(sobj2)
gc()
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- logNormCounts(sce)
norm <- as.Seurat(sce)

qsave(norm, "/mnt/disks/cellranger/Accelerator_qs_files_oct24/norm2.qs")
norm<-qread("/mnt/disks/cellranger/Accelerator_qs_files_oct24/norm2.qs")

##Reduce dimensions of scran-normalized merged object
print("Finding variable features...")
norm <- FindVariableFeatures(object=norm, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
print("Scaling...")
norm <- ScaleData(object=norm, verbose = FALSE) 
print("Running PCA...")
norm <- RunPCA(object=norm, npcs = 30, verbose = FALSE)


norm <- RunUMAP(norm, reduction = "pca", dims=1:30)
norm <- FindNeighbors(norm, reduction = "pca", dims=1:30)
norm <- FindClusters(norm, resolution = 1)

output_dir<-"~/Immune_project/analysis_oct24_fullDB/"

norm[['ident']]<-Idents(norm)
##Plot UMAP embedding of corrected dataset & cell-type-specific markers to evaluate integration performance and annotate clusters
print("Plotting UMAPs...")
pdf(paste0(output_dir,"/PreCorrection_UMAP.pdf"),width = 10)
print(DimPlot(object = norm, reduction="umap", pt.size = .5, label=TRUE))
dev.off()

mxm<-unique(grep("MXMERZ",norm@meta.data$Cohort, value=TRUE))
norm@meta.data$Cohort[norm@meta.data$Cohort%in%mxm]<-"MXM"

##Plot PCA scatterplots colored by covariates post-correction
covariates<-c("Sample_ID","Pool","Tissue","Cohort")
pdf(paste0(output_dir,"/PreCorrection_Covariate_UMAP.pdf"),width=30,height=20)
for(i in 1:length(covariates)){
  print(DimPlot(object = norm, reduction="umap", pt.size =1, group.by = covariates[i]))
}
dev.off()

#EloRD-BM batch
norm<-RenameIdents(norm, c('3'="1",'14'="7"))
norm[['ident']]<-Idents(norm)
DimPlot(norm, label=T)

##Batch_markers
pdf(paste0(output_dir,"/CD3D_Feature_precorrection.pdf"))
FeaturePlot(norm, "CD3D")
dev.off()
batch_markers<-FindMarkers(norm, ident.1=7, ident.2 = 1, logfc.threshold = 1, min.pct = 0.1)
batch_markers$genes<-rownames(batch_markers)
##Integrate with Harmony
integrated <- RunHarmony(norm, group.by.vars = "Pool")

##Cluster and embed corrected dataset
print("Clustering corrected dataset...")
integrated <- RunUMAP(integrated, reduction = "harmony", dims=1:30)
integrated <- FindNeighbors(integrated, reduction = "harmony", dims=1:30)
integrated <- FindClusters(integrated, resolution = 1)

##Plot UMAP embedding of corrected dataset & cell-type-specific markers to evaluate integration performance and annotate clusters
print("Plotting UMAPs...")
pdf(paste0(output_dir,"/PostCorrection_UMAP.pdf"),width = 20)
print(DimPlot(object = integrated, reduction="umap", pt.size = .1, label=TRUE))
dev.off()

pdf(paste0(output_dir,"/PostCorrection_Covariate_UMAP.pdf"),width = 20)
for(i in 1:length(covariates)){
  print(DimPlot(object = integrated, reduction="umap", pt.size = .1, split.by = covariates[i], label=TRUE))
}
dev.off()

###Mark cell with fewer than 500 counts and log-transform nCount for plotting
integrated$nCount_500 <- integrated$nCount_RNA < 500
integrated$lognCount <- log10(integrated$nCount_RNA)
features <- c("doublet_scores","DoubletScore","Doublet","scDblFinder.score","hybrid_score","SoupX_rho","lognCount", "nCount_500", "percent.mt","HBB","PPBP","MS4A2","CD34","DNTT","CDK6","TOP2A","MKI67","TUBB","CXCL12","LEPR","KITLG","NES","CSPG4","BGLAP","CDH5","ACTA2","S100A4","ACAN","COL2A1","ENG","VIM","VCAM1","ISG15","IFI44L","IFI6","CCL3","CCL4","CD3D","CCR7","SELL","IL7R","CD4","FOXP3","TBX21","GATA3","RORA","RORC","KLRB1","CXCR5","BCL6","CXCR3","CXCR6","CCR6","CCL5","LGALS1","S100A9","S100A11","CD8A","GZMK","TOX","TIGIT","EOMES","PDCD1","GZMH","HLA-DRA","HLA-DRB1","CD74","GZMB","GNLY","PRF1","NKG7","XCL2","B3GAT1","TRAC","TRGC1","TRGC2","TRDC","LYZ","CD14","FCGR3A","FCER1A","MPO","AZU1","ELANE","MZB1","CLEC4C","MS4A1","CD19","CD24","CD27","MYB","TNFRSF17","SDC1","SLAMF7","CD38","IGKC","IGLC1","IGLC2","IGHG1","IGHG2","IGHG3","IGHA1","IGHA2","IGHD","IGHM","IGHE","CCND1","CCND2","CCND3","NSD2","FGFR3","MAF","MAFB","FRZB","DKK1","SPN","CD5","MME")
pdf(paste0(output_dir,"/PostCorrection_Feature_UMAP.pdf"),width = 12)
for(i in 1:length(features)){
  print(FeaturePlot(integrated, features=features[i]))
}
dev.off()

##Write integrated object as .rds object
print("Saving Seurat Object...")
saveRDS(integrated, paste0(output_dir,"/integrated2.rds"))

##Get DE markers for all clusters
print("Performing Differential Expression for all clusters...")
markers <- FindAllMarkers(integrated)
write_csv(markers,paste0(output_dir,"/AllMarkers.csv"))





