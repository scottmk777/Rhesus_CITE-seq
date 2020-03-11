library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(data.table)
library(gridExtra)
library(parallel)
library(ggplot2)
library(umap)
library(ComplexHeatmap)
set.seed(7)


samp_info <- fread("~/data/cell_information.csv")

### ADT (Protein)
adt_mat <- readMM("~/data/ADT_Matrix.mtx")
adt_rn <- fread("~/data/features_ADTmat.tsv", header = F)
adt_barcodes <- fread("~/data/barcodes_ADTmat.tsv", header = F)

colnames(adt_mat) = adt_barcodes$V1
rownames(adt_mat) = adt_rn$V1


### Gene Expr
gene_mat <- readMM(file="~/data/Gene_Matrix.mtx")
gene_nms <- fread("~/data/features_genemat.tsv", header = F)
gene_barcodes <-  fread("~/data/barcodes_genemat.tsv", header = F)

gene_nms <- ifelse(gene_nms$V1 == "", "Unknown", gene_nms$V1)

gene_mat <- as.matrix(gene_mat)
colnames(gene_mat) = gene_barcodes$V1
rownames(gene_mat) = gene_nms


# check order
all(gene_barcodes == adt_barcodes)

# create Seurat object 
srt <- CreateSeuratObject(gene_mat)
srt <- NormalizeData(srt, display.progress = FALSE)
srt <- ScaleData(srt, display.progress = FALSE)

# add and scale ADT data
srt[["ADT"]] <- CreateAssayObject(counts = adt_mat)
srt <- NormalizeData(srt, assay = "ADT",  normalization.method = "CLR")

# find variable features - only use genes here
srt <- FindVariableFeatures(srt, do.plot = F, display.progress = FALSE)
srt <- RunPCA(object = srt, verbose = FALSE)
srt <- FindNeighbors(object = srt, dims = 1:7)
srt <- FindClusters(object = srt, resolution = 0.2)
srt <- RunTSNE(object = srt, dims = 1:7, method = "FIt-SNE", check_duplicates = FALSE)
TSNEPlot(srt,label = TRUE)

saveRDS(file = "~/data/srt_aftertzzle.RDS", srt)


# this tsne plot should be nearly identical to Figure 4c, but will be a little different
# given the stochastic nature of tSNE. The exact embeddings for each cell that will perfectly
# reproduce 4c can be found in samp_info

# look at distribution ADTs
FeaturePlot(object = srt, features = c("CD8A", "CD4", "ADT-CD8", "ADT-CD4", "GZMB", "ADT-HLA-DR", "ADT-CD20"))
FeaturePlot(object = srt, features = c("ADT-CD11c", "ADT-CD8", "ADT-CD4" , "ADT-CD16", "ADT-CD14", "ADT-TCRa-b", 
                                       "ADT-CD69", "ADT-CD28", "ADT-CD95", "ADT-CD20", "ADT-HLA-DR", "ADT-CD123",
                                       "ADT-BDCA1"))


# identify markers between each cluster
srt.markers <- FindAllMarkers(object = srt, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

fwrite(file = "~/csvs/forgit_srt_markers.csv", srt.markers)

# big heatmaps
# Seurat can also do this for you
# Im toggling to the samp_info information as all the cell types are labeled
barcode_clust_df <- subset(samp_info, select = c("Barcode", "Cluster"))
# pick the first 100 cells
all_cell_barcodes <- split(as.data.table(barcode_clust_df), by = "Cluster")
subset_cell_barcodes <- lapply(all_cell_barcodes, function(x) x[1:100])
used_barcode_subset <- rbindlist(subset_cell_barcodes)

# now lets drop MT and Unknown genes for intepretability
srt.markers_clean <- srt.markers[-grep("^MT", srt.markers$gene),]
srt.markers_clean <- srt.markers_clean[-grep("Unknown", srt.markers_clean$gene),]


# calculate the difference between the % expr of a cluster
# vs other clusters
srt.markers_clean$pc_diff <- srt.markers_clean$pct.1-srt.markers_clean$pct.2

srt.markers_clean <- srt.markers_clean[order(srt.markers_clean$pc_diff, decreasing = TRUE),]
srt.markers_clean <- srt.markers_clean[!duplicated(srt.markers_clean$gene),]

clust_def_genes_ls <- split(as.data.table(srt.markers_clean), by = "cluster")
clust_def_genes_ls_sub <- lapply(clust_def_genes_ls, function(x) x[1:5,])
clust_def_genes_sub <-rbindlist(clust_def_genes_ls_sub)
clust_def_genes_sub$gene

proj_clust_genes <- gene_mat[which(rownames(gene_mat) %in% clust_def_genes_sub$gene),]
proj_clust_genes_sub <- proj_clust_genes[,which(colnames(gene_mat) %in% used_barcode_subset$Barcode)]

# order both rows and columns
# order columns
df_for_column_order = data.frame(id = colnames(proj_clust_genes_sub), 
                                 cluster = samp_info$Cluster[which(samp_info$Barcode %in% colnames(proj_clust_genes_sub))], 
                                 group = samp_info$Group[which(samp_info$Barcode %in% colnames(proj_clust_genes_sub))])
df_for_column_order_ordered <- df_for_column_order[order(as.numeric(as.character((df_for_column_order$cluster))), df_for_column_order$group),]
df_for_column_order_ordered$id[1:5]
df_for_column_order_ordered$cluster[1:25]
matched_cols <- match(as.character(df_for_column_order_ordered$id), colnames(clust_def_genes_sub))
column_order <- df_for_column_order_ordered$id

### order rows
row_ordered <- clust_def_genes_sub[order(as.numeric(clust_def_genes_sub$cluster)),]
row_ordered <- row_ordered[which(row_ordered$gene %in% rownames(proj_clust_genes_sub)),]
row_ordered_order <- row_ordered$gene

# actually do the reordering
proj_clust_genes_sub <- proj_clust_genes_sub[match(row_ordered_order, rownames(proj_clust_genes_sub)),]
proj_clust_genes_sub <-  proj_clust_genes_sub[,match(df_for_column_order_ordered$id, colnames(proj_clust_genes_sub))]

proj_clust_genes_sub[1:5,1:5]



# create the default ggplot 9 colors
library(scales)
my_color_palette <- hue_pal()(9)
cluster_colors <- my_color_palette
names(cluster_colors) <- as.character(seq(0:8)-1)

# and custom colors for the heatmap
library(circlize)
col_fun = colorRamp2(c(0, 2, 4), c("#330066", "black", "yellow"))


column_ha = HeatmapAnnotation(Cluster = df_for_column_order_ordered$cluster, 
                              col = list(Cluster = cluster_colors))


Heatmap(proj_clust_genes_sub, name = "Expr", top_annotation = column_ha,
        cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, 
        show_row_dend  = FALSE, cluster_rows = FALSE, row_split = row_ordered$cluster,
        col = col_fun)


######### Pathways etc
# for GSEA we need all genes
custom_wilcox <- function(x){
  return(tryCatch(wilcox.test(x ~ mycounts_has_grp_nona$Group, paired = FALSE), error=function(e) data.frame(p.value = "FAIL")))
}

### MAST Analysis
library(MAST)
library(GSEABase)
aclust = 2
mycounts_has_grp_ind <- which(samp_info$Group %in% c("DMSO", "Gag") & samp_info$Cluster == aclust)

my_mait <- gene_mat
my_cdat <- subset(samp_info, select = c("Barcode", "Group", "Cluster", "sequencing_round", "monkey"))
rownames(my_cdat) <- my_cdat$Barcode
scaRaw_2 <- FromMatrix(my_mait, my_cdat)


FCTHRESHOLD <- .05

cdr2 <-colSums(assay(scaRaw_2)>0)
colData(scaRaw_2)$cngeneson <- scale(cdr2)
cond<-factor(colData(scaRaw_2)$Group)
cond<-relevel(cond,"DMSO")
colData(scaRaw_2)$condition <- cond
zlmCond <- zlm(~condition + cngeneson, scaRaw_2)
summaryCond <- summary(zlmCond, doLRT='conditionGag') 

#print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='conditionGag' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='conditionGag' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(scaRaw_2)), by='primerid')
setorder(fcHurdleSig, fdr)
colnames(fcHurdleSig)[3] <- "Log2_FoldChange"

# run the bootstrap - this takes awhile
boots <- bootVcov1(zlmCond, R=10)

# use the Blood Transcription Modules (BTM)
module <- "BTM"
min_gene_in_module <- 5
packageExt <- system.file("extdata", package='MAST')
module_file <- list.files(packageExt, pattern = module, full.names = TRUE)
gene_set <- getGmt(module_file)
gene_ids <- geneIds(gene_set)
sets_indices <- limma::ids2indices(gene_ids, mcols(scaRaw_2)$primerid)
sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]

# run GSEA
gsea <- gseaAfterBoot(zlmCond, boots, sets_indices, CoefficientHypothesis("conditionGag")) 
z_stat_comb <- summary(gsea, testType='normal')

# pick sig modules that are informative
sigModules <- z_stat_comb[combined_adj<.05]
sigModules_noTBA <- sigModules[-grep("TBA", sigModules$set),]
gseaTable <- melt(sigModules_noTBA[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')


# which genes are in which set?
gene_ids_1 <- names(gene_ids)
gene_ids_l <- lapply(gene_ids, function(x){
  df_x <- as.data.frame(paste(x, collapse = " ,"))
  df_x$name <- names(x)
  return(df_x)
} )

gene_ids_df <- rbindlist(gene_ids_l)
gene_ids_df$path_name <- gene_ids_1
colnames(gene_ids_df)[1] <- "pathway"
gene_ids_df <- gene_ids_df[-grep("TBA",gene_ids_df$path_name),]






