#################################################################################
# Analysis of cancerous T lymphocytes 
########################################

library(Seurat)
library(scater)
library(dplyr)
library(ggplot2)
colors = data.frame(
  "cell_type_color" = c(
    "#765aa4ff",
    "#ad2070ff",
    "#fdd62eff",
    "#96c9e6ff",
    "#f48121ff",
    "#68c3a5ff",
    "#ef3e2cff",
    "#0d522aff",
    "#42b649ff",
    "#660a17ff",
    "#3f7fc1ff",
    "#189cd5ff",
    "#0e8342ff",
    "#f9ae1dff",
    "#552e8aff",
    "#8b87c0ff",
    "#984d9dff",
    "#fec964ff",
    "#126badff",
    "#a51d23ff",
    "#e5569fff",
    "#eae71cff"
  )
)
maindir = "/home/localadmin/Documents/Darrier/"


cell_markers = c("PTPRC", "CD3E","CD4","CD8A", "FOXP3", "CD8B","TIGIT","IL2RA","KLRB1","GNLY","TRAC","TRBC1",
                 "LAMP3","CD40",  "CD14", "CD16", "ITGAM", "HLA-DRA","HLA-DRQA1","HLA-DRQB1", "HLA-DRB1", "CD163", "CD1A", "CD207",
                 "MS4A2", "MS4A1",
                 "CDH5","PECAM1","CD34","KRT5","KRT10","KRT14","KRT1", "FABP5", "SPPR2G", "CDSN", 
                 "ITGA6", "ITGB1", "GRHL1", 
                 "CD200", "SOX9", "KRT19", "APOC1", "ACSL5", "ABCC3",
                 "LUM", "PDGFRB", "COL1A1",
                 "SOX10", "S100B", "DCT", "TYRP1", "MLANA")


# Directories -------------------------------------------------------------
output = file.path(maindir, "output");  if(!dir.exists(output)){dir.create(output)}
QCdir = file.path(output, "QC");  if(!dir.exists(QCdir)){dir.create(QCdir)}

# Replace GEO path with the path to the GEO directory downloaded on your computer (scRNA MM468)
input = file.path(maindir, "input") 

# Load Data ---------------------------------------------------------------
BCmodel <- "Darrier"

#################################################################################
# Loading & QC
####################

samples = paste0("DAR", c(2,3,4,7))
metadata_all = lapply(samples, function(i) read.csv(list.files(file.path(input, i), pattern = "meta", full.names = T), sep = "\t"))
metadata_all = do.call("rbind", metadata_all)

scRNA = list()
files = setNames(unlist(lapply(samples, function(i) file.path(input, i, i, "outs", "filtered_feature_bc_matrix"))), samples)
for(i in seq_along(files)){
  name = names(files)[i]
  scRNA[[name]] = Seurat::Read10X(data.dir = files[i])
  gc()
}
scRNA = do.call("cbind", scRNA)
qs::qsave(scRNA, file.path(output, "mat.qs"))

rownames(metadata_all) = paste0(metadata_all$orig.ident, "_", metadata_all$Cell.Barcode)
Seu = CreateSeuratObject(scRNA, meta.data = metadata_all)
qs::qsave(Seu, file.path(output, "Seu.qs"))

TCR_genes = grep("^TRBV|^TRAV", rownames(Seu), value = T)

# QC, normalization  ------------------------------------------------------
keep_feat <- (rowSums(Seu)>0)>0 #remove genes that are not expressed in any cells
Seu <- Seu[keep_feat,]
Seu = Seu[, Seu$Percent_mtRNA < 25]

#################################################################################
# Seurat Analysis
####################
library(IDclust)
seurat_dir = file.path(output, "Seurat"); if(!dir.exists(seurat_dir)) dir.create(seurat_dir)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Idents(Seu) = Seu$orig.ident
Seu <- NormalizeData(Seu,  verbose = FALSE)
Seu <- FindVariableFeatures(Seu, selection.method = "vst")
Seu <- ScaleData(Seu, features = rownames(Seu))
Seu <- RunPCA(Seu, features = VariableFeatures(Seu), ndims.print = 1:10, nfeatures.print = 10)
FeaturePlot(Seu, reduction = "pca", features = "nFeature_RNA", cols = rev(viridis::magma(100)))
DimHeatmap(Seu, dims = 1:3, cells = 1000, balanced = TRUE)
ElbowPlot(Seu)
VlnPlot(Seu, c("IL17A","IL17B","IL17C","IL17D", "IL17E", "IL17F", "IL17G"), ncol = 4 )

Seu <- FindNeighbors(Seu, dims = 1:20)
Seu <- RunUMAP(Seu, dims = 1:20)
Seu <- Seurat::FindClusters(Seu)
Seu <- CellCycleScoring(Seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Seu = IDclust::iterative_differential_clustering(Seu, assay = "RNA")
IDC_DA = read.csv("IDC_DA.csv")
IDC_summary = read.csv("IDC_summary.csv")

pdf(file.path(seurat_dir,paste0("scRNA.pdf")))
print(DimPlot(Seu, reduction = "umap", group.by = "seurat_clusters"))
print(DimPlot(Seu, reduction = "umap", group.by = "IDcluster"))
plot_cluster_network(Seu, IDC_summary = IDC_summary, color_by = "seurat_clusters", assay = "originalexp")
print(DimPlot(Seu, reduction = "umap", group.by = "orig.ident"))
print(DimPlot(Seu, reduction = "umap", group.by = "Phase"))
print(FeaturePlot(Seu,  features = "nCount_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
print(FeaturePlot(Seu,  features = "nFeature_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
dev.off()


png(file.path(seurat_dir,paste0("UMAP_sample.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(Seu, reduction = "umap", pt.size = 1.5, cols = colors$cell_type_color, shuffle = TRUE, group.by = "orig.ident") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha =  0.7 
tplot
dev.off()








topmarkers = IDclust::top_differential_markers(
  IDC_DA,
  top = 10,
  gene_col = "gene",
  logFC_col = "avg_log2FC",
  qvalue_col = "p_val_adj",
  order_by = "logFC_col",
  pseudogene_pattern = NULL
)
write.csv(IDC_DA, file = file.path(seurat_dir, paste0("DA_df.csv")), quote = TRUE, row.names = FALSE)

# Concatenate top 5 markers per cluster/cluster_of_origin
topmarkers = topmarkers %>% dplyr::group_by(cluster_of_origin, cluster) %>%
  dplyr::summarise(Term = paste(gene, collapse = " "))

library(enrichR)
pathway_df = top_enriched_pathways(
  IDC_DA,
  top = 25,
  gene_col = "gene",
  qval.th = 0.1, 
  max_genes_enriched = 1000
)

write.csv(pathway_df, file = file.path(seurat_dir, paste0("pathway_df.csv")), quote = TRUE, row.names = FALSE)

pdf(file.path(seurat_dir, paste0("IDC_tree_colored_by_IDClust_wolegend.pdf")), height = 10, width = 10, pointsize = 24)
set.seed(47)
print(plot_cluster_network(Seu,
                           IDC_summary = IDC_summary,
                           color_by = "IDcluster",
                           colors = sample(colors$cell_type_color, length(unique(Seu$IDcluster))),
                           legend = F)
)
dev.off()

pathway_df_filtered = pathway_df
pdf(file.path(seurat_dir, paste0("IDC_tree_colored_by_Annotation_wpaths_filtered.pdf")), height = 10, width = 10, pointsize = 24)
set.seed(45)
print(
  plot_cluster_network(Seu, 
                       IDC_summary = IDC_summary,
                       color_by = "IDcluster",
                       colors = sample(colors$cell_type_color, length(unique(Seu$IDcluster))),
                       legend = F,
                       edge_df = pathway_df_filtered,
                       edge.label.cex = 0.5)
)
dev.off() 

library(dplyr)
# Cell type  markers
for(i in c("All")){
  cell_marker_db = read.csv("/media/localadmin/T7//InstitutCurie/Documents/Data/Annotation/cell_marker_db.csv", header = TRUE)
  
  cell_marker_db = cell_marker_db[grep(c("Skin|Dermis|Epithelium|Blood|Immune|Lymph"), ignore.case = T, cell_marker_db$organ),]
  genes = intersect(toupper(IDC_DA$gene), cell_marker_db$official.gene.symbol)
  
  cell_marker_db = cell_marker_db[which(cell_marker_db$official.gene.symbol %in% genes),]
  table(cell_marker_db$cell.type)
  
  IDC_DA$cell_type = cell_marker_db$cell.type[match(toupper(IDC_DA$gene), cell_marker_db$official.gene.symbol)]
  IDC_DA = IDC_DA %>% dplyr::mutate(
    rank_group_activation = order(pct.1),
    rank_logFC = order(avg_log2FC),
    rank_qval = order(p_val_adj),
  )
  cell_type_markers = IDC_DA %>% dplyr::filter(!is.na(cell_type)) %>% 
    group_by(cluster_of_origin, cluster, cell_type) %>%
    summarise(total = n(),
              avg_rank_group_activation = mean(rank_group_activation),
              avg_rank_logFC = mean(rank_logFC),
              avg_rank_qval = mean(rank_qval),
              Gene = gene[1]
    ) %>% 
    mutate(
      combined_score_rank = nrow(IDC_DA) / ((avg_rank_group_activation + avg_rank_logFC + avg_rank_qval) / 3)
    ) %>% slice_max(total, with_ties = TRUE) %>% slice_max(combined_score_rank)
  cell_type_markers = cell_type_markers %>% mutate(Term = paste0(cell_type, " | ", total))
  cell_type_markers = cell_type_markers[,c("cluster_of_origin", "cluster", "Term")]
  
  pdf(file.path(seurat_dir, paste0("IDC_tree_colored_by_Annotation_wCellType_", i, ".pdf")), height = 10,
      width = 10, pointsize = 12)
  print(plot_cluster_network(Seu, 
                             IDC_summary = IDC_summary,
                             color_by = "IDcluster",
                             legend = F,
                             edge_df = cell_type_markers,
                             colors = sample(colors$cell_type_color, length(unique(Seu$IDcluster))),
                             edge.label.cex = 0.6)
  )
  dev.off() 
  
}

pdf(file.path(seurat_dir, paste0("IDC_tree_colored_by_Annotation_wCellType_Markers.pdf")), height = 10,
    width = 10, pointsize = 12)
cell_marker_db = read.csv("/media/localadmin/T7//InstitutCurie/Documents/Data/Annotation/cell_marker_db.csv", header = TRUE)

cell_marker_db = cell_marker_db[grep(c("Skin|Dermis|Epithelium|Blood|Immune|Lymph"), cell_marker_db$organ),]
genes = intersect(toupper(IDC_DA$gene), cell_marker_db$official.gene.symbol)

cell_marker_db = cell_marker_db[which(cell_marker_db$official.gene.symbol %in% genes),]
table(cell_marker_db$cell.type)

IDC_DA$cell_type = cell_marker_db$cell.type[match(toupper(IDC_DA$gene), cell_marker_db$official.gene.symbol)]
IDC_DA = IDC_DA %>% mutate(
  rank_group_activation = order(pct.1),
  rank_logFC = order(avg_log2FC),
  rank_qval = order(p_val_adj),
)
cell_type_markers = IDC_DA %>% dplyr::filter(!is.na(cell_type)) %>% 
  group_by(cluster_of_origin, cluster, cell_type) %>%
  summarise(total = n(),
            avg_rank_group_activation = mean(rank_group_activation),
            avg_rank_logFC = mean(rank_logFC),
            avg_rank_qval = mean(rank_qval),
            Gene = paste(gene, collapse = ",")
  ) %>% 
  mutate(
    combined_score_rank = nrow(IDC_DA) / ((avg_rank_group_activation + avg_rank_logFC + avg_rank_qval) / 3)
  ) %>% slice_max(total, with_ties = TRUE) %>% slice_max(combined_score_rank)
cell_type_markers = cell_type_markers %>% tidyr::separate_rows(Gene, sep = ",")

for(i in unique(cell_type_markers$Gene)){
  
  print(plot_cluster_network(Seu, 
                             IDC_summary = IDC_summary,
                             color_by = i,
                             legend = F, 
                             threshold_to_define_feature_active = 1
  )
  )
}
dev.off() 

pdf(file.path(seurat_dir, paste0("IDC_tree_colored_by_Annotation_wTop_Markers_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(unlist(strsplit(topmarkers$Term, split = " ")))){
  plot_cluster_network(Seu,  
                       IDC_summary = IDC_summary, 
                       color_by = i, 
                       legend = F, 
                       threshold_to_define_feature_active = 1 
  )
}
dev.off()

pdf(file.path(seurat_dir, paste0("Markers_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in unique(unlist(strsplit(topmarkers$Term, split = " ")))){
  print( FeaturePlot(Seu,  features = i))
}
dev.off()

pdf(file.path(seurat_dir, paste0("Cell_type_Markers_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in c("KRT1","CD4","CD8A","IL17A","IL17B","IL17D","IL17F", "CD3E")){
  print( FeaturePlot(Seu,  features = i))
}
dev.off()


pdf(file.path(seurat_dir, paste0("IL17_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in c("CD4","CD8A","IL17A","IL17B","IL17D","IL17F", "CD3E","KRT14")){
  print( FeaturePlot(Seu,  features = i))
}
dev.off()

qs::qsave(Seu, file.path(output, "Seu.qs"))



library(IDclust)
list_res = list()
for(i in c( "FrK",  "RUW", "SchH", "WaG" )){
  Seu. = Seu[,which(Seu$patient == i)]
  if(length(Seu$cell_id[Seu$Type == "Healthy_Tcell"]) >= 5){
    res = IDclust::differential_edgeR_pseudobulk_LRT(Seu.,
                                                     by = "Type",
                                                     assay = "originalexp",
                                                     logFC.th = 0.5, qval.th = 0.1, min.pct = 0.2)
    list_res[[i]] = res[which(res$cluster == "Tumor_Tcell"),]
  }  
}

WriteXLS::WriteXLS(list_res, ExcelFileName = file.path(seurat_dir, "..", "Tumor_Tcell_vs_Heatlhy_T_cell_by_patient.xlsx"),
                   SheetNames = names(list_res))

########################################
# Integrating with psoriasis
########################################
library(harmony)
Seu_psoriasis = qs::qread(file.path(maindir, "../Darrier_Disease_sc/output/Seu_1.qs"))
common_genes = intersect(rownames(Seu_psoriasis), rownames(Seu))
Seu_combined =  merge(Seu, y = Seu_psoriasis, add.cell.ids = c("Darrier", "Psoriasis"), project = "RNA")
Seu_combined$batch = gsub("_.*","", colnames(Seu_combined))
Seu_combined$batch[Seu_combined$batch == "Psoriasis"] = "Psoriasis/Healthy"

cell_markers = intersect(cell_markers, rownames(Seu_combined))

Idents(Seu_combined) = Seu_combined$orig.ident
Seu_combined <- NormalizeData(Seu_combined,  verbose = FALSE)
Seu_combined <- FindVariableFeatures(Seu_combined, selection.method = "vst")
Seu_combined <- ScaleData(Seu_combined, features = rownames(Seu_combined))
Seu_combined <- RunPCA(Seu_combined, features = VariableFeatures(Seu_combined), ndims.print = 1:10, nfeatures.print = 10)

Seu_combined <-  RunHarmony(Seu_combined, group.by.vars = "batch")
Seu_combined <- FindNeighbors(Seu_combined, reduction = "harmony", dims = 1:50)
Seu_combined = RunUMAP(Seu_combined, reduction = "harmony", dims = 1:50)
Seu_combined <- FindClusters(Seu_combined)

Seu_combined_T_cells = Seu_combined[,Seu_combined$seurat_clusters == 6]
Seu_combined_T_cells <- FindNeighbors(Seu_combined_T_cells, reduction = "harmony", dims = 1:50)
Seu_combined_T_cells <- FindClusters(Seu_combined_T_cells, resolution = 0.1)
Seu_combined_T_cells$seurat_clusters = as.numeric(Seu_combined_T_cells$seurat_clusters)
Seu_combined_T_cells$seurat_clusters[Seu_combined_T_cells$seurat_clusters == 1] = 29
Seu_combined_T_cells$seurat_clusters[Seu_combined_T_cells$seurat_clusters == 2] = 30
DimPlot(Seu_combined_T_cells, group.by = "seurat_clusters")

Seu_combined$seurat_clusters = as.numeric(Seu_combined$seurat_clusters)
Seu_combined$seurat_clusters[colnames(Seu_combined) %in% colnames(Seu_combined_T_cells)] = Seu_combined_T_cells$seurat_clusters
Seu_combined$seurat_clusters = as.factor(Seu_combined$seurat_clusters )
table(Seu_combined$seurat_clusters)
DimPlot(Seu_combined, group.by = "seurat_clusters", cols = sample(rainbow(30)))

seurat_dir_combined = file.path(output, "seurat_dir_combined")
if(!dir.exists(seurat_dir_combined)) dir.create(seurat_dir_combined)

qs::qsave(Seu_combined, file.path(seurat_dir_combined, "Seu_combined.qs"))

pdf(file.path(seurat_dir_combined,paste0("scRNA_combined.pdf")))
print(DimPlot(Seu_combined, reduction = "umap", group.by = "seurat_clusters"))
print(DimPlot(Seu_combined, reduction = "umap", group.by = "IDcluster"))
print(DimPlot(Seu_combined, reduction = "umap", group.by = "orig.ident"))
print(DimPlot(Seu_combined, reduction = "umap", group.by = "batch"))
print(DimPlot(Seu_combined, reduction = "umap", group.by = "Phase"))
print(FeaturePlot(Seu_combined,  features = "nCount_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
print(FeaturePlot(Seu_combined,  features = "nFeature_RNA",reduction = "umap", cols = rev(viridis::magma(100))))
dev.off()

png(file.path(seurat_dir_combined,paste0("UMAP_sample_All.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(Seu_combined, reduction = "umap", pt.size = 1.5, cols = colors$cell_type_color[c(8:9,6, 10,5,7, 4,3,2,1)],
                group.by = "orig.ident", order = order(Seu_combined$orig.ident, decreasing = F)) + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

for(group in c("Ctrl", "Pso", "DAR")){
  Seu_combined. = Seu_combined[,grep(group, Seu_combined$orig.ident)]
  if(group == "Ctrl") color = colors$cell_type_color[c(6,8,9)]
  if(group == "Pso") color = colors$cell_type_color[c(5,7,10)]
  if(group == "DAR") color = colors$cell_type_color[c(1:4)]
  png(file.path(seurat_dir_combined,paste0("UMAP_sample_", group,".png")), width = 1400, height = 1200, res = 250)
  tplot = DimPlot(Seu_combined., reduction = "umap", pt.size = 1.5, cols = color,
                  group.by = "orig.ident") + ggtitle("")
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
  dev.off()
}

png(file.path(seurat_dir_combined,paste0("UMAP_batch.png")), width = 1400, height = 1200, res = 250)
tplot = DimPlot(Seu_combined, reduction = "umap", pt.size = 1.5, cols = c("grey", "royalblue"),
                group.by = "batch", order = order(Seu_combined$batch, decreasing = F)) + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

pdf(file.path(seurat_dir_combined, paste0("Cell_type_Markers_expression.pdf")), height = 10,
    width = 10, pointsize = 12)
for(i in cell_markers){
  print( FeaturePlot(Seu_combined,  features = i, cols = c("grey", "red")))
}
dev.off()
THEMIS
print( FeaturePlot(Seu_combined,  features = "PTPRC", cols = c("grey", "red")))

png(file.path(seurat_dir_combined,paste0("Heatmap_cell_markers.png")), width = 1600, height = 1200, res = 250)
Seurat::DoHeatmap(Seu_combined, features = cell_markers)
dev.off()

map_cluster_celltype = unlist(list(
  "Epidermis_" = c(1, 7, 8, 10, 13, 14, 15, 21, 22) + 1 ,
  "Follicular epidermis_" = c(20, 26, 27) + 1,
  "T helper cells_" = c(29),
  "Cytotoxic T cells_" = c(30),
  "Macrophages and monocytes_" = c(12),
  "Mast cells_" = c(17),
  "Fibroblasts_" = c(0, 4, 9, 17, 18)+ 1,
  "Sebaceous glands_" = c(3),
  "Melanocytes_" = c(24),
  "Schwann cells_" = c(26),
  "Endothelial cells_" = c(3, 5, 12)+ 1,
  "Lymphatic vessel cells_" = c(19)+ 1,
  "Smooth muscle cells_" = c(24)+ 1
))

Seu_combined$celltype = gsub("_.*", "", names(map_cluster_celltype))[match(Seu_combined$seurat_clusters, map_cluster_celltype)]

col_to_celltype = read.csv("output/col_to_celltype.csv")

Seu_combined$celltype_color = col_to_celltype[match(Seu_combined$celltype, names(col_to_celltype))]
png(file.path(seurat_dir_combined,paste0("UMAP_cell_type_All.png")), width = 1800, height = 1200, res = 250)
tplot = DimPlot(Seu_combined, reduction = "umap", pt.size = 0.75,
                cols = unique(Seu_combined$celltype_color[order(Seu_combined$celltype)]),
                group.by = "celltype") + ggtitle("")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
tplot
dev.off()

for(group in c("Ctrl", "Pso", "DAR")){
  Seu_combined. = Seu_combined[,grep(group, Seu_combined$orig.ident)]
  png(file.path(seurat_dir_combined,paste0("UMAP_celltype_", group,".png")), width = 1800, height = 1200, res = 250)
  tplot = DimPlot(Seu_combined., reduction = "umap", pt.size = 0.75,
                  cols =  unique(Seu_combined$celltype_color[order(Seu_combined$celltype)]),
                  group.by = "celltype") + ggtitle("")
  tplot[[1]]$layers[[1]]$aes_params$alpha = 0.7
  print(tplot)
  dev.off()
}

# colors$cell_type_color[c(F1:8,15,12,14,10,13)]

png(file.path(seurat_dir_combined,paste0("Heatmap_cell_markers_cell_types.png")), width = 3200, height = 1600, res = 200)
Seurat::DoHeatmap(Seu_combined, features = cell_markers, group.by = "celltype")
dev.off()

png(file.path(seurat_dir_combined,paste0("DotPlot_cell_markers_cell_types.png")), width = 3200, height = 1400, res = 250)
Seurat::DotPlot(Seu_combined, features = cell_markers, group.by = "celltype", cols = c("grey", "gold")) + theme(axis.text.x =  element_text(angle = 90))
dev.off()

pdf(file.path(seurat_dir_combined, paste0("Cell_type_Markers_expression_violin.pdf")), height = 10,
    width = 15, pointsize = 12)
for(i in cell_markers){
  print( VlnPlot(Seu_combined,  features = i,  pt.size = 1.5, group.by = "seurat_clusters"))
}
dev.off()

Seu_combined$IL17A.F = colSums(Seu_combined[c("IL17A","IL17F"),])
png(file.path(seurat_dir_combined,paste0("IL17A_F_All.png")), width = 2800, height = 1400, res = 250)
print(VlnPlot(Seu_combined,  features = "IL17A.F",  pt.size = 2,  group.by = "celltype")) + 
  ylab("Sum of IL17A-IL17F normalized expression") + xlab("") + ggtitle("") +  NoLegend()
dev.off()

for(group in c("Ctrl", "Pso", "DAR")){
  Seu_combined. = Seu_combined[,grep(group, Seu_combined$orig.ident)]
  Seu_combined.$IL17A.F = colSums(Seu_combined.[c("IL17A","IL17F"),])
  png(file.path(seurat_dir_combined, paste0("IL17A_F_", group, ".png")), width = 2800, height = 1400, res = 250)
  print(VlnPlot(Seu_combined.,  features = "IL17A.F", pt.size = 2, group.by = "celltype") + 
          ylab("Sum of IL17A-IL17F normalized expression") + xlab("") + ggtitle("") + 
          NoLegend())
  dev.off()
}

png(file.path(seurat_dir_combined,paste0("IL17A_expression.png")), width = 3200, height = 1400, res = 250)
print(VlnPlot(Seu_combined,  features = "IL17A", group.by = "celltype")) + NoLegend()
dev.off()

png(file.path(seurat_dir_combined,paste0("IL17F_expression.png")), width = 3200, height = 1400, res = 250)
print(VlnPlot(Seu_combined,  features = "IL17F", group.by = "celltype")) + NoLegend()
dev.off()

png(file.path(seurat_dir_combined,paste0("IL17A_expression_DarrierOnly.png")), width = 3200, height = 1400, res = 250)
print(VlnPlot(Seu_combined[,Seu_combined$batch == "Darrier"],  features = "IL17A", group.by = "celltype")) + NoLegend()
dev.off()

png(file.path(seurat_dir_combined,paste0("IL17F_expression_DarrierOnly.png")), width = 3200, height = 1400, res = 250)
print(VlnPlot(Seu_combined[,Seu_combined$batch == "Darrier"],  features = "IL17F", group.by = "celltype")) + NoLegend()
dev.off()

library(IDclust)
res = IDclust::differential_edgeR_pseudobulk_LRT(Seu_combined,
                                                 by = "seurat_clusters",
                                                 assay = "RNA",
                                                 logFC.th = 0.5, qval.th = 0.1, min.pct = 0.2)

WriteXLS::WriteXLS(res, ExcelFileName = file.path(seurat_dir_combined, "DA_Markers.xlsx"))






# Composition analysis ---------------------------------------------------------
library(dittoSeq)
pdf(file.path(output_dir, "Composition_analysis.pdf"))
dittoBarPlot(Seu_combined, "batch", group.by = "orig.ident", color.panel = colors$cell_type_color, main = "Cell Count", scale = "count")  
dittoBarPlot(Seu_combined, "celltype", group.by = "orig.ident", color.panel = colors$cell_type_color, main = "Condition x Cell Type")  
dittoBarPlot(Seu_combined, "condition", group.by = "condition", scale = "count",  color.panel = unique(Seu_combined$condition_color[order(Seu_combined$condition)]))
dev.off()

g <- ggplot_build(tplot)
col_to_celltype = unique(cbind(g$plot$data$celltype, g$data[[1]]["colour"]))
col_to_celltype = setNames(col_to_celltype$colour, col_to_celltype$`g$plot$data$celltype`)

Seu_combined$condition = gsub("[0-9]*", "",Seu_combined$orig.ident)
d = dittoBarPlot(Seu_combined, "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  
png(file.path(seurat_dir_combined, "Sample x CellType x Condition.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(Seu_combined, "celltype", group.by = "orig.ident", color.panel = col_to_celltype[levels(d$data$label)], main = "Condition x Cell Type",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()

Seu_combined_immune = Seu_combined[, (Seu_combined$celltype %in% c("T helper cells", "Cytotoxic T cells", "Mast cells", "Macrophages and monocytes"))]
d = dittoBarPlot(Seu_combined_immune, "celltype", group.by = "orig.ident", color.panel = col_to_celltype, main = "Condition x Cell Type")  
png(file.path(seurat_dir_combined, "Sample x CellType x Condition Immune.png"), width = 3000, height = 1600, res = 300)
dittoBarPlot(Seu_combined_immune, "celltype", group.by = "orig.ident", color.panel = col_to_celltype[levels(d$data$label)], main = "Condition x Cell Type",
             split.adjust = list(scales = 'free'), split.by = "condition", split.nrow = 1)  + xlab("")
dev.off()



