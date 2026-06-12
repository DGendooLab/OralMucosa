#### libs

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(optparse)
library(ComplexHeatmap)
library(clustree)

###functs
optimize_kparam <- function(seurat_obj, 
                            dims = 1:20, 
                            k_values = c(10, 20, 30, 50), 
                            resolution = 0.1, 
                            assay = "RNA", 
                            slot = "data",
                            do_plot = TRUE) {
  
  require(Seurat)
  require(cluster)
  require(dplyr)
  require(ggplot2)
  
  results <- list()
  metrics <- data.frame(k.param = integer(), 
                        n_clusters = integer(), 
                        avg_silhouette = numeric())
  
  for (k in k_values) {
    message("Testing k.param = ", k)
    
    # Copy to avoid modifying original object
    obj_copy <- seurat_obj
    
    # Run neighborhood graph and clustering
    obj_copy <- FindNeighbors(obj_copy, dims = dims, k.param = k, verbose = FALSE)
    obj_copy <- FindClusters(obj_copy, resolution = resolution, verbose = FALSE)
    
    clusters <- obj_copy$seurat_clusters
    results[[paste0("k", k)]] <- clusters
    
    # Calculate silhouette width
    pca <- Embeddings(obj_copy, "pca")[, dims]
    sil <- cluster::silhouette(as.integer(clusters), dist(pca))
    avg_sil <- mean(sil[, 3])
    
    # Store metrics
    metrics <- rbind(metrics, data.frame(k.param = k, 
                                         n_clusters = length(unique(clusters)), 
                                         avg_silhouette = avg_sil))
  }
  
  if (do_plot) {
    p1 <- ggplot(metrics, aes(x = k.param, y = n_clusters)) + 
      geom_line() + geom_point() +
      labs(title = "Number of Clusters vs k.param", y = "Number of Clusters", x = "k.param")
    
    p2 <- ggplot(metrics, aes(x = k.param, y = avg_silhouette)) + 
      geom_line() + geom_point() +
      labs(title = "Average Silhouette vs k.param", y = "Avg. Silhouette Width", x = "k.param")
    
    print(p1)
    print(p2)
  }
  
  return(list(clusterings = results, metrics = metrics))
}





optimize_resolution_by_silhouette <- function(seurat_obj, 
                                              resolutions = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
                                              dims = 1:min_pc,
                                              verbose = TRUE) {
  require(Seurat)
  require(cluster)
  require(dplyr)
  require(ggplot2)
  
  # Make a copy to avoid overwriting the original object
  sobj <- seurat_obj
  
  # Build SNN graph
  sobj <- FindNeighbors(sobj, dims = dims, verbose = verbose)
  
  # Run clustering at multiple resolutions
  for (res in resolutions) {
    sobj <- FindClusters(
      sobj, resolution = res, verbose = verbose,
      algorithm = 1,               # optional, can be changed
      group.singletons = FALSE    # optional
    )
    # Rename the active identity column
    cluster_col <- paste0("res_", res)
    sobj[[cluster_col]] <- sobj@meta.data$seurat_clusters
  }
  
  # Function to compute silhouette from clustering column
  compute_silhouette <- function(obj, cluster_col) {
    pca_mat <- Embeddings(obj, "pca")[, dims]
    clusters <- obj[[cluster_col]][,1]
    if (length(unique(clusters)) < 2) return(NA)
    sil <- cluster::silhouette(as.integer(as.factor(clusters)), dist(pca_mat))
    mean(sil[, 3])
  }
  
  # Compute silhouette score for each resolution
  sil_scores <- sapply(resolutions, function(res) {
    col_name <- paste0("RNA_snn_res.", res)
    compute_silhouette(sobj, col_name)
  })
  
  sil_df <- data.frame(resolution = resolutions, silhouette = sil_scores)
  
  # Generate plot
  p <- ggplot(sil_df, aes(x = resolution, y = silhouette)) +
    geom_line() + geom_point(size = 2) +
    theme_minimal() +
    labs(title = "Silhouette Score vs Resolution",
         x = "Resolution", y = "Avg Silhouette Width")
  
  print(p)
  
  # Choose best resolution
  best_res <- sil_df$resolution[which.max(sil_df$silhouette)]
  if (verbose) {
    cat("✅ Best resolution by silhouette:", best_res,
        "(score =", round(max(sil_df$silhouette), 3), ")\n")
  }
  
  list(
    seurat = sobj,
    silhouette_df = sil_df,
    best_resolution = best_res,
    silhouette_plot = p
  )
}



get_optimal_pc_count <- function(seurat_obj, variance_cutoff = 90, percent_drop = 0.1) {
  stdv <- seurat_obj[["pca"]]@stdev
  sum_stdv <- sum(stdv)
  percent_stdv <- (stdv / sum_stdv) * 100
  cumulative <- cumsum(percent_stdv)
  
  # First condition: capture >90% variance but avoid very flat PCs
  co1 <- which(cumulative > variance_cutoff & percent_stdv < 5)[1]
  
  # Second condition: dropoff in variance explained > 0.1%
  diffs <- percent_stdv[1:(length(percent_stdv) - 1)] - percent_stdv[2:length(percent_stdv)]
  co2 <- sort(which(diffs > percent_drop), decreasing = TRUE)[1] + 1
  
  # Fallback if one of them is NA
  if (is.na(co1)) co1 <- Inf
  if (is.na(co2)) co2 <- Inf
  
  min_pc <- min(co1, co2)
  
  return(min_pc)
}


option_list <- list(
  make_option(c("-s", "--species"),
              type = "character",
              default = "human",
              help = "species to use [default = %default]. Options: human or mouse",
              metavar = "character"),
  make_option(c("-i", "--input"),
              type = "character",
              default = ,
              help = "input RDS dir to load",
              metavar = "character"),
  make_option(c("-n", "--n_cores"),
              type = "integer",
              default = 1,
              help = "Number of cores to use to run  [default = %default]",
              metavar = "character"),
  make_option(c("-c", "--cachedir"),
              type = "character",
              default = 1,
              help = "Path to save intermediate RDS caches [default = %default]",
              metavar = "character"),
  make_option(c("-o", "--outdir"),
              type = "character",
              default = 1,
              help = "Path to save results [default = %default]",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)

n_cores <- opt$n_cores
cache_dir <- opt$cachedir
outdir <- opt$outdir
species <- opt$species
input <- opt$input

# Store selected species
species <- tolower(opt$species)

species <- "mouse"
# Validate species input
if (!species %in% c("human", "mouse")) {
  stop("Invalid species. Use 'human' or 'mouse'.")
}

out_prefix_dir <- file.path(outdir, "Seperate_samples")

proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
analysis_dir <- file.path(proj_dir, "scRNAseqAnalysis/")
git_dir <- file.path(analysis_dir, "OralMucosa/HGK_HTK_VU40T")
out_prefix_dir <- file.path(git_dir, "Seperate_samples")
cache_dir <- file.path(proj_dir, "rds_cache")


###main
if (is.null(input)){
  seurat_singlets_list <- readRDS(file.path(cache_dir, paste0("VU40T_HTK_HGK_singlets_only_sep_samples_",species,".RDS")))
  } else {
    seurat_singlets_list <- input
  ### sanity check if using opt$input arg to make sure species matches seurat list RDS
  if (any(sapply(seurat_singlets_list, function(obj) any(grepl("^ENSMUS", rownames(obj[["RNA"]])))))) {
    species <- "mouse"
  } else {
    species <- "human"
  }
}



plots_dir <- file.path(out_prefix_dir, paste0(stringr::str_to_title(species), "/Plots"))
res_dir <- file.path(out_prefix_dir, stringr::str_to_title(species))

chk_dir_list <- list(plots_dir, res_dir)

for (path in chk_dir_list){
  if(!(dir.exists(path))){
    dir.create(path, recursive = T)
  }
}

species_df <- read.csv(file.path(git_dir, "/VU40T_species_assignment_by_reads.csv"))[,-c(1)]

# get rid of samples that didn't work for species

if (species == "mouse"){
  samples_to_exclude <- c("MW_HGK")
} else {
  samples_to_exclude <- c("MW_HTK2_HPV16", "MW_VU40T")
}


# 3. Remove them from your list
seurat_singlets_list[samples_to_exclude] <- NULL

message("Dropped samples: ", paste(samples_to_exclude, collapse = ", "))


for (sample in names(seurat_singlets_list)) {
  
  # Get the Seurat object for the current sample
  seurat_obj <- seurat_singlets_list[[sample]]
  
  # Filter species_df for this sample
  df_sample <- subset(species_df, sample == sample)
  
  # Make sure barcodes match
  df_sample <- df_sample[!duplicated(df_sample$barcode), ]
  rownames(df_sample) <- df_sample$barcode
  
  # Match the species to the Seurat object's barcodes
  matched_species <- df_sample[colnames(seurat_obj), "species_label"]
  
  # Add to metadata
  seurat_obj$species <- matched_species
  
  # Store back in list
  seurat_singlets_list[[sample]] <- seurat_obj
}


if(species == "mouse"){
  
  for (i in seq_along(seurat_singlets_list)) {
    
    seurat_obj <- seurat_singlets_list[[i]]
    
    # Keep only cells from dominant species
    cells_to_keep <- colnames(seurat_obj)[seurat_obj$species == "Mouse"]
    
    # Subset Seurat object
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    
    # Store it back
    seurat_singlets_list[[i]] <- seurat_obj
  }
} else {
  for (i in seq_along(seurat_singlets_list)) {
    
    seurat_obj <- seurat_singlets_list[[i]]
    
    # Keep only cells from dominant species
    cells_to_keep <- colnames(seurat_obj)[seurat_obj$species == "Human"]
    
    # Subset Seurat object
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    
    # Store it back
    seurat_singlets_list[[i]] <- seurat_obj
  }
}

pc_summary <- data.frame(sample = character(), optimal_pcs = integer(), stringsAsFactors = FALSE)

for (sample_name in names(seurat_singlets_list)) {
  seurat_obj <- seurat_singlets_list[[sample_name]]
  optimal_pc <- get_optimal_pc_count(seurat_obj)
  
  pc_summary <- rbind(pc_summary, data.frame(
    sample = sample_name,
    optimal_pcs = optimal_pc
  ))
}

for (id in names(seurat_singlets_list)){
  print(names(seurat_singlets_list[id]))
  print(ncol(seurat_singlets_list[[id]]))
}

for (id in names(seurat_singlets_list)){
  seurat_obj <- seurat_singlets_list[[id]]
  print(ElbowPlot(seurat_obj, ndims = 50))
  
  
  # Run JackStraw (slow - use num.replicate=100 as a perm subsample)
  seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
  
  # Score PCs:
  
  seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
  
  seurat_singlets_list[[id]] <- seurat_obj
}
for (id in names(seurat_singlets_list)){
  seurat_obj <- seurat_singlets_list[[id]]
  sample_id <- names(seurat_singlets_list[id])
  
  print(sample_id)
  print(ElbowPlot(seurat_obj, ndims = 50))
  png(file = file.path(plots_dir, paste0(sample_id, "_Jackstraw_integrated.png")), 
      width = 16, height = 12, units = "in", res = 300)

  print(JackStrawPlot(seurat_obj, dims = 1:20))
  dev.off()
}

## optimum PCs: 20, 20, 12 (HTK_HPV18)

optimumPC_df <- data.frame(
  sample = c("MW_HTK_HPV18", "MW_HTK2_HPV16", "MW_VU40T"),
  nPCs = c(12,20,20))

# for (i in seq_along(seurat_singlets_list)){
#   sample_name <- names(seurat_singlets_list)[i]
#   min_pc <- pc_summary$optimal_pcs[pc_summary$sample == sample_name]
#   k_param_res[[i]] <- optimize_kparam_by_silhouette(seurat_singlets_list[[i]], k_values = c(0.1, 0.25, 0.5, 0.75, 1.0), dims = 1:min_pc)
# }

#}
# `

## rerun clustering with optimised dims

for (i in seq_along(seurat_singlets_list)){
  current_sample <- names(seurat_singlets_list)[i]
  
  # 2. Extract the specific nPCs value cleanly as a numeric scalar
  min_pc <- optimumPC_df$nPCs[optimumPC_df$sample == current_sample]
  
  # Safety check: Ensure a matching sample was actually found in your dataframe
  if (length(min_pc) == 1) {
    # 3. Run FindNeighbors with the dynamic PC cutoff
    seurat_singlets_list[[i]] <- FindNeighbors(seurat_singlets_list[[i]], dims = 1:min_pc)
    # 4. Run FindClusters (uncommented, ensuring 'resolution' is defined globally)
    # seurat_singlets_list[[i]] <- FindClusters(seurat_singlets_list[[i]], resolution = resolution)
  } else {
    warning(paste("Sample", current_sample, "not found or ambiguous in optimumPC_df. Skipping..."))
  }
  
  resolutions <- seq(0.1, 1.5, by = 0.1)
    
  for (res in resolutions) {
    seurat_singlets_list[[i]] <- FindClusters(
      seurat_singlets_list[[i]],
      resolution = res,
      verbose = FALSE, 
      random.seed = 666
    )

}


  
  # Save clustering results manually
  seurat_singlets_list[[i]][[paste0("RNA_snn_res.", res)]] <- Idents(seurat_singlets_list[[i]])
  
}

for (id in names(seurat_singlets_list)){
  
  seurat_obj <- seurat_singlets_list[[id]]
  
  sample_id <- names(seurat_singlets_list[id])
  
  # 2. Extract the specific nPCs value cleanly as a numeric scalar
  min_pc <- optimumPC_df$nPCs[optimumPC_df$sample == sample_id]
  
  png(file = file.path(plots_dir, paste0(sample_id,"_",min_pc, "_PCs_ClustTree_VU40T_stabilityscores_integrated.png")), width = 16, height = 12, units = "in", res = 300)
    p1 <- clustree(seurat_obj, prefix = "RNA_snn_res.", node_colour = "sc3_stability") ## best resolution lies between 0.6-0.8
    print(p1)
    dev.off()
    
  png(file = file.path(plots_dir, paste0(sample_id,"_",min_pc, "_PCs_ClustTree_VU40T_integrated.png")), width = 16, height = 12, units = "in", res = 300)
    p1 <- clustree(seurat_obj, prefix = "RNA_snn_res.") ## best resolution lies between 0.6-0.8
    print(p1)
    dev.off()
}

#optimum res hpv18 possibly 0.3,  hpv16 = 0.4, VU40T = 0.

for (i in seq_along(seurat_singlets_list)){
  sample_name <- names(seurat_singlets_list)[i]
  
  plot1 <- DimPlot(seurat_singlets_list[[i]], reduction = "umap", label = T)
  
  
  png(file = file.path(plots_dir, paste0 ("SingletsOnly_UMA`P_and_t-sne_", sample_name, ".png")), width = 6, height = 5, units = "in", res = 300, )
  
  print(
    plot1 +
      patchwork::plot_annotation(title = paste0("Doublet-Filtered ", species, " only UMAP and t-SNE:\n", sample_name))
  )
  
  dev.off()
}

## run manual annotation of all samples as before

## genelist 3 -> mouse and human epithelial, fibroblast and EMT markers
genelist <- as.data.frame(read.csv("~/Markers_for_dotplots_2_ep_2nd_set.csv", header = T, skip = 1))
genelist <- as.data.frame(lapply(genelist, toupper))
genelist <- lapply(genelist, function(x) gsub("CDH2", "", x))

if (species == "mouse") {
  # 1. Subset the list to keep only the 4th and 5th elements
  genelist <- genelist[c(4, 5)] 
  
  names(genelist) <- gsub("_[HM]", "", names(genelist))
  
  # 3. Clean the gene strings inside the vectors to be mouse-standard Title Case
  genelist <- lapply(genelist, stringr::str_to_title)
}else{
  genelist <- genelist[, c(1:2)] ## both human and mouse markers, no need for EMT here
  colnames(genelist) <- gsub("_[HM]", "", colnames(genelist))
}
names(genelist) <- gsub("_", " ", names(genelist))
names(genelist) <- gsub("markers", "Markers", names(genelist))
## convert to named vectors
genevectors <- lapply(genelist, function(x) unique(x[x != "" & !is.na(x)]))

if(species == "mouse"){
  genevectors <- lapply(genevectors, stringr::str_to_title)
}

for (sample in names(seurat_singlets_list)){
  seurat_obj <- seurat_singlets_list[[sample]]
  
  png(file = file.path(
    plots_dir,
    paste0("CombinedDotplot_epi_fibro_Markers_clusterResolution",species,"_",sample,".png")),
    width = 8, height = 4, units = "in", res = 300)
  
  p <- DotPlot(
    seurat_obj,
    features = genevectors,
    group.by = "seurat_clusters",
  ) +
    RotatedAxis() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = "italic")  # Optional: italics for gene names
    ) + labs(y = "Clusters")  # 👈 Relabel y-axis
  print(p)  
  dev.off()
}

# can be certain all mouse clusters are fibroblast

# now running caf list

for (sample in names(seurat_singlets_list)){
  seurat_obj <- seurat_singlets_list[[sample]]
  if (all(seurat_obj$species == "Mouse")){
  
  
  ## make stacked barplot of num fibroblasts per clust
  # clust_cell_df <- as.data.frame(table(seurat_obj$seurat_clusters))
  # colnames(clust_cell_df) <- c("Cluster", "No. Cells")
  # write.csv(clust_cell_df, file = file.path(git_dir, "Integrated/Mouse/MouseOnly_fibroblast_cells_per_clust.csv"))
  # clust_cell_df$`Cell Type` <- "Fibroblast"
  # ## plot barplot
  # 
  # png(file = file.path(
  #   plots_dir,
  #   "Integrated/Mouse/Plots/MouseONLY_Stacked_Bar_numCells_per_cluster_Resolution",optimum_res,"_VU40T_combined.png"),
  #   width = 4, height = 6, units = "in", res = 300)
  # 
  # p <- clust_cell_df %>% 
  #   ggplot(aes(x = `Cell Type`, y = `No. Cells`, fill = Cluster)) + 
  #   geom_bar(position="stack", stat="identity") + 
  #   theme_minimal() + xlab("")
  # print(p)
  # dev.off()
  
  # png(file = file.path(
  #   git_dir,
  #   "Integrated/Mouse/Plots/MouseONLY_Stacked_Bar_percent_Cells_per_cluster_Resolution0.6_VU40T_combined.png"),
  #   width = 4, height = 6, units = "in", res = 300)
  # p <- clust_cell_df %>% 
  #   ggplot(aes(x = `Cell Type`, y = `No. Cells`, fill = Cluster)) + 
  #   geom_bar(position="fill", stat="identity") + 
  #   theme_minimal() +
  #   ylab("Proportion of Cells") + xlab("") + 
  #   scale_y_continuous(labels = scales::percent)
  # print(p)
  # dev.off()
  # genelist <- read.csv("~/Fibroblasts_dotplots_markers.csv", header = F)
  # names(genelist) <- "Fibroblast Markers"
  # 
  # png(file = file.path(
  #   git_dir,
  #   "Integrated/Mouse/Plots/MouseONLY_MarkerDotplots/Dotplot_Fibro_Markers_clusterResolution0.6_VU40T_combined.png"),
  #   width = 12, height = 4, units = "in", res = 300)
  # 
  # p <- DotPlot(
  #   seurat_obj,
  #   features = genelist,
  #   group.by = "seurat_clusters",
  # ) +
  #   scale_x_discrete(labels = function(x) stringr::str_to_title(x)) +
  #   RotatedAxis() +
  #   scale_color_gradient(low = "lightgrey", high = "red") +
  #   theme(
  #     axis.text.x = element_text(angle = 45, hjust = 1),
  #     axis.text.y = element_text(face = "italic")  # Optional: italics for gene names
  #   ) + labs(y = "Clusters")  # 👈 Relabel y-axis
  # print(p)
  # dev.off()
  
  ## fibroblast heatmap
  # heatmap_markers <- read.csv("~/fibroblast_heatmap.csv", header = T)
  fibro2_marker_ls <- list()
  for (i in seq(1,6)){
    fibro2_marker_ls[[i]] <- readxl::read_xlsx(file.path(proj_dir ,"/scRNAseqAnalysis/Fibroblasts_dotplots.xlsx"), sheet = i)
  }
  heatmap_markers <- fibro2_marker_ls[[2]]
  ### for heatmap
  # heatmap_markers <- genelist[, c(1:5)]
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "_", " ")
  ## for fibro list 2
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "Oxid phosph", "Oxidative phosphorylation")
  colnames(heatmap_markers) <- stringr::str_to_title(colnames(heatmap_markers))
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "Ecm", "ECM")
  colnames(heatmap_markers) <- stringr::str_replace(colnames(heatmap_markers), "To", "to")
  heatmap_markers <- as.data.frame(apply(heatmap_markers, 2 , stringr::str_to_title))
  marker_long <- heatmap_markers %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Group", values_to = "Gene") %>%
    dplyr::filter(!is.na(Gene) & Gene != "")
  # Filter for genes that are present
  marker_long <- marker_long %>% dplyr::filter(Gene %in% rownames(seurat_obj))
  
  ### line wrapping func.
  marker_long$Group <- sapply(marker_long$Group, function(x) {
    if (nchar(x) > 15 && grepl(" ", x)) {
      gsub(" ", "\n", x)  # replace . or _ with a line break
    } else {
      x  # keep as is
    }
  })
  # Compute average expression
  avg_expr <- AverageExpression(seurat_obj, features = unique(marker_long$Gene), return.seurat = FALSE, assays = "RNA")$RNA
  avg_expr@Dimnames[[2]] <- stringr::str_replace(avg_expr@Dimnames[[2]], "g", "")
  # Order genes as in the original df
  marker_long <- marker_long %>%
    distinct(Gene, .keep_all = TRUE)  # drop duplicate entries
  
  # Reorder avg_expr to match marker order
  avg_expr <- avg_expr[marker_long$Gene, ]
  
  scaled_expr <- t(scale(t(avg_expr)))
  # Clip values to range [-2, 2] for visual clarity
  scaled_expr[scaled_expr > 2] <- 2
  
  grp <- as.character(marker_long$Group)
  lev <- sort(unique(grp))
  set.seed(123)
  pal <- setNames(
    circlize::rand_color(length(lev)),
    lev
  )
  
  
  row_annot <- rowAnnotation(
    Group = grp,                        
    annotation_legend_param = list(
      Group = list(title = "Group") # legend title
    ))
  
  # # Annotation (for row group labels)
  
  # Select 3-color palette from RColorBrewer
  color_palette <- RColorBrewer::brewer.pal(n = 3, name = "YlOrRd")
  
  # Build color function
  col_fun <- circlize::colorRamp2(
    breaks = c(-2, 0, 2),
    colors = color_palette
  )
  png(file = file.path(
    plots_dir,
    paste0("Heatmap_FibroMarkers_set2_",sample,".png")),
    width = 8, height = 14, units = "in", res = 300)
  p <- Heatmap(
    matrix = as.matrix(scaled_expr),
    name = "Z-scaled Expr",
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = T,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 0,
    row_split = marker_long$Group,      # still split by group
    row_title = NULL,                   # <- hide group slice labels
    row_title_gp = gpar(fontsize = 0),  # <- belt-and-braces: ensure invisible
    left_annotation = row_annot,
    col = col_fun,
    heatmap_legend_param = list(title = "Z-scaled\nExpression", fontsize = 8),
    use_raster = TRUE,
  )
  ## because it's mouse let's modify gene names so they fit convention: sentence case
  p@row_names_param[["labels"]] <- stringr::str_to_sentence(p@row_names_param[["labels"]])
  p@row_names_param[["anno"]]@var_env[["value"]] <- stringr::str_to_sentence(p@row_names_param[["anno"]]@var_env[["value"]])
  
  draw(p,   padding = unit(c(10, 20, 10, 10), "mm"))  # prevent clipping
  dev.off()
  
  ### violinplot
  
  ViolinGenes <- fibro2_marker_ls[[4]]
  # make sure header gets included in vector
  colname <- colnames(ViolinGenes)[1]
  ViolinGenes <- c(colname, as.character(ViolinGenes[[1]]))
  ## add extra markers
  ViolinGenes <- append(ViolinGenes, 
                        values = c( # "ACTA2", <- wasn't that informative, all show minimal expression
                          "P4HA1", "CTSD")
  )
  ViolinGenes <- stringr::str_to_title(ViolinGenes)
  png(file = file.path(
    plots_dir,
    paste0("ViolinPlot_FibroMarkers_set2_",sample,".png")),
    width = 6, height = 10, units = "in", res = 300)
  p <- VlnPlot(
    seurat_obj,
    features = ViolinGenes,
    alpha = 1, pt.size = 0,
    combine = FALSE
  ) 
  
  # remove legends + tighten margins
  p <- lapply(p, function(pp) {
    pp + theme(
      legend.position = "none",
      plot.margin = margin(2, 0.5, 2, 1),
      axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5),
    )
  })
  
  for (i in seq_along(p)){
    p[[i]][["labels"]][["x"]] <- ""
    p[[i]][["labels"]][["y"]] <- ""
  } 
  
  p_stack <- wrap_plots(p, ncol = 2, nrow = 5,  guides = "collect")
  
  print(p_stack)
  dev.off()
  
  ### v2 violinplots
  
  
  ViolinGenes <- as.data.frame(fibro2_marker_ls[[5]])
  ViolinGenes <- as.data.frame(apply(ViolinGenes, 2, stringr::str_to_title))
  
  # desired height per subplot in inches
  per_gene_height <- 1   
  # 1) get a clean character vector of genes for this panel/column
  for (panel in colnames(ViolinGenes)) {
    # 1) get a clean character vector of genes for this panel/column
    genes <- ViolinGenes[[panel]]
    genes <- unique(na.omit(as.character(genes)))
    # keep only genes present in the object
    genes <- intersect(genes, rownames(seurat_obj))
    if (length(genes) == 0) next
    
    plots <- lapply(seq_along(genes), function(i) {
      g <- genes[i]
      gp <- VlnPlot(
        seurat_obj,
        features = g,
        group.by = "seurat_clusters",
        pt.size  = 0,
        combine  = TRUE
      ) + 
        NoLegend() +
        labs(title = "", x = "",
             y = stringr::str_to_sentence(tolower(g))) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.6, size = 14),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 14),
          axis.text.y = element_text(size = 8, angle = 0),
          plot.margin = margin(2, 5, 2, 5)
        )
      
      # Remove x-axis elements for all but the last plot
      if (i < length(genes)) {
        gp <- gp +
          theme(
            axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank()
          )
      } else {
        gp <- gp +
          labs(x = "") +
          theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
      }
      gp
    })
    
    # 3) stack vertically and add a global title
    panel_chr <- as.character(panel)
    p_stack <- wrap_plots(plots, ncol = 1, guides = "collect") +
      plot_annotation(
        title = panel_chr,
        theme = theme(plot.title = element_text(hjust = 0.6, size = 16, face = "bold"))
      )
    
    # output height = per-gene height * number of genes + extra for title
    h_in <- per_gene_height * length(plots) + 1
    # 4) save with png()/dev.off()
    png(
      file = file.path(
        plots_dir,
        paste0(
          "V2_ViolinPlot_FibroMarkers_",
          gsub("[^A-Za-z0-9._-]+", "_", panel_chr),
          "_", sample, ".png"
        )
      ),
      width = 5, height = h_in, units = "in", res = 300
    )
    print(p_stack)
    dev.off()
  }

  
  ### overlay umaps
  
  umap_markers <- fibro2_marker_ls[[6]]
  colname <- colnames(umap_markers)[1]
  umap_markers <- c(colname, as.character(umap_markers[[1]]))
  umap_markers<- append(umap_markers, "ZBTB7B")
  umap_markers <- stringr::str_to_title(umap_markers)
  
  png(file = file.path(
    plots_dir,
    paste0("Marker_Overlay_UMAPs_FibroMarkers_set2_",sample,".png")),
    width = 15, height = 25, units = "in", res = 300)
  p <- FeaturePlot(seurat_obj, 
                   features = umap_markers, 
                   label = F,# label.size = 3, repel = T,
  ) & 
    scale_color_viridis_c()
  
  print(p)
  dev.off()
  }
}

# optional quick annotation
bp <- BiocParallel::MulticoreParam(workers = 8)
ref <- readRDS("./R/HPCA_reference.rds")

library(SingleR)
for (i in seq_along(seurat_singlets_list)) {
  seurat_obj <- seurat_singlets_list[[i]]
  sample_name <- names(seurat_singlets_list)[i]

  # outfile <- file.path(git_dir, "Seperate_samples/SingleR_annotations", paste0(sample_name, "_SingleR_asSepCells.rds"))

  # if (file.exists(outfile)) {
  #   message(sample_name, " already annotated. Skipping.")
  #   next
  # }

  message("Annotating ", sample_name)

  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj)

  # Get cluster labels
  cluster_labels <- seurat_obj$seurat_clusters

  # Run SingleR annotation
  pred <- SingleR::SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.fine,
    clusters =  cluster_labels,
    #assay.type.test = 1,
    BPPARAM=bp
  )

  # Save prediction
  # saveRDS(pred, outfile)

  # Apply cluster annotations to each cell
  seurat_obj[["SingleR_cluster_label"]] <- pred$labels[
    match(seurat_obj$seurat_clusters, rownames(pred))
  ]

  ## per cell annotations

  # # Check alignment (rownames(pred) should match colnames(seurat_obj))
  # if (!all(rownames(pred) == colnames(seurat_obj))) {
  #   stop("Cell names in SingleR result do not match Seurat object.")
  # }
  #
  #
  # # (Optional) assign pruned labels or other metadata
  # seurat_obj$SingleR_pruned <- pred$pruned.labels
  #
  # # # Apply cluster annotations to each cell
  # seurat_obj[["SingleR_label"]] <- pred$labels


  # Store back into list
  seurat_singlets_list[[i]] <- seurat_obj
}

##checkpoint
saveRDS(seurat_singlets_list, "annotated_VU40T_sepsamples_singlets_only_mouse.rds")

for (i in seq_along(seurat_singlets_list)) {
  sample_name <- names(seurat_singlets_list)[i]
  seurat_obj <- seurat_singlets_list[[i]]

  p1 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE)

  p2 <- DimPlot(seurat_obj, group.by = "SingleR_cluster_label", label = F, repel = TRUE)
  png(file = file.path(plots_dir, paste0("Seurat_clusters_vs_SingleR_UMAP_ascluster_", sample_name, ".png")), 
      width = 16, height = 5, units = "in", res = 300)
  print(p1 + p2 + patchwork::plot_annotation(title = sample_name))
  dev.off()
}
# output_dir <- file.path(git_dir, "Seperate_samples/Plots")

SeuratList_markers <- list()
for (i in seq_along(seurat_singlets_list)){
  SeuratList_markers[[i]] <- FindAllMarkers(seurat_singlets_list[[i]],
                                            only.pos = T,
                                            min.pct = 0.25,
                                            logfc.threshold = 0.25)
}

# Set names to match SeuratList for easy referencing


names(SeuratList_markers) <- names(seurat_singlets_list)


for (i in seq_along(SeuratList_markers)) {
  # Extract sample name and corresponding marker table
  sample_name <- names(SeuratList_markers)[i]
  markers <- SeuratList_markers[[i]]
  seurat_obj <- seurat_singlets_list[[sample_name]]
  
  # Split markers by cluster and get top 3 per cluster
  topMarkers <- split(markers, markers$cluster)
  top10 <- lapply(topMarkers, head, n = 10)
  top10 <- do.call("rbind", top10)
  
  # Display or save FeaturePlot
  print(paste("Top 10 markers for", sample_name))
  # print(top10)
  
  
  # Generate FeaturePlots
  plots <- lapply(top10$gene, function(gene) {
    FeaturePlot(seurat_obj, features = gene) + ggtitle(gene)
  })
  
  # Combine with patchwork, arrange in 10 columns
  combined_plot <- patchwork::wrap_plots(plots, ncol = 10)
  
  # Save to PNG
  output_file <- file.path(plots_dir, paste0("top10_markers_", sample_name, ".png"))
  
  # Adjust image dimensions
  n_rows <- ceiling(length(plots) / 10)
  width_px <- 10 * 300   # ~300 px per plot
  height_px <- n_rows * 350
  
  png(filename = output_file,
      width = width_px, height = height_px, res = 150)
  print(combined_plot)
  dev.off()
}

for (i in seq_along(seurat_singlets_list)) {
  sample_name <- names(SeuratList_markers)[i]
  markers <- SeuratList_markers[[i]]
  seurat_obj <- seurat_singlets_list[[sample_name]]
  
  
  # Split markers by cluster and get top 3 per cluster
  topMarkers <- split(markers, markers$cluster)
  top10 <- lapply(topMarkers, head, n = 10)
  top10 <- do.call("rbind", top10)
  
  # # Display or save FeaturePlot
  # print(paste("Top 10 markers for", comb.names[i]))
  # # print(top10)
  # 
  
  # Generate FeaturePlots
  plots <- lapply(top10$gene, function(gene) {
    FeaturePlot(seurat_obj, features = gene) + ggplot2::ggtitle(gene)
  })
  
  # Combine with patchwork, arrange in 10 columns
  combined_plot <- patchwork::wrap_plots(plots, ncol = 10)
  
  # Save to PNG
  output_file <- file.path(git_dir, paste0("/Seperate_samples/Mouse/Plots/top10Markers_", sample_name, ".png"))
  
  # Adjust image dimensions
  n_rows <- ceiling(length(plots) / 10)
  width_px <- 10 * 300   # ~300 px per plot
  height_px <- n_rows * 350
  
  png(filename = output_file,
      width = width_px, height = height_px, res = 150)
  print(combined_plot)
  dev.off()
}



for (i in seq_along(SeuratList_markers)) {
  sample_name <- names(SeuratList_markers)[i]
  markers <- SeuratList_markers[[i]][SeuratList_markers[[i]]$p_val_adj <= 0.05, ]
  
  write.csv(markers, file = file.path(git_dir, paste0("/Seperate_samples/Mouse/ClusterMarkers_", sample_name, ".csv")), row.names = FALSE)
}


