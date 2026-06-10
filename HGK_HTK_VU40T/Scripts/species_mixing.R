# Load packages
library(Seurat)
library(Matrix)
library(tidyverse)
library(optparse)

# Helper function: compute species % per barcode
compute_species_df <- function(human_counts, mouse_counts) {
  all_barcodes <- union(names(human_counts), names(mouse_counts))
  
  df <- data.frame(
    barcode = all_barcodes,
    human_counts = human_counts[all_barcodes],
    mouse_counts = mouse_counts[all_barcodes]
  )
  
  df[is.na(df)] <- 0
  df$total_counts <- df$human_counts + df$mouse_counts
  df$percent_human <- df$human_counts / df$total_counts * 100
  df$percent_mouse <- df$mouse_counts / df$total_counts * 100
  
  return(df)
}

convert_ensembl_to_symbol <- function(seurat_obj, species = c("human", "mouse"), assay = "RNA") {
  # Load biomaRt if needed
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("biomaRt")
  }
  library(biomaRt)
  
  # Match species and set Ensembl dataset
  species <- match.arg(species)
  dataset <- switch(
    species,
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl"
  )
  symbol_attr <- switch(
    species,
    human = "hgnc_symbol",
    mouse = "mgi_symbol"
  )
  
  # Connect to Ensembl
  message("Connecting to Ensembl for ", species, "...")
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Extract Ensembl IDs
  ens_ids <- rownames(seurat_obj[[assay]])
  
  # Get gene symbol mapping
  gene_map <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", symbol_attr),
    values = ens_ids,
    mart = ensembl
  )
  
  # Clean map
  gene_map <- gene_map[gene_map[[symbol_attr]] != "", ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  common_ids <- intersect(rownames(seurat_obj[[assay]]), gene_map$ensembl_gene_id)
  
  # Replace Ensembl IDs with symbols
  new_names <- toupper(gene_map[[symbol_attr]][match(common_ids, gene_map$ensembl_gene_id)])
  rownames(seurat_obj[[assay]])[match(common_ids, rownames(seurat_obj[[assay]]))] <- new_names
  
  message("Replaced ", length(common_ids), " Ensembl IDs with gene symbols.")
  return(seurat_obj)
}


## func to get aggregate rows w/ duplicated mouse gene names
aggregate_by_gene_fast <- function(mat) {
  # Convert sparse matrix to triplet format (i,j,x)
  m_coo <- summary(mat)
  genes <- rownames(mat)
  barcodes <- colnames(mat)
  
  # Add gene names to triplet
  m_coo$gene <- genes[m_coo$i]
  m_coo$barcode <- barcodes[m_coo$j]
  
  # Convert to data.table
  dt <- as.data.table(m_coo)
  
  # Sum duplicated gene-barcode combinations
  dt_agg <- dt[, .(x = sum(x)), by = .(gene, barcode)]
  
  # Convert back to sparse matrix
  gene_levels <- unique(dt_agg$gene)
  barcode_levels <- unique(dt_agg$barcode)
  
  sparse_mat <- Matrix::sparseMatrix(
    i = match(dt_agg$gene, gene_levels),
    j = match(dt_agg$barcode, barcode_levels),
    x = dt_agg$x,
    dimnames = list(gene_levels, barcode_levels)
  )
  
  return(sparse_mat)
}

option_list <- list(
  make_option(c("-h", "--humandir"),
              type = "character",
              default = "human",
              help = "species to use [default = %default]. Options: human or mouse",
              metavar = "character"),
  make_option(c("-m", "--mousedir"),
              type = "character",
              default = ,
              help = "input RDS dir to load",
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
# print(opt)

# For each sample:
species_list <- list()
seurat_list_human <- list()
seurat_list_mouse <- list()

humandir <- opt$humandir
mousedir <- opt$mousedir
outdir <- opt$outdir
proj_dir <- "/rds/projects/g/gendood-3dmucosa/"
outdir <- file.path(proj_dir, "scRNAseqAnalysis/OralMucosa/HGK_HTK_VU40T/")

plots_dir <- file.path(outdir, "Mixed_species_plots/")
if (!(dir.exists(plots_dir))) {
  dir.create(plots_dir)
}
## wd
mtx_dir_mouse <- file.path(proj_dir, "Globus/hgk_htk_VU40T_QC_and_counts_Mouse/cellranger/mtx_conversions")
mtx_dir_human <- file.path(proj_dir, "Globus/hgk_htk_VU40T_QC_and_counts_human/cellranger/mtx_conversions")
Samples <- list.files(path = mtx_dir_mouse, 
                      pattern = "^GB",
                      recursive = F
) ##check we have all samples

## get only cellbender seurat_obj paths
seurat_objs_human <- list.files(path = mtx_dir_human, 
                          pattern = "^GB.*cellbender.*.seurat.rds$",
                          full.names = T,
                          recursive = T
)

seurat_objs_mouse <- list.files(path = mtx_dir_mouse, 
                                pattern = "^GB.*cellbender.*.seurat.rds$",
                                full.names = T,
                                recursive = T
)

# remove HGK_rep2 from list as it has been qc excluded

#sample_to_exclude <- "GB240624-MW_Rep2_HGK"

#seurat_objs_human <- seurat_objs_human[!grepl(sample_to_exclude, seurat_objs_human)]
#seurat_objs_mouse <- seurat_objs_human[!grepl(sample_to_exclude, seurat_objs_mouse)]

for (i in Samples) {
  # Load CellBender matrices
  seurat_human <- readRDS(paste0(mtx_dir_human, "/", i,"/" , i, "_cellbender_filter_matrix.seurat.rds"))
  seurat_mouse <- readRDS(paste0(mtx_dir_mouse, "/", i, "/", i, "_cellbender_filter_matrix.seurat.rds"))
  
  ##get total counts per barcode
  human_counts_total <- Matrix::colSums(GetAssayData(seurat_human, slot = "counts"))
  mouse_counts_total <- Matrix::colSums(GetAssayData(seurat_mouse, slot = "counts"))
  
  # Compute species score df
  species_df <- compute_species_df(human_counts_total, mouse_counts_total)
  species_df$sample <- i
  species_df$species_label <- case_when(
    species_df$percent_human >= 70 ~ "Human",
    species_df$percent_mouse >= 70 ~ "Mouse",
    TRUE ~ "Mixed"
  )
  # Store species df and Seurat object for later
  species_list[[i]] <- species_df
  
  # For visualization, pick one object to hold metadata (e.g. human Seurat object)
  barcodes_in_seurat <- colnames(seurat_human)
  species_df_sub <- species_df[species_df$barcode %in% barcodes_in_seurat, ]
  species_df_sub <- species_df_sub[match(barcodes_in_seurat, species_df_sub$barcode), ]
  
  ## set quite a liberal threshold for species mixing 
  ## due to expected ~85% sequence homology between mus musculus and homo sapiens
  seurat_human$percent_human <- species_df_sub$percent_human
  seurat_human$percent_mouse <- species_df_sub$percent_mouse
  seurat_human$species_label <- case_when(
    seurat_human$percent_human >= 70 ~ "Human",
    seurat_human$percent_mouse >= 70 ~ "Mouse",
    TRUE ~ "Mixed"
  )
  
  seurat_mouse$percent_human <- species_df_sub$percent_human
  seurat_mouse$percent_mouse <- species_df_sub$percent_mouse
  seurat_mouse$species_label <- case_when(
    seurat_mouse$percent_human >= 70 ~ "Human",
    seurat_mouse$percent_mouse >= 70 ~ "Mouse",
    TRUE ~ "Mixed"
  )
  
  # Add sample ID
  seurat_human$sample <- paste0("Sample", i)
  seurat_mouse$sample <- paste0("Sample", i)
  
  seurat_list_human[[i]] <- seurat_human
  seurat_list_mouse[[i]] <- seurat_mouse
}

# ## inspect species distribution
plot_list <- list()

for (sample_name in names(species_list)) {
  
  df <- species_list[[sample_name]]
  
  # create percent_human once
  df$percent_human <- 100 - df$percent_mouse
  
  p <- ggplot(df, aes(x = percent_human)) +
    geom_histogram(bins = 50, fill = "#3c8dbc", color = "black") +
    geom_vline(xintercept = c(30, 70), linetype = "dashed") +
    theme_bw() +
    ggtitle(sample_name) +
    ylab("Cell Count") +
    scale_x_continuous(
      name = "Percent Human",
      limits = c(0, 100),
      
      # only the TOP axis is reversed
      sec.axis = sec_axis(
        trans = ~ 100 - .,
        name = "Percent Mouse"
      )
    ) +
    annotate("text", x = 12, y = Inf, label = "Mouse",
             vjust = 1.5, size = 4, colour = "red") +
    annotate("text", x = 50, y = Inf, label = "Mixed",
             vjust = 1.5, size = 4, colour = "red") +
    annotate("text", x = 88, y = Inf, label = "Human",
             vjust = 1.5, size = 4, colour = "red")
  
  plot_list[[sample_name]] <- p
}

# Combine all into one figure
combined_plot <- patchwork::wrap_plots(plot_list, ncol = (length(species_list))/2, nrow = 3)  # change columns as you like

#human only data in HGK, limited in HGK rep 2 (sample failed QC, insufficient cells)
# mouse only data: VU40T (human peak 100 cells), HTK-HPV16 (HTK_Rep1_A)
# both human and mouse worked: HTK-HPV18 (HTK), VU40T (human-ish with peak at 100 cells)

## samples to exclude from mouse analysis: both HGKs

ggsave(
  file.path(plots_dir,
    "percent_species_distribution_histograms.png"
  )
  ,
  combined_plot,
  width = 14, height = 10, dpi = 300
)

## output this for further analysis on seperate alignments

species_df_all <- do.call(rbind, species_list)

species_df_all$sample <- stringr::str_remove(species_df_all$sample, "GB240624-")

species_df_all$sample <- stringr::str_remove(species_df_all$sample, "_S[0-9]")

mapping <- c(
  "MW_HGK"        = "MW_HGK",
  "MW_HTK"        = "MW_HTK_HPV18",
  "MW_Rep2_VU40T" = "MW_VU40T", 
  "MW_Rep1_A"     = "MW_HTK2_HPV16",
  "MW_Rep2_HGK"   = "MW_HGK_Rep2" 
)


species_df_all$sample <- mapping[species_df_all$sample]

write.csv(species_df_all, 
          file = file.path(
            outdir, 
            "VU40T_species_assignment_by_reads.csv")
          )

all_species <- list(
  human = seurat_list_human,
  mouse = seurat_list_mouse
)

for (species in names(all_species)) {
  alignment_list <- all_species[[species]]
  message("--- Processing alignment species: ", species, " ---")
  
  # Get the names of the individual samples in this alignment list
  sample_names <- names(alignment_list)
  
  for (sample_id in sample_names) {
    message("  Generating diagnostic plots for sample: ", sample_id)
    
    # Extract the single sample Seurat object
    seu_obj <- alignment_list[[sample_id]]
    
    # Run standard processing on the single sample to see unvarnished local coordinates
    seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
    seu_obj <- FindVariableFeatures(seu_obj, verbose = FALSE)
    seu_obj <- ScaleData(seu_obj, verbose = FALSE)
    seu_obj <- RunPCA(seu_obj, verbose = FALSE)
    seu_obj <- RunUMAP(seu_obj, dims = 1:20, verbose = FALSE)
    
    # -------------------------------------------------------------
    # Plot 1: Continuous Species Scores (percent_human vs percent_mouse)
    # -------------------------------------------------------------
    fn_scores <- file.path(
      plots_dir, 
      paste0(species, "_alignment_", sample_id, "_percent_species_UMAP.png")
    )
    
    png(filename = fn_scores, width = 12, height = 6, units = "in", res = 300)
    p1 <- FeaturePlot(seu_obj, features = c("percent_human", "percent_mouse"), ncol = 2) +
      plot_annotation(title = paste0(sample_id, " (", species, " alignment): Species Transcript Percentages"))
    print(p1)
    dev.off()
    
    # -------------------------------------------------------------
    # Plot 2: Discrete Species Assignment Label
    # -------------------------------------------------------------
    fn_labels <- file.path(
      plots_dir, 
      paste0(species, "_alignment_", sample_id, "_species_assignment_UMAP.png")
    )
    
    png(filename = fn_labels, width = 8, height = 6, units = "in", res = 300)
    p2 <- DimPlot(seu_obj, group.by = "species_label", label = TRUE, repel = TRUE) + 
      ggtitle(paste0(sample_id, " (", species, " alignment):\nSpecies Assignment")) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    print(p2)
    dev.off()
  }
}
