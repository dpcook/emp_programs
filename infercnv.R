library(infercnv)

file_name <- commandArgs(trailingOnly = TRUE) #only 1 arg

# Load count matrix from RDS
mat <- readRDS(paste0("./data/", file_name, "_counts.rds"))
meta <- read.csv(paste0("./data/", file_name, "_meta.csv"),
                 row.names=1, header=T)

# Set up reference names as all cell types other than malignant cells in a good order
refs <- meta[-grep("_Epi", unique(meta$CellType)),]
refs <- refs[-grep("_Epi", refs)]
refs <- factor(refs, levels = c("Unknown", "Platelets", "DC", 
                                "NK_cell", "B_cell", "T_cells",
                                "Monocyte", "Macrophage", 
                                "Endothelial_cells", "Smooth_muscle_cells",
                                "Fibroblasts"))
refs <- unique(as.character(refs[order(refs)]))

# Create infercnv object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = mat,
                                     annotations_file = meta,
                                     gene_order_file = "./data/gencode_v19_gene_pos_infercnv.txt",
                                     ref_group_names=refs,
                                     max_cells_per_group = 500)
# Run infercnv

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir=paste0("./output/CNV_", file_name),
                             cluster_by_groups = T,
                             cluster_references = F,
                             denoise = T,
                             HMM = F,
                             num_threads=8,
                             output_format = "pdf",
                             useRaster = T)

#Can get subclusters from infercnv_obj@tumor_subclusters$subclusters
