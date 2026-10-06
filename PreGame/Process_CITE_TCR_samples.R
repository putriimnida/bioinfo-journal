
# --------------------------------------------
# Libraries and directories
# --------------------------------------------
source("/path/to/utils.R")
data_dir =  "/path/to/cellranger/outputs/"
out_dir = "/path/to/output/"
pdf_dir = "/path/to/save/pdfs/"
# --------------------------------------------

# --------------------------------------------
# Read in cellranger aligment data for sample x
# --------------------------------------------
sample = "sampleX"

# Gene and FB (protein) expression in demultiplexed cells assigned to the sample 
GEX_FB_sample = Seurat::Read10X_h5(paste0(data_dir, "step3_cellranger_multi/", sample, "_final/outs/per_sample_outs/", sample, "_final/count/sample_filtered_feature_bc_matrix.h5"))

# Raw gene and FB expression 
GEX_FB_raw = Seurat::Read10X_h5(paste0(data_dir, "step3_cellranger_multi/", sample, "_final/outs/multi/count/raw_feature_bc_matrix.h5"))

# List of predicted cells and their sample assignments. 
cr_cell_barcodes = read.csv(paste0(data_dir, "step1_demultiplexing/demultiplexed_samples/outs/multi/multiplexing_analysis/assignment_confidence_table.csv"), row.names = 1) 
# --------------------------------------------

# --------------------------------------------
# Determine background expression for dsb
# --------------------------------------------

FB_raw = GEX_FB_raw$`Antibody Capture`
GEX_raw = GEX_FB_raw$`Gene Expression`

rna_size = log10(Matrix::colSums(GEX_raw)+1)
prot_size = log10(Matrix::colSums(FB_raw)+1)
ngene = Matrix::colSums(GEX_raw > 0)
mtgene = grep(pattern = "^MT-", rownames(GEX_raw), value = TRUE)
propmt = Matrix::colSums(GEX_raw[mtgene, ]) / Matrix::colSums(GEX_raw)
md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
md$bc = rownames(md)

pdf(paste0(pdf_dir, 'test.pdf'), width = 16, height = 8)
p1 = ggplot(md[md$rna_size > 0, ], aes(x = rna_size)) + geom_histogram(fill = "dodgerblue") + ggtitle("RNA library size \n distribution")
p2 = ggplot(md[md$prot_size> 0, ], aes(x = prot_size)) + geom_density(fill = "firebrick2") + ggtitle("Protein library size \n distribution")
cowplot::plot_grid(p1, p2, nrow = 1)
dev.off()

# If samples were multiplxed: 
# Add CR cell assignment to QC plots 
cr_cell_assign = cr_cell_barcodes[,c("Barcode", "Assignment")]
colnames(cr_cell_assign) = c("bc", "Assignment")
t = left_join(md, cr_cell_assign, by = "bc")
t$is_cell = ifelse(is.na(t$Assignment) == T, "Background", "Cell")
t$Assignment = ifelse(is.na(t$Assignment) == T, "Background", as.character(t$Assignment))

# Same QC plots as above but split for CR cell assignment
pdf(paste0(pdf_dir, 'test.pdf'), width = 16, height = 8)
p1 = ggplot(t[t$rna_size > 0, ], aes(x = rna_size, fill = Assignment)) +
  geom_histogram() +
  ggtitle("RNA library size \n distribution")
p2 = ggplot(t[t$prot_size> 0, ], aes(x = prot_size, fill = Assignment)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_vline(xintercept = 2, linetype = "dashed") +
  ggtitle("Protein library size \n distribution")
cowplot::plot_grid(p1, p2, nrow = 1)
dev.off()


# Same QC plots as above but split for CR cell assignment and filtered for sample of analysis
t = subset(t, t$Assignment %ni% "C0252") #change based on multiplexed tag used

pdf(paste0(pdf_dir, 'test.pdf'), width = 16, height = 8)
p1 = ggplot(t[t$rna_size > 0, ], aes(x = rna_size, fill = Assignment)) +
  geom_histogram() +
  ggtitle("RNA library size \n distribution")
p2 = ggplot(t[t$prot_size> 0, ], aes(x = prot_size, fill = Assignment)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 2.75, linetype = "dashed") +
  geom_vline(xintercept = 2, linetype = "dashed") +
  ggtitle("Protein library size \n distribution")
cowplot::plot_grid(p1, p2, nrow = 1)
dev.off()


# Select thresholds for high quality cells per QC plots above
if (sample == "sampleX") { 
  background_drops = t[t$prot_size > 1.6 & t$prot_size < 3 & t$ngene < 80 & t$is_cell == "Background", ]$bc 
  positive_cells = t[t$prot_size >= 3 & t$ngene > 200 & t$propmt < 0.15 & t$Assignment == "C0251", ]$bc
  
}


negative_mtx_rawprot = FB_raw[ , background_drops] %>% as.matrix()
cells_mtx_rawprot = GEX_FB_sample$`Antibody Capture`[ , positive_cells[which(positive_cells %in% colnames(GEX_FB_sample$`Antibody Capture`))]] %>% as.matrix()

t$Bg_thresh = ifelse(t$bc %in% background_drops, "Background", 
                     ifelse(t$bc %in% positive_cells, "Cells",
                            "Excluded"))
table(t$Bg_thresh, useNA="ifany")
table(t$Bg_thresh, t$Assignment,useNA="ifany")

# Inspect final QC plot:
pdf(paste0(pdf_dir, 'test.pdf'), width = 16, height = 8)
p1 = ggplot(t[t$rna_size > 0, ], aes(x = rna_size, fill = Bg_thresh)) +
  geom_histogram() +
  ggtitle("RNA library size \n distribution")
p2 = ggplot(t[t$prot_size> 0, ], aes(x = prot_size, fill = Bg_thresh)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  geom_vline(xintercept = 1.75, linetype = "dashed") +
  ggtitle("Protein library size \n distribution")
cowplot::plot_grid(p1, p2, nrow = 1)
dev.off()


# Run DSB FB normalization 
dsb_norm_prot = DSBNormalizeProtein(
  cell_protein_matrix = cells_mtx_rawprot, # cell containing droplets
  empty_drop_matrix = negative_mtx_rawprot, # estimate ambient noise with the background drops 
  denoise.counts = TRUE, # model and remove each cell's technical component
  use.isotype.control = TRUE, # use isotype controls to define the technical component
  isotype.control.name.vec = rownames(cells_mtx_rawprot)[grep("Ctrl", rownames(cells_mtx_rawprot))] # names of isotype control abs
)

# --------------------------------------------


# --------------------------------------------
# Create Seurat object using GEX 
# --------------------------------------------

GEX.seurat =
  GEX_FB_sample$`Gene Expression`[ , positive_cells
                                   [which(positive_cells %in% colnames(GEX_FB_sample$`Gene Expression`))]] %>% 
  as.matrix() %>%
  CreateSeuratObject(., project = paste0(sample), min.cells = 3, min.features = 200) %>% 
  PercentageFeatureSet(., pattern = "^MT-", col.name = "percent.mt") %>%
  CellCycleScoring(., s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F) %>%
  SCTransform(.,                               
              vars.to.regress = "percent.mt", 
              method="glmGamPoi",
              vst.flavor = "v2",
              variable.features.n = 3000,       
              return.only.var.genes = F,
              verbose = T) %>% 
    RunPCA(., 
         features = VariableFeatures(object = .),
         verbose = FALSE) 
  
# Choose optimal number of PCs for downstream clustering (e.g. by elbow plot, etc.)
dims = c(1:13) 
GEX.seurat =
  GEX.seurat %>%
  RunUMAP(., dims = dims) %>%
  FindNeighbors(., dims = dims) %>%
  FindClusters(., resolution = 0.5)

# --------------------------------------------

# --------------------------------------------
# Add FB data from above
# --------------------------------------------

GEX.seurat[["FB"]] = CreateAssayObject(data = dsb_norm_prot[,rownames(GEX.seurat@meta.data)])
GEX.seurat@assays$FB@counts = cells_mtx_rawprot[,rownames(GEX.seurat@meta.data)]
GEX.seurat = ScaleData(GEX.seurat, assay = "FB")

# Cluster on FB data
DefaultAssay(GEX.seurat) = "FB"

# Choose optimal number of PCs for downstream clustering (e.g. by elbow plot, etc.)
fb_dims = c(1:12)

GEX.seurat =
  GEX.seurat %>% 
  RunPCA(.,
         features = rownames(.),
         reduction.name = "pca_fb",
         reduction.key = "PC.FB_",
         verbose = F) %>%
  RunUMAP(.,
          dims = fb_dims, 
          reduction = "pca_fb",
          reduction.name = "umap_fb",
          reduction.key = "UMAP.FB_") %>%
  FindNeighbors(.,
                dims = fb_dims,
                reduction = "pca_fb") %>%
  FindClusters(.,
               graph.name = "FB_snn",
               resolution = 0.5)

DefaultAssay(GEX.seurat) = "SCT"
Idents(GEX.seurat) = GEX.seurat$SCT_snn_res.0.5
# --------------------------------------------


# --------------------------------------------
# Get and Format BCR, TCR, and gdTCR data
# --------------------------------------------

# BCR paths 
VDJ_B = read.table(paste0(data_dir, "step3_cellranger_multi/", sample, "_final/outs/per_sample_outs/", sample, "_final/vdj_b/filtered_contig_annotations.csv"), 
                   header = T, sep = ",", stringsAsFactors = F)
VDJ_B.clntyp = read.table(paste0(data_dir, "step3_cellranger_multi/", sample, "_final/outs/per_sample_outs/", sample, "_final/vdj_b/clonotypes.csv"), 
                          header = T, sep = ",", stringsAsFactors = F)
# TCR paths 
VDJ_T = read.table(paste0(data_dir, "step3_cellranger_multi/", sample, "_final/outs/per_sample_outs/", sample, "_final/vdj_t/filtered_contig_annotations.csv"), 
                   header = T, sep = ",", stringsAsFactors = F)
VDJ_T.clntyp = read.table(paste0(data_dir, "step3_cellranger_multi/", sample, "_final/outs/per_sample_outs/", sample, "_final/vdj_t/clonotypes.csv"), 
                          header = T, sep = ",", stringsAsFactors = F)

# Merge BCR file data and select needed columns 
VDJ_B.fulljoin = 
  VDJ_B %>% dplyr::select(barcode, raw_clonotype_id) %>% distinct(barcode, .keep_all = TRUE) %>% 
  full_join(., VDJ_B.clntyp, by=c("raw_clonotype_id"="clonotype_id")) %>% 
  arrange(desc(frequency)) %>% 
  mutate(topClone = if_else(raw_clonotype_id == "clonotype1", "topClone", "minorClones"))

# Merge TCR file data and select needed columns 
VDJ_T.fulljoin =
  VDJ_T %>% dplyr::select(barcode, raw_clonotype_id) %>% distinct(barcode, .keep_all = TRUE) %>% 
  full_join(., VDJ_T.clntyp, by=c("raw_clonotype_id"="clonotype_id")) %>% 
  arrange(desc(frequency)) # %>% dplyr::select(-c(inkt_evidence, mait_evidence))

# Merge gd TCR file data from each alignment, save re-named and unique gd T cell clonotypes and keep relevant barcode level columns
VDJ_T_GD.fulljoin = parse_gd_tcrs(sample, data_dir)

# --------------------------------------------


# --------------------------------------------
# Add VDJ data to Seurat
# --------------------------------------------

GEX.seurat.meta = GEX.seurat@meta.data

dim(VDJ_B.fulljoin)
VDJ_B.fulljoin = subset(VDJ_B.fulljoin, VDJ_B.fulljoin$barcode %in% Intersect(list(VDJ_B.fulljoin$barcode, rownames(GEX.seurat.meta))))
dim(VDJ_B.fulljoin)

dim(VDJ_T.fulljoin)
VDJ_T.fulljoin =  subset(VDJ_T.fulljoin, VDJ_T.fulljoin$barcode %in% Intersect(list(VDJ_T.fulljoin$barcode, rownames(GEX.seurat.meta))))
dim(VDJ_T.fulljoin)

dim(VDJ_T_GD.fulljoin)
VDJ_T_GD.fulljoin =  subset(VDJ_T_GD.fulljoin, VDJ_T_GD.fulljoin$barcode %in% Intersect(list(VDJ_T_GD.fulljoin$barcode, rownames(GEX.seurat.meta))))
dim(VDJ_T_GD.fulljoin)

# Recalculate TCR / BCR proportions based on cells retained 
VDJ_B.fulljoin = VDJ_B.fulljoin %>%
  dplyr::group_by(raw_clonotype_id) %>%
  dplyr::mutate(frequency = length(raw_clonotype_id)) %>%
  ungroup() %>%
  dplyr::mutate(proportion = frequency / nrow(.))

# Merge all rows of AB and GD TCRs, make unique, and recalculate combined TCR frequency 
VDJ_TCR.fulljoin = VDJ_T.fulljoin %>%
  dplyr::select(-frequency, -proportion, -inkt_evidence, -mait_evidence) %>%
  dplyr::mutate(raw_clonotype_id = paste0("AB_",raw_clonotype_id)) %>%
  rbind(VDJ_T_GD.fulljoin %>% 
          dplyr::select(-comment, -clonotype_key) %>%
          dplyr::mutate(raw_clonotype_id = paste0("GD_",raw_clonotype_id))) %>%
  dplyr::group_by(raw_clonotype_id) %>%
  dplyr::mutate(frequency = length(raw_clonotype_id)) %>%
  ungroup() %>%
  dplyr::mutate(proportion = frequency / length(unique(barcode)))

# Use above merge dataframe to join total TCR frequencies back to AB TCRs
VDJ_T.fulljoin = VDJ_T.fulljoin %>%
  dplyr::mutate(raw_clonotype_id = paste0("AB_",raw_clonotype_id)) %>%
  dplyr::select(-frequency,-proportion) %>%
  left_join(., VDJ_TCR.fulljoin %>%
              dplyr::select(raw_clonotype_id,frequency,proportion),
            by = "raw_clonotype_id",
            multiple = "first") %>%
  dplyr::mutate(raw_clonotype_id = gsub("AB_", "", raw_clonotype_id)) %>%
  dplyr::select(barcode, raw_clonotype_id, frequency, proportion, cdr3s_aa, cdr3s_nt,
                inkt_evidence, mait_evidence)

# Same as above comment for GD TCRs. The frequency and proportion for GD TCRs is now demultiplexed and scaled for total AB + GD
VDJ_T_GD.fulljoin = VDJ_T_GD.fulljoin %>%
  dplyr::mutate(raw_clonotype_id = paste0("GD_",raw_clonotype_id)) %>%
  left_join(., VDJ_TCR.fulljoin %>%
              dplyr::select(raw_clonotype_id,frequency,proportion),
            by = "raw_clonotype_id",
            multiple = "first") %>%
  dplyr::mutate(raw_clonotype_id = gsub("GD_", "", raw_clonotype_id)) %>%
  dplyr::select(barcode, raw_clonotype_id, frequency, proportion, cdr3s_aa, cdr3s_nt, clonotype_key, comment)


# Add TCR / BCR columns to metadata
GEX.seurat.meta2 = 
  GEX.seurat.meta %>% 
  rownames_to_column(var="barcode") %>% 
  
  # join VDJ_B
  full_join(., VDJ_B.fulljoin %>% 
              distinct(barcode, .keep_all = T) %>% 
              dplyr::select(barcode, raw_clonotype_id, frequency, proportion, cdr3s_aa, cdr3s_nt) %>% 
              rename_with(.cols = names(.), .fn = ~ gsub("^", "VDJ_B_", .x)), 
            by=c("barcode"="VDJ_B_barcode"), keep=T) %>% 
  # join VDJ_T
  full_join(., VDJ_T.fulljoin %>% 
              distinct(barcode, .keep_all = T) %>% 
              rename_with(.cols = names(.), .fn = ~gsub("^", "VDJ_T_", .x)),
            by=c("barcode"="VDJ_T_barcode"), keep=T) %>% 
  # join GD-VDJ_T
  full_join(., VDJ_T_GD.fulljoin %>% 
              distinct(barcode, .keep_all = T) %>% 
              rename_with(.cols = names(.), .fn = ~gsub("^", "VDJ_T_GD_", .x)),
            by=c("barcode"="VDJ_T_GD_barcode"), keep=T) %>% 
  
  ## cell_type_v1 based on BCR/TCR detected
  mutate(cell_type_v1 =  case_when(!is.na(VDJ_B_barcode) & !is.na(VDJ_T_barcode) & is.na(VDJ_T_GD_barcode) ~ "abTBcell",
                                   !is.na(VDJ_B_barcode) & !is.na(VDJ_T_barcode) & !is.na(VDJ_T_GD_barcode) ~ "abgdTBcell",
                                   !is.na(VDJ_B_barcode) & is.na(VDJ_T_barcode) & !is.na(VDJ_T_GD_barcode) ~ "gdTBcell",
                                   !is.na(VDJ_B_barcode) & is.na(VDJ_T_barcode) & is.na(VDJ_T_GD_barcode) ~ "Bcell",
                                   is.na(VDJ_B_barcode) & !is.na(VDJ_T_barcode) & is.na(VDJ_T_GD_barcode) ~ "abTcell",
                                   is.na(VDJ_B_barcode) & is.na(VDJ_T_barcode) & !is.na(VDJ_T_GD_barcode) ~ "gdTcell",
                                   is.na(VDJ_B_barcode) & !is.na(VDJ_T_barcode) & !is.na(VDJ_T_GD_barcode) ~ "abgdTcell",
                                   is.na(VDJ_B_barcode) & is.na(VDJ_T_barcode) & is.na(VDJ_T_GD_barcode) ~ "Other")) %>%
  
  ## Bcell_clonal
  mutate(Bcell_clonal = if_else(VDJ_B_raw_clonotype_id == "clonotype1", "monoclonal", 
                                if_else(is.na(VDJ_B_raw_clonotype_id), "NA", "polyclonal"))) %>% 
  
  ## write out rownames
  column_to_rownames(var="barcode") %>%
  dplyr::select(-VDJ_B_barcode, -VDJ_T_barcode, -VDJ_T_GD_barcode)

identical(rownames(GEX.seurat.meta2), rownames(GEX.seurat@meta.data)) #check TRUE
GEX.seurat@meta.data = GEX.seurat.meta2

# --------------------------------------------


# --------------------------------------------
# WNN Clustering 
# --------------------------------------------

DefaultAssay(GEX.seurat) = "SCT"
GEX.seurat = FindMultiModalNeighbors(GEX.seurat,
                             reduction.list = list("pca", "pca_fb"), 
                             dims.list = list(dims, fb_dims), 
                             modality.weight.name = "SCT.weight", 
                             verbose = T) %>%
  RunUMAP(.,
          nn.name = "weighted.nn", 
          reduction.name = "wnn.umap", 
          reduction.key = "wnnUMAP_") %>%
  FindClusters(.,
               graph.name = "wsnn",
               algorithm = 3,
               resolution = 1, 
               verbose = T, 
               random.seed = 1990)

# Done, save per sample object for cluster annotations / downstream integration
saveRDS(GEX.seurat, paste0(out_dir, "Multiome.seurat.",sample, ".final.Rds"))



# --------------------------------------------
