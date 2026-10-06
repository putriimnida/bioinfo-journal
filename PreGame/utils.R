require(Seurat)
require(harmony)
require(dsb)
require(ggplot2)
require(dplyr)
require(tidyverse) 
require(RColorBrewer) 
require(ggsignif)
require(DESeq2)
require(glmGamPoi)
require(sctransform)
source("/path/to/palettes.R")

'%ni%' <- Negate('%in%')

pdff = function( ...,  width = 8, height = 7) {
  plots = list(...)
  print(paste0("saving ", length(plots), " plot(s)"))
  pdf(paste0(out_dir, "test.pdf"), width, height)
  for (i in 1:length(plots)) {
    print(plots[[i]])
    print(paste0("saved plot ",i))
  }
  dev.off()
}

fb_names = function(seurat, grep_patt = NA, assay = "FB") {
  # Useful because i always forget protein names 
  if (is.na(grep_patt)) {
    print(rownames(GetAssay(seurat, assay)))
  } else {
    print(rownames(GetAssay(seurat, assay))[grep(grep_patt, rownames(GetAssay(seurat, assay)))])
  }
}


find_gdTCR_results_dir = function(directory) {
  # arguments
  # directory - character string, parent directory to search for 'gdTCR_all' directory in
  
  files = list.files(directory)
  if ("gdTCR_all" %in% files) {
    print(paste0("Found gd TCR directory: ", paste0(directory, "gdTCR_all/" )))
    return(paste0(directory, "gdTCR_all/"))
  } else {
    daughter_dirs = list.dirs(directory, recursive = F)
    for (i in 1:length(daughter_dirs)) {
      files = list.files(daughter_dirs[i])
      if ("gdTCR_all" %in% files) {
        print(paste0("Found gd TCR directory: ", paste0(daughter_dirs[i], "/gdTCR_all/" )))
        return(paste0(daughter_dirs[i], "/gdTCR_all/"))
      }
    }
  }
  
  return("failed to locate /gdTCR_all/ directory in first 2 levels of provided directory, please check such a directory exists or if the provided directory is valid...")
  
}

find_analysis_dir = function(directory) {
  # arguments
  # directory - character string, parent directory to search for 'gdTCR_all' directory in
  
  files = list.files(directory)
  if ("seq_batch_gd_clonotypes" %in% files) {
    return(paste0(directory, "seq_batch_gd_clonotypes/"))
  } else {
    daughter_dirs = list.dirs(directory, recursive = F)
    for (i in 1:length(daughter_dirs)) {
      files = list.files(daughter_dirs[i])
      if ("seq_batch_gd_clonotypes" %in% files) {
        return(paste0(daughter_dirs[i], "/seq_batch_gd_clonotypes/"))
      }
    }
  }
  
  return("failed to locate /seq_batch_gd_clonotypes/ directory in first 2 levels of provided directory, please check such a directory exists or if the provided directory is valid...")
  
}

parse_gd_tcrs = function(sample, sample_dir, save_only = F, save_dir = NULL, verbose = F) { 
  
  # Arguments: 
  # sample - character string, sample to parse gdTCRs, only used for output file name 
  # sample_dir - character string, alignment directory for the given sample that contains a 'gdTCR_all' directory somewhere
  # save_only - logical, whether to only save the parsed gdTCRs (non-barcode flattened), default will save and return parsed gdTCRs for addition to big_mak@meta.data
  # save_dir - character string, optional, if not provided, parsed gdTCRs will be saved to "seq_batch_gd_clonotypes" in either the myeloma or AML analysis directories
  #                                       if provided, will save to the provided directory.
  # verbose - logical, whether to print more intermediate messaging. 
  
  
  # Notes on sample_dir structure: 
  # In this study we used two different gamma delta TCR sequencing primers, generating two distinct libraries / sets of fastqs.
  # The primer sets are derived from Mimitou et al., and Gherardin et al., (see Methods), referred to as "M" and "G" downstream.
  # Three alignments were performed:
  #     1. fastqs from the M library alone (M)
  #     2. fastqs from the G library alone (G)
  #     3. fastqs from both M and G libraries (GM)
  
  # The outputs are as they would be for typical abTCR alignment by cellranger. 
  # A directory called 'gdTCR_all' needs to be created within the parent directory specified as 'sample_dir' in arguments.
  # The "clonotypes.csv" and "filtered_contig_annotations.csv" output files for each of the 3 alignments, needs to be linked (or copied, though only tested with links, eg. ln -s [...])
  # into this directory, as follows: 
  
  # {sample_dir}
  #              /gdTCR_all/
  #                           GM_clonotypes
  #                           GM_contig
  #                           G_clonotypes
  #                           G_contig
  #                           M_clonotypes
  #                           M_contig
  
  
  # This above structure and correct files will be checked for first: 
  # Determine proper directory to search for the gdTCR alignment outputs
  print("checking for gdTCR_all directory in the provided sample dir...")
  
  if (grepl("gdTCR_all", sample_dir)) {
    gd_data_dir = sample_dir
  } else {
    gd_data_dir = find_gdTCR_results_dir(sample_dir)
  }
  
  # Check to ensure proper files are present
  gd_contigs = list.files(gd_data_dir, pattern = ".*_contig")
  stopifnot(length(gd_contigs) == 3)
  gd_clonotypes = list.files(gd_data_dir, pattern = ".*_clonotypes")
  stopifnot(length(gd_clonotypes) == 3)
  
  print("successfully located raw alignment files...")

  # Loop through the files from each alignment and merge contig and clonotype files as is done with VDJ_T and VDJ_B files 
  VDJ_T_GD.fulljoin.list = list()
  try_errors = c()
  try_error_on_primer_set = c("G", "GM", "M")
  for (i in 1:length(gd_contigs)) {
    VDJ_T_GD = try(read.table(paste0(gd_data_dir, gd_contigs[i]),
                          header = T, sep = ",", stringsAsFactors = F))
    VDJ_T_GD.clntyp = try(read.table(paste0(gd_data_dir, gd_clonotypes[i]),
                                 header = T, sep = ",", stringsAsFactors = F))
    
    if (class(VDJ_T_GD) == "try-error" | class(VDJ_T_GD.clntyp) == "try-error") {
      print(paste0("in reading in 10X vdj outputs, found empty file from ", try_error_on_primer_set[i], " primer set..."))
      try_errors = c(try_errors, try_error_on_primer_set[i])
    } else {
      # subset alpha/gamma chains, convert TRA to TRG, combine relevant columns to keep for later, 
      # and collapse any non-distinct information into 1 row per clonotype.
      # For clonotypes with multiple CDR3s per chain, collapse the CDRs by ";" as cellranger does, but keep the full_sequence 
      # info separated by a "?" 
      # For clonotypes with multiple full_sequence options (usually small changes in fwr2 or c_gene), collapse that information 
      # in full_sequence by "/"
      VDJ_T_GD_alpha_gamma_chain = as.data.frame(VDJ_T_GD %>%
                                                   filter(chain %in% c("TRA", "TRG")) %>%
                                                   group_by(cdr3) %>%
                                                   mutate(cdr3_freq = length(cdr3)) %>%
                                                   ungroup() %>%
                                                   mutate(chain = "TRG",
                                                          v_gene = gsub("TRA", "TRG", v_gene),
                                                          d_gene = gsub("TRA", "TRG", d_gene),
                                                          j_gene = gsub("TRA", "TRG", j_gene),
                                                          c_gene = gsub("TRA", "TRG", c_gene)) %>%
                                                   mutate(clonotype_key = paste(v_gene,d_gene,j_gene,cdr3_nt, sep = ";")) %>%
                                                   mutate(full_sequence = paste0(chain,";", v_gene, ";", d_gene, ";", j_gene, ";", c_gene,
                                                                                 "|fwr1:", fwr1_nt, "|cdr1:",cdr1_nt,
                                                                                 "|fwr2:", fwr2_nt, "|cdr2:",cdr2_nt,
                                                                                 "|fwr3:", fwr3_nt, "|cdr3:",cdr3_nt, "|fwr4:", fwr4_nt)) %>%
                                                   mutate(cdr3 = paste0("TRG:",cdr3),
                                                          cdr3_nt = paste0("TRG:",cdr3_nt),) %>%
                                                   dplyr::select(raw_clonotype_id, raw_consensus_id, cdr3, cdr3_nt, clonotype_key, full_sequence, cdr3_freq) %>% # barcode
                                                   distinct(., .keep_all = T) %>%
                                                   group_by(raw_clonotype_id, raw_consensus_id) %>%
                                                   mutate(full_sequence = paste(full_sequence, collapse = "/")) %>%
                                                   distinct(., .keep_all = T) %>%
                                                   ungroup() %>%
                                                   dplyr::select(-raw_consensus_id) %>%
                                                   group_by(raw_clonotype_id) %>%
                                                   # arrange(str_length(cdr3)) %>% 
                                                   arrange(desc(cdr3_freq)) %>%
                                                   dplyr::select(-cdr3_freq) %>%
                                                   mutate(cdr3 = paste(cdr3, collapse = ";"),
                                                          cdr3_nt = paste(cdr3_nt, collapse = ";"),
                                                          clonotype_key = paste(clonotype_key, collapse = ";"),
                                                          full_sequence = paste(full_sequence, collapse = "?")) %>%
                                                   distinct(., .keep_all = T))
      
      # Do the same for beta/delta chains
      VDJ_T_GD_beta_delta_chain = as.data.frame(VDJ_T_GD %>%
                                                  filter(chain %in% c("TRB", "TRD")) %>%
                                                  group_by(cdr3) %>%
                                                  mutate(cdr3_freq = length(cdr3)) %>%
                                                  ungroup() %>%
                                                  mutate(chain = "TRD",
                                                         v_gene = gsub("TRB", "TRD", v_gene),
                                                         d_gene = gsub("TRB", "TRD", d_gene),
                                                         j_gene = gsub("TRB", "TRD", j_gene),
                                                         c_gene = "TRDC") %>%
                                                  mutate(clonotype_key = paste(v_gene,d_gene,j_gene,cdr3_nt, sep = ";")) %>%
                                                  mutate(full_sequence = paste0(chain,";", v_gene, ";", d_gene, ";", j_gene, ";", c_gene,
                                                                                "|fwr1:", fwr1_nt, "|cdr1:",cdr1_nt,
                                                                                "|fwr2:", fwr2_nt, "|cdr2:",cdr2_nt,
                                                                                "|fwr3:", fwr3_nt, "|cdr3:",cdr3_nt, "|fwr4:", fwr4_nt)) %>%
                                                  mutate(cdr3 = paste0("TRD:",cdr3),
                                                         cdr3_nt = paste0("TRD:",cdr3_nt),) %>%
                                                  dplyr::select(raw_clonotype_id, raw_consensus_id, cdr3, cdr3_nt, clonotype_key, full_sequence, cdr3_freq) %>% #barcode 
                                                  distinct(., .keep_all = T) %>%
                                                  group_by(raw_clonotype_id, raw_consensus_id) %>%
                                                  mutate(full_sequence = paste(full_sequence, collapse = "/")) %>%
                                                  distinct(., .keep_all = T) %>%
                                                  ungroup() %>%
                                                  dplyr::select(-raw_consensus_id) %>%
                                                  group_by(raw_clonotype_id) %>%
                                                  # arrange(str_length(cdr3)) %>% 
                                                  arrange(desc(cdr3_freq)) %>%
                                                  dplyr::select(-cdr3_freq) %>%
                                                  mutate(cdr3 = paste(cdr3, collapse = ";"),
                                                         cdr3_nt = paste(cdr3_nt, collapse = ";"),
                                                         clonotype_key = paste(clonotype_key, collapse = ";"),
                                                         full_sequence = paste(full_sequence, collapse = "?")) %>%
                                                  distinct(., .keep_all = T))
      
      
      # Now add the new clonotype columns to VDJ_T_GD.clntyp 
      VDJ_T_GD_manual_clntypes = VDJ_T_GD.clntyp %>%
        dplyr::select(clonotype_id, frequency, proportion, cdr3s_aa) %>% #remove cdr3s_aa later
        distinct(., .keep_all = T) %>%
        left_join(., y = VDJ_T_GD_beta_delta_chain, 
                  by=c("clonotype_id" = "raw_clonotype_id")) %>%
        left_join(., y = VDJ_T_GD_alpha_gamma_chain, 
                  by=c("clonotype_id" = "raw_clonotype_id"),
                  suffix = c(".delta", ".gamma")) %>%
        mutate(cdr3 = case_when(!is.na(cdr3.delta) & !is.na(cdr3.gamma) ~ paste(cdr3.delta, cdr3.gamma, sep = ";"),
                                !is.na(cdr3.delta) & is.na(cdr3.gamma) ~ cdr3.delta,
                                is.na(cdr3.delta) & !is.na(cdr3.gamma) ~ cdr3.gamma),
               cdr3_nt = case_when(!is.na(cdr3_nt.delta) & !is.na(cdr3_nt.gamma) ~ paste(cdr3_nt.delta, cdr3_nt.gamma, sep = ";"),
                                   !is.na(cdr3_nt.delta) & is.na(cdr3_nt.gamma) ~ cdr3_nt.delta,
                                   is.na(cdr3_nt.delta) & !is.na(cdr3_nt.gamma) ~ cdr3_nt.gamma),
               clonotype_key = case_when(!is.na(clonotype_key.delta) & !is.na(clonotype_key.gamma) ~ paste(clonotype_key.delta, clonotype_key.gamma, sep = "|"),
                                         !is.na(clonotype_key.delta) & is.na(clonotype_key.gamma) ~ clonotype_key.delta,
                                         is.na(clonotype_key.delta) & !is.na(clonotype_key.gamma) ~ clonotype_key.gamma),
               full_sequence = case_when(!is.na(full_sequence.delta) & !is.na(full_sequence.gamma) ~ paste(full_sequence.delta, full_sequence.gamma, sep = ">"),
                                         !is.na(full_sequence.delta) & is.na(full_sequence.gamma) ~ full_sequence.delta,
                                         is.na(full_sequence.delta) & !is.na(full_sequence.gamma) ~ full_sequence.gamma)
        ) %>%
        dplyr::select(clonotype_id, frequency, proportion, cdr3, cdr3_nt, clonotype_key, full_sequence) #add cdr3_aa for checks against cellranger outputs
      
      VDJ_T_GD.fulljoin =
        VDJ_T_GD %>% 
        dplyr::select(barcode, raw_clonotype_id) %>% 
        dplyr::mutate(Primer_set = gsub("_.*", "", gd_contigs[i])) %>%
        full_join(., VDJ_T_GD_manual_clntypes, by=c("raw_clonotype_id"="clonotype_id")) %>% 
        # distinct(., .keep_all = TRUE) %>%
        arrange(desc(frequency)) # %>% dplyr::select(-c(inkt_evidence, mait_evidence))
      
      # Now we have a per barcode list of clonotypes with all associated sequencing data, from the given primer set, for cloning
      VDJ_T_GD.fulljoin.list[[i]] = VDJ_T_GD.fulljoin
    }
    
  }

  if (length(try_errors) == 3) {
    return("No results found for each primer set alignment...")
  }
  
  print(paste0("successfully processed contig files into clonotypes for ", 
               paste0(try_error_on_primer_set[which(try_error_on_primer_set %ni% try_errors)], collapse = ", "),
               " alignments. "))

  # Merge results from each of the 3 alignments together
  VDJ_T_GD.fulljoin = rbind(VDJ_T_GD.fulljoin.list[[1]], VDJ_T_GD.fulljoin.list[[2]], VDJ_T_GD.fulljoin.list[[3]])
  print("Dimensions of full rbind between M, G, and GM, clonotypes:")
  print(dim(VDJ_T_GD.fulljoin))
  
  # Re-level / re-name colonotypes based on cdr3 nucleotide sequence, since in each alignment the exact clonotype name 
  # can be different despite being the same cdr3 i.e. the same "ground truth clonotype" 
  
  print("determining unique clonotypes between the 3 alignments...")
  VDJ_T_GD.fulljoin.unique = as.data.frame(VDJ_T_GD.fulljoin %>%
                                             left_join(., y = as.data.frame(table(VDJ_T_GD.fulljoin$clonotype_key)) %>%
                                                         arrange(desc(Freq)) %>%
                                                         dplyr::mutate(clonotype_id = paste0("clonotype", rownames(.))) %>%
                                                         dplyr::select(-Freq), 
                                                       by = c("clonotype_key" = "Var1") ) %>%
                                             dplyr::select(barcode, clonotype_id, cdr3, cdr3_nt, clonotype_key, full_sequence) %>% # , inkt_evidence, mait_evidence
                                             distinct(.keep_all = T) %>% # remove any 100% duplicate rows, don't do by barcode
                                             group_by(barcode) %>%
                                             mutate(n_diff_clonotypes = length(barcode)) 
  )
  
  VDJ_T_GD.fulljoin.unique.save = VDJ_T_GD.fulljoin.unique %>%
    mutate(source = case_when(clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "All",
                              clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "G",
                              clonotype_key %ni% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "M",
                              clonotype_key %ni% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "C",
                              clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "MG_notC",
                              clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "GC_notM",
                              clonotype_key %ni% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "MC_notG")
  )
    
  
  print("sizes of VDJ_T_GD.fulljoin.unique before full_sequence removal:")
  print(dim(VDJ_T_GD.fulljoin.unique))
  print(table(VDJ_T_GD.fulljoin.unique$n_diff_clonotypes))
  
  if (missing(save_dir)) {
    write.table(VDJ_T_GD.fulljoin.unique.save,
                paste0("/path/to/save/dir/", 
                       sample, "_clonotypes_merged_alignments.txt"), sep = "\t", quote = F)
    print(paste0("writing all clonotypes to file ... /path/to/save/dir/",
                 sample, "_clonotypes_merged_alignments.txt"))
  } else {
    write.table(VDJ_T_GD.fulljoin.unique.save,
                paste0(save_dir, sample, "_clonotypes_merged_alignments.txt"), sep = "\t", quote = F)
    print(paste0("writing all clonotypes to file ... ", save_dir, 
                 sample, "_clonotypes_merged_alignments.txt"))
    
  }
  
  if (save_only == T) {
    return(paste0("Finished saving primer-merged gd clonotypes for sample: ", sample, ". Returning."))
  }
  
  
  # Now remove full_sequence column since it doesnt define clonotype
  VDJ_T_GD.fulljoin.unique = VDJ_T_GD.fulljoin.unique %>%
    dplyr::select(barcode, clonotype_id, cdr3, cdr3_nt, clonotype_key) %>% 
    distinct(.keep_all = T) %>% 
    group_by(barcode) %>%
    mutate(n_diff_clonotypes = length(barcode)) 
  
  print("sizes of VDJ_T_GD.fulljoin.unique after full_sequence removal:")
  print(dim(VDJ_T_GD.fulljoin.unique))
  print(table(VDJ_T_GD.fulljoin.unique$n_diff_clonotypes))
  
  # Need to loop through the above dataframe and in each case where there are multiple cdrs3/clonotype names for a barcode 
  # (i.e. n_diff_clonotypes > 1) we need to choose 1 (for integration with Seurat) and add a comment. 
  
  if (length(unique(VDJ_T_GD.fulljoin.unique$n_diff_clonotypes)) > 1) {
    VDJ_T_GD.fulljoin.unique.multibc = VDJ_T_GD.fulljoin.unique %>%
      filter(n_diff_clonotypes > 1) %>%
      dplyr::select(-n_diff_clonotypes)
    collapsed_VDJ_T_GD.fulljoin.unique.multibc = c()
    
    unique_barcodes = unique(VDJ_T_GD.fulljoin.unique.multibc$barcode)
    
    for (i in 1:length(unique_barcodes)) { 
      barcode_to_group = VDJ_T_GD.fulljoin.unique.multibc %>%
        filter(barcode == unique_barcodes[i]) %>%
        mutate(chains = case_when( grepl("TRD", cdr3_nt) &  grepl("TRG", cdr3_nt) ~ "TRD:TRG",
                                   grepl("TRD", cdr3_nt) &  !grepl("TRG", cdr3_nt) ~ "TRD",
                                   !grepl("TRD", cdr3_nt) &  grepl("TRG", cdr3_nt) ~ "TRG")
        ) %>%
        mutate(chains_priority = case_when(chains == "TRD:TRG" & str_count(cdr3_nt, "TRG|TRD") > 2 ~ 20,
                                           chains == "TRD:TRG" & str_count(cdr3_nt, "TRG|TRD") == 2 ~ 10,
                                           chains == "TRD" ~ 0,
                                           chains == "TRG" ~ 0)
        ) %>%
        mutate(source = case_when(clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "All",
                                  clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "G",
                                  clonotype_key %ni% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "M",
                                  clonotype_key %ni% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "C",

                                  clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "MG_notC",
                                  clonotype_key %in% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %ni% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "GC_notM",
                                  clonotype_key %ni% VDJ_T_GD.fulljoin.list[[1]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[2]]$clonotype_key & clonotype_key %in% VDJ_T_GD.fulljoin.list[[3]]$clonotype_key ~ "MC_notG")
        ) %>%
        mutate(source_priority = case_when(source == "All" ~ chains_priority + 4,
                                           source %in% c("MG_notC") ~ chains_priority + 3,  
                                           source %in% c("GC_notM", "MC_notG") ~ chains_priority + 2, 
                                           source %in% c("G", "M") ~ chains_priority + 1,
                                           source == "C" ~ chains_priority + 0
        )
        ) %>%
        mutate(source_priority = case_when(as.numeric(gsub("clonotype", "", clonotype_id)) == max(as.numeric(gsub("clonotype", "", clonotype_id))) ~ source_priority + 1,
                                           as.numeric(gsub("clonotype", "", clonotype_id)) != max(as.numeric(gsub("clonotype", "", clonotype_id))) ~ source_priority)
        ) %>%
        mutate(clonotypes = paste(unique(clonotype_id), collapse = ";"),
               sources = paste(unique(source), collapse = ";"),
               source_prioritys = paste((source_priority), collapse = ";")) %>%
        filter(source_priority == max(source_priority)) %>%
        mutate(comment = case_when(length(source_priority) == 1 ~ paste0("clonotypes:",clonotypes,
                                                                         "|scores:",source_prioritys,
                                                                         "|sources:",sources),
                                   length(source_priority) > 1 ~ paste0("clonotypes:",clonotypes,
                                                                        "|scores:",source_prioritys,
                                                                        "|sources:",sources,
                                                                        "|otherCdr3sTied:", paste(clonotype_key, collapse = ";")))
        ) %>%
        dplyr::select(-chains, -chains_priority, -source, -source_priority, -clonotypes, -sources, -source_prioritys)
      
      if (verbose) {
        print(paste0("finished collapsing barcode: ", unique_barcodes[i]))
        print("attaching:")
        print(as.data.frame(barcode_to_group))
      }
            
      collapsed_VDJ_T_GD.fulljoin.unique.multibc = rbind(collapsed_VDJ_T_GD.fulljoin.unique.multibc, barcode_to_group[1,]) # adding [1,] to the end ensures only 1 entry is added in case of priority ties
    }
    
    
    print("finished collapsing multi-clonotype barcodes...")
    print(paste0("# barcodes = ", length(collapsed_VDJ_T_GD.fulljoin.unique.multibc$barcode)))
    print(paste0("# unique barcodes = ", length(unique(collapsed_VDJ_T_GD.fulljoin.unique.multibc$barcode))))
    
    # Now combine all the results back together
    VDJ_T_GD.fulljoin.unique = rbind(VDJ_T_GD.fulljoin.unique %>%
                                       filter(n_diff_clonotypes == 1) %>%
                                       dplyr::select(-n_diff_clonotypes) %>%
                                       mutate(comment = ""),
                                     collapsed_VDJ_T_GD.fulljoin.unique.multibc)
    
    print("finished merging collapsed barcodes and unique barcodes...")
    print(paste0("# barcodes = ", length(VDJ_T_GD.fulljoin.unique$barcode)))
    print(paste0("# unique barcodes = ", length(unique(VDJ_T_GD.fulljoin.unique$barcode))))
    
  } else {
    VDJ_T_GD.fulljoin.unique = VDJ_T_GD.fulljoin.unique %>%
      dplyr::select(-n_diff_clonotypes) %>%
      mutate(comment = "")
  }
  
  stopifnot(length(VDJ_T_GD.fulljoin.unique$barcode) == length(unique(VDJ_T_GD.fulljoin.unique$barcode)))
  #stopifnot(length(VDJ_T_GD.fulljoin$barcode) == length(unique(VDJ_T_GD.fulljoin$barcode)))
  print("successfully flattened list of clonotypes by barcode, returning barcode level clonotype results for seurat...")
  
  VDJ_T_GD.fulljoin.unique = as.data.frame(VDJ_T_GD.fulljoin.unique)
  colnames(VDJ_T_GD.fulljoin.unique) = c("barcode", "raw_clonotype_id", "cdr3s_aa", "cdr3s_nt", "clonotype_key", "comment")
  return(VDJ_T_GD.fulljoin.unique)
  
}


split_gdTCR = function(seq, sep) {
  return(split_dual_TCR(seq, sep, chain = "gd"))
}


split_dual_TCR = function(seq, sep, chain = "gd") {
  # quick function that takes a string input of multi gdTCRs and returns a dataframe of possible two-way combinations  
  # seq example = TRD:CALGEHPYHLYWGFTTDKLIF;TRG:CATWDGYYYKKLF;TRG:CATWDRFYYKKLF
  # sep = character to separate the gdTCR sequences by, e.g. ";"
  split_TCRs = str_split(seq, sep)
  
  if (length(unlist(split_TCRs)) == 1) {
    possible_TCRs = data.frame(TCRs = unlist(split_TCRs))
    return(possible_TCRs)
  }
  
  if (chain == "gd") {
    possible_TCRs = as.data.frame(paste0(expand.grid(unlist(split_TCRs), unlist(split_TCRs))$Var1, ";", expand.grid(unlist(split_TCRs), unlist(split_TCRs))$Var2)) %>%
      dplyr::rename(., TCRs = names(.)[1]) %>%
      filter(grepl("^TRD:", TCRs) == T & grepl(";TRG:", TCRs) == T) %>%
      distinct()
  } else if (chain == "ab") {
    possible_TCRs = as.data.frame(paste0(expand.grid(unlist(split_TCRs), unlist(split_TCRs))$Var1, ";", expand.grid(unlist(split_TCRs), unlist(split_TCRs))$Var2)) %>%
      dplyr::rename(., TCRs = names(.)[1]) %>%
      filter(grepl("^TRB:", TCRs) == T & grepl(";TRA:", TCRs) == T) %>%
      distinct()
  } else {
    return("please specify a valid chain")
  }
  return(possible_TCRs)
}

split_by_chain = function(df1, chain, clonoColumn) {
  
  if (missing(chain) | missing(clonoColumn)) {
    print("In split_by_chain, chain and/or clonoColumn not specified, returning unchanged df.")
    return(df1)
  }
  
  if (chain %in% c("gd", "ab")) {
    return(df1)
  }
  
  chain_split = case_when(chain == "g" ~ "TRG:", chain == "d" ~ "TRD:", chain == "a" ~ "TRA:", chain == "b" ~ "TRB:")
  for (i in 1:nrow(df1)) {
    if (grepl(chain_split, df1[i,clonoColumn])) {
      split_clono = str_split(df1[i,clonoColumn], ";")[[1]]
      split_clono = paste0(split_clono[grepl(chain_split, split_clono)], collapse = ";")
    } else {
      split_clono = NA
    }
    df1[i,clonoColumn] = split_clono
  }
  
  df1 = subset(df1, !is.na(df1[,clonoColumn]))
  return(df1)
  
}

match_dual_tcrs = function(seurat_metadata) {
  # This function will take a seurat metadata object containing VDJ_T_GD columns and check whether any clonotypes can 
  # be combined in the case where dual TCRs are ordered differently, e.g. TRDX:TRGY:TRGZ and TRDX:TRGZ:TRGY should be 
  # counted as the same clonotype, so long as they have the same gene usage in the clonotype key column. 
  
  starting_rownames = rownames(seurat_metadata)
  starting_rows = nrow(seurat_metadata)
  starting_cols = ncol(seurat_metadata)
  
  # Identify which metadata rows have dual TCRs
  dual_tcr_rows = c()
  for (i in 1:nrow(seurat_metadata)) {
    if (!is.na(seurat_metadata$VDJ_T_GD_cdr3s_aa[i])) {
      if (str_count( seurat_metadata$VDJ_T_GD_cdr3s_aa[i], "TRD:") > 1 | str_count(seurat_metadata$VDJ_T_GD_cdr3s_aa[i], "TRG:") > 1) {
        dual_tcr_rows = c(dual_tcr_rows, i)
      }
    }
  }
  
  if (is.null(dual_tcr_rows)) {
    print("no dual TCRs were found, returning metadata unchanged.")
    return(seurat_metadata)
  } else {
    print(paste0("found ", length(dual_tcr_rows), " dual TCRs"))
  }
  
  # Identify which dual TCRs could be merged if they were re-ordered
  unique_dual_tcrs = unique(seurat_metadata[dual_tcr_rows,"VDJ_T_GD_cdr3s_aa"])
  
  # Testing data:
  # unique_dual_tcrs = c("TRD:CALGELTFLRLLYWGIDSRPLIF;TRG:CATWDVGYSNYYKKLF;TRG:CATWDTNYYKKLF;TRG:CATWDRQRARKKLF" ,
  #                      "TRD:CALGELTFLRLLYWGIDSRPLIF;TRG:CATWDRQRARKKLF;TRG:CATWDTNYYKKLF;TRG:CATWDVGYSNYYKKLF",
  #                      "TRD:CACDLRLPFGKWLEPYWGIATDKLIF;TRG:CATWDGRCDYKKLF;TRG:CATWDRQRARKKLF",
  #                      "TRD:CACDLRLPFGKWLEPYWGIATDKLIF;TRG:CATWDRQRARKKLF;TRG:CATWDGRCDYKKLF",
  #                      "TRD:CALGELTFLRLLYWGIDSRPLIF;TRG:CATWDRQRARKKLF;TRG:CATWDVGYSNYYKKLF;TRG:CATWDTNYYKKLF",
  #                      "TRD:CALGETALSLIRFWGISVFLGPLIF;TRG:CALWELYYYKKLF;TRG:CATWDVVKLF")
  
  unique_dual_tcrs_remaining = unique_dual_tcrs
  matches = c()
  
  for (i in 1:length(unique_dual_tcrs)) {
    split_TCR_query = split_dual_TCR(unique_dual_tcrs[i], ";")
    
    if (unique_dual_tcrs[i] %in% unique_dual_tcrs_remaining) {
      unique_dual_tcrs_remaining = unique_dual_tcrs_remaining[-which(unique_dual_tcrs_remaining == unique_dual_tcrs[i])]
    }
    
    match_indices = c()
    if (length(unique_dual_tcrs_remaining) > 0) {
      for (j in 1:length(unique_dual_tcrs_remaining)) {
        split_TCR_other =  split_dual_TCR(unique_dual_tcrs_remaining[j], ";")
        if (all(split_TCR_query$TCRs %in% split_TCR_other$TCRs) & all(split_TCR_other$TCRs %in% split_TCR_query$TCRs)) {
          matches = rbind(matches, data.frame(Query = unique_dual_tcrs[i],
                                              Match = unique_dual_tcrs_remaining[j]))
          match_indices = c(match_indices, j)
          print(paste0("made a match i = ", i, " j = ", j))
        }
      }
    }
    if (!is.null(match_indices)) {
      unique_dual_tcrs_remaining = unique_dual_tcrs_remaining[-match_indices]
    }
  }
  
  if (!is.null(matches)) {
    # Combine all matches for a given query
    matches = matches %>%
      group_by(Query) %>% 
      mutate(Match = paste0(Match, collapse = "|")) %>% 
      distinct() %>% 
      as.data.frame()
    
    # Now for each match, go through and update metadata
    for (i in 1:nrow(matches)) {
      update_cols = c("orig.ident", "VDJ_T_GD_raw_clonotype_id",
                      "VDJ_T_GD_frequency", "VDJ_T_GD_proportion",
                      "VDJ_T_GD_cdr3s_nt", "VDJ_T_GD_clonotype_key")
      seurat_metadata$Matched = NA
      query_values = c()
      match_values = c()
      # Change CDR3aa column to unified match CDR3aa and flag a change was made
      for (j in 1:nrow(seurat_metadata)) {
        if (any(seurat_metadata$VDJ_T_GD_cdr3s_aa[j] %in% as.character(str_split_fixed(matches[i,"Match"], "\\|", Inf)))) {
          print(paste0("found match at j = ", j, " for i = ", i))
          seurat_metadata$VDJ_T_GD_cdr3s_aa[j] = matches[i,"Query"]
          seurat_metadata$Matched[j] = i
          match_values = rbind(match_values, seurat_metadata[j, update_cols])
        } else if (seurat_metadata$VDJ_T_GD_cdr3s_aa[j] %in% matches[i,"Query"]) {
          seurat_metadata$Matched[j] = i
          query_values = rbind(query_values, seurat_metadata[j, update_cols])
          rownames(query_values) = NULL
        }
      }
      
      # Update other related columns based on flags 
      match_values = match_values %>% distinct()
      query_values = query_values %>% distinct()
      
      for (j in 1:nrow(query_values)) {
        if (query_values$orig.ident[j] %in% match_values$orig.ident) {
          query_values$VDJ_T_GD_frequency[j] = query_values$VDJ_T_GD_frequency[j] + match_values[which(match_values$orig.ident == query_values$orig.ident[j]),"VDJ_T_GD_frequency"]
          query_values$VDJ_T_GD_proportion[j] = query_values$VDJ_T_GD_proportion[j] + match_values[which(match_values$orig.ident == query_values$orig.ident[j]),"VDJ_T_GD_proportion"]
        }
      }
      
      for (j in 1:nrow(seurat_metadata)) {
        if (!is.na(seurat_metadata$Matched[j])) {
          if (seurat_metadata$orig.ident[j] %in% query_values$orig.ident) {
            seurat_metadata[j, update_cols] = query_values[which(query_values$orig.ident == seurat_metadata$orig.ident[j]),]
          } else {
            stopifnot("query_values CDR3nt sequence not unique" = length(unique(query_values$VDJ_T_GD_cdr3s_nt)) == 1)
            stopifnot("query_values clonotype_key sequence not unique" = length(unique(query_values$VDJ_T_GD_clonotype_key)) == 1)
            
            seurat_metadata[j, "VDJ_T_GD_cdr3s_nt"] = unique(query_values$VDJ_T_GD_cdr3s_nt)
            seurat_metadata[j, "VDJ_T_GD_clonotype_key"] = unique(query_values$VDJ_T_GD_clonotype_key)
            
          }
        }
      }
      
    }
    
    seurat_metadata = seurat_metadata %>% 
      dplyr::select(-Matched) %>% 
      as.data.frame()
    
    stopifnot("Before final return, rownames did not match" = rownames(seurat_metadata) == starting_rownames)
    stopifnot("Before final return, nrow did not match" = nrow(seurat_metadata) == starting_rows)
    stopifnot("Before final return, ncol did not match" = ncol(seurat_metadata) == starting_cols)
    return(seurat_metadata)
    
  } else {
    print("No dual TCRs were matched, returning metadata unchanged.")
    return(seurat_metadata)
  }
  
}

Clone_distribution <- function(clone_df, cols = NA, cols_overide = NA, radius_factor = 0.8, cell_size = 2) {
  
  require(ggplot2)
  require(ggthemes)
  require(packcircles) 
  require(dplyr)
  set.seed(123)
  
  if (ncol(clone_df) == 2) {
    color_column = F
  } else if (ncol(clone_df) == 3) {
    color_column = T
  } else {
    return("please supply a dataframe of either 2 or 3 columns corresponding to clone id, cell id, and (optionally) color")
  }
  
  if (color_column == T) {
    color_df = clone_df[,c(2:3)]
    colnames(color_df) <- c('cell_id', 'color')
    
    # Palette: 
    if (is.numeric(color_df$color)) {
      cols_pal = scale_color_gradient2(low = "grey", high = "blue")
    } else {
      color_df$color = as.character(color_df$color)
      if (!missing(cols)) {
        cols_pal = scale_color_manual(values = cols)
      } else {
        cols_pal = scale_color_manual(values = gg_color_hue(length(unique(color_df$color))))
      }
    }
  }
  
  if (!missing(cols_overide)) {
    cols_pal = cols_overide
  }
  
  clone_df = clone_df[,c(1:2)]
  colnames(clone_df) <- c('clone_id', 'cell_id')
  
  clone_sizes <- clone_df %>% dplyr::count(clone_id) %>% arrange(desc(n))
  
  packing <- circleProgressiveLayout(clone_sizes$n, sizetype = 'area')
  packing$radius <- packing$radius * seq(1, 0.8, length.out = nrow(packing))
  packing$clone_id <- clone_sizes$clone_id
  vertices <- circleLayoutVertices(packing) %>%
    mutate(clone_id = rep(packing$clone_id, each = nrow(.) / length(packing$clone_id)))
  
  clone_positions <- clone_sizes %>%
    mutate(x_center = packing$x, y_center = packing$y, radius = packing$radius)
  
  cell_positions <- list()
  
  for(i in 1:nrow(clone_positions)) {
    Clone_id <- clone_positions$clone_id[i]
    clone_size <- clone_positions$n[i]
    x_center <- clone_positions$x_center[i]
    y_center <- clone_positions$y_center[i]
    radius <- clone_positions$radius[i]
    
    cell_radius <- sqrt((radius^2) / clone_size) * radius_factor 
    
    cell_coords <- data.frame()
    current_radius <- cell_radius
    cells_added <- 0
    
    while (cells_added < clone_size) {
      num_cells_in_layer <- ceiling(2 * pi * current_radius / (2 * cell_radius))
      angles <- seq(0, 2 * pi, length.out = num_cells_in_layer + 1)[-1]
      
      layer_coords <- data.frame(
        x = x_center + current_radius * cos(angles),
        y = y_center + current_radius * sin(angles),
        clone_id = Clone_id,
        cell_radius = cell_radius
      )
      
      if (cells_added + nrow(layer_coords) > clone_size) {
        layer_coords <- layer_coords[1:(clone_size - cells_added), ]
      }
      
      cell_coords <- rbind(cell_coords, layer_coords)
      cells_added <- nrow(cell_coords)
      current_radius <- current_radius + 2 * cell_radius  
    }
    
    cell_coords$cell_id <- filter(clone_df, clone_id == Clone_id) %>% pull(cell_id)
    cell_positions[[i]] <- cell_coords
  }
  
  all_cells <- bind_rows(cell_positions) 
  
  if (nrow(all_cells) != nrow(clone_df)) {
    stop(" check here")
  }
  
  
  if (color_column == T) {
    all_cells = all_cells %>% 
      left_join(color_df, by = "cell_id")
    
    p <- ggplot() +
      geom_point(data = all_cells, aes(x, y, color = color), size = all_cells$cell_radius * cell_size) +
      coord_equal() + theme_void() + cols_pal +
      theme(legend.position = "bottom") + labs(color = "")
    
  } else if (color_column == F) {
    p <- ggplot() +
      geom_point(data = all_cells, aes(x, y, color = clone_id), size = all_cells$cell_radius * cell_size) +
      # geom_polygon(data = subset(vertices, id<11), aes(x, y, group = id), fill=NA, color = "black") +
      coord_equal() + theme_void() +
      theme(legend.position = "none") 
  } 
  
  # print(p)
  return(list(plot = p, coords=all_cells, circle=vertices))
}

