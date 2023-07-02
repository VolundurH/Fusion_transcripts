# Process raw FusionCatcher output for machine learning for both TCGA and SCAN-B data sets.
# FPKM for each fusion gene is added w/ expression matrices obtained from both cohorts.
fusioncatcher <- read_tsv("fusioncatcher_all_cancer_types_filtered_with_fpkm.txt")


# Repeatmasker information -------------------------------------------------------------------------
# obtained from http://repeatmasker.org/ - RepeatMasker open-4.0.5 
rmsk <- read_tsv("rmsk.txt", col_names = c("bin","swScore","milliDiv","milliDel", "milliIns","genoName","genoStart","genoEnd", "genoLeft","strand","repName","repClass","repFamily","repStart", "repEnd","repLeft","id")) %>% 
  mutate(genoName = str_remove(genoName, "chr")) %>% 
  filter(nchar(genoName) < 3) 

rmsk_ranges <- GenomicRanges::GRanges(seqnames = rmsk$genoName, ranges = IRanges::IRanges(start = rmsk$genoStart, end =rmsk$genoEnd))


# Look at overlap with fusioncatcher junctions (hg38)

fiveprime_junction_rmsk_overlap <- GenomicRanges::GRanges(seqnames =  fusioncatcher$fiveprime_chr, ranges = IRanges::IRanges(start = fusioncatcher$fiveprime_junction_hg38, width = 1)) %>% 
  GenomicRanges::findOverlaps(rmsk_ranges) %>% 
  as_tibble() %>% 
  rename(fusioncatcher_row_id = queryHits, rmsk_row_id = subjectHits)

threeprime_junction_rmsk_overlap <- GenomicRanges::GRanges(seqnames =  fusioncatcher$threeprime_chr, ranges = IRanges::IRanges(start = fusioncatcher$threeprime_junction_hg38, width = 1)) %>% 
  GenomicRanges::findOverlaps(rmsk_ranges) %>% 
  as_tibble() %>% 
  rename(fusioncatcher_row_id = queryHits, rmsk_row_id = subjectHits)

test <- fusioncatcher %>% 
  mutate(row_id = row_number()) %>% 
  mutate(barcode = str_sub(barcode, 1, 15)) %>% 
  left_join(fiveprime_junction_rmsk_overlap, by = c("row_id" = "fusioncatcher_row_id")) %>% 
  left_join(
    rmsk %>%
      mutate(row_id = row_number()) %>%
      select(row_id, repName, repClass, repFamily),
    by = c("rmsk_row_id" = "row_id")) %>% 
  rename(fiveprime_repname = repName, fiveprime_repclass = repClass, fiveprime_repfamily = repFamily) %>% 
  select(-rmsk_row_id) %>% 
  left_join(threeprime_junction_rmsk_overlap, by = c("row_id" = "fusioncatcher_row_id")) %>% 
  left_join(
    rmsk %>%
      mutate(row_id = row_number()) %>%
      select(row_id, repName, repClass, repFamily),
    by = c("rmsk_row_id" = "row_id")) %>% 
  rename(threeprime_repname = repName, threeprime_repclass = repClass, threeprime_repfamily = repFamily) %>% 
  select(-rmsk_row_id) 

fusion_junction_repeat_families <- test %>% 
  select(fusion_id, fiveprime_repfamily, threeprime_repfamily) %>% 
  filter(!(is.na(fiveprime_repfamily) & is.na(threeprime_repfamily))) %>% 
  group_by(fusion_id) %>% 
  summarise(fiveprime_repfamily = paste0(fiveprime_repfamily, collapse = ","), 
            threeprime_repfamily = paste0(threeprime_repfamily, collapse = ","))

fusion_junction_repeat_families <- fusion_junction_repeat_families %>% 
  mutate(across(c(fiveprime_repfamily, threeprime_repfamily), ~ifelse(.x == "NA", NA, .x)))

fusioncatcher <- fusioncatcher %>% 
  left_join(fusion_junction_repeat_families, by = "fusion_id")

fusioncatcher <- fusioncatcher %>% 
  mutate(across(.cols = c(fiveprime_repfamily, threeprime_repfamily), ~replace_na(.x, "None"))) 


# add miRNA and snoRNA host status -----------------------------------------------------------------------
# https://mirbase.org/ftp.shtml
# 

miRNA_hosts <- read_tsv("miRNA_hosts.txt") %>% 
  distinct(host_ID) %>% 
  pull(host_ID)

snoRNA_hosts <- read_tsv("snoRNA_host_genes.txt") %>% 
  distinct(host_ID) %>% 
  as_vector()

fusioncatcher <- fusioncatcher %>% 
  mutate(
    fiveprime_is_mirna_host = fiveprimepartner%in%miRNA_hosts,
    threeprime_is_mirna_host = threeprimepartner%in%miRNA_hosts,
    fiveprime_is_snorna_host = fiveprimepartner%in%snoRNA_hosts,
    threeprime_is_snorna_host = threeprimepartner%in%snoRNA_hosts) 


# add multi-junction information ---------------------------------------------------------------------
# if we find multiple FC events for the same gene combo in the same sample it is classified 
# as a multi-junction event

fusion_multi_junction_status <- fusioncatcher %>% 
  select(fusion_id, barcode, fiveprimepartner, threeprimepartner, fiveprime_junction_hg38, threeprime_junction_hg38) %>% 
  group_by(barcode, fiveprimepartner, threeprimepartner) %>% 
  nest() %>% 
  mutate(multi_juncion = map_lgl(data, ~nrow(.x)>1)) %>% 
  unnest() %>% 
  select(fusion_id, multi_juncion)

fusioncatcher <- fusioncatcher %>% 
  left_join(fusion_multi_junction_status %>% ungroup %>%  distinct(fusion_id, multi_juncion))


# Add reciprocal fusion status --------------------------------------------
# only label fusion as reciprocal if there are two non-blacklisted fusions present. 

fusion_reciprocal_status <- fusioncatcher %>% 
  select(fusion_id, barcode, fiveprimepartner, threeprimepartner) %>% 
  mutate(key = paste0(fiveprimepartner, threeprimepartner),
         key_rev = paste0(threeprimepartner, fiveprimepartner)) %>% 
  group_by(barcode) %>% 
  nest() %>% 
  mutate(data = map(data, ~mutate(.x, is_reciprocal = key%in%key_rev ))) %>% 
  unnest(data) %>% 
  select(fusion_id, is_reciprocal)

fusioncatcher <- fusioncatcher %>% 
  left_join(fusion_reciprocal_status %>% ungroup %>% distinct(fusion_id, is_reciprocal))

fusioncatcher <- fusioncatcher %>% 
  mutate(is_normal = str_sub(barcode, 14, 14) == 1)


# label a fusion partner if it is common (top 10) in normal tissue for that tissue type ------------------------

top10_normal <- fusioncatcher %>% 
  select(fusion_id, cancer_type, fiveprimepartner, threeprimepartner, is_normal) %>% 
  filter(is_normal) %>% 
  pivot_longer(cols = c(fiveprimepartner, threeprimepartner), names_to = "partner", values_to = "gene") %>% 
  group_by(cancer_type, gene) %>% 
  tally(sort = T) %>% 
  group_by(cancer_type) %>% 
  slice_max(order_by = n, n = 10) %>% 
  ungroup() %>% 
  select(-n) %>% 
  mutate(is_top10_normal = T)

fusioncatcher <- fusioncatcher %>% 
  left_join(top10_normal, by = c("cancer_type", "fiveprimepartner" = "gene")) %>% 
  rename(fiveprime_is_top10_normal = is_top10_normal) %>% 
  left_join(top10_normal, by = c("cancer_type", "threeprimepartner" = "gene")) %>% 
  rename(threeprime_is_top10_normal = is_top10_normal) %>% 
  mutate(across(c(fiveprime_is_top10_normal, threeprime_is_top10_normal), ~ifelse(is.na(.x), F, .x)))

no_normals <- fusioncatcher %>%
  count(cancer_type, is_normal) %>%
  pivot_wider(names_from = "is_normal", values_from = "n") %>%
  filter(is.na(`TRUE`)) %>% pull(cancer_type)

fusioncatcher <- fusioncatcher %>% 
  mutate(across(c(fiveprime_is_top10_normal, threeprime_is_top10_normal),
                ~ifelse(cancer_type %in% no_normals, NA, .x))) 

# label a fusion partner if it is common (top 10) in tumor samples for that tissue type ------------------------

top10_tumor <- fusioncatcher %>% 
  select(fusion_id, cancer_type, fiveprimepartner, threeprimepartner, is_normal) %>% 
  filter(!is_normal) %>% 
  pivot_longer(cols = c(fiveprimepartner, threeprimepartner), names_to = "partner", values_to = "gene") %>% 
  group_by(cancer_type, gene) %>% 
  tally(sort = T) %>% 
  group_by(cancer_type) %>% 
  slice_max(order_by = n, n = 10) %>% 
  ungroup() %>% 
  select(-n) %>% 
  mutate(is_top10_tumor = T)

fusioncatcher <- fusioncatcher %>% 
  left_join(top10_tumor, by = c("cancer_type", "fiveprimepartner" = "gene")) %>% 
  rename(fiveprime_is_top10_tumor = is_top10_tumor) %>% 
  left_join(top10_tumor, by = c("cancer_type", "threeprimepartner" = "gene")) %>% 
  rename(threeprime_is_top10_tumor = is_top10_tumor) %>% 
  mutate(across(c(fiveprime_is_top10_tumor, threeprime_is_top10_tumor), ~ifelse(is.na(.x), F, .x)))


# Add same strand infomation, same chr information ----------------------------------------------

fusioncatcher <- fusioncatcher%>% 
  mutate(is_same_strand = fiveprime_strand == threeprime_strand, 
         is_same_chr = fiveprime_chr == threeprime_chr)


# Add cosmic tags ---------------------------------------------------------

cosmic <- read_csv("cancer_gene_census.csv") %>% 
  select(hallmark = Hallmark, role = `Role in Cancer`, synonyms = Synonyms) %>% 
  mutate(gene_id = str_extract(synonyms, "ENSG\\d+")) %>% 
  select(-synonyms) 

fusioncatcher <- fusioncatcher %>% 
  left_join(cosmic, by = c("fiveprimepartner" = "gene_id")) %>% 
  rename(fiveprime_cosmic_hallmark = hallmark, 
         fiveprime_cosmic_role = role) %>% 
  left_join(cosmic,  by = c("threeprimepartner" = "gene_id")) %>% 
  rename(threeprime_cosmic_hallmark = hallmark, 
         threeprime_cosmic_role = role) %>% 
  mutate(across(.cols = c(fiveprime_cosmic_hallmark, threeprime_cosmic_hallmark), ~replace_na(.x, "No"))) %>% 
  mutate(across(.cols = c(fiveprime_cosmic_role, threeprime_cosmic_role), ~replace_na(.x, "None"))) 

fusioncatcher <- fusioncatcher %>% 
  mutate(across(c(fiveprime_cosmic_hallmark, threeprime_cosmic_hallmark), ~ifelse(.x == "Yes", T, F))) 

# add spatial information -------------------------------------------------

# hg19 fusion junctions
# hg38 junction distance (if same chr)
# hg19 junction distance (if same chr)

fusioncatcher %>% 
  select(fiveprime_chr, fiveprime_junction_hg38, fusion_id) %>% 
  mutate(fiveprime_chr = paste0("chr", fiveprime_chr)) %>% 
  mutate(end = fiveprime_junction_hg38) %>% 
  select(fiveprime_chr, fiveprime_junction_hg38, end, fusion_id) %>% 
  write_tsv("fusioncatcher_fiveprime_junction_liftover_hg38.bed", col_names = F)

fusioncatcher %>% 
  select(threeprime_chr, threeprime_junction_hg38, fusion_id) %>% 
  mutate(threeprime_chr = paste0("chr", threeprime_chr)) %>% 
  mutate(end = threeprime_junction_hg38) %>% 
  select(threeprime_chr, threeprime_junction_hg38, end, fusion_id) %>% 
  write_tsv("fusioncatcher_threeprime_junction_liftover_hg38.bed", col_names = F)

fiveprime_junctions_hg19 <- read_tsv("fusioncatcher_fiveprime_junction_liftover_hg19.bed", 
                                     col_names = c("chr", "start", "coord", "fusion_id", "misc")) %>% 
  select(fusion_id, chr, coord) %>% 
  mutate(chr = str_remove(chr, "chr"))

threeprime_junctions_hg19 <- read_tsv("fusioncatcher_threeprime_junction_liftover_hg19.bed", 
                                     col_names = c("chr", "start", "coord", "fusion_id", "misc")) %>% 
  select(fusion_id, chr, coord) %>%
  mutate(chr = str_remove(chr, "chr"))


fusioncatcher <- fusioncatcher %>% 
  left_join(fiveprime_junctions_hg19 %>% distinct(), by = c("fusion_id", "fiveprime_chr"="chr")) %>% 
  rename(fiveprime_junction_hg19 = coord) %>% 
  left_join(threeprime_junctions_hg19 %>% distinct(), by = c("fusion_id", "threeprime_chr"="chr")) %>% 
  rename(threeprime_junction_hg19 = coord) %>% 
  mutate(junction_dist_hg38 = ifelse(
    fiveprime_chr == threeprime_chr,
    log2(abs(threeprime_junction_hg38 - fiveprime_junction_hg38)),
    NA)) %>% 
  mutate(junction_dist_hg19 = ifelse(
    fiveprime_chr == threeprime_chr,
    log2(abs(threeprime_junction_hg19 - fiveprime_junction_hg19)),
    NA)) %>% 
  mutate(junction_dist_diff = junction_dist_hg38 - junction_dist_hg19)

# 5' gene length # 3' gene length 

gencode <- read_tsv("gencode_genes_v27_hg38.txt") # https://www.gencodegenes.org/human/release_27.html

fusioncatcher <- fusioncatcher %>% 
  left_join(gencode %>% 
              select(gene_id, fiveprime_chr = chr, fiveprime_start = start, fiveprime_end = end) %>% 
              mutate(fiveprime_chr = str_remove(fiveprime_chr, "chr")), 
            by = c("fiveprimepartner" = "gene_id", "fiveprime_chr")) %>% 
  left_join(gencode %>% 
              select(gene_id, threeprime_chr = chr, threeprime_start = start, threeprime_end = end) %>% 
              mutate(threeprime_chr = str_remove(threeprime_chr, "chr")), 
            by = c("threeprimepartner" = "gene_id", "threeprime_chr")) %>% 
  mutate(fiveprie_gene_length = fiveprime_end - fiveprime_start,
         threeprime_gene_length = threeprime_end - threeprime_start)

# gene location 
chr_sizes <- read_tsv("chromosome_sizes.txt") %>% 
  mutate(chr = str_remove(chr, "chr")) %>% 
  rename(size = len)

fusioncatcher <- fusioncatcher %>% 
  left_join(chr_sizes, by = c("fiveprime_chr" = "chr")) %>% 
  mutate(fiveprime_rel_gene_location = fiveprime_start /  size  ) %>% 
  select(-size) %>% 
  left_join(chr_sizes, by = c("threeprime_chr" = "chr")) %>% 
  mutate(threeprime_rel_gene_location = threeprime_start /  size  ) %>% 
  select(-size)


# FPKM difference (log2)
fusioncatcher <- fusioncatcher %>% 
  mutate(fpkm_difference_log2 = log2(fiveprime_fpkm/threeprime_fpkm)) 


# N fusions in sample 

fusioncatcher <- fusioncatcher %>% 
  group_by(barcode) %>% 
  mutate(n_fusions_in_sample = n()) %>% 
  mutate(n_fusions_in_sample_bins = cut(n_fusions_in_sample,  c(0,10,50,100,500,Inf)))

# Junction location as percentage into gene -------------------------------


fusioncatcher <- fusioncatcher %>% 
  left_join(gencode %>% 
              select(gene_id, chr, start, end) %>% 
              mutate(chr = str_remove(chr, "chr")), by = c("fiveprimepartner" = "gene_id", "fiveprime_chr" = "chr")) %>% 
  mutate(fiveprime_junction_gene_location = round((fiveprime_junction_hg38 - start) / (end - start), 2 )) %>% 
  select(-c(start, end)) %>% 
  left_join(gencode %>% 
              select(gene_id, chr, start, end) %>% 
              mutate(chr = str_remove(chr, "chr")), by = c("threeprimepartner" = "gene_id", "threeprime_chr" = "chr")) %>%
  mutate(threeprime_junction_gene_location = round((threeprime_junction_hg38 - start) / (end - start), 2 )) %>% 
  select(-c(start, end)) 

# replace junction distance NAs with -1

fusioncatcher <- fusioncatcher %>% 
  mutate(across(c(junction_dist_diff, junction_dist_hg38, junction_dist_hg19), ~replace_na(.x, -1))) 

# add validation status ---------------------------------------------------

fusioncatcher <- fusioncatcher %>% 
  mutate(support_disc_reads = "unknown")
  
fusioncatcher <- fusioncatcher %>% 
  select(fusion_id, support_disc_reads, everything())


# Split Predicted effect --------------------------------------------------

fusioncatcher <- fusioncatcher %>% 
  mutate(Predicted_effect = ifelse(str_detect(Predicted_effect, "/"), Predicted_effect, paste0(Predicted_effect, "/", Predicted_effect))) %>% 
  separate(Predicted_effect, into=c("fiveprime_predicted_effect", "threeprime_predicted_effect"), sep = "/", remove = F) 



# move mX descriptor to own column ----------------------------------------


fusioncatcher <- fusioncatcher %>% 
  mutate(mX = str_extract(Fusion_description, "m\\d+") %>% str_remove("m") %>% as.integer()) %>%  # create int column mX, where values are X
  mutate(Fusion_description = str_remove_all(Fusion_description, "m\\d+,*")) %>%                  # rm all mX from original col
  mutate(Fusion_description = if_else(Fusion_description == "", "none", Fusion_description)) %>%  # add "none" if description empty
  mutate(Fusion_description = str_remove_all(Fusion_description, "^,|,$")) %>%                    # rm leading or trailing commas
  mutate(Fusion_description = str_replace(Fusion_description, ",{2,}", ",")) %>%                  # rm duplicated commas
  mutate(mX = replace_na(mX, -1))                                                                 # replace NA with -1


# Add gene type -----------------------------------------------------------

fusioncatcher <- fusioncatcher %>% 
  left_join(gencode %>% distinct(gene_id, gene_type),
            by = c("fiveprimepartner" = "gene_id")) %>% 
  rename(fiveprime_gene_type = gene_type) %>% 
  left_join(gencode %>% distinct(gene_id, gene_type),
            by = c("threeprimepartner" = "gene_id")) %>% 
  rename(threeprime_gene_type = gene_type) 

# remove unneeded columns

fusioncatcher <- fusioncatcher %>% 
  ungroup() %>% 
  select(-c(id_1, new_id, category, fiveprime_symbol, threeprime_symbol)) %>% 
  rename(sample = barcode) %>% 
  mutate(cohort = "TCGA")

# Add quantile FPKM values ------------------------------------------------------------------


quantiles <- read_tsv("~TCGA_expression_quantiles.txt") %>% 
  mutate(cancer_type = str_replace(cancer_type, "_", "-" )) %>% 
  mutate(gene_id = str_sub(gene_id, 1, 15))

quantiles <- quantiles %>% 
  distinct() %>% 
  group_by(cancer_type, gene_id) %>% 
  filter(quant_100 == max(quant_100))

submitter_ids <- fusioncatcher %>% 
  mutate(submitter_id = str_sub(sample, 1, 12)) %>% 
  distinct(submitter_id) %>% 
  as_vector()

project_ids <- GenomicDataCommons::cases() %>% 
  GenomicDataCommons::filter(submitter_id == submitter_ids) %>% 
  GenomicDataCommons::select(c("submitter_id", "project.project_id")) %>% 
  GenomicDataCommons::results_all() %>% 
  as_tibble() %>% 
  unnest(project) %>% 
  select(submitter_id, project_id)

fusioncatcher <- fusioncatcher %>% 
  mutate(submitter_id = str_sub(sample, 1, 12)) %>% 
  left_join(project_ids, by = "submitter_id") %>% 
  left_join(quantiles, by = c("project_id" = "cancer_type", "fiveprimepartner"="gene_id")) %>% 
  rename_with(~paste0("fiveprime_cohort_expression_", .x), .cols = starts_with("quant")) %>% 
  rename(fiveprime_cohort_expression_sd = sd) %>% 
  left_join(quantiles, by = c("project_id" = "cancer_type", "threeprimepartner"="gene_id")) %>% 
  rename_with(~paste0("threeprime_cohort_expression_", .x), .cols = starts_with("quant")) %>% 
  rename(threeprime_cohort_expression_sd = sd) %>% 
  select(-c(submitter_id, project_id)) 

rm(project_ids, submitter_ids)

# write file output 
fusioncatcher %>% 
  filter(cohort == "TCGA") |> 
  write_tsv("221124_tcga_output_processed.txt")

fusioncatcher |> 
  filter(cohort == "SCANB") |> 
  write_tsv("221124_scanb_tnbc_output_processed.txt")



