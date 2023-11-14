library(tidyverse)
arriba_brcA_luad_fusions <- read_tsv("arriba_brca_luad_discordant_read_pairs.txt")


# process disc reads further ----------------------------------------------

 
arriba_disc_read_clusters <- arriba_brcA_luad_fusions %>% 
  mutate(width = nchar(SEQ), 
    start = as.numeric(POS)) %>% 
  nest(.by = c(fusion_id, QNAME)) |> 
  mutate(QNAME_ID = row_number()) %>%  
  unnest(data) |> 
  select(fusion_id, QNAME_ID, query, width, start, fiveprime_strand, threeprime_strand) %>% 
  mutate(end = start + width) %>% 
  mutate(strand = ifelse(query == "fiveprime_reads", fiveprime_strand, threeprime_strand)) %>% 
  select(-c(fiveprime_strand, threeprime_strand)) %>% 
  mutate(start = ifelse( (query == "fiveprime_reads" & strand == "+") | (query == "threeprime_reads" & strand == "-") , start-50, start-500)) %>% 
  mutate(end = ifelse( (query == "fiveprime_reads" & strand == "+") | (query == "threeprime_reads" & strand == "-") , end+500, end+50)) %>% 
  group_by(fusion_id, query) %>% 
  nest() %>% 
  mutate(iranges = map(data, function(df){
    IRanges::IRanges(df$start, end = df$end) %>% IRanges::reduce() %>% IRanges::as.data.frame() %>% as_tibble()
  })) %>% 
  unnest(iranges) %>% 
  rename(istart = start, iend = end, iwidth = width) %>% 
  unnest(data) %>% 
  ungroup()  %>% 
  mutate(filt = between(start, istart, iend)) %>% 
  filter(filt) %>% 
  select(-c(start, end, width, filt)) %>% 
  rename(start = istart, end = iend, width = iwidth) %>% 
  group_by(fusion_id, query, strand, start, end, width) %>% 
  summarise(n_disc_reads_in_cluster = n(), disc_reads_in_cluster = paste0(QNAME_ID, collapse = ";")) %>% 
  ungroup()

arriba_disc_read_clusters <- arriba_disc_read_clusters |> 
  group_by(fusion_id) %>% 
  nest() %>% 
  mutate(data = map(data, function(df){
    df %>% 
      mutate(cluster_id = 1:nrow(.))
  })) %>% 
  unnest(data) %>% 
  mutate(disc_read_cluster_id = paste0(fusion_id, "_", cluster_id)) %>% 
  select(-cluster_id) %>% 
  ungroup()


arriba_disc_read_cluster_key <- arriba_disc_read_clusters %>% 
  select(fusion_id, query,disc_reads_in_cluster, disc_read_cluster_id) %>% 
  mutate(disc_reads_in_cluster = str_split(disc_reads_in_cluster, ";")) %>% 
  unnest(disc_reads_in_cluster) %>% 
  left_join(
    arriba_disc_read_clusters %>% 
      select(fusion_id, query,disc_reads_in_cluster, disc_read_cluster_id) %>% 
      mutate(disc_reads_in_cluster = str_split(disc_reads_in_cluster, ";")) %>% 
      unnest(disc_reads_in_cluster) %>% 
      rename(associated_disc_read_cluster = disc_read_cluster_id,
        associated_query = query), 
    join_by(fusion_id, disc_reads_in_cluster)
  ) %>% 
  filter(query != associated_query) %>% 
  group_by(fusion_id, query, disc_read_cluster_id) %>% 
  summarise(
    disc_reads_in_cluster = paste0(disc_reads_in_cluster, collapse = ";"),
    associated_disc_read_clusters = paste0(unique(associated_disc_read_cluster), collapse = ";")) %>% 
  ungroup()


disc_read_clusters <- disc_read_clusters %>% 
  left_join(input %>% 
      select(fusion_id, path, fiveprime_chr, threeprime_chr) %>% 
      distinct()) %>% 
  mutate(chr = ifelse(query == "fiveprime_reads", fiveprime_chr, threeprime_chr)) %>% 
  select(-c(fiveprime_chr, threeprime_chr))


# -------------------------------------------------------------------------

arriba_brcA_luad_fusions |> 
  distinct(QNAME, fusion_id) |> 
  count(fusion_id) |> 
  ggplot(aes(n)) + 
  geom_histogram() +
  scale_x_log10()

# which fusions can we validate? 

search_table <- read_tsv("arriba_brca_luad_fusion_search_table_ran_on_bianca_as_input.txt")
fusion_table_full <- read_tsv("arriba_brca_luad_fusions_search_table.txt")

# filter out fusions in fusion_table_full that were removed from search_table
# these fusions may or may not be real but we do not have the files on bianca to check. 
search_table |> select(fusion_id)

fusion_table_full <- fusion_table_full |> 
  filter(fusion_id %in% search_table$fusion_id)

fusion_table_full <- fusion_table_full |> distinct()

# some numbers:
# Fusion table full is 16,136 rows. 
# These represent 13860 unique fusion transcripts detected by arriba. 
# The fusion transcripts are in 9458 unique gene/sample combinations 

fusion_table_full |>
  distinct(fusion_id, gene_id1, gene_id2) |> 
  count(gene_id1, gene_id2, sort = T) 

fusion_table_full |> 
  distinct(fusion_id) |> 
  left_join(arriba_brcA_luad_fusions |> 
      count(fusion_id) |> 
      mutate(n = n / 2)) |> 
  mutate(support = ifelse(is.na(n), F, T)) |> 
  count(support)

# We found support for 4657 fusions, and no support for 9203 fusions
# roughly a 33% validation ratio

# for comparison, we got a 25% validation ratio in our methods paper using high confidence fusions 
# in BRCA and GBM 
# 

# set up a simple ML model to see if we can improve this at all. 
# we need to be careful with how we set up the input table - since the fusion table full has more rows than fusion IDs. 
# Standard features - FPKM, gene locations within chromosome, distance between partners, cosmic status, repeatmasker information, oncofuse label? 
# Arriba specific features - fusion strand and gene strand (might make a big difference), info from the fusion table full
# Logical retained protein coding domains?
# We can split the filters column into a wide format and use as features too
# peptide sequence length, can do length for each partner too. 

# Take out duplicated WGS paths and we have a nice start to our table.
fusion_ml_table <- fusion_table_full |> 
  select(-path) |> 
  distinct()

fusion_ml_table <- fusion_ml_table |> 
  mutate(support = fusion_id %in% arriba_brcA_luad_fusions$fusion_id) |> 
  select(fusion_id, support, everything()) |> 
  select(-value)

fusion_ml_table <- fusion_ml_table |> 
  select(-c(closest_genomic_breakpoint1, closest_genomic_breakpoint2, read_identifiers, transcript_id1, transcript_id2, pass_5, pass_3))

# we need to separate the transcript and peptide sequence into gene1 and gene2. problem is some events have homolgy which is assigned to both partners. the string would then look like ATGAT|AAT|TGTAAT. Because of this, and we dont really need the full sequence for either gene, i think it is a better idea to create the new columns with a mutate instead of separate. 
# We are probably most interested in the dinucleotide sequence at the fusion boundary. 
# all peptide homology is a maximum of 1 amino acid long, so we will just have the sequence itself as a feature

fusion_ml_table <- fusion_ml_table |> 
  separate(`strand1(gene/fusion)`, into = c("strand1_gene", "strand1_fusion"), sep = "/") |> 
  separate(`strand2(gene/fusion)`, into = c("strand2_gene", "strand2_fusion"), sep = "/") |> 
  mutate(
    has_transcript_homology = str_count(fusion_transcript, "\\|") >= 2, 
    has_peptide_homology = str_count(peptide_sequence, "\\|") >= 2,
    ) |> 
  mutate(
    peptide_sequence_gene1_length = str_extract(peptide_sequence, "^[^\\|]+") |> str_length(),
    peptide_sequence_gene1_junction = str_extract(peptide_sequence, "[^\\|](?=\\|)"),
    peptide_sequence_gene2_length = str_extract(peptide_sequence, "[^\\|]+$") |> str_length(), 
    peptide_sequence_gene2_junction = str_extract(peptide_sequence, "[^\\|]+$") |> str_extract("[^\\.]{1}"),

    peptide_homology_sequence = str_extract(peptide_sequence, "(?<=\\|)[^\\|]+(?=\\|)")
  ) |>
  mutate(
    transcript_sequence_gene1_length = str_extract(fusion_transcript, "^[^\\|]+") |> str_length(),
    transcript_sequence_gene1_junction_2nt = str_extract(fusion_transcript, "[^\\|]{2}(?=\\|)"),
    transcript_sequence_gene2_length = str_extract(fusion_transcript, "[^\\|]+$") |> str_length(),
    transcript_sequence_gene2_junction_2nt = str_extract(fusion_transcript, "[^\\|]+$") |> str_extract("[^\\.]{2}"),
    transcript_homology_sequence = str_extract(fusion_transcript, "(?<=\\|)[^\\|]+(?=\\|)"), 
    transcript_homology_length = str_length(transcript_homology_sequence)
  ) 

# now we can remove transctipt and peptide sequence, they are not useful by themselves as features

fusion_ml_table <- fusion_ml_table |> 
  select(-c(fusion_transcript, peptide_sequence)) |> 
# the filters column should be separated into multiple columns  
  mutate(filters = str_split(filters, ",")) |> 
  unnest(filters) |> 
  separate(filters, into = c("filter_name", "filter_val"), sep = "\\(") |> 
  mutate(filter_val = str_remove(filter_val, "\\)") |> as.numeric()) |>
  mutate(filter_name = paste0("filter_", filter_name)) |> 
  pivot_wider(names_from = filter_name, values_from = filter_val, values_fill = 0) |>  
  rename(has_filters = filter_.) |> 
  mutate(has_filters = ifelse(is.na(has_filters), F, T)) 

# Add features:
# Cosmic status

cosmic_genes <- read_csv("cosmic_cancer_gene_census_220805.csv") |> 
  select(Synonyms, cosmic_role = "Role in Cancer", is_cosmic_hallmark = "Hallmark") |> 
  mutate(name = str_extract(Synonyms, "ENSG\\d+")) |> 
  select(-Synonyms) |> 
  mutate(is_cosmic_hallmark = ifelse(is.na(is_cosmic_hallmark), F, T)) |> 
  mutate(cosmic_role = str_remove_all(cosmic_role, "\\s")) |> 
  mutate(cosmic_role = str_split(cosmic_role, ",")) |> 
  unnest() |> 
  mutate(role = T) |> 
  drop_na(name) |> 
  pivot_wider(id_cols = c(name, is_cosmic_hallmark), names_from = cosmic_role, values_from = role, values_fill = F, names_prefix = "cosmic_") |> 
  select(-cosmic_NA)

fusion_ml_table <- fusion_ml_table |> 
  mutate(
    fiveprimepartner = str_remove(gene_id1, "\\.\\d+$"),
    threeprimepartner = str_remove(gene_id2, "\\.\\d+$"),
    ) |> 
  left_join(cosmic_genes, join_by(fiveprimepartner == name)) |> 
  mutate(across(.cols =  c(is_cosmic_hallmark, cosmic_oncogene, cosmic_TSG, cosmic_fusion), ~replace_na(.x, F) )) |> 
  rename_with(.cols = c(is_cosmic_hallmark, cosmic_oncogene, cosmic_TSG, cosmic_fusion), ~paste0("fiveprime_", .x)) |> 
  left_join(cosmic_genes, join_by(threeprimepartner == name)) |> 
  mutate(across(.cols =  c(is_cosmic_hallmark, cosmic_oncogene, cosmic_TSG, cosmic_fusion), ~replace_na(.x, F) )) |> 
  rename_with(.cols = c(is_cosmic_hallmark, cosmic_oncogene, cosmic_TSG, cosmic_fusion), ~paste0("threeprime_", .x)) 

# mirna host status

mirna_hosts <- read_tsv("miRNA_hosts.txt") |> 
  distinct(host_ID) |> 
  as_vector()

fusion_ml_table <- fusion_ml_table |> 
  mutate(
    fiveprime_is_mirna_host = fiveprimepartner %in% mirna_hosts,
    threeprime_is_mirna_host = threeprimepartner %in% mirna_hosts
  )

# snorna host status

snorna_hosts <- read_tsv("snoRNA_host_genes.txt") |> 
  distinct(host_ID) |> 
  drop_na() |> 
  as_vector()

fusion_ml_table <- fusion_ml_table |> 
  mutate(
    fiveprime_is_snorna_host = fiveprimepartner %in% snorna_hosts,
    threeprime_is_snorna_host = threeprimepartner %in% snorna_hosts
  ) 

# RMSK info


rmsk <- read_tsv("rmsk.txt", col_names = c("bin","swScore","milliDiv","milliDel", "milliIns","genoName","genoStart","genoEnd",
  "genoLeft","strand","repName","repClass","repFamily","repStart",
  "repEnd","repLeft","id")) |> 
  filter(nchar(genoName) < 6) 

rmsk_ranges <- GenomicRanges::GRanges(seqnames = rmsk$genoName, ranges = IRanges::IRanges(start = rmsk$genoStart, end =rmsk$genoEnd))

fiveprime_junction_rmsk_overlap <- GenomicRanges::GRanges(seqnames =  fusion_ml_table$fiveprime_chr, ranges = IRanges::IRanges(start = fusion_ml_table$fiveprime_junction_hg38, width = 1)) |>  
  GenomicRanges::findOverlaps(rmsk_ranges) |>  
  as_tibble() |>  
  rename(fusion_table_row_id = queryHits, rmsk_row_id = subjectHits) |> 
  left_join(rmsk |> select(repClass) |> mutate(rmsk_row_id = row_number())) |> 
  group_by(fusion_table_row_id) |> 
  summarise(repClass = paste0(repClass, collapse = ","))

threeprime_junction_rmsk_overlap <- GenomicRanges::GRanges(seqnames =  fusion_ml_table$threeprime_chr, ranges = IRanges::IRanges(start = fusion_ml_table$threeprime_junction_hg38, width = 1)) |>  
  GenomicRanges::findOverlaps(rmsk_ranges) |>  
  as_tibble() |>  
  rename(fusion_table_row_id = queryHits, rmsk_row_id = subjectHits) |> 
  left_join(rmsk |> select(repClass) |> mutate(rmsk_row_id = row_number())) |> 
  group_by(fusion_table_row_id) |> 
  summarise(repClass = paste0(repClass, collapse = ","))

fusion_ml_table <- fusion_ml_table |> 
  mutate(fusion_table_row_id = row_number()) |> 
  left_join(fiveprime_junction_rmsk_overlap) |> 
  rename(fiveprime_rmsk_repClass = repClass) |> 
  left_join(threeprime_junction_rmsk_overlap) |> 
  rename(threeprime_rmsk_repClass = repClass) |> 
  mutate(
    fiveprime_junction_overlaps_rmsk = !is.na(fiveprime_rmsk_repClass),
    threeprime_junction_overlaps_rmsk = !is.na(threeprime_rmsk_repClass),
    )

rm(rmsk, rmsk_ranges)
# multi-junction fusion

multi_fusion_junctions <- fusion_ml_table |> 
  distinct(sample_barcode, fiveprimepartner, threeprimepartner, fiveprime_junction_hg38, threeprime_junction_hg38) |> 
  count(sample_barcode, fiveprimepartner, threeprimepartner, sort =T ) |> 
  rename(n_fusion_junctions = n)

fusion_ml_table <- fusion_ml_table |> 
  left_join(multi_fusion_junctions) 

# reciprocal fusions detected?

fusion_reciprocal_status <- fusion_ml_table %>% 
  select(fusion_id, sample_barcode, fiveprimepartner, threeprimepartner) %>% 
  mutate(key = paste0(fiveprimepartner, threeprimepartner),
    key_rev = paste0(threeprimepartner, fiveprimepartner)) %>% 
  group_by(sample_barcode) %>% 
  nest() %>% 
  mutate(data = map(data, ~mutate(.x, is_reciprocal = key%in%key_rev ))) %>% 
  unnest(data) %>% 
  select(fusion_id, is_reciprocal) |> 
  ungroup() |> 
  select(-sample_barcode)

fusion_ml_table <- fusion_ml_table |> 
  left_join(fusion_reciprocal_status)

rm(fusion_reciprocal_status)
# chromosome locations as in the FC ML

fusion_ml_table <- fusion_ml_table |> 
  mutate(is_same_chr = fiveprime_chr == threeprime_chr)

# junction distance hg38 and hg19, if on same chr

fusion_ml_table <- fusion_ml_table |> 
  mutate(
    junction_distance_hg38_log2 = ifelse(is_same_chr, log2(abs(fiveprime_junction_hg38 - threeprime_junction_hg38)), NA),
    junction_distance_hg19_log2 = ifelse(is_same_chr, log2(abs(fiveprime_junction_hg19 - threeprime_junction_hg19)), NA),
    junction_distance_hg38_hg19_difference = junction_distance_hg38_log2 - junction_distance_hg19_log2
  )

# gene length

gencode <- read_tsv("gencode.v38.annotation.gtf.gz", comment = "#", col_names = F) |> 
  filter(X3 == "gene") |> 
  mutate(
    gene_id = str_extract(X9, '(?<=gene_id ")[^"]+'),
    gene_type = str_extract(X9, '(?<=gene_type ")[^"]+'),
    gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')
  ) |> 
  select(gene_id, gene_type, gene_name, start = X4, end = X5, strand = X7) |> 
  mutate(width = end-start)

fusion_ml_table <- fusion_ml_table |> 
  left_join(gencode |> select(gene_id1 = gene_id, fiveprime_width = width)) |> 
  left_join(gencode |> select(gene_id2 = gene_id, threeprime_width = width)) 

# gene location, relative to chr length (% into chromosome)

chr_sizes <- read_tsv("chromosome_sizes_hg38") |> 
  mutate(chr = paste0("chr", chr))

fusion_ml_table <- fusion_ml_table |> 
  left_join(chr_sizes, by = c("fiveprime_chr" = "chr")) |> 
  mutate(fiveprime_rel_gene_location = round(fiveprime_junction_hg38 /  size, 2)  ) |> 
  select(-size) |>  
  left_join(chr_sizes, by = c("threeprime_chr" = "chr")) |> 
  mutate(threeprime_rel_gene_location = round(threeprime_junction_hg38 /  size, 2)  ) |> 
  select(-size) 

#fpkm

brca_expression <- read_tsv("~/projects/tcga_data/expression_matrices/TCGA_BRCA_Primary_Tumor_expression_matrix_fpkm.txt")
luad_expression <- read_tsv("~/projects/tcga_data/expression_matrices/TCGA_LUAD_Primary_Tumor_expression_matrix_fpkm.txt")

brca_expression <- brca_expression |> 
  mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |> 
  pivot_longer(cols = -gene_id, names_to = "sample", values_to = "FPKM")

luad_expression <- luad_expression |> 
  mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |> 
  pivot_longer(cols = -gene_id, names_to = "sample", values_to = "FPKM")

fusion_ml_table |>  write_tsv("arriba_ml_table.txt")
# add TCGA barcode to the expression

expression_file_key <- brca_expression |> 
  distinct(sample) |> 
  as_vector() |> 
  TCGAutils::UUIDtoBarcode("file_id") |> 
  as_tibble() |> 
  mutate(sample_barcode = str_sub(associated_entities.entity_submitter_id, 1, 15)) |> 
  select(file_id, sample_barcode) |> 
  bind_rows(
    luad_expression |> 
      distinct(sample) |> 
      as_vector() |> 
      TCGAutils::UUIDtoBarcode("file_id") |> 
      as_tibble() |> 
      mutate(sample_barcode = str_sub(associated_entities.entity_submitter_id, 1, 15)) |> 
      select(file_id, sample_barcode)
  )

brca_expression <- brca_expression |> 
  left_join(expression_file_key, join_by(sample == file_id)) |> 
  select(-sample)

brca_expression <- brca_expression |> 
  group_by(sample_barcode) |> 
  nest() |> 
  mutate(rows = map_dbl(data, ~nrow(.x))) |> 
  arrange(-rows) |> 
  mutate(data =
      ifelse(
        rows == 60660, 
        data, 
        map(
          data, 
          ~group_by(.x, gene_id) |> 
            filter(FPKM == max(FPKM)) |> 
            distinct()
        ))) |> 
  select(-rows) |> 
  unnest(data)
  
luad_expression <- luad_expression |> 
  left_join(expression_file_key, join_by(sample == file_id)) |> 
  select(-sample)

luad_expression <- luad_expression |> 
  group_by(sample_barcode) |> 
  nest() |> 
  mutate(rows = map_dbl(data, ~nrow(.x))) |> 
  arrange(-rows) |> 
  mutate(data =
      ifelse(
        rows == 60660, 
        data, 
        map(
          data, 
          ~group_by(.x, gene_id) |> 
            filter(FPKM == max(FPKM)) |> 
            distinct()
        ))) |> 
  select(-rows) |> 
  unnest(data)



all_expression <- brca_expression |> 
  mutate(project_id = "BRCA") |> 
  bind_rows(luad_expression |> 
      mutate(project_id = "LUAD"))

rm(brca_expression, luad_expression)


fusion_ml_table <- fusion_ml_table |> 
  left_join(all_expression, join_by(sample_barcode, fiveprimepartner == gene_id)) |> 
  rename(fiveprime_FPKM = FPKM) |> 
  left_join(all_expression, join_by(sample_barcode, threeprimepartner == gene_id, project_id)) |> 
  rename(threeprime_FPKM = FPKM)
  
# some samples do not have expression information. most of them are normal or have a barcode ending in 06
fusion_ml_table |>
  select(sample_barcode, fiveprimepartner, 101:103) |> 
  filter(is.na(fiveprime_FPKM)) |> 
  count(sample_barcode) |> 
  arrange(-n)


# n fusions in sample

fusion_ml_table |> 
  mutate(n_fusion_transcripts_in_sample = n(), .by = c(sample_barcode))

# finally put ID columns first along with the outcome

fusion_ml_table <- fusion_ml_table |> 
  select(support, fusion_id, fiveprimepartner, threeprimepartner, gene_id1, gene_id2, gene1, gene2, retained_protein_domains, sample_barcode, RNA_file_ID, everything())

fusion_ml_table |>  write_tsv("arriba_ml_table.txt")

# ML classifier -----------------------------------------------------------

library(tidymodels)

fusion_ml_table <- read_tsv("arriba_ml_table.txt")

fusion_ml_table <- fusion_ml_table |> 
  mutate(support = as_factor(support) |> fct_rev())

# lgbm_vars <- names(fusion_ml_table)
# lgbm_roles = c("outcome", rep("ID", 10), rep("predictor", 92))

set.seed(123)
splits <- initial_split(fusion_ml_table, strata = support)
training <- training(splits)
crossvalidation <- vfold_cv(training, strata = support, v = 10)

lgbm_recipe <- recipe(training, formula = support ~ .) |> 
  step_rm("fusion_id", "fiveprimepartner", "threeprimepartner", "gene_id1","gene_id2", "gene1", "gene2", "retained_protein_domains", "sample_barcode", "RNA_file_ID", "fusion_table_row_id") |> 
  step_dummy_extract(type, site1, site2, sep = "/" ) |> 
  step_nzv() |> 
  prep()

library(bonsai)

lgbm_model <- boost_tree(
  mode = "classification",
  mtry = tune(),
  trees = tune(),
  min_n = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  loss_reduction = tune(),
  sample_size = tune(),
  stop_iter = tune()) |>
  set_engine("lightgbm") 

# workflow

fusion_wf <- workflow() |> 
  add_recipe(lgbm_recipe) |> 
  add_model(lgbm_model)

# tuning grid
set.seed(123)
lgbm_grid <- fusion_wf |> 
  extract_parameter_set_dials() |> 
  finalize(juice(lgbm_recipe)) |> 
  grid_latin_hypercube(size = 250)

# tuning 
set.seed(123)
lgbm_tuning <- fusion_wf |> 
  finetune::tune_race_anova(
    resamples = crossvalidation,
    grid = lgbm_grid,
    control = finetune::control_race(verbose = T, verbose_elim = T),
    metrics = metric_set(brier_class, pr_auc)
  )

finetune::plot_race(lgbm_tuning) +
  theme_classic() + 
  scale_x_continuous(breaks = c(1:10))

ggsave("ML_results_tuning_race.pdf", height = 80, width = 80, units = "mm")

lgbm_tuning |> show_best("pr_auc")
lgbm_tuning |> show_best("brier_class")


# mtry trees min_n tree_depth learn_rate loss_reduction sample_size stop_iter .metric .estimator  mean     n std_err .config               
# <int> <int> <int>      <int>      <dbl>          <dbl>       <dbl>     <int> <chr>   <chr>      <dbl> <int>   <dbl> <chr>                 
# 1    98  1404    15         13     0.0294        0.00359       0.632        14 pr_auc  binary     0.824    10 0.00791 Preprocessor1_Model211
# 2    16  1899    17         14     0.0186        0.637         0.674         8 pr_auc  binary     0.816    10 0.00786 Preprocessor1_Model033
# 3     9   914    28         15     0.0366        0.211         0.658        16 pr_auc  binary     0.814    10 0.00794 Preprocessor1_Model019


# mtry trees min_n tree_depth learn_rate loss_reduction sample_size stop_iter .metric     .estimator  mean     n std_err .config               
# <int> <int> <int>      <int>      <dbl>          <dbl>       <dbl>     <int> <chr>       <chr>      <dbl> <int>   <dbl> <chr>                 
# 1    98  1404    15         13     0.0294        0.00359       0.632        14 brier_class binary     0.122    10 0.00340 Preprocessor1_Model211
# 2    16  1899    17         14     0.0186        0.637         0.674         8 brier_class binary     0.122    10 0.00272 Preprocessor1_Model033
# 3     9   914    28         15     0.0366        0.211         0.658        16 brier_class binary     0.123    10 0.00273 Preprocessor1_Model019

# model 211 performs well in pr_auc and brier_class

lgbm_tuning |>  saveRDS("arriba_fusion_ml_lgbm_tuning_object.Rds")

lgbm_best <- lgbm_tuning |> select_best("brier_class")


# final model

final_model <- fusion_wf |> 
  finalize_workflow(lgbm_best)

final_model |> write_rds("arriba_fusion_ml_lgbm_final_model.Rds")

fusion_last_fit <- final_model |> 
  last_fit(splits)

# Metrics: 
fusion_last_fit |> 
  collect_predictions() |> 
  calc_all_metrics(support, .pred_class, .pred_TRUE)

# Accuracy 0.836
# Kappa 0.621
# Multinomal log loss 0.388
# roc auc 0.898
fusion_last_fit |> 
  collect_predictions() |> 
  conf_mat(support, estimate = .pred_class)

fusion_last_fit |> 
  collect_predictions() |>
  mutate(.pred_class = as_factor(.pred_TRUE > 0.23) |> fct_rev() ) |> 
  metrics(support, estimate = .pred_class, .pred_TRUE)

# pr auc 0.828
fusion_last_fit |> 
  collect_predictions() |> 
  pr_auc(support, .pred_TRUE)

# precision 0.794
fusion_last_fit |> 
  collect_predictions() |> 
  precision(support, .pred_class)

# recall/sensitivity 0.692
fusion_last_fit |> 
  collect_predictions() |> 
  recall(support, .pred_class)

# brier class 0.12
fusion_last_fit |> 
  collect_predictions() |> 
  brier_class(support, .pred_TRUE)

# F1 score 0.739
fusion_last_fit |> 
  collect_predictions() |> 
  f_meas(support, .pred_class)

# specificity 0.909
fusion_last_fit |> 
  collect_predictions() |> 
  specificity(support, .pred_class)

# Youden's J
fusion_last_fit |> 
  collect_predictions() |> 
  roc_curve(support, .pred_TRUE) |> 
  mutate(J = sensitivity + specificity - 1) |> 
  arrange(-J) |> 
  ggplot(aes(x = .threshold, y = J)) + 
  geom_line() + 
  theme_classic() + 
  geom_vline(xintercept = 0.23, linetype = 2) + 
  annotate("text", x = .4, 0.25,label = "max=0.23") + 
  labs(x = "Threshold", y = "Youden's J")


# Classification probability score vs truth
fusion_last_fit |> 
  collect_predictions() |> 
  mutate(pred = as_factor(.pred_TRUE > 0.23) |> fct_rev()) |>
  mutate(truth = ifelse(support == TRUE, "WGS support", "No WGS support") |> fct_rev()) |> 
  ggplot(aes(.pred_TRUE, fill = pred)) +
  geom_histogram() +
  facet_wrap(~truth, ncol = 1, scales = "free_y") + 
  labs(fill = "Predicted class", x = "Predicted probability", y = "Count") + 
  theme_classic() +
  scale_fill_grey() +
  theme(legend.position = "top")

ggsave("ML_results_predicted_probability_vs_truth.pdf", width = 80, height = 80, units = "mm" )

# PR curve
fusion_last_fit |> 
  collect_predictions() |> 
  pr_curve(support, .pred_TRUE, ) |> 
  bind_rows(tibble(.threshold = 0, recall = 1, precision = 0)) |> 
  ggplot(aes(x = recall, y = precision)) + 
  geom_line() + 
  coord_cartesian(xlim = c(0,1)) + 
  # coord_equal() + 
  theme_classic() + 
  labs(x = "Recall", y = "Precision") +
  annotate(geom = "text", label = "AUC=0.828", x = 0.35, y = 0.65)

ggsave("ML_results_pr_curve.pdf", width = 80, height = 80, units = "mm" )


# ROC curve 
fusion_last_fit |> 
  collect_predictions() |> 
  roc_curve(support, .pred_TRUE, ) |> 
  ggplot(aes(x = 1-specificity, y = sensitivity)) + 
  geom_line() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  theme_classic()+
  labs(x = "1 - Specificity", y = "Sensitivity") + 
  annotate(geom = "text", label = "AUC=0.899", x = 0.35, y = 0.65)

ggsave("ML_results_roc_curve.pdf", width = 80, height = 80, units = "mm" )

  
# Feature importance
feature_importance <- fusion_last_fit |> 
  extract_fit_parsnip() |> 
  vip::vi()

fusion_last_fit |> 
  extract_fit_parsnip() |> 
  vip::vip(num_features = 20) + 
  theme_classic() 

ggsave("ML_results_feature_importance.pdf", width = 100, height = 80, units = "mm" )

# confidence matrix

# Truth
# Prediction TRUE FALSE
# TRUE   980   472
# FALSE  185  1829

fusion_last_fit |> 
  collect_predictions() |> 
  mutate(pred = as_factor(.pred_TRUE > 0.23) |> fct_rev()) |>
  conf_mat(support, pred) |> 
  autoplot() 

ggsave("ML_results_confidence_matrix.pdf", width = 80, height = 80, units = "mm" )


# How does the ML compare to "normal" filtering of the fusions?
# If we filter them to retain only high / medium confidence, contain only X or more discordant reads, filter out common mapping reads etc? 


fusion_ml_table |> 
  mutate(test_filter = confidence %in% c("high", "medium") |> as_factor() |> fct_rev()) |># |> as.numeric()) |> #  
  # count(test_filter, support)
  f_meas(support, test_filter)


fusion_ml_table |> 
  mutate(test_filter = as_factor(discordant_mates > 2 & confidence %in% c("high", "medium")) |> fct_rev()) |>
  # mutate(test_filter = as_factor(discordant_mates > 2) |> fct_rev()) |>
  f_meas(support, test_filter)

fusion_ml_table |> 
  mutate(test_filter = as.numeric(discordant_mates > 2 & confidence %in% c("high", "medium"))) |>
  # mutate(test_filter = as.numeric(discordant_mates > 2)) |>
  pr_auc(support, test_filter)


fusion_ml_table |> 
  mutate(test_filter = as_factor(reading_frame == "in-frame") |> fct_rev()) |> 
  select(support, test_filter, reading_frame) |>
  f_meas(support, test_filter)


# Can we simulate something similar to a roc curve if we filter over a range of discordant read pair values?

fusion_ml_table |> 
  select(support, discordant_mates) |> 
  mutate(support = as.logical(support)) |> 
  group_by(discordant_mates) |> 
  nest() |> 
  arrange(discordant_mates) |> 
  mutate(n_pos = map_int(data, ~pull(.x, support) |> which() |> length())) |> 
  mutate(test = map2_dbl(n_pos, data, ~.x/nrow(.y)) ) |> 
  arrange(-discordant_mates) |> 
  ggplot(aes(x = discordant_mates, y = test, size = n_pos)) + 
  geom_point()


# Evaluate model on SCAN-B TNBC arriba fusions -----------------------------------------------

# ML table of Arriba derived SCAN-B fusions was created in process_arriba_scanb_tnbc_fusions.R script

arriba_scanb_ml_table <- read_tsv("~/projects/wgs_fusions/reviewer_comments_code/arriba_fusions/arriba_scanb_ml_table.txt")

arriba_scanb_ml_table <- arriba_scanb_ml_table |> 
  mutate(support = as_factor(support) |> fct_rev())

scanb_predictions <- predict(extract_workflow(fusion_last_fit_test), arriba_scanb_ml_table, type = "prob")
scanb_predictions

arriba_scanb_ml_table <- arriba_scanb_ml_table |> 
  bind_cols(scanb_predictions |> mutate(pred = .pred_TRUE > 0.5)) |> 
  select(fusion_id, support, pred, .pred_TRUE, .pred_FALSE, everything())

arriba_scanb_ml_table <- arriba_scanb_ml_table |>
  mutate(pred = as_factor(pred) |> fct_rev()) 

calc_all_metrics <- function(df, truth, pred_class, estimate){
  df |> precision({{truth}}, {{pred_class}}) |> 
    bind_rows(df |> recall({{truth}}, {{pred_class}})) |> 
    bind_rows(df |> specificity({{truth}}, {{pred_class}})) |> 
    bind_rows(df |> accuracy({{truth}}, {{pred_class}})) |> 
    bind_rows(df |> f_meas({{truth}}, {{pred_class}})) |> 
    bind_rows(df |> kap({{truth}}, {{pred_class}})) |> 
    bind_rows(df |> mn_log_loss({{truth}}, {{estimate}})) |> 
    bind_rows(df |> roc_auc({{truth}}, {{estimate}})) |> 
    bind_rows(df |> pr_auc({{truth}}, {{estimate}})) |> 
    bind_rows(df |> brier_class({{truth}}, {{estimate}})) 
}

arriba_scanb_ml_table |> 
  mutate(pred = as_factor(.pred_TRUE > 0.23) |> fct_rev()) |>
  calc_all_metrics(support, pred, .pred_TRUE)

arriba_scanb_ml_table |> 
  mutate(pred = as_factor(.pred_TRUE > 0.23) |> fct_rev()) |>
  conf_mat(support, pred)

arriba_scanb_ml_table |> 
  write_tsv("~/projects/wgs_fusions/reviewer_comments_code/arriba_fusions/arriba_scanb_ml_table_w_predictions.txt")

fusion_last_fit |> 
  collect_predictions() |> 
  write_tsv("~/projects/wgs_fusions/reviewer_comments_code/arriba_fusions/arriba_tcga_fusions_testing_data_predictions.txt")


# scan-b filtering metrics ------------------------------------------------


arriba_scanb_ml_table |> 
  mutate(test_filter = as_factor(confidence == "high") |> fct_rev()) |> conf_mat(support, test_filter)
  calc_all_metrics(support, test_filter, .pred_TRUE)

arriba_scanb_ml_table |> 
  mutate(test_filter = as_factor(confidence %in% c("high", "medium")) |> fct_rev()) |> 
  calc_all_metrics(support, test_filter, .pred_TRUE)

arriba_scanb_ml_table |> 
  mutate(test_filter = as_factor(discordant_mates > 2) |> fct_rev()) |> 
  calc_all_metrics(support, test_filter, .pred_TRUE)

arriba_scanb_ml_table |> 
  mutate(test_filter = as_factor(discordant_mates > 2 & confidence %in% c("high", "medium") ) |> fct_rev()) |> 
  calc_all_metrics(support, test_filter, .pred_TRUE)

arriba_scanb_ml_table |> 
  mutate(test_filter = as_factor(reading_frame == "in-frame") |> fct_rev()) |> 
  calc_all_metrics(support, test_filter, .pred_TRUE)


# Kinase fusions  ---------------------------------------------------------

kinases <- read_tsv("~/projects/wgs_fusions/human_kinases.txt")
kinases

kinase_vec <- kinases |> pull(`HGNC Name`)

arriba_scanb_ml_table |> 
  mutate(kinase_fusion = gene1 %in% kinase_vec | gene2 %in% kinase_vec) |> 
  filter(kinase_fusion & support == F) |> View()
  # count(kinase_fusion, support)

arriba_scanb_ml_table |> 
  mutate(kinase_fusion = gene1 %in% kinase_vec | gene2 %in% kinase_vec) |> 
  filter(kinase_fusion) |> mutate(pred = as_factor(.pred_TRUE > 0.23) |> fct_rev()) |> calc_all_metrics(support, pred, .pred_TRUE)

arriba_scanb_ml_table |> 
  mutate(kinase_fusion = gene1 %in% kinase_vec | gene2 %in% kinase_vec) |> 
  filter(kinase_fusion) |>
  mutate(pred = as_factor(confidence %in% c("high", "medium")) |> fct_rev()) |>
  calc_all_metrics(support, pred, .pred_TRUE)
