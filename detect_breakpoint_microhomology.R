library(tidyverse)


# functions ---------------------------------------------------------------
complement <- function(dna){
  # Create a character vector of complements
  complements <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # Use the `chartr` function to replace each base with its complement
  dna_compl <- chartr(paste0(names(complements), collapse = ""), paste0(complements, collapse = ""), toupper(dna))
  return(dna_compl)
}

reverse_complement <- function(dna) {
  
  # Reverse the complemented DNA string and return it
  return(stringi::stri_reverse(complement(dna)))
}

detect_microhomology <- function(seq1, seq2, start_pos1, start_pos2, seq2_rev_strand = F, softclip_len, expected) {
  microhomology <- ""
  
  if(expected == "right" & seq2_rev_strand) {
    start_pos2 = start_pos2 + as.numeric(softclip_len)
  }
  
  seq1_short = paste0(str_sub(tolower(seq1), start_pos1 - 6, start_pos1-1), "*", str_sub(tolower(seq1), start_pos1, start_pos1 + 5))
  seq2_short = paste0(str_sub(tolower(seq2), start_pos2 - 6, start_pos2-1), "*", str_sub(tolower(seq2), start_pos2, start_pos2 + 5))
  
  # if flag&&16 & side == left we reverse complement
  # if !flag&&16 & side == left we reverse
  if (seq2_rev_strand & expected == "left") {
    seq2_short = reverse_complement(seq2_short) %>% tolower()
  } else if (!seq2_rev_strand & expected == "left") {
    seq2_short = stringi::stri_reverse(seq2_short) %>% tolower()
  } else if (seq2_rev_strand & expected == "right") {
    seq2_short = reverse_complement(seq2_short) %>% tolower()
  }
  
  for (i in 6:1) {
    if (substr(seq1_short, i, i) == substr(seq2_short, i, i)) {
      microhomology <- paste0(microhomology, substr(seq1_short, i, i))
    } else {
      break
    }
  }
  microhomology <- stringi::stri_reverse(microhomology)
  microhomology <- paste0(microhomology, "*")
  
# microhomology should only occur on the right
  # for (i in 8:13) {
  #   if (substr(seq1_short, i, i) == substr(seq2_short, i, i)) {
  #     microhomology <- paste0(microhomology, substr(seq1_short, i, i))
  #   } else {
  #     break
  #   }
  # }
  is_homology = microhomology != "*"
  
  return(tibble(seq1_break = seq1_short, seq2_break = seq2_short, sequence_microhomology = microhomology))
}

extract_random_char <- function(string) {
  index <- sample(1:nchar(string), 1)
  return(substr(string, index, index))
}

extract_random_substring <- function(string, n) {
  if(n > nchar(string)) stop("substring length should be less than or equal to the string length")
  indexes <- sample(1:(nchar(string) - n + 1), 1)
  return(substr(string, indexes, indexes + n - 1))
}

compare_homology <- function(str1, str2){
  homology = ""
  str1 = tolower(str1)
  str2 = tolower(str2)
  for (i in min(nchar(str1), nchar(str2)):1)  {
    
    if (substr(str1, i, i) == substr(str2, i, i)) {
      homology <- paste0(homology, substr(str1, i, i))
    } else {
      break
    }
    
  }
  return(homology)
}

# main  -------------------------------------------------------------------
sam_fields <- c("QNAME", "FLAG","RNAME",	"POS","MAPQ","CIGAR","MRNM","MPOS","ISIZE","SEQ",	"QUAL",	"OPT")
cancer_types <- c("gbm", "luad", "cesc_esca", "lgg")
early_cancer_types <- c("brca", "blca", "ov", "kirp", "kirc", "kich")

softclips <- cancer_types %>% 
  as_tibble() %>% 
  rename(cancer_type = value) %>% 
  mutate(softclips = map(cancer_type, ~read_tsv(paste0(.x, "/output/HQ_softclips.txt"), col_types = "cdcddccddcccccccccccc"))) %>% 
  unnest(softclips) %>% 
  bind_rows( 
    early_cancer_types %>% 
      as_tibble() %>% 
      rename(cancer_type = value) %>% 
      mutate(softclips = map(cancer_type, ~read_tsv(paste0(.x, "/", toupper(.x), "_HQ_softclips.txt"), col_types = "cdcddccddcccccccccccc"))) %>% 
      unnest(softclips)
  )


alignment <- cancer_types %>% 
  as_tibble() %>% 
  rename(cancer_type = value) %>% 
  mutate(alignment = map(cancer_type, ~read_tsv(paste0(.x, "/output/softclip_alignment.txt"), comment = "@", col_names = sam_fields))) %>% 
  unnest(alignment) %>% 
  filter(RNAME != "*") %>% 
  bind_rows(
    early_cancer_types %>% 
      as_tibble() %>% 
      rename(cancer_type = value) %>% 
      mutate(alignment = map(cancer_type, ~read_tsv(paste0(.x, "/", toupper(.x), "_softclip_alignment.txt"), comment = "@", col_names = sam_fields))) %>% 
      unnest(alignment) %>% 
      filter(RNAME != "*")
  )


cluster_index_key <-  cancer_types %>% 
  as_tibble() %>% 
  rename(cancer_type = value) %>% 
  mutate(key = map(cancer_type, ~read_tsv(paste0(.x, "/output/index_region_key.txt"), col_types = "ccddc"))) %>% 
  unnest(key) %>% 
  bind_rows(
    early_cancer_types %>% 
      as_tibble() %>% 
      rename(cancer_type = value) %>% 
      mutate(key = map(cancer_type, ~read_tsv(paste0(.x, "/", toupper(.x), "_index_region_key.txt"), col_types = "ccddc"))) %>% 
      unnest(key) 
  )


cluster_alignment_key <-cancer_types %>% 
  as_tibble() %>% 
  rename(cancer_type = value) %>% 
  mutate(key = map(cancer_type, ~read_tsv(paste0(.x, "/output/cluster_alignment_key.txt")))) %>% 
  unnest(key) %>% 
  bind_rows(
    early_cancer_types %>% 
      as_tibble() %>% 
      rename(cancer_type = value) %>% 
      mutate(key = map(cancer_type, ~read_tsv(paste0(.x, "/", toupper(.x),  "_cluster_alignment_key.txt")))) %>% 
      unnest(key)
  )

index_seqs <- list.files(paste0(cancer_types,  "/output/index_seqs"), full.names = T) %>% 
  as_tibble() %>% 
  bind_rows(
    list.files(paste0(early_cancer_types, "/output/index_seqs"), full.names = T) %>% 
      as_tibble()
  ) %>% 
  filter(str_detect(value, ".fa$")) %>% 
  mutate(id = str_extract(value, "index_reg.*$")) %>% 
  mutate(file = map_chr(value, ~read_tsv(.x, skip = 1, col_names = "seq", show_col_types = F) %>% as_vector() %>%  paste0(collapse = ""))) %>%
  unnest(file)

index_seqs <- index_seqs %>% 
  mutate(id = str_remove(id, ".fa$")) %>% 
  mutate(cancer_type = str_extract(value, "\\w+")) %>% 
  select(cancer_type, id, file)

homology <- alignment %>% 
  # mutate(across(c(QNAME, RNAME), ~str_remove(.x, "^[:alpha:]+_")))%>%
  mutate(seq_id = str_extract(QNAME, "^[:alpha:]*_?\\d+_\\d+_\\d+"),
    cluster = str_extract(QNAME, "^[:alpha:]*_?\\d+_\\d+")) %>% 
  left_join(cluster_index_key %>% select(cancer_type, cluster, index_region_id), 
    by = c("cluster", "cancer_type")) %>% 
  select(cancer_type, seq_id, cluster, index_region_id, RNAME, aligned_FLAG = FLAG, aligned_POS = POS, aligned_cigar = CIGAR, aligned_seq= SEQ) %>% 
  left_join(softclips %>% select(cancer_type, seq_id, FLAG, POS, CIGAR, SEQ, expected, primary_alignment_partner = fusion_gene_partner), by = c("cancer_type", "seq_id")) %>% 
  left_join(
    index_seqs %>%
      mutate(id =  ifelse(
        cancer_type %in% early_cancer_types, 
        paste0(toupper(cancer_type), "_", id),
        id
        )),
    by = c("index_region_id" = "id", "cancer_type")) %>% 
  rename(index_reg1_seq = file) %>% 
  left_join(index_seqs %>%
      mutate(id =  ifelse(
        cancer_type %in% early_cancer_types, 
        paste0(toupper(cancer_type), "_", id),
        id
      )),
    by = c("RNAME"= "id", "cancer_type")) %>% 
  rename(index_reg2_seq = file) %>% 
  mutate(tmp = ifelse(expected == "right", str_extract(CIGAR, "\\d+(?=S$)"), str_extract(CIGAR, "^\\d+"))) %>% 
  mutate(softclip = ifelse(
    expected == "right",
    str_sub(SEQ, nchar(SEQ) - as.numeric(tmp) + 1,  nchar(SEQ)),
    str_sub(SEQ, 1, tmp))) %>% 
  mutate(M_start = ifelse(
    str_detect(CIGAR, "^\\d+S"),
    as.numeric(str_extract(CIGAR, "^\\d+(?=S)")) + 1, 
    1),
    M_end = ifelse(
      str_detect(CIGAR, "\\d+S$"),
      nchar(SEQ) - as.numeric(str_extract(CIGAR, "\\d+(?=S$)")),
      nchar(SEQ) 
    )) %>% 
  mutate(primary_m_seq = str_sub(SEQ, M_start, M_end)) 

homology <- homology %>% 
  mutate(seq_to_align = ifelse(
    expected == "right" & str_count(CIGAR, "M") > 1,
    str_sub(primary_m_seq, nchar(primary_m_seq) - as.numeric(str_extract(CIGAR, "\\d+(?=M\\d+S$)")) + 1 ,nchar(primary_m_seq)), 
    str_sub(primary_m_seq, 1, as.numeric(str_extract(CIGAR, "\\d+(?=M)")))
  )) %>% 
  mutate(locate = map2(index_reg1_seq, seq_to_align, ~str_locate(tolower(.x), tolower(.y)))) %>% 
  rowwise() %>% 
  mutate(breakpoint_ref1 = ifelse(
    expected == "right", 
    locate[[2]] + 1,
    locate[[1]]
  ))

homology <- homology %>% 
  drop_na(breakpoint_ref1) %>% 
  mutate(aligned_rev_strand = bitwAnd(aligned_FLAG, 16) != 0) %>% 
  mutate(test = pmap(list(index_reg1_seq, index_reg2_seq, breakpoint_ref1, aligned_POS, aligned_rev_strand, tmp, expected),
    ~detect_microhomology(seq1 = ..1, seq2 = ..2, start_pos1 = ..3, start_pos2 = ..4, seq2_rev_strand = ..5, softclip_len = ..6, expected = ..7 ))) %>%
  select(test, everything()) %>% 
  unnest(test)

# add random breaks as a comparison 
set.seed(1234)
homology <- homology %>%
  mutate(
    rand_first = map_chr(index_reg1_seq, ~extract_random_substring(.x, 5)),
    rand_second = map_chr(index_reg2_seq, ~extract_random_substring(.x, 5)),
    rand_homology = map2_chr(rand_first, rand_second, ~compare_homology(.x, .y)))

cesc <- read_tsv("cesc_esca/output/output_master_table.txt") %>% 
  select(fusion_id, project_id_short) %>% 
  filter(project_id_short == "CESC") %>% 
  pull(fusion_id)

esca <- read_tsv("cesc_esca/output/output_master_table.txt") %>% 
  select(fusion_id, project_id_short) %>% 
  filter(project_id_short == "ESCA") %>% 
  pull(fusion_id)

homology <- homology %>% 
  mutate(fusion_id = str_extract(seq_id, "\\d+")) %>% 
  mutate(cancer_type = case_when(
    fusion_id %in% cesc ~ "cesc", 
    fusion_id %in% esca ~ "esca", 
    T ~ cancer_type)) %>% 
  filter(cancer_type != "cesc_esca")

homology <- homology %>% 
  mutate(cancer_type = ifelse(cancer_type %in% c("kirc", "kich", "kirp"), "kidney", cancer_type))


# plot observed microhomology vs random simulated breaks ------------------
full_cancer_names <- homology %>% 
  distinct(cancer_type) %>% 
  mutate(full_cancer_name = c("Glioblastoma", "Lung", "Cervix", "Esophagus", "Lower grade glioma", "Breast", "Bladder", "Ovary", "Kidney"))

homology %>% 
  mutate(left = str_detect(sequence_microhomology, "[^\\*]+\\*")) %>% 
  filter(left) %>% 
  count(cancer_type, homology = nchar(sequence_microhomology) - 1) %>% 
  mutate(category = case_when(
    homology == 0 ~ "No\nhomology",
    homology > 0 & homology < 6 ~ "Observed\nmicrohomology",
    homology > 5 ~ "6+ nt\nhomology",
    # T ~ paste0(as.character(homology), " nt")
  )) %>% 
  group_by(cancer_type,category) %>% 
  summarise(n = sum(n)) %>% 
  bind_rows(homology %>% 
      count(cancer_type, category = ifelse(nchar(rand_homology) > 0, "Random\nmicrohomology", "rand_no_homology")) %>% 
      filter(category != "rand_no_homology")
  ) %>% 
  mutate(category = fct_relevel(category, c("Random\nmicrohomology", "Observed\nmicrohomology", "6+ nt\nhomology" ))) %>% 
  left_join(full_cancer_names) %>% 
  ggplot(aes(category, n)) +
  geom_col() +
  facet_wrap(~full_cancer_name, scales = "free_y") + 
  theme_classic() + 
  scale_y_continuous(labels = ~format(.x, scientific = F)) +
  labs(x=NULL, y = "Breakpoint events")

ggsave("../breakpoint_analysis/230119_microhomology_vs_random_breaks.pdf", device = "pdf",
  width = 260, 
  height = 140,
  units = "mm")

homology %>% 
  count(cancer_type, homology = nchar(sequence_microhomology) - 1) %>%
  mutate(category = case_when(
    homology == 0 ~ "No homology",
    homology > 0 & homology < 6 ~ "Microhomology",
    homology > 5 ~ "Larger homology\n(6+nt)",
    # T ~ paste0(as.character(homology), " nt")
  )) %>% 
  group_by(cancer_type,category) %>% 
  summarise(n = sum(n)) %>% 
  ggplot(aes(fct_rev(category), n)) + 
  geom_col() +
  facet_wrap(~cancer_type, scales = "free_y") +
  theme_classic() + 
  labs(y = "Breakpoint events", x = NULL, title = NULL)

# breakdown of microhomology lengths 


homology %>% 
  count(cancer_type, homology = nchar(sequence_microhomology) - 1) %>%
  mutate(cat = "Observed") %>% 
  bind_rows(homology %>% 
      count(cancer_type, homology = nchar(rand_homology)) %>% 
      mutate(cat = "Simulated")) %>% 
  mutate(category = case_when(
    homology == 0 ~ "No\nhomology",
    # homology > 0 & homology < 6 ~ "Microhomology",
    homology > 5 ~ "6+ nt",
    T ~ paste0(as.character(homology), " nt")
  )) %>% 
  # group_by(cancer_type, category) %>% 
  # summarise(n = sum(n)) %>% 
  filter(category != "No\nhomology") %>% 
  ggplot(aes(category, n, fill = cat)) + 
  geom_col(position = "dodge") +
  facet_wrap(~cancer_type, scales = "free_y") +
  theme_classic() + 
  labs(y = "Breakpoint events", x = NULL, fill = NULL) +
  scale_fill_grey()



ggsave("../breakpoint_analysis/230119_microhomology_length.pdf",
  device = "pdf", width = 250, height = 140, units = "mm")


# stacked bar plot 

homology %>% 
  count(cancer_type, homology = nchar(sequence_microhomology) - 1) %>%
  mutate(cat = "Observed") %>% 
  bind_rows(homology %>% 
      count(cancer_type, homology = nchar(rand_homology)) %>% 
      mutate(cat = "Simulated")) %>% 
  mutate(category = case_when(
    homology == 0 ~ "No\nhomology",
    # homology > 0 & homology < 6 ~ "Microhomology",
    homology > 5 ~ "6+ nt",
    T ~ paste0(as.character(homology), " nt")
  )) %>% 
  # group_by(cancer_type, category) %>% 
  # summarise(n = sum(n)) %>% 
  filter(category != "No\nhomology") %>% 
  drop_na() %>% 
  left_join(full_cancer_names) %>% 
  ggplot(aes(cat, n, fill = fct_rev(category))) + 
  geom_col(position = "stack") +
  facet_wrap(~full_cancer_name, scales = "free_y") +
  theme_classic() + 
  labs(y = "Breakpoint events", x = NULL, fill = "Level of\nhomology") +
  scale_fill_grey()

ggsave(filename = "../breakpoint_analysis/230119_microhomology_length_stacked.pdf",
  device = "pdf",
  width = 200,
  height = 140,
  units = "mm")

# using distinct fusion breakpoints 
# we randomly select a breakpoint and simulated breakpoint for each fusion event

tmp <- validation %>% 
  filter(support_breakpoint) %>% 
  distinct(fusion_id, cancer_type) %>% 
  mutate(hom = 0, rand_hom = 0) %>% 
  filter(!fusion_id %in% homology$fusion_id) %>% 
  left_join(fusioncatcher %>% select(cancer_type, project_id, fusion_id)) %>% 
  select(-cancer_type) %>% 
  left_join(cancer_names_full %>% distinct(project_id, cancer)) %>% 
  select(-project_id) 

homology_plot <- homology %>% 
  mutate(hom = nchar(sequence_microhomology)-1, rand_hom = nchar(rand_homology) ) %>% 
  mutate(rand_hom = replace_na(rand_hom, 0)) %>%
  distinct(cancer_type, fusion_id, .keep_all = T) %>% 
  select(cancer_type, fusion_id, hom, rand_hom) %>% 
  mutate(cancer_type = toupper(cancer_type)) %>% 
  left_join(cancer_names_full %>% distinct(project_id, cancer), by = c("cancer_type" = "project_id")) %>% 
  mutate(cancer = ifelse(cancer_type == "KIDNEY", "Kidney", cancer)) %>% 
  select(-cancer_type) %>% 
  bind_rows(tmp) %>%  # some fusion events with breakpoints get removed from the process. they will need to be added manually in and given values of 0 for both homologies. 
  pivot_longer(cols = c(hom, rand_hom), names_to = "cat", values_to = "n") %>% 
  group_by(cancer, cat, n) %>% 
  tally() %>% 
  mutate(n = as_factor(n)) %>% 
  # filter(n!=0) %>%
  mutate(cat = ifelse(cat == "hom", "Observed", "Simulated"),
    n = ifelse(n == 6, "6+ nt", paste0(n, " nt"))) %>% 
  drop_na(cancer)

homology_plot %>% 
  filter(n!="0 nt") %>%
  ggplot(aes(cat, nn, fill = fct_rev(n))) + 
  geom_col(position = "stack") +
  scale_fill_grey() +
  theme_classic() +
  labs(x = NULL, y = "Breakpoint events", fill = "Level of\nhomology") + 
  facet_wrap(~cancer, scales = "free_y")

ggsave(filename = "~/projects/wgs_fusions/fusion_plots/breakpoint_analysis/230207_microhomology_length_stacked_all_cancers.pdf",
  device = "pdf",
  width = 170,
  height = 140,
  units = "mm")

homology_plot %>% 
  filter(n!="0 nt") %>%
  group_by(cat, n) %>% 
  summarise(nn = sum(nn)) %>% 
  ggplot(aes(cat, nn, fill = fct_rev(n))) + 
  geom_col(position = "stack") +
  scale_fill_grey() +
  theme_classic() +
  labs(x = NULL, y = "Breakpoint events", fill = "Level of\nhomology")

ggsave(filename = "~/projects/wgs_fusions/fusion_plots/breakpoint_analysis/230207_microhomology_length_stacked_all_cancers_combined.pdf",
  device = "pdf",
  width = 85,
  height = 55,
  units = "mm")

# statistical test of detected microhomologies vs random breakpoints 

# Mann-Whitney U test, accounts for length of microhomologies
homology %>% 
  select(cancer_type, sequence_microhomology, rand_homology) %>% 
  mutate(n_homology = nchar(sequence_microhomology) - 1,
    n_rand_homology = nchar(rand_homology)) %>% 
  group_by(cancer_type) %>% 
  nest() %>% 
  mutate(test = map_dfc(data, ~wilcox.test(.x$n_homology, .x$n_rand_homology) %>% broom::tidy())) %>% 
  select(-data) %>% 
  unnest(test) %>% 
  select(-statistic) %>% 
  write_tsv("../breakpoint_analysis/230119_microhomology_mann_whitney_u_test_vs_random_breaks_table.txt")

# McNemar test, paired chi-square test. Length of microhomology not accounted for, just presence/absence
homology %>% 
  select(cancer_type, sequence_microhomology, rand_homology) %>% 
  mutate(a = nchar(sequence_microhomology)-1>0,
    b = nchar(rand_homology)>0) %>% 
  group_by(cancer_type) %>% 
  nest() %>% 
  mutate(test = map(data, ~mcnemar.test(.x$a, .x$b) %>% broom::tidy())) %>% 
  unnest(test)

homology %>% 
  write_tsv("../breakpoint_analysis/230119_fusion_microhomology_data.txt")

# breakpoint analysis: location of breakpoints ----------------------------------------------------
homology <- read_tsv("~/projects/wgs_fusions/breakpoint_analysis/230119_fusion_microhomology_data.txt")
gencode_full <- read_tsv("~/projects/wgs_fusions/gencode/gencode.v27lift37.annotation.gff3", comment = "#", col_names = F)

# are the detected breakpoints located mostly in introns, exons, UTR? make a plot 
# compare with gencode annotation? 

cancer_types <- c("gbm", "luad", "cesc_esca", "lgg")
early_cancer_types <- c("brca", "blca", "ov", "kirp", "kirc", "kich")

# get chr coordinates of fusion events

clusters <- cancer_types %>% 
  as_tibble() %>% 
  mutate(path = paste0("~/projects/wgs_fusions/FusionCatcher_wgs_validation/", value, "/output/wgs_spanning_read_clusters.txt")) %>% 
  rbind(
    early_cancer_types %>% 
      as_tibble() %>% 
      mutate(path = paste0("~/projects/wgs_fusions/FusionCatcher_wgs_validation/", value, "/", toupper(value), "_wgs_spanning_read_clusters.txt"))
  ) %>% 
  mutate(df = map(path, ~read_tsv(.x, show_col_types = F, col_types = "dcccddcdcc"))) %>% 
  rename(cancer_type = value)

clusters <- clusters %>% 
  unnest(df)

homology
breakpoint_coords <- homology %>% 
  # filter(cancer_type=="blca") %>% 
  mutate(breakpoint = breakpoint_ref1 + POS) %>%
  select(cancer_type, primary_alignment_partner, fusion_id, breakpoint) %>%
  left_join(clusters %>% distinct(cancer_type, fusion_id, fusion_gene_read, gene, chr),
    by = c("cancer_type", "fusion_id", "primary_alignment_partner"="fusion_gene_read")
    ) %>% 
  mutate(chr = paste0("chr", chr)) %>% 
  drop_na(gene)

gencode_full <- gencode_full %>% 
  filter(X3 != "transcript") %>% 
  mutate(gene = str_extract(X9, "ENSG\\d{11}")) %>% 
  select(gene, chr = X1, start = 4, end = 5, type = 3)

breakpoint_coords <- breakpoint_coords %>% 
  distinct()

test <- GenomicRanges::findOverlaps(
  query = GenomicRanges::GRanges(seqnames = breakpoint_coords$gene, ranges = IRanges::IRanges(start = breakpoint_coords$breakpoint, width = 1)),
  subject = GenomicRanges::GRanges(seqnames = gencode_full$gene, ranges = IRanges::IRanges(start = gencode_full$start, end = gencode_full$end))
) %>% as_tibble() %>% 
  left_join(gencode_full %>%
      mutate(row = row_number()) %>% 
      select(type, row),
    by = c("subjectHits" = "row"))

test %>% 
  group_by(queryHits) %>% 
  mutate(intron_overlap = all(type == "gene")) %>% 
  mutate(overlap = ifelse(intron_overlap, "intron", type)) %>% 
  ungroup() %>% 
  count(overlap)

breakpoint_coords %>% 
  mutate(queryHits = row_number()) %>% 
  left_join(test) %>% 
  group_by(queryHits) %>% 
  mutate(intron_overlap = all(type == "gene")) %>% 
  mutate(overlap = ifelse(intron_overlap, "intron", type)) %>% 
  ungroup() %>% 
  count(cancer_type, overlap) %>% 
  ggplot(aes(x = overlap, y = n)) + 
  geom_col() + 
  theme_classic() + 
  coord_flip()+
  facet_wrap(~cancer_type, scales = "free_x")

# can we comepare this to the average lengths of each gene type? 
gencode_full %>% 
  mutate(width = end - start + 1) %>%
  group_by(type) %>%
  summarise(mean_width = mean(width)) %>% 
  filter(!type %in% c("gene", "UTR", "exon")) %>% 
  bind_rows(tibble(type = "intron", mean_width=3300))

breakpoint_coords %>% 
  mutate(queryHits = row_number()) %>% 
  left_join(test) %>% 
  group_by(queryHits) %>% 
  mutate(intron_overlap = all(type == "gene")) %>% 
  mutate(overlap = ifelse(intron_overlap, "intron", type)) %>% 
  ungroup() %>% 
  count(overlap) %>% 
  filter(!overlap %in% c("gene", "exon"))

matrix(
  c(50000,1000,337000,29800, 
    150, 130, 3300, 550), 
  ncol = 4, byrow = T, dimnames = list(c("observed", "expected"), c ("CDS", "fiveprime_utr", "intron", "threeprime_utr"))) %>%
  chisq.test() -> gene_chisq


# Are the breakpoints overlapping a repetetive element? check rmsk. 
rmsk <- read_tsv("~/projects/wgs_fusions/rmsk.txt",col_names = F)

test_rmsk <- GenomicRanges::findOverlaps(
  query = GenomicRanges::GRanges(seqnames = breakpoint_coords$chr, ranges = IRanges::IRanges(start = breakpoint_coords$breakpoint, width = 1)),
  subject = GenomicRanges::GRanges(seqnames = rmsk$X6, ranges = IRanges::IRanges(start = rmsk$X7, end = rmsk$X8))
) %>% as_tibble()  %>% 
  left_join(rmsk %>% 
      mutate(row = row_number()) %>% 
      select(row, X11, X12, X13),
    by = c("subjectHits"="row")
    )

breakpoint_coords %>% 
  mutate(queryHits = row_number()) %>% 
  left_join(test_rmsk) %>% 
  filter(X13%in%c("L1", "L2")) %>% 
  count(X12, sort=T) 
