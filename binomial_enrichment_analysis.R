# A possible way to do an enrichment analysis with gene repeats  

# Example:

# A pathway has 4 genes, g1, g3, g4, g5
# We have a list of genes: g1, g2, g4, g4, (notice g4 is repeated)
# Our universe is all expressed genes: g1, g2, g3, g4, g5, g6. 

# In this case, three genes in our list belong to the pathway (g1, g4, g4).
# The probablility of seeing three or more genes that belong to the pathway is:

1 - pbinom(
  q = 2,       # number of genes in the list belonging to the pathway, minus one. 
  size = 4,    # number of genes in our list
  prob = 4/6)  # number of genes in pathway / number of genes in the universe

# This can be simulated to confirm:

# simulate to sanity check
bag <-  c(T, T, T, T, F, F) 
tibble(iter = 1:100000) %>% 
  mutate(n = 4, # total number of genes in our list
    draw5 = map_dbl(n, ~sample(bag, .x, replace = T) %>% which() %>% length()), # sample n times from the bag with replacement
    pass = draw5 >= 3) %>% # iteration passes if it drew a T 3 or more times. 
  count(pass) |> 
  mutate(n = n/sum(n))

# implement in a function:

binom_enrich <- function(list_of_genes, universe, gene_sets){
  
  results <- gene_sets %>% 
    filter(gene%in%universe$gene_id) %>% 
    group_by(term) %>% 
    nest() %>% 
    mutate(
      q = map_dbl(data, ~which(list_of_genes%in%.x$gene) %>% length())-1) %>% 
    filter(q > 0) %>%
    mutate(
      size = length(list_of_genes), 
      prob = map_dbl(data, ~nrow(.x)/nrow(universe))
    ) %>% 
    arrange(-q) %>% 
    ungroup() %>% 
    mutate(one_plus_log_pval = pbinom(q, size, prob, log.p = T)) %>% 
    mutate(data = map(data, ~filter(.x, gene %in% list_of_genes)))
  
  return(results)
}


ml_predictions <- read_tsv("221216_all_tcga_fusions_classifier_predictions.txt")


val_fusions %>% 
  filter(cancer_type == "breast") %>% 
  select(-enrichment_results) %>% 
  mutate(binom_enrich = map2(data, universe, ~binom_enrich(list_of_genes = as_vector(.x), universe =  .y, gene_sets = msigdb_kegg_reactome))) %>% 
  select(binom_enrich) %>% 
  unnest()

test <- val_fusions %>% 
  mutate(partner = str_remove(partner, "partner")) %>%
  mutate(binom_enrich = map2(data, universe, ~binom_enrich(list_of_genes = as_vector(.x), universe =  .y, gene_sets = msigdb_kegg_reactome)))


all_binom_enrichment <- val_fusions %>% 
  mutate(group = "Validated") %>% 
  mutate(partner = str_remove(partner, "partner")) %>% 
  bind_rows(pred_fusions %>% mutate(group = "Predicted")) %>% 
  bind_rows(all_pred_fusions %>% mutate(group = "All predicted")) %>% 
  bind_rows(non_val_fusions %>% mutate(group = "Non-validated")) %>% 
  select(-enrichment_results) %>% 
  mutate(binom_enrich = map2(data, universe, ~binom_enrich(list_of_genes = as_vector(.x), universe =  .y, gene_sets = msigdb_kegg_reactome)))


plot_binom_enrich <-  function(type = "breast", top_n = 10) {
  top10 = all_binom_enrichment %>% 
    filter(cancer_type == type) %>% 
    select(1,2,5,6) %>% 
    unnest(binom_enrich) %>% 
    group_by(group, cancer_type, partner) %>% 
    slice_min(order_by = padjust, n = top_n, with_ties = F) %>% 
    distinct(term) %>% 
    pull(term)
  
  enrichment_plot = all_binom_enrichment %>% 
    filter(cancer_type == type) %>% 
    select(1,2,5,6) %>% 
    unnest(binom_enrich) %>% 
    filter(term %in% top10) %>% 
    mutate(sig = padjust < 0.05) %>% 
    mutate(term = tolower(term)) %>% 
    ggplot(aes(group %>% 
        fct_relevel(c("Validated","Predicted", "All predicted", "Non-validated" )), 
      term %>% fct_inorder(), 
      size = -log10(padjust),
      col = sig)) + 
    geom_point() + 
    theme_minimal() +
    facet_wrap(~partner, ncol = 2) +
    scale_color_viridis_d() +
    labs(x = NULL, y = NULL)
  
  return(enrichment_plot)
  
}

test <- test %>% 
  ungroup %>% 
  distinct(cancer_type) %>% 
  mutate(binom_enrich_plots = map(cancer_type, ~plot_binom_enrich(type = .x, top_n = 10)))


test$binom_enrich_plots[7]

test %>% 
  pwalk(~ggsave(filename = paste0("binom_enirchment_plot_", ..1, "_expressed_universe.pdf"),
    plot = ..2,
    width = 15,
    height = 8))




# -validated fusions ------------------------------------------------------


binom_enrich_results <- fusioncatcher %>% 
  left_join(validation %>% select(fusion_id, support_disc_reads)) %>% 
  filter(support_disc_reads) %>% 
  distinct(cancer_type, sample, fiveprimepartner, threeprimepartner) %>% 
  pivot_longer(cols = c(fiveprimepartner, threeprimepartner), names_to = "partner", values_to = "gene") %>% 
  select(-sample) %>% 
  group_by(cancer_type,partner) %>% 
  nest() %>% 
  left_join(universes %>% rename(universe = data) %>% mutate(cancer_type = str_to_sentence(cancer_type)), by = "cancer_type") %>% 
  mutate(binom_enrich = map2(data, universe, ~binom_enrich(list_of_genes = as_vector(.x), universe =  .y, gene_sets = msigdb_kegg_reactome)))

binom_enrich_results <- binom_enrich_results %>% 
  select(1,2,binom_enrich) %>% 
  unnest(binom_enrich) %>% 
  ungroup() %>% 
  mutate(genes = map_chr(data, ~paste0(as_vector(.x), collapse = "|"))) %>% 
  arrange(-one_plus_log_pval) %>% 
  mutate(rank = row_number(), 
    BH_adjust = (rank/nrow(.)) * 1.0 ) %>% 
  mutate(padjust = 1 - exp(one_plus_log_pval * nrow(.))) %>% 
  mutate(bh_sig = BH_adjust < 0.05, 
    bf_sig = padjust < 0.25) %>% 
  filter(bh_sig)
# count(cancer_type, partner)
# View()




# predicted fusions in RNA+DNA --------------------------------------------

binom_enrich_results_predicted <- validation %>%
  select(-cancer_type) %>% 
  left_join(ml_predictions %>% select(fusion_id, .pred_TRUE, is_banned)) %>% 
  left_join(fusioncatcher %>% select(fusion_id, cancer_type, fiveprimepartner, threeprimepartner)) %>% 
  drop_na(is_banned) %>% 
  filter(.pred_TRUE >= 0.2) %>% 
  distinct(cancer_type, wgs_barcode, fiveprimepartner, threeprimepartner) %>% 
  pivot_longer(cols = c(fiveprimepartner, threeprimepartner), names_to = "partner", values_to = "gene") %>% 
  select(-wgs_barcode) %>% 
  group_by(cancer_type,partner) %>% 
  nest() %>% 
  left_join(universes %>% rename(universe = data) %>% mutate(cancer_type = str_to_sentence(cancer_type)), by = "cancer_type") %>% 
  mutate(binom_enrich = map2(data, universe, ~binom_enrich(list_of_genes = as_vector(.x), universe =  .y, gene_sets = msigdb_kegg_reactome)))

binom_enrich_results_predicted <- binom_enrich_results_predicted %>% 
  select(1,2,binom_enrich) %>% 
  unnest(binom_enrich) %>% 
  select(-data) %>% 
  ungroup() %>% 
  arrange(-one_plus_log_pval) %>% 
  mutate(rank = row_number(), 
    BH_adjust = (rank/nrow(.)) * 1.0 ) %>% 
  mutate(padjust = 1 - exp(one_plus_log_pval * nrow(.))) %>% 
  mutate(bh_sig = BH_adjust < 0.05, 
    bf_sig = padjust < 0.25) %>% 
  filter(bh_sig)

# -all predicted fusions -------------------------------------------------------

binom_enrich_results_all_predicted <- fusioncatcher %>% 
  left_join(ml_predictions %>% 
      mutate(pred = .pred_TRUE >= 0.2) %>% 
      select(fusion_id, pred)
  ) %>% 
  filter(pred) %>% 
  distinct(cancer_type, sample, fiveprimepartner, threeprimepartner) %>% 
  pivot_longer(cols = c(fiveprimepartner, threeprimepartner), names_to = "partner", values_to = "gene") %>% 
  select(-sample) %>% 
  group_by(cancer_type,partner) %>% 
  nest() %>% 
  left_join(universes %>% rename(universe = data) %>% mutate(cancer_type = str_to_sentence(cancer_type)), by = "cancer_type") %>% 
  mutate(binom_enrich = map2(data, universe, ~binom_enrich(list_of_genes = as_vector(.x), universe =  .y, gene_sets = msigdb_kegg_reactome)))

binom_enrich_results_all_predicted <- binom_enrich_results_all_predicted %>% 
  select(1,2,binom_enrich) %>% 
  unnest(binom_enrich) %>% 
  select(-data) %>% 
  ungroup() %>% 
  arrange(-one_plus_log_pval) %>% 
  mutate(rank = row_number(), 
    BH_adjust = (rank/nrow(.)) * 1.0 ) %>% 
  mutate(padjust = 1 - exp(one_plus_log_pval * nrow(.))) %>% 
  mutate(bh_sig = BH_adjust < 0.05, 
    bf_sig = padjust < 0.25) %>% 
  filter(bh_sig)


# -whole fusioncatcher output --------------------------------------------------


binom_enrich_results_all_fusioncatcher <- fusioncatcher %>% 
  mutate(sample = str_sub(barcode, 1, 15 )) %>% 
  distinct(cancer_type, sample, fiveprimepartner, threeprimepartner) %>% 
  pivot_longer(cols = c(fiveprimepartner, threeprimepartner), names_to = "partner", values_to = "gene") %>% 
  select(-sample) %>% 
  group_by(cancer_type,partner) %>% 
  nest() %>% 
  left_join(universes %>% rename(universe = data) %>% mutate(cancer_type = str_to_sentence(cancer_type)), by = "cancer_type") %>% 
  mutate(binom_enrich = map2(data, universe, ~binom_enrich(list_of_genes = as_vector(.x), universe =  .y, gene_sets = msigdb_kegg_reactome)))

binom_enrich_results_all_fusioncatcher <- binom_enrich_results_all_fusioncatcher %>% 
  select(1,2,binom_enrich) %>% 
  unnest(binom_enrich) %>% 
  select(-data) %>% 
  ungroup() %>% 
  arrange(-one_plus_log_pval) %>% 
  mutate(rank = rank(-one_plus_log_pval, ties.method = "min"), 
    BH_adjust = (rank/nrow(.)) * 1.0 ) %>% 
  mutate(padjust = 1 - exp(one_plus_log_pval * nrow(.))) %>% 
  mutate(bh_sig = BH_adjust < 0.05, 
    bf_sig = padjust < 0.25) %>% 
  filter(bh_sig)

# save 

binom_enrich_results_all <- binom_enrich_results %>% 
  mutate(group = "Validated") %>% 
  bind_rows(binom_enrich_results_predicted %>% mutate(group = "Predicted")) %>% 
  bind_rows(binom_enrich_results_all_predicted %>% mutate(group = "All predicted")) %>% 
  bind_rows(binom_enrich_results_all_fusioncatcher %>% mutate(group = "All fusioncatcher"))

binom_enrich_results_all %>% 
  write_tsv("enrichment_analysis/230202_binomial_enrichment_results_all.txt")
rm(binom_enrich_results, binom_enrich_results_all_fusioncatcher, binom_enrich_results_all_predicted, binom_enrich_results_predicted)

binom_enrich_results_all <-read_tsv("enrichment_analysis/230202_binomial_enrichment_results_all.txt")


binom_enrich_results_all %>% 
  filter(cancer_type == "Glioblastoma", partner == "fiveprimepartner") %>% 
  select(group, term, BH_adjust) %>% 
  pivot_wider(names_from = group, values_from = BH_adjust) %>% 
  View()

results_wide <- binom_enrich_results_all %>% 
  select(cancer_type, partner, group, term, BH_adjust) %>% 
  group_by(cancer_type, partner) %>% 
  nest()  %>% 
  mutate(data = map(data, ~pivot_wider(.x, names_from = group, values_from = BH_adjust))) %>% 
  arrange(cancer_type)


# PCA ---------------------------------------------------------------------


library(tidymodels)

results_wide$data[[2]]

pca_plot <- function(df, cancer_type, partner, w_legend = F){
  
  pcaplot <- df %>% 
    mutate(across(-c(term), ~replace_na(.x, 1))) %>% 
    pivot_longer(-term, names_to = "group", values_to = "p_val") %>% 
    pivot_wider(names_from = term, values_from = p_val) %>% 
    recipe(group~.) %>% 
    step_normalize() %>% 
    step_pca(all_numeric(), num_comp = 3) %>% 
    prep() %>% 
    juice() %>% 
    ggplot(aes(x = PC1, y  = PC2, col = group)) + 
    geom_jitter(width = 0.5, height = 0.5) + 
    theme_classic() +
    scale_color_viridis_d() +
    labs(title = paste0(cancer_type, "\n", partner)) +
    theme(legend.position = ifelse(w_legend, "right", "none"),
      title = element_text(size = 8))
  
  return(pcaplot)
}

results_wide <- results_wide %>% 
  bind_rows(results_wide %>% ungroup %>% slice_tail()) %>% 
  ungroup %>% 
  mutate(show_legend = row_number() == 19) %>% 
  slice_head(n = 19) %>% 
  mutate(row = row_number()) %>% 
  mutate(partner = ifelse(row%%2, "5' genes", "3' genes")) %>% 
  mutate(pca = pmap(list(data, cancer_type, partner, show_legend), ~pca_plot(..1, ..2, ..3, ..4)))

results_wide$data[[6]]
results_wide$pca[[6]] 

cowplot::plot_grid(plotlist = results_wide$pca[1:19], nrow = 5, ncol = 4) 

ggsave("~/projects/wgs_fusions/enrichment_analysis/230206_enrichment_pca_plots.pdf", 
  device = "pdf", 
  width = 210,
  height = 200, 
  units = "mm")

results_wide %>% 
  saveRDS("~/projects/wgs_fusions/enrichment_analysis/binomial_enrichment_results_w_pca.rds")

results_wide %>% 
  ungroup() %>% 
  slice_head() %>% 
  select(-c(4,5)) %>% 
  unnest(data) %>% 
  View

tidy_pca <- function(df){
  result <- df %>% 
    prcomp(center = T, scale. = T )
  
  return(result[["rotation"]] %>% 
      as_tibble(rownames = "ID"))
}

results_wide %>% 
  ungroup %>% 
  select(-c(4,5)) %>% 
  unnest(data) %>% 
  mutate(across(-c(term, cancer_type, partner), ~replace_na(.x, 1))) %>% 
  select(4,5,6,7) %>% 
  tidy_pca() %>% 
  ggplot(aes(PC1, PC2, col = ID)) +
  geom_point() + 
  theme_classic() +
  scale_color_viridis_d() +
  labs(title = "All enrichment results") 

