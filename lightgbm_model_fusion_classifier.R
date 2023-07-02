library(tidyverse)
library(tidymodels)

fusioncatcher <- read_tsv("221124_tcga_output_processed.txt") %>% 
  rbind(read_tsv("221124_scanb_tnbc_output_processed.txt"))

ml_features_lgbm <- fusioncatcher %>% 
  mutate(support_disc_reads = as.factor(support_disc_reads) %>% fct_rev()) %>% 
  select(fusion_id, support_disc_reads, sample, cancer_type, cohort, everything())

lgbm_vars <- ml_features_lgbm %>% names()
lgbm_roles <-  c("ID", "outcome", "ID",  "predictor", "ID", rep("predictor", length(lgbm_vars)-5))

# split data based on cohort
set.seed(123)
splits_lgbm <- group_initial_split(ml_features_lgbm %>% select(-c(fiveprimepartner, threeprimepartner)), group = cohort)
training_lgbm <- training(splits_lgbm) 
cv_lgbm <- vfold_cv(training_lgbm, strata = support_disc_reads)

lgbm_recipe <- training_lgbm %>% 
  recipe(support_disc_reads ~ .) %>% 
  update_role(fusion_id, sample, cohort, new_role = "ID") %>%
  step_dummy_extract(Fusion_description, sep = ",") %>% 
  step_dummy_extract(Fusion_finding_method, sep = ";") %>% 
  step_zv()

# There are lots of aliases for the hyperparameters, the comments show the names from the lightgbm documentation
library(bonsai)
lgbm_model <- 
  boost_tree(
    trees = tune(),  # n_estimators, high number, will be determined by early stopping
    tree_depth = tune(), # max_depth() , range between 4-64
    mtry = tune(),  # feature_fraction_bynode
    min_n = tune(), # min_data_in_leaf
    learn_rate = tune() # learning rate
  ) %>% 
  set_engine(
    "lightgbm", 
    num_leaves = tune(), # 8-4096
    scale_pos_weight = tune(),  # 2-100
  ) %>% 
  set_mode("classification")

lgbm_workflow <- 
  workflow() %>% 
  add_recipe(lgbm_recipe) %>% 
  add_model(lgbm_model)

set.seed(123)
lgbm_grid <- 
  grid_max_entropy(
    trees(),
    finalize(mtry(), training_lgbm),
    min_n(),
    tree_depth(),
    num_leaves(),
    learn_rate(),
    scale_pos_weight(c(1,5)),
    size = 50
  )

set.seed(123)
lgbm_results <- lgbm_workflow %>% 
  finetune::tune_race_anova(
    resamples = cv_lgbm, 
    grid = lgbm_grid,
    metrics = metric_set(pr_auc),
    control = finetune::control_race(
      verbose = T, 
      verbose_elim = T, 
      save_pred = T)
  )

finetune::plot_race(lgbm_results) +
  theme_classic()

lgbm_best <-lgbm_results %>% 
  collect_metrics() %>% 
  filter(mean == max(mean)) %>% 
  slice_head(n= 1)

# > lgbm_best
# # A tibble: 1 × 13
# mtry trees min_n tree_depth learn_rate num_leaves scale_pos_weight .metric .estimator  mean     n std_err .config              
# <int> <int> <int>      <int>      <dbl>      <int>            <dbl> <chr>   <chr>      <dbl> <int>   <dbl> <chr>                
# 37   290    36          9     0.0380         32             2.82 pr_auc  binary     0.801    10 0.00718 Preprocessor1_Model36

lgbm_results %>% 
  collect_predictions(parameters = lgbm_best) %>% 
  pr_curve(support_disc_reads, .pred_TRUE) %>% 
  filter(recall > 0.02) %>% 
  ggplot(aes(precision, recall)) +   
  geom_line() +  
  theme_classic()
  # pr_curve(support_disc_reads, .pred_TRUE) %>% 
  # mutate(f1 = 2*(precision*recall)/(precision+recall)) %>% 
  # arrange(-f1) %>% 
  # filter(precision > recall)

lgbm_results %>% 
  collect_predictions(parameters = lgbm_best) %>% 
  pr_curve(support_disc_reads, .pred_TRUE) %>% 
  mutate(f1 = 2*(precision*recall)/(precision+recall)) %>% 
  arrange(-f1) %>% 
  filter(precision > recall)
  
lgbm_results %>% 
  collect_predictions(parameters = lgbm_best) %>%   
  mutate(pred = factor(.pred_TRUE > 0.20) %>% fct_rev) %>% 
  precision(support_disc_reads, pred)
  conf_mat(support_disc_reads, pred)
  recall(support_disc_reads, pred)
  f_meas(support_disc_reads, pred)
  
  pr_auc(support_disc_reads, .pred_TRUE)
  count(support_disc_reads, pred)
  mcc(support_disc_reads, pred)

pr_curve_plot <- lgbm_results %>% 
  collect_predictions(parameters = lgbm_best) %>% 
  pr_curve(support_disc_reads, .pred_TRUE) %>% 
  filter(recall > 0.02) %>% 
  mutate(split = "Training (TCGA)")
  

### last fit -----

lgbm_last_fit <- lgbm_workflow %>% 
  finalize_workflow(lgbm_best)

set.seed(123)
lgbm_last_fit <- lgbm_last_fit %>% 
  last_fit(splits_lgbm)

lgbm_last_fit %>% 
  collect_predictions() %>% 
  mutate(pred = factor(.pred_TRUE > 0.20) %>% fct_rev) %>% 
  pr_auc(support_disc_reads, .pred_TRUE)
  precision(support_disc_reads, pred)
  kap(support_disc_reads, pred)
  recall(support_disc_reads, pred)
  conf_mat(support_disc_reads, pred)
  f_meas(support_disc_reads, pred)

lgbm_last_fit %>% 
  collect_predictions() %>% 
  pr_curve(support_disc_reads, .pred_TRUE) %>% 
  mutate(f1 = 2*(precision*recall)/(precision+recall)) %>% 
  arrange(-f1) %>% 
  filter(precision > recall) %>% 
  View()

lgbm_last_fit %>% 
  collect_predictions() %>% 
  pr_curve(support_disc_reads, .pred_TRUE) %>% 
  ggplot(aes(precision, recall)) +
  geom_line() +  
  theme_classic() +
  labs(title = "PR AUC = 0.78")
  
lgbm_last_fit %>% 
  collect_predictions()

lgbm_model_importance <- lgbm_last_fit %>% 
  extract_fit_engine() %>% 
  lgb.importance() %>% 
  as_tibble()

lgbm_model_importance %>% 
  write_tsv("final_lgbm_model_feature_importance.txt")

lgbm_model_importance %>% 
  slice_head(n = 20) %>% 
  mutate(Category = case_when(
    str_detect(Feature, "Spanning|description|finding|anchor|in_sample|Predicted_effect") ~ "Chimeric mRNA features",
    str_detect(Feature, "fpkm") ~ "Expression information", 
    str_detect(Feature, "junction|length|chr|cancer_type") ~ "Genetic features",
    T ~ "Other"
  )) %>% 
  ggplot(aes(Gain, fct_reorder(Feature, Gain), fill = Category)) + 
  geom_col() + 
  theme_classic() + 
  labs(y = NULL) +
  scale_fill_brewer(palette = "Dark2")

ggsave("machine_learning/plots/final_lgbm_model_feature_importance_top20.pdf",
  device = "pdf",
  width = 140,
  height = 70,
  units = "mm")

lgbm_last_fit %>% 
  extract_workflow() %>% 
  saveRDS("")

lgbm_last_fit %>% 
  saveRDS("machine_learning/221207_final_lgbm_model.RDS")
lgbm_last_fit

lgbm_last_fit %>% 
  extract_fit_engine() %>% 
  saveRDS.lgb.Booster("machine_learning/221207_final_lgbm_model_booster.RDS")

pr_curve_plot %>% 
  rbind(lgbm_last_fit %>% 
      collect_predictions() %>% 
      pr_curve(support_disc_reads, .pred_TRUE) %>% 
      mutate(split = "Testing (SCAN-B)")) %>% 
  ggplot(aes(precision, recall, col = split)) +
  geom_line() + 
  theme_classic()

fusioncatcher_predict <- read_tsv("fusioncatcher_all_filtered_tcga_fusions_processed.txt")

predictions <- lgbm_last_fit %>% 
  extract_workflow() %>% 
  predict(fusioncatcher_predict, type = "prob") 

predictions %>%
  bind_cols(fusioncatcher_predict) %>% 
  write_tsv("machine_learning/fusioncatcher_all_filitered_tcga_fusions_processed.txt")

fusioncatcher_predict <- predictions %>%
  bind_cols(fusioncatcher_predict) 

fusioncatcher_predict %>% 
  mutate(pred = factor(.pred_TRUE > 0.2) %>% fct_rev) %>% 
  count(pred, cancer_type) %>% 
  View

lgbm_last_fit %>% 
  extract_workflow() %>% 
  predict(head(fusioncatcher_predict), type = "prob")


# -------------------------------------------------------------------------

model <- read_rds("machine_learning/221207_final_lgbm_model.RDS")
model %>% 
  collect_predictions() %>%
  mutate(pred = as_factor(.pred_TRUE > 0.11) %>%
      fct_rev()) %>% 
  pr_curve(support_disc_reads, .pred_TRUE) %>% 
  ggplot(aes(precision, recall)) + 
  geom_line() +
  coord_equal() +
  theme_classic() +
  labs(x = "Precision", y = "Recall")

ggsave("machine_learning/plots/final_fit_pr_curve.pdf", 
  device = "pdf", 
  width = 70,
  height = 70,
  units = "mm")

model |> 
  collect_predictions() |> 
  mutate(pred = as_factor(.pred_TRUE > 0.2) %>%
      fct_rev()) |> 
  accuracy(support_disc_reads, pred)
  conf_mat(support_disc_reads, pred)
  f_meas(support_disc_reads, pred)

model %>% 
  extract_fit_engine() %>% 
  lightgbm::lgb.importance()

lgbm_booster <- readRDS.lgb.Booster("machine_learning/221207_final_lgbm_model_booster.RDS")
lgbm_booster %>% 
  lgb.importance() %>% 
  as_tibble() %>% 
  View()
