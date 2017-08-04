# from cytominer audit.R
median_pairwise_correlation <- function(df, variables, group_by, method = "correlation", K = 20, sigma = 0.5) {
  if (method == "correlation") {
    #print("in median_pairwise_correlation: using correlation")
    df %>%
      dplyr::group_by_(.dots = group_by) %>%
      do(tibble::data_frame(correlation = median(cor(t(as.matrix(.[variables])))))) 
    
  } else  if (method == "affinity") {
    #print("in median_pairwise_correlation: using affinity")
    df %>%
      dplyr::group_by_(.dots = group_by) %>%
      do(tibble::data_frame(correlation = median( SNFtool::affinityMatrix(as.matrix(dist(.[variables])),K = K, sigma = sigma))))
  } else {
    warning("method not found, using method = 'correlation' ")
    df %>%
      dplyr::group_by_(.dots = group_by) %>%
      do(tibble::data_frame(correlation = median(cor(t(as.matrix(.[variables]))))))
  }

}


fraction_strong <- function(population, temp_feature, group_by, method = "correlation"){ 
  #print("correlation")
  correlations <-  population %>%
    median_pairwise_correlation(temp_feature, group_by, method = method, K = 20, sigma = 0.5)
  
  #print("null thrshld")
  null_threshold <- population %>%
    tidyr::unite_("group_by", group_by) %>%
    mutate(group_by = sample(group_by)) %>%
    median_pairwise_correlation(temp_feature, "group_by", method = method, K = 20, sigma = 0.5) %>%
    magrittr::extract2("correlation") %>%
    quantile(0.95, na.rm = TRUE)
  
  result <-
    tibble::data_frame(
      null_threshold = null_threshold,
      fraction_strong = (sum(correlations$correlation > null_threshold) / nrow(correlations)),
      mean_correlation = correlations %>% 
        ungroup() %>% 
        summarise(mean_correlation = mean(correlation)) %>% 
        extract2("mean_correlation"), 
      method = method
      )
  
  return(result)
}

# select all feature and area features 
na_per_column <- function(input_df){ 
  temp_feature <-
    colnames(input_df) %>%
    str_subset("^Area|^Texture|^Track")
  
  temp_df <- input_df %>% 
    select_(.dots = temp_feature) 
  
  mean_na <- colMeans(is.na(temp_df))
  result <- tibble(feature_names = temp_feature, na_per_column = mean_na)
  return(result)
}


fraction_strong_affinity <- function(W, metadata_rows, group_by) {
  Wt <-  W %>% cbind(.,metadata_rows) %>%
    select(-Metadata_id) %>%
    group_by_(.dots = group_by) %>%
    summarise_each(funs(median)) %>%
    ungroup() %>%
    select(-Metadata_condition, -Metadata_dose, -Metadata_matrix, -Metadata_date) %>%
    t 
  
  M <- cbind(metadata_rows, Wt) %>%
    select(-Metadata_id) %>%
    group_by(Metadata_matrix, Metadata_condition, Metadata_dose) %>%
    summarise_each(funs(median)) %>%
    ungroup() %>%
    select(-Metadata_condition, -Metadata_dose, -Metadata_matrix, -Metadata_date) %>%
    as.matrix()
  
  threshold <- M[lower.tri(M)] %>%
    quantile(0.95)
  
  tibble(threshold = threshold, fraction_strong = sum(diag(M) > threshold) / nrow(M), mean_affinity = mean(diag(M)))
}

