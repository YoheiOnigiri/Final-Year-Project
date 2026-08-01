library(tidyverse)
library(grid)
library(here)


plot_path = here('output')
rds_path = here('data', 'processed')

model1_list = c(
  'BM_test3_20260218_124508',
  'BM_test3_20260218_132156',
  'BM_test3_20260218_135829',
  'BM_test3_20260218_143534',
  'BM_test3_20260218_151240',
  'BM_test3_20260218_154935',
  'BM_longrun1_20260303_064355',
  'BM_chain2_20260408_205107',
  'BM_chain3_20260408_210640'
)

model1_res = model1_list %>%
  map(~ {
    filename = paste0(.x, '.rds')
    readRDS(here(rds_path, filename))
  })


# 置換用のヘルパー関数を作成
rename_sigma_param <- function(res) {
  
  # 1. final_params の名前
  if (!is.null(names(res$final_params))) {
    names(res$final_params)[names(res$final_params) == "sigma"] <- "sigma_re"
  }
  
  # 2. chain の列名
  if (!is.null(colnames(res$chain))) {
    colnames(res$chain)[colnames(res$chain) == "sigma"] <- "sigma_re"
  }
  
  # 3. config$init_params の名前
  if (!is.null(names(res$config$init_params))) {
    names(res$config$init_params)[names(res$config$init_params) == "sigma"] <- "sigma_re"
  }
  
  # 4. config$init_cov の行名と列名
  if (!is.null(rownames(res$config$init_cov))) {
    rownames(res$config$init_cov)[rownames(res$config$init_cov) == "sigma"] <- "sigma_re"
    colnames(res$config$init_cov)[colnames(res$config$init_cov) == "sigma"] <- "sigma_re"
  }
  
  return(res)
}

# mapを使って全データに適用
model1_res_updated <- model1_res %>%
  map(rename_sigma_param)

# purrr::walk2を使って、更新されたリストを元のファイル名・場所で上書き保存
walk2(model1_res_updated, model1_list, ~{
  filename = paste0(.y, '.rds')
  saveRDS(.x, here(rds_path, filename))
})

# 完了メッセージを出しておくと安心です
message("All RDS files have been successfully updated and overwritten!")
