# workflowr::wflow_publish("analysis/eda.Rmd", all = TRUE, update = TRUE, republish = TRUE, verbose = TRUE)
# workflowr::wflow_publish("analysis/eda.Rmd", all = TRUE, verbose = TRUE)
workflowr::wflow_build(c("analysis/01A-eda-whole_dataset-fpr_0.001.Rmd", "analysis/01B-eda-whole_dataset-fpr_0.01.Rmd", "analysis/01C-eda-whole_dataset-fpr_0.05.Rmd", "analysis/01D-eda-whole_dataset-fpr_0.1.Rmd"), verbose = TRUE, log_dir = here::here("logs_workflowr"))
