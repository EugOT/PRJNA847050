library(tidyverse)
metrics_summary_SRR19578591 <- read_csv("/data/PRJNA847050/cellranger/SRR19578591/outs/metrics_summary.csv") %>% mutate(Run = "SRR19578591")
metrics_summary_SRR19578590 <- read_csv("/data/PRJNA847050/cellranger/SRR19578590/outs/metrics_summary.csv") %>% mutate(Run = "SRR19578590")
metrics_summary_SRR19578589 <- read_csv("/data/PRJNA847050/cellranger/SRR19578589/outs/metrics_summary.csv") %>% mutate(Run = "SRR19578589")
metrics_summary <-
  bind_rows(
    metrics_summary_SRR19578591,
    metrics_summary_SRR19578590,
    metrics_summary_SRR19578589)

metrics_summary |>
  select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, here::here("metrics_summary.tsv"))

