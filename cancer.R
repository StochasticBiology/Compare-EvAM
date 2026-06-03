library(readxl)
library(dplyr)
library(tidyr)
library(hyperinf)

# https://www.nature.com/articles/s41586-019-1907-7/figures/1
df = read_excel("41586_2019_1907_MOESM4_ESM.xlsx", sheet = 2)

df_wide <- df %>%
  select(samplename, chrom_arm, histology_abbreviation) %>%
  distinct() %>%                       # remove duplicate rows if any
  mutate(present = 1) %>%             # indicator for presence
  pivot_wider(
    names_from = chrom_arm,
    values_from = present,
    values_fill = 0
  )
ncol(df_wide)
unique(df_wide$histology_abbreviation)

mat = as.matrix(df_wide[df_wide$histology_abbreviation == "Kidney-RCC.clearcell",3:ncol(df_wide)])
mat <- mat[, order(as.numeric(gsub("[^0-9].*", "", colnames(mat))))]
mat <- mat[order(rowSums(mat)), ]
heatmap(mat,
        Rowv = NA,
        Colv = NA,
        col = c("white", "black"),
        scale = "none")

this.types = c("Ovary-AdenoCA", "Kidney-RCC.clearcell")
#this.types = c("Panc-Endocrine", "Kidney-RCC.clearcell")
#this.types = c("Panc-Endocrine", "ColoRect-AdenoCA")
#this.types = c("Panc-Endocrine", "CNS-GBM")

this.df = df_wide[df_wide$histology_abbreviation %in% this.types,]

top.together = TRUE
if(top.together) {
n <- 8  # number of columns you want to keep

df_top <- this.df %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "col", values_to = "sum") %>%
  arrange(desc(sum)) %>%
  slice_head(n = n) %>%
  pull(col) -> top_cols

df_new <- this.df %>%
  select(samplename, histology_abbreviation, all_of(top_cols))
} else {
  n <- 5  # number of top columns per group
  
  top_cols <- this.df %>%
    select(-samplename) %>%
    pivot_longer(-histology_abbreviation, names_to = "col", values_to = "val") %>%
    group_by(histology_abbreviation, col) %>%
    summarise(sum = sum(val, na.rm = TRUE), .groups = "drop") %>%
    group_by(histology_abbreviation) %>%
    slice_max(sum, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(col) %>%          # <- union across groups
    pull(col)
  
  df_new <- this.df %>%
    select(samplename, histology_abbreviation, all_of(top_cols))
}

cancer.fits = data.mat = list()
for(this.type in this.types) {
  data.mat[[this.type]] = as.matrix(df_new[df_new$histology_abbreviation==this.type, 3:ncol(df_new)])
  cancer.fits[[this.type]] = hyperinf(data.mat[[this.type]], boot.parallel = 10)
}

#plotHypercube.lik.trace(cancer.fits[[2]])

plot_hyperinf_compare_orderings(cancer.fits[[1]], cancer.fits[[2]],
                                expt.names = this.types,
                                thetastep = 3, threshold = 0.33)

