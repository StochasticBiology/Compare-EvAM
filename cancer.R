library(readxl)
library(dplyr)
library(tidyr)
library(hyperinf)
# ^ to get this, run remotes::install_github("StochasticBiology/hyperinf")

sf = 2

run.code = FALSE
# this datafile is from
# https://www.nature.com/articles/s41586-019-1907-7/figures/1
# and contains chromosomal aberrations found in multiple cancer samples
df = read_excel("41586_2019_1907_MOESM4_ESM.xlsx", sheet = 2)

# pull this into wide format, keeping ID and cancer type
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
cancer.types = unique(df_wide$histology_abbreviation)

# take a look at an example slice through the dat (kidney) to compare to the article
#mat = as.matrix(df_wide[df_wide$histology_abbreviation == "Kidney-RCC.clearcell",3:ncol(df_wide)])
#mat = as.matrix(df_wide[df_wide$histology_abbreviation == "Myeloid-MPN",3:ncol(df_wide)])
mat = as.matrix(df_wide[df_wide$histology_abbreviation == "CNS-GBM",3:ncol(df_wide)])
mat <- mat[, order(as.numeric(gsub("[^0-9].*", "", colnames(mat))))]
mat <- mat[order(rowSums(mat)), ]
heatmap(mat,
        Rowv = NA,
        Colv = NA,
        col = c("white", "black"),
        scale = "none")

#### 

# pick some pairs of cancer types to compare
this.types = c("Ovary-AdenoCA", "Kidney-RCC.clearcell")
#this.types = c("Panc-Endocrine", "Kidney-RCC.clearcell")
#this.types = c("Panc-Endocrine", "ColoRect-AdenoCA")
#this.types = c("Panc-Endocrine", "CNS-GBM")

if(run.code == TRUE) {
  res.list = list()
  
  for(p1 in 1:(length(cancer.types)-1)) {
    for(p2 in (p1+1):length(cancer.types)) {
      this.types = c(cancer.types[p1], cancer.types[p2])
      
      # subset out this pair of cancer types
      this.df = df_wide[df_wide$histology_abbreviation %in% this.types,]
      
      # we have two choices here -- pick N features most represented in the union of the two types, or the 2 N/2 features most represented in each type individually
      top.together = TRUE
      if(top.together) {
        # the first choice
        n <- 8  
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
        # the second choice
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
      
      # df_new now stores data on the top N features for our pair
      cancer.fits = data.mat = list()
      
      if(nrow(df_new[df_new$histology_abbreviation==this.types[1],]) > 1 &
         nrow(df_new[df_new$histology_abbreviation==this.types[2],]) > 1) {
        
        # fit bootstrapped HyperHMM models to both cancer types
        for(this.type in this.types) {
          data.mat[[this.type]] = as.matrix(df_new[df_new$histology_abbreviation==this.type, 3:ncol(df_new)])
          cancer.fits[[this.type]] = hyperinf(data.mat[[this.type]], boot.parallel = 100)
        }
        
        pair.comp.rel = compare_orderings(cancer.fits[[this.types[1]]], cancer.fits[[this.types[2]]], type = "relative")
        pair.comp.abs = compare_orderings(cancer.fits[[this.types[1]]], cancer.fits[[this.types[2]]], type = "absolute")
        
        #plot_hyperinf_compare_orderings(cancer.fits[[this.types[1]]], cancer.fits[[this.types[2]]])
        pair.list = list(types = this.types,
                         comp.rel = pair.comp.rel,
                         comp.abs = pair.comp.abs)
        res.list[[length(res.list)+1]] = pair.list
      }
    }
  }
  
  save(res.list, file = "res-list-100.Rdata")
} else {
  load("res-list.Rdata")
}

hits = data.frame()
for(i in 1:length(res.list)) {
  #nhit = nrow(res.list[[i]]$comp.rel) + 
  nhit = nrow(res.list[[i]]$comp.abs) 
  if(nhit > 0) {
    hits = rbind(hits, data.frame(ref=i, Count = nhit, 
                                  type.1 = gsub("-", "\n", res.list[[i]]$types[1]),
                                  type.2 = gsub("-", "\n", res.list[[i]]$types[2])))
  }
}

library(ggraph)
library(igraph)

hits.g = graph_from_data_frame(hits[hits$Count>2,c(3,4,2)])
png("diff-graph.png", width=400*sf, height=300*sf, res=72*sf)
ggraph(hits.g, layout = "kk") + geom_edge_link(aes(width=Count), alpha=0.2) + 
  geom_node_text(aes(label=name), size=2.3) + theme_void()
dev.off()

prune.g = delete_vertices(hits.g, V(hits.g)[name == "Myeloid\nMPN"])
png("diff-prune-graph.png", width=400*sf, height=300*sf, res=72*sf)
ggraph(prune.g, layout = "kk") + geom_edge_link(aes(width=Count), alpha=0.2) + 
  geom_node_text(aes(label=name), size=2.3) + theme_void()
dev.off()

ref = 353
this.types = res.list[[ref]]$types
#this.types = c("CNS-PiloAstro", "CNS-GBM")
#this.types = c("Kidney-RCC.papillary", "Kidney-RCC.clearcell")
#this.types = c("CNS-Oligo", "CNS-PiloAstro")
#this.types = c("Breast-AdenoCA", "Thy-AdenoCA")
#this.types = c("Ovary-AdenoCA", "Uterus-AdenoCA")
# subset out this pair of cancer types
this.df = df_wide[df_wide$histology_abbreviation %in% this.types,]

# we have two choices here -- pick N features most represented in the union of the two types, or the 2 N/2 features most represented in each type individually
top.together = TRUE
if(top.together) {
  # the first choice
  n <- 8  
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
  # the second choice
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

# df_new now stores data on the top N features for our pair
cancer.fits = data.mat = list()

recolor2 = scale_fill_manual(values=c("steelblue", "tomato", "white"))
nominorx = scale_x_continuous(breaks = scales::breaks_width(1), minor_breaks = NULL)
nominory = scale_y_continuous(breaks = scales::breaks_width(1), minor_breaks = NULL)

# fit bootstrapped HyperHMM models to both cancer types
for(this.type in this.types) {
  data.mat[[this.type]] = as.matrix(df_new[df_new$histology_abbreviation==this.type, 3:ncol(df_new)])
  cancer.fits[[this.type]] = hyperinf(data.mat[[this.type]], boot.parallel = 100)
}

cancer.co = plot_hyperinf_compare_orderings(cancer.fits[[1]], cancer.fits[[2]],
                                            expt.names = this.types,
                                            thetastep = 2) + scale_x_continuous(breaks=1:8)

# visualise the comparison of inferred orderings
png("cancer-example.png", width=400*sf, height=300*sf, res=72*sf)
print(cancer.co + recolor2)
dev.off()


####### KpAMR example

# pulled from Kp paper
load("example-countries.Rdata")
country.list = country.sub.list

# here we see differences
fit.x = multiple_fits_to_booted_fit(country.list$Gambia)
fit.y = multiple_fits_to_booted_fit(country.list$South_Korea)

kp.amr.plot = plot_hyperinf_compare_orderings(fit.x, fit.y, 
                                              sqrt.trans=TRUE, thetastep=3,
                                              expt.names=c("Gambia", "South Korea"), 
                                              feature.names = fit.x$boots[[1]]$featurenames)

kp.amr.plot + recolor2
png("part-kp-plot.png", width=600*sf, height=400*sf, res=72*sf)
print(kp.amr.plot + recolor2)
dev.off()

####

graph.summary = ggraph(hits.g, layout = "kk") + geom_edge_link(aes(width=Count), alpha=0.2) + 
  geom_node_text(aes(label=name), size=2.3) + theme_void()

png("examples-plot.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(kp.amr.plot + recolor2,
          ggarrange(cancer.co + recolor2, graph.summary, labels=c("B", "C")),
          labels=c("A", ""), nrow = 2, heights=c(1.5,1))
dev.off()

