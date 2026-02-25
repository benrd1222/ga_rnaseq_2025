# Load Libraries and data -----
library(conflicted)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library("org.Mm.eg.db")
#library("biomaRt")
library(circlize)
library(gridtext)
library(clusterProfiler)
library(BiocParallel)
#library(VennDiagram)

dat <- read_delim("./data/gaston_counts.txt")
meta <- readxl::read_xlsx("./data/gaston_metadata.xlsx")

set.seed(12345)

# data cleaning ----

# need to remove the samples that failed intial QC and post sequencing QC

# data organization ----

directories <- c(
  "./results",
  "./results/data",
  "./results/volcano",
  "./results/GSEA"
)

if (all(dir.exists(directories)) != TRUE) {
  for (i in 1:length(directories)) {
    if (!dir.exists(directories[i])) {
      dir.create(directories[i])
    }
  }
}


# Run DESeq2 ------

smallestGroupSize <- 5
dds <- dds[rowSums(counts(dds) >= 10) >= smallestGroupSize, ]

dds <- DESeq(dds)

saveRDS(dds,"./results/dds.rds")


# Extracting results ----


results_dds_names <- resultsNames(dds)[-(resultsNames(dds) == "Intercept")]

res_list <- map(results_dds_names, ~ results(dds, name = .x))

names(res_list) <- results_dds_names

enes_annot <- as.data.frame(genes_annot)
rownames(genes_annot) <- genes_annot$mgi_symbol

# raw
paths <- paste0("./results/data/", results_dds_names, ".csv")
tmp <- map(
  res_list,
  ~ merge(as.data.frame(.x), genes_annot, by = 'row.names', all = FALSE)
)
walk2(tmp, paths, write.csv)
test <- merge(
  as.data.frame(res_list[[1]]),
  genes_annot,
  by = 'row.names',
  all = FALSE
)

#filtered
res_list_filtered <- map(res_list, ~ .x[.x$padj < 0.05 & !is.na(.x$padj), ])
tmp <- map(
  res_list_filtered,
  ~ merge(as.data.frame(.x), genes_annot, by = 'row.names', all = FALSE)
)

paths_filtered <- paste0(
  "./results/data/",
  results_dds_names,
  "_filtered_padj05.csv"
)

walk2(tmp, paths_filtered, write.csv)


# summary of filtered

summarize_dds_res <- function(x, y) {
  sink(sprintf("results/data/%s_filtered_summary.txt", y))
  print(summary(x))
  sink()
}

walk2(res_list_filtered, results_dds_names, summarize_dds_res)

# Making the results into filtered, ordered, and annotated dataframes
cutoff <- 1
res_workable <- list()
for (i in 1:length(res_list_filtered)) {
  tmp <- as.data.frame(res_list_filtered[[i]])
  tmp <- merge(tmp, genes_annot, by = 'row.names', all = FALSE)
  
  tmp <- tmp |>
    mutate(
      dir = case_when(
        log2FoldChange > cutoff ~ "UP",
        log2FoldChange < -1 * cutoff ~ "DOWN",
        T ~ "NS"
      )
    ) |>
    arrange(desc(log2FoldChange))
  
  res_workable[[i]] <- tmp
}

names(res_workable) <- results_dds_names

#Volcano Plots:

volcano_annot <- genes_annot[, -2, drop = FALSE]

volcano_dds <- function(x, y, ...) {
  tmp = merge(as.data.frame(x), volcano_annot, by = 'row.names', all = FALSE)
  
  png(
    sprintf("results/Volcano/%s_.png", y),
    width = 800,
    height = 800,
    units = "px",
    pointsize = 12
  )
  print(EnhancedVolcano(
    tmp,
    lab = tmp[, 1],
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05
  ))
  dev.off()
}

walk2(res_list, results_dds_names, ~ volcano_dds(.x, .y, volcano_annot))

# Clustering -------

# for PCAs we need normalized counts and not LFC
rld_all <- rlog(dds, blind = FALSE)
rld_all_df <- as.data.frame(assay(rld_all))

intgroup = "treatment" 

PCA_plot_all <- DESeq2::plotPCA(
  rld_all,
  intgroup = intgroup,
  returnData = TRUE,
  ntop = 1000
)
percentVar <- round(100 * attr(PCA_plot_all, "percentVar"))

PCA <- ggplot(
  PCA_plot_all,
  aes_string(x = "PC1", y = "PC2", color = intgroup, label = intgroup),
) +
  geom_point(size = 5) +
  scale_x_continuous(
    name = paste0("PC1: ", percentVar[1], "% variance"),
    limits = c(-20, 20)
  ) +
  scale_y_continuous(
    name = paste0("PC2: ", percentVar[2], "% variance"),
    limits = c(-20, 20)
  ) +
  coord_fixed() +
  ggforce::geom_mark_ellipse(aes_string(label = "NULL", color = intgroup)) +
  geom_text_repel(
    size = 7,
    vjust = "inward",
    hjust = "inward",
    show.legend = FALSE,
    point.padding = 10
  ) +
  ggtitle("PCA") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 20)
  )

ggsave('results/basic_PCA.jpg', plot = PCA, width = 15, height = 15)

# GSEAs -----------

# Ensure that the everything is in MGI

# Generalize to produce a GSEA for each res_list entry

# In:
# - res_list, res_list names
# - ontology type, default to BP
# Out:
# - ontology rds
# - ontology csv
# - dotplot

gsea_util <- function(result, names, ont = "BP", x = 1) {
  tmp <- as.data.frame(result)
  tmp <- tmp |>
    arrange(desc(log2FoldChange))
  
  list <- tmp$log2FoldChange
  names(list) <- rownames(tmp)
  
  bp_param <- SnowParam(
    workers = x,
    type = "SOCK",
    RNGseed = 42,
  )
  
  gse <- gseGO(
    geneList = list,
    ont = ont,
    keyType = "SYMBOL",
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = "org.Mm.eg.db",
    pAdjustMethod = "BH",
    seed = TRUE,
    by = "fgsea",
    BPPARAM = bp_param
  )
  
  saveRDS(gse, str_glue("./results/GSEA/{name}.rds"))
  
  write.csv(gse@result, "./Results/GSEA/{name}.csv")
  
  png(
    "Results/GSEA/dotplot_{name}.png",
    width = 800,
    height = 800,
    units = "px",
    pointsize = 12
  )
  print(dotplot(gse, showCategory = 15, x = "NES"))
  dev.off()
}

walk2(res_list, res_dds_list, ~ gsea_util(.x, .y, x = 12))
