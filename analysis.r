# Load Libraries and data -----
library(conflicted)
library(readxl)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library("org.Mm.eg.db")
library(circlize)
library(gridtext)
library(clusterProfiler)
library(BiocParallel)
library(writexl)
library(enrichplot)
library("pheatmap")
library(RColorBrewer)
library(gt)
library(sva)

dat <- read_delim("./data/gaston_counts.txt")
meta <- readxl::read_xlsx("./data/gaston_metadata.xlsx")

set.seed(12345)

# data cleaning ----

meta <- meta[, 1:3]
colnames(meta) <- c("group", "drug", "run")
meta <- meta[, c("run", "group", "drug")]

meta <- as.data.frame(na.omit(meta))

rownames(meta) <- meta[, 1]

# if we are going to do this all with mgi need to ensure there are no duplicates
counts <- dat[, -c(1, 2, 4)]

counts <- counts %>%
  group_by(external_gene_name) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup()

sum(duplicated(counts$external_gene_name))

counts <- as.data.frame(counts)
rownames(counts) <- counts$external_gene_name
counts <- counts[, -1]

#filter colnames to match metadata
colnames(counts) <- sub("^.*?-", "", colnames(counts))

## Dealing with the sample renames: ----------
# the sample id changes by the core listed
#
# GA-2 -> GA-29
# GA-4 -> GA-30
# GA-6 -> GA-31
# GA-9 -> GA-32
# GA-17 -> GA-33

# I will reverse the naming conventions applied to the counts data so that it
# matches the original metdata
current <- NULL
for (i in 29:33) {
  current <- c(current, sprintf("GA-%i", i))
}
revert <- c("GA-2", "GA-4", "GA-6", "GA-9", "GA-17")

colnames(counts)[colnames(counts) %in% current] <- revert

# we need to order the rownames and colnames in the same order

counts <- counts[, order(as.numeric(sub("GA-", "", colnames(counts))))]

meta <- meta[order(as.numeric(sub("GA-", "", rownames(counts)))), ]
meta <- na.omit(meta)

#now let's rerun differences and see where we are at

colset <- colnames(counts)
rowset <- rownames(meta)

diff1 <- base::setdiff(colset, rowset) # all elements in counts can be found in meta
diff2 <- base::setdiff(rowset, colset) # there is extra metadata

# Combine the differences
all_differences <- base::union(diff1, diff2)
all_differences # should be 3,13,15,24 which all need to be removed from the metadata

meta <- meta[!rownames(meta) %in% all_differences, ] # removes extra metadata for non-existant samples

meta_grouped <- meta

meta_grouped$treatment <- base::as.factor(paste(
  meta_grouped$group,
  meta_grouped$drug,
  sep = "_"
))

meta_grouped <- meta_grouped |>
  dplyr::filter(drug != "XIB4035") |>
  dplyr::select(treatment)

counts <- counts[, colnames(counts) %in% rownames(meta_grouped)]

counts <- as.matrix(counts)

meta_grouped_con <- meta_grouped
meta_grouped_cis <- meta_grouped

meta_grouped_con$treatment <- relevel(
  meta_grouped$treatment,
  ref = "CONTROL_Vehicle"
)
meta_grouped_cis$treatment <- relevel(
  meta_grouped$treatment,
  "CISPLATIN_Vehicle"
)

# directory organization ----

directories <- c(
  "./results",
  "./results/data",
  "./results/volcano",
  "./results/GSEA",
  "./results/oras",
  "./results/oras/emap",
  "./results/oras/cnet"
)

if (all(dir.exists(directories)) != TRUE) {
  for (i in 1:length(directories)) {
    if (!dir.exists(directories[i])) {
      dir.create(directories[i])
    }
  }
}

# Run DESeq2 ------

dds_con <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_grouped_con,
  design = ~treatment
)
smallestGroupSize <- 4
dds_con <- dds_con[rowSums(counts(dds_con) >= 10) >= smallestGroupSize, ]

# surrogate variable analysis based correction to tissue contamination
dds_con <- estimateSizeFactors(dds_con)

dat_sur <- counts(dds_con, normalized = TRUE)
idx <- rowMeans(dat_sur) > 1
dat_sur <- dat_sur[idx, ]

mod <- model.matrix(~treatment, colData(dds_con))
mod0 <- model.matrix(~1, colData(dds_con))

# Identify and run SVA
# n.sv  <- num.sv(dat_sur, mod, method = "leek") identification of significant
# sv's through this method resulted in 12 which is too high
n.sv <- 3
svobj <- svaseq(dat_sur, mod, mod0, n.sv = n.sv)

# Myh4, Ckm, Tnnt3 are all highly correlated with SV1 around 0.7
# therefore I believe SV1 will account for the noise of muscle contamination well
# we have 4 samples per treatment so I lean towards only including the correction
# for this expected addition of noise

# Update design to include the SVs
dds_con[["SV1"]] <- svobj$sv[, 1]
# dds_con[["SV2"]] <- svobj$sv[, 2]
cor(dds_con$SV1, as.numeric(counts(dds_con["Tnnt3", ], normalized = TRUE)))
# Myh4, Ckm, Tnnt3 are all highly correlated with SV1 around 0.7
# therefore I believe SV1 will account for the noise of muscle contamination well
# we have 4 samples per treatment so I lean towards only including the correction
# for this expected addition of noise

# TODO: consider running with SV2 as well to correct for batch effect from
# core sample reruns. SV1 is clearly associated with muscle contam so we
# can see if SV2 will correct for any batch effect

design(dds_con) <- formula("~ treatment + SV1")

dds_con <- DESeq(dds_con)

#saveRDS(dds_con, "./results/dds_con.rds")

# with cisplatin as reference level
dds_cis <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta_grouped_cis,
  design = ~treatment
)
dds_cis <- dds_cis[rowSums(counts(dds_cis) >= 10) >= smallestGroupSize, ]

dds_cis <- estimateSizeFactors(dds_cis)

dat_sur <- counts(dds_cis, normalized = TRUE)
idx <- rowMeans(dat_sur) > 1
dat_sur <- dat_sur[idx, ]

mod <- model.matrix(~treatment, colData(dds_cis))
mod0 <- model.matrix(~1, colData(dds_cis))

# Identify and run SVA
# n.sv  <- num.sv(dat_sur, mod, method = "leek") identification of significant
# sv's through this method resulted in 12 which is too high
n.sv <- 3
svobj <- svaseq(dat_sur, mod, mod0, n.sv = n.sv)

# Update design to include the SVs
dds_cis[["SV1"]] <- svobj$sv[, 1]
cor(dds_cis$SV1, as.numeric(counts(dds_cis["Tnnt3", ], normalized = TRUE)))
# should be the same just a different reference level

design(dds_cis) <- formula("~ treatment + SV1")

dds_cis <- DESeq(dds_cis)

# saveRDS(dds_cis, "./results/dds_cis.rds")

# Need to essentially run everything below twice...
list_dds <- list(dds_con, dds_cis)

### functions

summarize_dds_res <- function(x, y) {
  sink(sprintf("results/data/%s_filtered_summary.txt", y))
  print(summary(x))
  sink()
}

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

gsea_util <- function(result, name, ont = "BP", x = 1) {
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

  write.csv(gse@result, str_glue("./Results/GSEA/{name}.csv"))

  png(
    str_glue("Results/GSEA/dotplot_{name}.png"),
    width = 800,
    height = 800,
    units = "px",
    pointsize = 12
  )
  print(dotplot(gse, showCategory = 25, x = "NES"))
  dev.off()
}

ora_util <- function(dat, map, cut = 1) {
  degs_up <- dat |>
    dplyr::filter(log2FoldChange >= cut)

  degs_down <- dat |>
    dplyr::filter(log2FoldChange <= -cut)

  ora_up <- enrichGO(
    gene = rownames(degs_up),
    universe = map,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    keyType = "SYMBOL",
    readable = TRUE,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

  ora_down <- enrichGO(
    gene = rownames(degs_down),
    universe = map,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    keyType = "SYMBOL",
    readable = TRUE,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

  out <- list(
    "up" = ora_up,
    "down" = ora_down,
    "degs_up" = degs_up,
    "degs_down" = degs_down
  )
  return(out)
}

ORA_bar <- function(ora, top = 20, color = "sandybrown") {
  # returns a bar plot of an ORA

  tmp <- as.data.frame(ora@result) %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      Description = forcats::fct_reorder(Description, neg_log10_padj) # Order by significance
    ) %>%
    arrange(desc(neg_log10_padj)) %>%
    head(top)

  p <- ggplot(tmp, aes(x = neg_log10_padj, y = Description)) +
    geom_col(fill = color) +
    labs(
      x = "-log10(Adjusted p-value)",
      y = "GO Term",
      title = "ORA: Enriched Biological Process Terms"
    ) +
    coord_cartesian(xlim = c(0, 20)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, face = "plain", color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "dotted", color = "grey"),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )

  return(p)
}

call_plot <- function(ora, name) {
  ORA_bar(ora$up)
  ggsave(str_glue("./results/oras/{name}_ora_up.png"), bg = "white")

  try(
    {
      ora_up_sim <- pairwise_termsim(ora$up)
      emapplot(ora_up_sim)
      ggsave(
        str_glue("./results/oras/emap/{name}_enrichment_map_up.png"),
        bg = "white"
      )
    },
    silent = TRUE
  )

  try(
    {
      cnetplot(
        ora$up,
        categorySizeBy = ~itemNum,
        foldChange = ora$degs_up$log2FoldChange
      )
      ggsave(
        str_glue("./results/oras/cnet/{name}_cnetplot_up.png"),
        bg = "white"
      )
    },
    silent = TRUE
  )

  ORA_bar(ora$down, color = "dodgerblue")
  ggsave(str_glue("./results/oras/{name}_ora_down.png"), bg = "white")

  try(
    {
      ora_down_sim <- pairwise_termsim(ora$down)
      emapplot(ora_down_sim)
      ggsave(
        str_glue("./results/oras/emap/{name}_enrichment_map_down.png"),
        bg = "white"
      )
    },
    silent = TRUE
  )

  try(
    {
      cnetplot(
        ora$down,
        categorySizeBy = ~itemNum,
        foldChange = ora$degs_down$log2FoldChange
      )
      ggsave(
        str_glue("./results/oras/cnet/{name}_cnetplot_down.png"),
        bg = "white"
      )
    },
    silent = TRUE
  )
}

convert_export_ora <- function(ora_list) {
  out <- list()

  for (i in 1:length(ora_list)) {
    new <- list("up" = c(), "down" = c())
    new$up <- as.data.frame(ora_list[[i]]$up@result)
    new$down <- as.data.frame(ora_list[[i]]$down@result)
    out[[i]] <- new
  }

  names(out) <- names(ora_list)
  return(out)
}

for (dds in list_dds) {
  # Extracting results ----

  results_dds_names <- resultsNames(dds)[-(resultsNames(dds) == "Intercept")]
  results_dds_names <- results_dds_names[!grepl("SV\\d", results_dds_names)]

  res_list <- map(results_dds_names, ~ results(dds, name = .x))

  names(res_list) <- results_dds_names

  genes_annot <- dat |>
    dplyr::select(external_gene_name, description) |>
    distinct(external_gene_name, .keep_all = TRUE)

  genes_annot <- as.data.frame(genes_annot)

  rownames(genes_annot) <- genes_annot$external_gene_name
  genes_annot <- genes_annot[, -1, drop = FALSE]

  # raw
  paths <- paste0("./results/data/", results_dds_names, ".csv")
  tmp <- map(
    res_list,
    ~ merge(as.data.frame(.x), genes_annot, by = 'row.names', all = FALSE)
  )
  walk2(tmp, paths, write.csv)

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

  #TODO: gsea_util crashes if there are no significant groupings
  # could throw a try catch around it
  #walk2(res_list, results_dds_names, ~ gsea_util(.x, .y, x = 12))

  # ORAs -------

  map <- rownames(counts)
  ora_list <- map(res_list_filtered, ~ ora_util(as.data.frame(.x), map))

  paths_oras <- paste0(
    "./results/oras/",
    results_dds_names,
    "_ora.xlsx"
  )

  # write_xlsx takes care of the writing both the up and down ora data to one
  # excel file one per each sheet
  # ahhhh but the issue is that I need the ora in full for the plot and then I
  # need the ora@result to write as an excel file

  ora_list_export <- convert_export_ora(ora_list)
  walk2(ora_list_export, paths_oras, write_xlsx)

  # then save the bar graphs for each ora
  walk2(ora_list, results_dds_names, call_plot)
}

# finally have to consider extra comparisons if needed

# Heatmaps to highlight major changes ----

rld_all <- rlog(dds_con, blind = FALSE)

select <- order(
  rowMeans(counts(dds_con, normalized = TRUE)),
  decreasing = TRUE
)[1:50]

df <- base::as.data.frame(colData(dds_con)[, c("treatment")])
df <- data.frame(treatment = colData(dds_con)$treatment)
rownames(df) <- colnames(dds_con)

pheatmap(
  assay(rld_all)[select, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  annotation_col = df
)


# this row means version is just showing the highest difference, in the genes with the most abundant counts

# We can pull a visualize the genes with the most change in overexpression
# and under expression

# over and under expression heat maps ----
# I am ordering top changed by the cisplatin_47 vs Cisplatin treatment
res <- read.csv(
  "./results/data/treatment_CISPLATIN_47_vs_CISPLATIN_Vehicle.csv"
)
select <- res |>
  arrange(log2FoldChange) |>
  dplyr::select(Row.names)

select_downregulated <- as.vector(select[1:50, ])
select_upregulated <- as.vector(select[nrow(select) - 50:nrow(df), ])


# the scale = row option allows us to z-score the heatmap to show the change in
# genes count by standard deviations away from the mean count of the gene.
# I think this better highlights the differences when zooming in on these

pheatmap(
  assay(rld_all)[select_downregulated, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  annotation_col = df
)

pheatmap(
  assay(rld_all)[select_upregulated, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  annotation_col = df
)

# Cisplatin vs Control ----

res <- read.csv(
  "./results/data/treatment_CISPLATIN_Vehicle_vs_CONTROL_Vehicle.csv"
)
select <- res |>
  arrange(log2FoldChange) |>
  dplyr::select(Row.names)

select_downregulated <- as.vector(select[1:50, ])
select_upregulated <- as.vector(select[nrow(select) - 50:nrow(df), ])

pheatmap(
  assay(rld_all)[select_downregulated, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  annotation_col = df
)

pheatmap(
  assay(rld_all)[select_upregulated, ],
  cluster_rows = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  annotation_col = df
)

# Lastly I want to organize a few heatmaps and tables from the dds results by the
# story outlined by the GO terms

# get data in
ora_1 <- read_xlsx(
  "./results/oras/treatment_CISPLATIN_47_vs_CISPLATIN_Vehicle_ora.xlsx",
  sheet = 1
)
ora_2 <- read_xlsx(
  "./results/oras/treatment_CISPLATIN_47_vs_CISPLATIN_Vehicle_ora.xlsx",
  sheet = 2
)

ora_3 <- read_xlsx(
  "./results/oras/treatment_CISPLATIN_Vehicle_vs_CONTROL_Vehicle_ora.xlsx",
  sheet = 1
)

ora_4 <- read_xlsx(
  "./results/oras/treatment_CISPLATIN_Vehicle_vs_CONTROL_Vehicle_ora.xlsx",
  sheet = 2
)

ora_5 <- read_xlsx(
  "./results/oras/treatment_CONTROL_47_vs_CONTROL_Vehicle_ora.xlsx",
  sheet = 1
)

ora_6 <- read_xlsx(
  "./results/oras/treatment_CONTROL_47_vs_CONTROL_Vehicle_ora.xlsx",
  sheet = 2
)

ora_list <- list(ora_6, ora_5, ora_4, ora_3, ora_2, ora_1)
# a better way would have been to isolate the filepaths names, then loop read them
# directly into the list

# grab the rows that match a grepl query

q1 <- "erythrocyte"
q2 <- "leukocyte|myeloid"
q3 <- "neutrophil|humoral immune"
q4 <- "defense|defense response|interferon|innate immune|killing"
q5 <- "metabolic|catabolic|lipid|thermo"
q6 <- "interleukin"
q7 <- "synapse"

query_vector <- c(q1, q2, q3, q4, q5, q6, q7)
# for each ora in my list, run through query for each query in list, add genes,
# to vector/. Save the unique entries of vector to list for same length as number
# of queries
gene_list <- vector(mode = "list", length = length(query_vector))

for (ora in ora_list) {
  for (i in 1:length(query_vector)) {
    tmp <- ora
    filtered_df <- tmp[
      grepl(query_vector[i], tmp$Description, ignore.case = TRUE),
    ]
    genes <- unlist(strsplit(filtered_df$geneID, "/"))

    if (length(genes) > 0) {
      gene_list[[i]] <- unique(c(gene_list[[i]], genes))
    }
  }
}

# now each gene list is associated with a query that we can use in a title
palette <- brewer.pal(4, "RdGy")
colors = list(
  treatment = c(
    CONTROL_Vehicle = palette[1],
    CISPLATIN_47 = palette[2],
    CISPLATIN_Vehicle = palette[3],
    CONTROL_47 = palette[4]
  )
)

# prepping data list for
# need to use the filtered datasets
data_names <- list.files("./results/data", pattern = ".csv", full.names = T)
data_names <- data_names[grepl(
  "vs_control.*padj05.csv$",
  data_names,
  ignore.case = T,
  perl = T
)]

data_list <- list()
for (i in 1:length(data_names)) {
  tmp <- read_csv(data_names[[i]])
  data_list[[i]] <- tmp[, -1]
}

# organization of heatmaps and tables per query
# I have replaced all NAs from the left join with 0's to represent no difference
# from the control state

tab_list <- list()

for (i in 1:length(query_vector)) {
  # make a heatmap and a table
  clean <- str_replace_all(str_glue("{query_vector[i]}"), "\\|", "_")
  clean <- gsub("[*.]", "", clean)
  name <- str_glue("results/Heatmap_{clean}.png")

  png(
    name,
    width = 800,
    height = 800,
    units = "px",
    pointsize = 12
  )
  print(pheatmap(
    assay(rld_all)[gene_list[[i]], ],
    cluster_rows = FALSE,
    show_rownames = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    annotation_col = df,
    main = str_glue("Heatmap of genes in go terms including {clean}"),
    annotation_colors = colors,
    color = colorRampPalette(c("blue", "white", "red"))(100)
  ))
  dev.off()

  # make a nice table for the Cis effects, cis.47 effects, and 47 effects
  # showing the genes from each query arrange by log2foldchange

  # get dataframe you want, filter by gene_list[[i]]
  # need to filter through each comparison we want to add to the table, then merge
  # to the list of genes with a full join

  tab <- data.frame(`Row.names` = gene_list[[i]])
  lfc_colnames <- c()
  for (j in 1:length(data_list)) {
    data <- data_list[[j]] |>
      dplyr::select(`Row.names`, log2FoldChange)
    tab <- left_join(tab, data, by = "Row.names")

    tmp <- str_extract(data_names[j], "(?<=_).*(?=_vs)")
    lfc_colnames[j] <- str_glue("L2FC_{tmp}")
  }
  colnames(tab) <- c("Gene", lfc_colnames)

  tab[is.na(tab)] <- 0

  tab_list[[i]] <- tab

  # this is good code for general tables including everything that comes up in a
  # query of the GO terms

  gt_table <- tab |>
    relocate(L2FC_CISPLATIN_Vehicle, .after = Gene) |>
    mutate(
      Group = case_when(
        L2FC_CISPLATIN_Vehicle >= 0.5 & L2FC_CISPLATIN_47 < 0.5 ~
          "Effects Cisplatin UP with Cisplatin and 47 opposite",
        L2FC_CISPLATIN_Vehicle <= -0.5 & L2FC_CISPLATIN_47 > -0.5 ~
          "Effects Cisplatin DOWN with Cisplatin and 47 opposite",
        L2FC_CISPLATIN_Vehicle < 0.5 & L2FC_CISPLATIN_47 >= 0.5 ~
          "Effects Cisplatin and 47 UP with Cisplatin opposite",
        L2FC_CISPLATIN_Vehicle > -0.5 & L2FC_CISPLATIN_47 <= -0.5 ~
          "Effects Cisplatin and 47 DOWN with Cisplatin opposite",
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(Group)) |>
    group_by(Group) |>
    arrange(Gene) |>
    gt() |>
    tab_header(
      title = md(str_glue(
        "**Differentially Expressed Genes in the GO Terms That Contain the terms {clean}**"
      )),
      subtitle = "Log2 Fold Change across treatment groups compared to the control"
    ) |>
    tab_spanner(
      label = "Treatment LFC",
      columns = starts_with("L2FC")
    ) |>
    cols_label(
      L2FC_CISPLATIN_47 = "Cisplatin and 47",
      L2FC_CISPLATIN_Vehicle = "Cisplatin",
      L2FC_CONTROL_47 = "47",
      Gene = "Gene Symbol"
    ) |>
    fmt_number(
      columns = starts_with("L2FC"),
      decimals = 2
    ) |>
    data_color(
      columns = starts_with("L2FC"),
      fn = scales::col_bin(
        palette = c("#3B6895", "#9FCAE6", "#f7f7f7", "#f4a582", "#d73027"),
        domain = c(-5, 5),
        bins = c(-Inf, -2, -0.5, 0.5, 2, Inf)
      )
    ) |>
    opt_stylize(color = "gray", style = 1)

  tab_fp <- str_glue("./results/table_GO_{clean}.html")
  gtsave(gt_table, tab_fp)
}

# This table format is pretty nice for noting interesting trends.
# TODO: think of more ways to cluster genes to look at
# * those related to pain
# * some specific to DRG's or cisplatin

# decided to try and use a biplot to see if I can match some of the top explanatory
# genes from the pca

# the biplot has not had the correction for the muscle contamination applied
# and so we can see the clear grouping of muscle genes in the third quadrant

library(PCAtools)

mat <- assay(rld_all)

p <- pca(mat, metadata = colData(rld_all), removeVar = 0.1)

PCAtools::biplot(
  p,
  showLoadings = TRUE,
  ntopLoadings = 10,
  colby = 'treatment',
  labSize = 3,
  legendPosition = 'right'
)

ggsave("./results/biplot.png")

# corrected pca
library(limma)

vsd <- vst(dds_con, blind = FALSE)

# 2. Extract the matrix of values
mat <- assay(vsd)

mat_corrected <- removeBatchEffect(
  mat,
  covariates = dds_con$SV1,
  design = model.matrix(~treatment, colData(dds_con))
)

vsd_corrected <- vsd
assay(vsd_corrected) <- mat_corrected

mat <- assay(vsd_corrected)

p <- pca(mat, metadata = colData(vsd_corrected), removeVar = 0.1)

PCAtools::biplot(
  p,
  showLoadings = TRUE,
  ntopLoadings = 10,
  colby = 'treatment',
  labSize = 3,
  legendPosition = 'right'
)

ggsave("./results/biplot.png")
