#!/usr/bin/R

################################
# Antibody/hashtag exploration #
################################

# This code will cuantify the overlap between two libraries of Gex and Antibody
# Capture + some other metrics

optlist <- list(
  optparse::make_option(
    opt_str = c("-e", "--edata"), type = "character",
    help = "10x gene expression library."
  ),
  optparse::make_option(
    opt_str = c("-c", "--capture"), type = "character",
    help = "Feature Barcodes."
  ),
  optparse::make_option(
    opt_str = c("-o", "--outdir"), type = "character",
    help = "Out put directory."
  ),
  optparse::make_option(
    opt_str = c("-m", "--min_count"), type = "numeric", default = 100,
    help = "Minimum count for top feature count (min(max(CB))) in cell barcode."
  ),
  optparse::make_option(
    opt_str = c("-r", "--ratio_second"), type = "numeric", default = 3,
    help = "Fold change between the first and second hashtag."
  ),
  optparse::make_option(
    opt_str = c("-p", "--prefix"), type = "character", default = NULL,
    help = "Prefix. Sample name matching the Gex and CITE."
  ),
  optparse::make_option(
    opt_str = c("-a", "--abodies"), type = "character", default = NULL,
    help = "Antibodies to take into account."
  ),
  optparse::make_option(
    opt_str = c("-d", "--separator"), type = "character", default = "-",
    help = paste0("Separator to identify metadata in hashtag names,\n\t\t",
    "e. g., in D1-ITU-C0256-6 is '-', which is the default.")
  )
)
optparse <- optparse::OptionParser(option_list = optlist)
opt <- optparse::parse_args(optparse)

# opt <- list(
#   edata        = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_biopsy/raw/cellranger/count/Biopsy1_Hu_45X_2P_Gex/outs/filtered_feature_bc_matrix",
#   capture      = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_biopsy/raw/cellranger/count/Biopsy1_Hu_45X_2P_CITE/outs/raw_feature_bc_matrix",
#   outdir       = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma_biopsy/results/ab_demux/biop1to7_100th",
#   min_count    = 100,
#   ratio_second = 3,
#   prefix       = "Biopsy1_Hu_45X_2P",
#   abodies      = "NIH",
#   separator    = "-"
# )

### Loading packages and functions ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R'); rm(.Last)
# Link is loading: dircheck, joindf, remove.factors
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
})
theme_set(theme_cowplot())

### Pre-processing parameters ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str(opt)
if(is.null(opt$abodies)) opt$abodies = ""
if(is.null(opt$prefix)){
  prefix <- basename(gsub("outs.*", "", unlist(opt[1:2])))
  prefix <- gsub("_Gex$|_CITE$", "", prefix)
  prefix <- paste0(unique(prefix), collapse = "_VS_")
}else{
  prefix <- opt$prefix
}

### Reading data ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_dir <- paste0(
  if(!grepl("scratch|beegfs", getwd(), ignore.case = TRUE)){
    cat("No scratch folder involved; careful about temp files...\n")
    dircheck(opt$outdir)
  }else{ "./" }, prefix
)
dir.create(output_dir, recursive = TRUE); setwd(output_dir)
cat("Working in:", getwd(), "\n")

# Load in the UMI matrix
writeLines(text = opt$capture, con = "_capture_path")
htos_count <- Seurat::Read10X(opt$capture)
str(htos_count)

# Checking if a subset of features will be analysed
if(grepl("c\\(", opt$abodies)) eval(parse(text = paste("abodies =", opt$abodies))) else abodies <- opt$abodies
abodies <- paste0(abodies, collapse = "|")
found_bodies <- grepl(abodies, rownames(htos_count))
if(any(found_bodies)){
  cat("Filtering by specified antibodies:", abodies, "\n")
  tvar <- rownames(htos_count)[!found_bodies]
  htos_count <- htos_count[found_bodies, ]
  if(length(tvar) > 0) warning("Discarding antibodies: ", commas(tvar))
}

### Initial visualisation ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("Range of UMI:", paste0(range(htos_count), collapse = " to "), "\n");
min_thr <- min(c(10, matrixStats::rowMaxs(as.matrix(htos_count)) - 1))
filtered_cb <- matrixStats::colMaxs(as.matrix(htos_count)) > min_thr
cat(sum(filtered_cb), "of", ncol(htos_count), "who's top feature has more than", min_thr, "counts.\n")
htos_count <- htos_count[, filtered_cb]
cat("Range after:", paste0(range(htos_count), collapse = " to "), "\n");

ddf <- data.frame(t(as.matrix(htos_count)))
ddf$Total <- rowSums(ddf)
ddf <- reshape2::melt(ddf)
ddf$value[ddf$value == 0] <- 1
p <- ggplot(ddf, aes(x = value, fill = variable)) +
  geom_density(alpha = .5) +
  # geom_density(aes(y = ..scaled..), alpha = .5) +
  geom_vline(xintercept = unique(c(opt$min_count, 500))) +
  facet_wrap(~variable, ncol = 1) +
  theme_minimal() + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + labs(fill = "Hashtag", x = "UMI") + scale_x_log10()
pdf(paste0("step_0_distribution_1pre.pdf"), width = 7, height = 10)
print(p)
graphics.off()

filtered_cb <- matrixStats::colMaxs(as.matrix(htos_count)) > opt$min_count
cat(sum(filtered_cb), "of", ncol(htos_count), "who's top feature has more than", opt$min_count, "counts\n")
if(sum(filtered_cb) == 0){
  stop("Nothing to do here... not enough counts per hashtag to perform the analysis")
}
htos_count <- htos_count[, filtered_cb]

ddf <- data.frame(t(as.matrix(htos_count)))
ddf$Total <- rowSums(ddf)
ddf <- reshape2::melt(ddf)
ddf$value[ddf$value == 0] <- 1
p <- ggplot(ddf, aes(x = value, fill = variable)) +
  geom_density(alpha = .5) +
  geom_vline(xintercept = unique(c(opt$min_count, 500))) +
  facet_wrap(~variable, ncol = 1) +
  theme_minimal() + theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + labs(fill = "Hashtag", x = "UMI") + scale_x_log10()
pdf(paste0("step_0_distribution_2post.pdf"), width = 7, height = 10)
print(p)
graphics.off()

ddf <- reshape2::melt(t(as.matrix(htos_count)))
ddf$value <- log2(ddf$value + 1)
p <- ggplot(data = ddf, aes_string(x = "Var2", y = "value", fill = "Var2")) +
  geom_violin(width = 0.8, alpha = 0.7, trim = TRUE, adjust = 1, scale = 'width') +
  geom_boxplot(width=0.1, fill = "white", alpha = 0.25, outlier.shape = NA, color = "black") +
  labs(x = NULL, y = "log2(UMI + 1)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  )
pdf(paste0("step_0_vlnplot_per_feature.pdf")); print(p); graphics.off()

### Checking overlap ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
tvar <- readLines(list.files(opt$edata, pattern = "barcode", full.names = TRUE))
gex_count <- t(data.frame(mock = 1:length(tvar), row.names = tvar))
joint.bcs <- intersect(colnames(gex_count), colnames(htos_count))
sum_statement <- paste0(
  "Gex: ", length(colnames(gex_count)),
  "\nHashtag: ", length(colnames(htos_count)),
  "\nIntersection: ", length(joint.bcs)
)
tvar <- data.frame(
  Gex = length(colnames(gex_count)),
  HT = length(colnames(htos_count)),
  Intersection = length(joint.bcs)
)
write.csv(tvar, file = paste0("step_0_intersection.csv"))
gexnames = colnames(gex_count); rm(gex_count)
cat(sum_statement, "\n")

### Setup Seurat object ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("Creating Seurat HTO object.\n")
ht_object <- CreateSeuratObject(
  counts = htos_count, assay = 'HTO',
  min.cells = ceiling(ncol(htos_count) * 0.001),
  names.delim = opt$separator # this is the separator I use for my hashtag names
);
preserve_names <- rownames(htos_count)
names(preserve_names) <- rownames(ht_object)
rm(htos_count); gc()

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ht_object <- NormalizeData(ht_object, assay = "HTO", normalization.method = "CLR")

## Demultiplex cells based on HTO enrichment
# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using the
# default settings
cat("Demultiplexing.\n")
void <- try(ht_object <- HTODemux(ht_object, assay = "HTO", positive.quantile = 0.99))
void <- try(ht_object <- MULTIseqDemux(ht_object, assay = "HTO", autoThresh = TRUE, maxiter = 10))
if(class(void) != 'try-error'){
  tvar <- as.character(ht_object$MULTI_ID) %in% c('Negative', 'Doublet')
  ht_object$MULTI_classification.global <- ifelse(tvar, as.character(ht_object$MULTI_ID), 'Singlet')
}; rm(void); gc() # PROBLEM! Duplicating object.

## Visualize demultiplexing results
# Global classification results
mode_classes <- c(HTO = "hash.ID", MULTI = "MULTI_ID")
mode_classes <- mode_classes[mode_classes %in% colnames(ht_object[[]])]
if(1){
  for(i in names(mode_classes)){
    cat(i, "classification:")
    print(unname(reshape2::melt(table(ht_object[[paste0(i, "_classification.global")]]))))
  }
  if(length(mode_classes) == 2){
    table(
      ht_object@meta.data[, paste0(names(mode_classes)[1], "_classification.global")],
      ht_object@meta.data[, paste0(names(mode_classes)[2], "_classification.global")]
    )
  }
}
ht_object$origlib <- basename(gsub("outs.*", "", opt$edata))
annot <- ht_object@meta.data[, grepl("MULTI|HTO|hash|origlib", colnames(ht_object[[]]))]
cat("Calculating fold-change: 1st / 2nd.\n")
ids_umis <- data.frame(data.table::rbindlist(lapply(1:nrow(annot), function(thiscell){
  assigned_id <- as.character(annot[thiscell, tail(mode_classes, 1)])
  all_ids <- ht_object@assays$HTO@counts[, rownames(annot[thiscell, ])]
  rest_ids <- all_ids[!names(all_ids) %in% assigned_id]
  assigned_umi <- c(all_ids[assigned_id]) # take hashtag ID (it's NA if Negative or Doublet)
  y <- data.frame(
    hash.ID = names(which.max(all_ids)), # this is the highest regardless of the assigned
    # first ID's UMI, take assigned; if not assigned, taking the first in rest
    firstUMI = unname(ifelse(is.na(assigned_umi), max(rest_ids), assigned_umi)),
    # second ID, take the top in rest; if not assigned, taking the second in rest
    secondID = unname(ifelse(is.na(assigned_umi), names(sort(rest_ids))[length(rest_ids) - 1], names(which.max(rest_ids)))),
    secondUMI = unname(ifelse(is.na(assigned_umi), sort(rest_ids)[length(rest_ids) - 1], max(rest_ids)))
  ); colnames(y) <- paste0(tail(names(mode_classes), 1), "_", colnames(y))
  y
})), row.names = rownames(annot))
annot <- remove.factors(joindf(annot, ids_umis))
# annott <- annot
for(i in colnames(annot)){
  if(any(names(preserve_names) %in% annot[, i])){
    annot[, i] <- ifelse(is.na(preserve_names[annot[, i]]), annot[, i], preserve_names[annot[, i]])
  }
}
# tvar <- reshape2::melt(table(annot[, tail(mode_classes, 1)], annott[, tail(mode_classes, 1)]))
# tvar <- remove.factors(tvar); colnames(tvar) <- c("Old", "New", "Number")
# tvar <- tvar[tvar$Number > 0, ]; tvar <- tvar[tvar$Old != tvar$New, ]; print(tvar)
cnames <- grep(tail(names(mode_classes), 1), colnames(annot), value = TRUE)
annot$HT_FoldChange <- annot[, grep("firstUMI$",cnames,v=T)] / annot[, grep("secondUMI$",cnames,v=T)]
annot$HT_FoldChange[is.infinite(annot$HT_FoldChange)] <- max(annot$HT_FoldChange[is.finite(annot$HT_FoldChange)])
# assigned feature does not have the highest UMI count
annot$HT_FoldChange_conflict <- annot[, grep("firstUMI$",cnames,v=T)] < annot[, grep("secondUMI$",cnames,v=T)]
print(summary(annot$HT_FoldChange))
table(annot[, tail(mode_classes, 1)])

# Resolving conflictive barcodes where the assigned doesn't have the highest number of UMI
annot$HT_ID <- as.character(annot[, tail(mode_classes, 1)])
tvar <- annot$HT_FoldChange_conflict & annot$HT_ID != 'Doublet'
if(any(tvar)){ # Swapping conflictive ones that are not doublets
  cat("No. of false highest assigned:", sum(tvar), "\n")
  print(table(annot[tvar, ]$HT_ID))
  annot[tvar, ]$HT_ID <- annot[tvar, paste0(tail(names(mode_classes), 1), "_secondID")]
  # annot[tvar, ]$HT_FoldChange <- max(annot$HT_FoldChange) + 3
  annot[tvar, ]$HT_FoldChange <- annot[tvar, grep("secondUMI$",cnames,v=T)] / annot[tvar, grep("firstUMI$",cnames,v=T)]
}
tvar <- annot$HT_FoldChange >= opt$ratio_second & annot$HT_ID == "Negative"
if(any(tvar)){ # Rescuing Negatives
  cat("Negative declassification:", sum(tvar), "\n")
  tmp <- as.character(annot[, paste0(tail(names(mode_classes), 1), "_hash.ID")])
  annot$HT_ID <- ifelse(tvar, tmp, annot$HT_ID)
}
# if less than the FC threshold and different from negative - # increasing doublets!
# Misses FC > 3 that are doublets already!!!
annot$HT_ID <- ifelse(annot$HT_FoldChange < opt$ratio_second & annot$HT_ID != "Negative", "Doublet", as.character(annot$HT_ID))
tvar <- as.character(annot$HT_ID) %in% c('Negative', 'Doublet') # mark the doublets and negative
annot$HT_classification.global <- ifelse(tvar, as.character(annot$HT_ID), 'Singlet') # newly classified re-named as singlets
annot$HT_classification <- annot$HT_classification.global
head(annot[annot$HT_ID == "Doublet", c(grep("secon|first", colnames(annot), v=T), "HT_ID", "HT_FoldChange")])

if(1){
  cat("Changes with re-classification/fixes")
  print(table(annot[, paste0(tail(names(mode_classes), 1), "_classification.global")], annot$HT_classification.global))
  tvar <- reshape2::melt(table(annot[, tail(mode_classes, 1)], annot$HT_ID));
  tvar <- remove.factors(tvar); colnames(tvar) <- c("Old", "New", "Number")
  tvar <- tvar[tvar$Number > 0, ]; tvar <- tvar[tvar$Old != tvar$New, ]; print(tvar)
  annot$in_gex <- rownames(annot) %in% gexnames
  print(table(annot$in_gex))
  print(reshape2::melt(table(annot$HT_ID)))
  cat("Didn't pass the ratio filter")
  print(table(annot[annot$HT_FoldChange < opt$ratio_second, 'HT_ID']))
}
captured_htdf0$HT_ID.global <- ifelse(
  captured_htdf0$HT_ID %in% c("Doublet", "Negative"), captured_htdf0$HT_ID, "Singlet"
)

### Saving results ### %%%%%%%%%%%%%
save(annot, file = paste0("step_0_annotation.rdata"))

# Calculate tSNE embeddings with a distance matrix
cat("Dimentional reduction.\n")
cat("PCA\n")
ht_object <- ScaleData(ht_object, assay = 'HTO')
if(nrow(ht_object[["HTO"]]) > 2){
  ht_object <- RunPCA(ht_object, assay = 'HTO', features = rownames(ht_object[["HTO"]]))
  cat("t-SNE PCA\n")
  ht_object <- RunTSNE(
    ht_object, assay = 'HTO', features = rownames(ht_object[["HTO"]]), perplexity = 75,
    reduction.name = "tsne_pca",
    tsne.method = "FIt-SNE", fast_tsne_path = '/mnt/BioHome/ciro/bin/FIt-SNE2/bin/fast_tsne'
  )
}
ht_object <- RunUMAP(ht_object, assay = 'HTO', features = rownames(ht_object[["HTO"]]))

ht_object@assays <- ht_object@assays['HTO']
ht_object@meta.data <- joindf(ht_object@meta.data, annot)
save(ht_object, file = paste0("step_6_object.rdata"))

# Group cells based on the max HTO signal
ht_object@meta.data$in_gex <- rownames(ht_object@meta.data) %in% gexnames
annot <- ht_object@meta.data
# mode_classes <- c(HTO = "HTO_maxID", MULTI = "MULTI_ID")
mode_classes <- c(HTO = "hash.ID", MULTI = "MULTI_ID", HT = "HT_ID")
mode_classes <- mode_classes[mode_classes %in% colnames(annot)]
if(length(mode_classes) == 0){
  stop("Nothing to do here. Everything failed. Go cry.")
}
for(i in tail(1:length(mode_classes), 1)){
  mode_class <- mode_classes[i]
  myclassification <- paste0(names(mode_class), "_classification")
  prefixt <- paste0("step_method", names(mode_class))
  cat("Using:", mode_class, "\n")

  cat("Class per hashtag\n")
  annot$tmp <- as.character(annot[, paste0(myclassification, ".global")])
  annot$tmp[!annot$in_gex] <- "Gex_missed"
  mytab <- t(table_pct(annot[, c("tmp", unname(mode_class))]))
  write.table(mytab, file = paste0(prefixt, "_0_table_gex.txt"), sep = "\t", quote = FALSE)
  ddf <- reshape2::melt(mytab[-nrow(mytab), -ncol(mytab)])
  ddf <- ddf[!(grepl("Doublet|Negative", ddf[, 1]) & grepl("Singlet", ddf[, 2])), ]
  ddf <- ddf[!(grepl("Doublet|Negative|Gex_missed", ddf[, 2]) & ddf$value == 0), ]
  tvar <- c("Doublet", "Singlet", "Gex_missed")
  tvar <- tvar[tvar %in% as.character(ddf$Var2)]
  ddf$Var2 <- factor(as.character(ddf$Var2), tvar)

  p <- ggplot(ddf, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_minimal() + theme(
      axis.text.x = element_text(size=10, angle = 45, hjust = 1, face = "bold")
    ) + labs(x = "Tag", y = "Number of cells", fill = "Class", subtitle = sum_statement) +
    scale_fill_brewer(palette = "Set1") + scale_y_log10()
  pdf(paste0(prefixt, "_0_table.pdf"))
  print(p)
  graphics.off()
  cat("Class per hashtag, percentage\n")
  # Seurat annoyingly changes "_" to "-"
  # tvar <- rownames(ht_object[["HTO"]])[rownames(ht_object[["HTO"]]) %in% gsub("_", "-", rownames(mytab))]
  tvar <- rownames(mytab)[gsub("_", "-", rownames(mytab)) %in% rownames(ht_object[["HTO"]])]
  ddf <- mytab[tvar, -ncol(mytab)] / mytab[tvar, ncol(mytab)]
  ddf <- reshape2::melt(ddf * 100)
  tvar <- c("Doublet", "Singlet", "Gex_missed")
  tvar <- tvar[tvar %in% as.character(ddf$Var2)]
  ddf$Var2 <- factor(as.character(ddf$Var2), tvar)
  p <- ggplot(ddf, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_minimal() + theme(
      axis.text.x = element_text(size=10, angle = 45, hjust = 1, face = "bold")
    ) + labs(x = "Tag", y = "Number of cells", fill = "Class", subtitle = sum_statement) +
    scale_fill_brewer(palette = "Set1")
  pdf(paste0(prefixt, "_0_table_percentage.pdf"))
  print(p)
  graphics.off()

  cat("Plotting densities.\n")
  Idents(ht_object) <- unname(mode_class)
  p <- RidgePlot(ht_object, assay = "HTO", features = rownames(ht_object[["HTO"]]), combine = FALSE)
  p <- lapply(p, function(x) x + labs(x = NULL, y = NULL) + theme(legend.position = "none") )
  p[[1]] <- p[[1]] + labs(caption = "Scaled can be deceiving")
  pdf(paste0(prefixt, "_1_ridge_per_tag.pdf"), width = 14, height = 12)
  print(cowplot::plot_grid(plotlist = p))
  graphics.off()

  # Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
  cat("Plotting some scatters per pair.\n")
  hto_pairs <- gtools::combinations(nrow(ht_object[["HTO"]]), r = 2, v = rownames(ht_object[["HTO"]]), set = TRUE, repeats.allowed = FALSE)
  pdf(paste0(prefixt, "_2_scatter_pairs.pdf"), width = 7, height = 7)
  for(i in head(1:nrow(hto_pairs), 10)){
    print(hto_pairs[i, ])
    htos <- paste0("hto_", hto_pairs[i, ])
    p <- FeatureScatter(ht_object, feature1 = htos[1], feature2 = htos[2])
    print(p)
  }; graphics.off()

  # Compare number of UMIs for singlets, doublets and negative cells
  Idents(ht_object) <- paste0(myclassification, ".global")
  cat("UMI per classification (Single, Doublet, Negative).\n")
  p <- VlnPlot(ht_object, features = grep("nCount_", colnames(annot), value = TRUE), pt.size = 0.1, log = TRUE)
  pdf(paste0(prefixt, "_3_vlnplot_classes.pdf"), width = 10, height = 7)
  print(p)
  graphics.off()
  tvar <- round(min(annot[annot$HT_classification.global == "Singlet", "HT_FoldChange"], na.rm = TRUE), 2)
  p <- VlnPlot(ht_object, features = "HT_FoldChange", pt.size = 0.1, log = TRUE) +
    labs(caption = paste("Min fold change in singlets:", tvar))
  pdf(paste0(prefixt, "_3_vlnplot_classes_foldchange.pdf"), width = 10, height = 7)
  print(p)
  graphics.off()

  mylevels = c(Doublet = "red", Singlet = "blue", Negative = "#BEBEBE")[levels(Idents(ht_object))]
  ht_object@meta.data$is_there_gex <- ifelse(ht_object@meta.data$in_gex, "In Gex", "No Gex")
  ht_object@meta.data$UMIs <- log2(annot[, grep("nCount_", colnames(annot), value = TRUE)[1]] + 1)
  # ht_object@meta.data$UMIs[ht_object@meta.data$UMIs<log2(500+1)] <- NA
  for(redu in names(ht_object@reductions)){
    p1 <- DimPlot(ht_object, reduction = redu, cols = mylevels)
    p2 <- DimPlot(ht_object, reduction = redu, group.by = unname(mode_class))
    p3 <- DimPlot(ht_object, reduction = redu, group.by = 'is_there_gex', cols = c("In Gex" = "blue", "No Gex" = "red"))
    p4 <- Seurat::FeaturePlot(ht_object, features = "UMIs", reduction = redu) +
      ggplot2::scale_fill_gradientn(colors = c("lightgrey", "blue"), breaks = pretty(ht_object@meta.data$UMIs))
    pdf(paste0(prefixt, "_4_", redu, "_classes.pdf"), width = 14, height = 12)
    print((p1 | p2) / (p3 | p4))
    graphics.off()
  }

  cat("Heatmap.\n")
  # To increase the efficiency of plotting, you can subsample cells using the num.cells argument
  tvar <- if(names(mode_class) == "HT"){
    ht_object@meta.data[, tvar] <- gsub("_", "-", ht_object@meta.data[, tvar]); "HT_ID"
  }else{ myclassification }
  p <- try(HTOHeatmap(
    object = ht_object, assay = "HTO", ncells = 3000,
    classification = tvar, global.classification = paste0(myclassification, ".global")
  ), silent = TRUE)
  if(class(p)[1] != "try-error"){
    pdf(paste0(prefixt, "_5_heatmap_classes.pdf"), width = 10, height = 7)
    print(p)
    graphics.off()
  }
}
cat("Finished demultiplexing.\n")
