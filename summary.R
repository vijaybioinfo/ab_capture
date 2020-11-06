#!/usr/bin/R

################################
# Antibody/hashtag exploration #
################################

# This code will cuantify the overlap between two libraries of Gex and Antibody
# Capture + some other metrics

library(optparse)

optlist <- list(
  make_option(
    opt_str = c("-c", "--captures"), type = "character",
    help = "ab_counts output."
  ),
  make_option(
    opt_str = c("-m", "--min_count"), type = "numeric", default = 100,
    help = "Minimum count for top feature count (min(max(CB))) in cell barcode."
  ),
  make_option(
    opt_str = c("-s", "--selected"), type = "character", default = "",
    help = paste0("Pattern in library names to include in the anaylisis.\n\t\t",
    "It can also take aggregation.csv from the aggr routine from Cell Ranger.")
  ),
  make_option(
    opt_str = c("-t", "--tag_str"), type = "character", default = NULL,
    help = paste0("Column names for each section in hashtag names.\n\t\t",
    "Eg.: D1-ITU-C0256-6 gives c('donor', 'ht_severity', 'htID', 'HTN')\n\t\t",
    "or separated by '~', donor~ht_severity~htID~HTN.\n\t\t",
    "Always include 'donor' as one of them.")
  ),
  make_option(
    opt_str = c("-p", "--prefix"), type = "character", default = "",
    help = "Prefix for file names."
  ),
  make_option(
    opt_str = c("-d", "--separator"), type = "character", default = "-",
    help = paste0("Separator to create the metadata, e. g., in \n\t\t",
    "D1-ITU-C0256-6 is '-', which is the default.")
  ),
  make_option(
    opt_str = c("-r", "--replace"), type = "character", default = NULL,
    help = paste0("List of variable names to replace.\n\t\t",
    "Example: list(c('var_name', 'new_name'), c('ugly', 'bonito')).")
  ),
  make_option(
    opt_str = c("-e", "--metadata"), type = "character", default = "none",
    help = paste0("Donors' metadata file. Indicate the donor names\n\t\t",
                  "column after a '~', e. g., file_name.csv~patient")
  ),
  make_option(
    opt_str = c("-v", "--meta_vars"), type = "character", default = "orig~|",
    help = paste0("Label columns. Default is 'orig' but you can indicate it\n\t\t",
                  "either at the beginning ('orig~|'), for a prefix,\n\t\t",
                  "or the end ('|~tag'), for suffix.")
  )
)
optparse <- OptionParser(option_list = optlist)
opt <- parse_args(optparse)

# opt <- list(
#   captures = '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/fungal_allergy/results/ab_demux/all_fgal_100th',
#   min_count = 100,
#   selected = "/mnt/BioAdHoc/Groups/vd-vijay/cramirez/fungal_allergy/raw/NV035/aggr/all_nv035/outs/aggregation.csv",
#   tag_str = "donor~hashtag_n~hashtag_id",
#   separator = "-",
#   metadata = "/home/ciro/fungal_allergy/info/metadata_donor.csv~donor",
#   meta_vars = "orig~|"
# )

str(opt)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

### Functions ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R'); rm(.Last)
# Link is loading: dircheck, theObjectSavedIn, remove.factors, translist, get_grouping_cols
mybarplot <- function(mytab, xax, yax, fax){
  ggplot(mytab, aes_string(x = xax, y = yax, fill = fax)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_minimal() + theme(
      strip.text.y = element_text(angle = 0),
      axis.text.x = element_text(size=10, angle = 45, hjust = 1, face = "bold")
    ) + labs(x = NULL, y = NULL, fill = "Class") +
    scale_fill_brewer(palette = "Set1")
}

setwd(opt$captures); cat("Working at", getwd(), "\n")
### Getting the subset if necessary ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(file.exists(opt$selected)){
  clust_annot <- read.csv(opt$selected, stringsAsFactors = FALSE)
  libnames <- gsub("^[0-9]{1,}_|_Gex", "", as.character(clust_annot[, 1])) # substitions for compatibility
  # This preserves the aggregation file order
  fnames <- paste0(dircheck(opt$captures), libnames, "/", libnames, "_0_annotation.rdata")
  opt$prefix <- paste0(opt$prefix, basename(gsub('.outs.*', '', opt$selected)))
  opt$selected <- paste0(libnames, collapse = "|")
}else{
  if(is.null(opt$prefix)) opt$prefix <- "all"
  fnames <- list.files(pattern = "_0_annotation.rdata", recursive = TRUE, full.names = TRUE)
  if(is.null(opt$selected)) opt$selected <- paste0(dirname(fnames), collapse = "|")
  fnames <- paste0(dircheck(opt$captures), fnames[grepl(opt$selected, fnames)])
}
cat("Selection pattern:", opt$selected, "\n")

cat("ID:", opt$prefix, "\n")
fname <- list.files(pattern = "_path", recursive = TRUE, full.names = TRUE)
cat("Total files:", length(fname), "\n")

### Summary data.frame ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# You just need to declare 'opt' and 'fnames'
cat("Summarising annotation:", fnames, sep = "\n")
captured_htl <- lapply(fnames, function(x){
  if(!file.exists(x)) return("NULL")
  y <- theObjectSavedIn(x)
  rownames(y) <- gsub("\\-.*", paste0("-", which(fnames %in% x)), rownames(y))
  y
})
tvar <- !sapply(sapply(captured_htl, nrow), is.null); tvar
captured_htl <- captured_htl[tvar]
cat("Suffixes in barcodes:\n");
tvar <- data.frame(t(sapply(captured_htl, function(x) gsub(".*\\-(.*)", "-\\1", head(rownames(x), 5)) )))
dimnames(tvar) <- list(basename(dirname(fnames)), paste0("BC", 1:5))
print(tvar)

## Operations
captured_htdf0 <- data.frame(
  data.table::rbindlist(captured_htl, fill = TRUE),
  row.names = unlist(lapply(captured_htl, rownames))
)
table(captured_htdf0$HT_ID)
# captured_htdf0 <- captured_htdf0[grepl("_", captured_htdf0$MULTI_hash.ID), ] # Seurat removes the "_"...
captured_htdf0$HT_ID.global <- ifelse(captured_htdf0$HT_ID %in% c("Doublet", "Negative"), captured_htdf0$HT_ID, "Singlet")
tvar <- captured_htdf0$HT_ID.global != captured_htdf0$MULTI_classification.global
if(1){
  cat("Global fold changes:\n")
  print(sapply(unique(captured_htdf0$HT_ID.global), function(x){
    summary(captured_htdf0[captured_htdf0$HT_ID.global == x, "HT_FoldChange"])
  }))
  cat("Discrepancy between our final classification and MULTI:", sum(tvar), "/", nrow(captured_htdf0))
  print(table(captured_htdf0$HT_ID.global, captured_htdf0$MULTI_classification.global))
}

if(!any(grepl(opt$separator, as.character(captured_htdf0$HT_ID)))){
  cat("Changing separator\n"); opt$separator <- "_"
}
cat("Separator to demultiplexing feature names: ", opt$separator, ".\n", sep = "")
tvar <- strsplit(as.character(captured_htdf0$HT_ID), opt$separator)
maxln <- max(sapply(tvar, length))
cat("Creating donor metadata with", maxln, "columns.\n")
meta_donor <- data.table::rbindlist(lapply(tvar, function(x){
  if(length(x) < maxln && length(x) == 1) x <- rep(x, length.out = maxln)
  if(length(x) < maxln) x <- rep(paste0(x, collapse = "-"), length.out = maxln)
  as.data.frame(t(x))
}))
meta_donor <- remove.factors(data.frame(meta_donor, row.names = rownames(captured_htdf0)))
if(!is.null(opt$tag_str)) colnames(meta_donor) <- translist(opt$tag_str)[[1]]
headtail(meta_donor)
cat("Checking created columns\n")
sapply(meta_donor, table, useNA = 'always')
rnname <- ifelse(grepl("~", opt$metadata), gsub(".*~", "", opt$metadata), 1)
opt$metadata <- gsub("~.*", "", opt$metadata)
if(file.exists(opt$metadata)){
  cat("Adding given metadata\n")
  extra_meta_donor <- read.csv(opt$metadata, stringsAsFactors = FALSE, row.names = rnname)
  extra_meta_donor[extra_meta_donor == ""] <- NA
  extra_meta_donor_e <- extra_meta_donor[meta_donor$donor, ]; rownames(extra_meta_donor_e) <- rownames(meta_donor)
  sapply(extra_meta_donor_e, table, useNA = 'always')
  meta_donor <- joindf(meta_donor, extra_meta_donor_e)
}
headtail(meta_donor)
captured_htdf0 <- joindf(captured_htdf0, meta_donor)

## Final table ## -------------
tvar <- !grepl("^HTO_|^MULTI_|hash", colnames(captured_htdf0))
tvar <- tvar | colnames(captured_htdf0) %in% "MULTI_ID"
captured_htdf <- captured_htdf0[, tvar]
# str(captured_htdf)
tvar <- grepl(opt$meta_vars, colnames(captured_htdf)) # check what to tag
tvar <- tvar & !grepl("HT_c|origlib|in_gex|MULTI|HTO", colnames(captured_htdf))
if(any(tvar)){
  tagged_vars <- get_grouping_cols(
    metadat = captured_htdf,
    onames = colnames(captured_htdf)[tvar],
    maxn = 50,
    v = TRUE
  )
  tvar <- colnames(captured_htdf) %in% c(tagged_vars, "donor")
  cat("Tagging names:", paste0(colnames(captured_htdf)[tvar], collapse = ", "), "\n")
  if(grepl("~\\|", opt$meta_vars)){
    cat("As prefix\n")
    colnames(captured_htdf)[tvar] <- paste0(gsub("~.*", ".", opt$meta_vars), colnames(captured_htdf)[tvar])
  }else{
    cat("As suffix\n")
    colnames(captured_htdf)[tvar] <- paste0(colnames(captured_htdf)[tvar], gsub(".*~", ".", opt$meta_vars))
  }
}
if(!is.null(opt$replace)){
  for(i in translist(opt$replace)){
    colnames(captured_htdf) <- gsub(i[[1]], i[[2]], colnames(captured_htdf))
  }
}
str(captured_htdf)
save(captured_htdf, file = paste0(opt$prefix, ".rdata"))
list.files(pattern = "rdata")
cat("Written to:", paste0(opt$prefix, ".rdata"), "\n\n")

captured_htdf$tmp <- ifelse(captured_htdf$in_gex == "Missed", "Missing_HT", captured_htdf$HT_classification.global)
captured_htdf$tmp <- gsub("_", " ", captured_htdf$tmp)
p <- ggplot(data = captured_htdf, aes(x = tmp, y = log2(nCount_HTO + 1), fill = tmp)) +
  geom_violin(width = 0.8, alpha = 0.7, trim = TRUE, adjust = 1, scale = 'width') +
  geom_boxplot(width=0.1, fill = "white", alpha = 0.25, outlier.shape = NA, color = "black") +
  labs(x = "Category", y = expression("Log"[2]*"(Total UMI + 1)")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
pdf(paste0(opt$prefix, "_0_vlnplot_per_class.pdf"), width = 14, height = 14);
print(p); print(p + facet_wrap(~ origlib));
graphics.off()

### Some plots... ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(TRUE){
  if(!is.null(opt$selected)) fname <- fname[grepl(opt$selected, fname)]
  cat("Checking distributions on", fname, sep = "\n")
  cat("Pre-filtering (>10 UMI)\n")
  summdfl <- lapply(fname, function(x){
    cat(basename(dirname(x)), "\n")
    htos_count <- Seurat::Read10X(readLines(x))
    filtered_cb <- matrixStats::colMaxs(as.matrix(htos_count)) > 10
    cat(sum(filtered_cb), "of", ncol(htos_count), "who's top feature has more than", 10, "counts\n")
    cat("Range:", paste0(range(htos_count), collapse = " to "), "\n");
    ddf <- data.frame(t(as.matrix(htos_count)))
    ddf$Total <- rowSums(ddf)
    ddf$Library <- basename(dirname(x))
    ddf[filtered_cb, ]
  })
  summdfl <- summdfl[sapply(summdfl, nrow) > 1]
  ddf <- reshape2::melt(data.table::rbindlist(summdfl, fill = TRUE))
  ddf <- ddf[!is.na(ddf$value), ]
  ddf$value[ddf$value == 0] <- 1

  cat("Hashtag per library\n")
  p <- ggplot(ddf, aes(x = value, fill = variable)) +
    geom_vline(xintercept = unique(c(opt$min_count, 500))) +
    facet_wrap(~ Library, ncol = ifelse(length(unique(ddf$Library)) > 10, 2, 1)) +
    theme_minimal() +
    labs(fill = "Hashtag", x = "UMI") + scale_x_log10()
  pdf(paste0(opt$prefix, "_summary_distribution_1pre_htxlib.pdf"), width = 12, height = 15)
  print(p + geom_density(alpha = .5) + labs(title = "Not scaled, go to the second page"))
  print(p + geom_density(aes(y = ..scaled..), alpha = .5) + labs(title = "Scaled"))
  graphics.off()

  cat("Library per hashtag\n")
  p <- ggplot(ddf, aes(x = value, fill = Library)) +
    geom_vline(xintercept = unique(c(opt$min_count, 500))) +
    facet_wrap(~ variable, ncol = 2) +
    theme_minimal() +
    labs(fill = "Hashtag", x = "UMI") + scale_x_log10()
  pdf(paste0(opt$prefix, "_summary_distribution_1pre_libxht.pdf"), width = 12, height = 15)
  print(p + geom_density(alpha = .5) + labs(title = "Not scaled, go to the second page"))
  print(p + geom_density(aes(y = ..scaled..), alpha = .5) + labs(title = "Scaled"))
  graphics.off()

  ### Closer look ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cat("Intersections\n")
  fname <- list.files(pattern = "intersection.csv", recursive = TRUE, full.names = TRUE)
  if(!is.null(opt$selected)) fname <- fname[grepl(opt$selected, fname)]
  ddfl <- lapply(fname, function(x){
    y <- read.csv(x, row.names = 1, stringsAsFactors = FALSE)
    y$set <- gsub("_intersection.*", "", basename(x))
    y
  })
  ddfl[[1]]
  ddf <- data.table::rbindlist(ddfl)
  ddf$set <- gsub("_0", "", ddf$set)
  ddf$Recovered_HT <- ddf$Intersection / ddf$Gex
  ddf$Missing_HT <- 1 - ddf$Recovered_HT
  ddf$set <- paste(ddf$set, "-", round(ddf$Recovered_HT * 100, 1))

  ddf1 <- reshape2::melt(ddf[, -c(1:3)])
  p0 <- ggplot(data = ddf1, aes(x=set, y=value, fill=variable)) +
    geom_bar(stat="identity") + coord_flip() +
    scale_fill_brewer(palette = "Set1", direction = -1) +
    theme_minimal() + labs(subtitle = "Percentage of Gex cells wich HT information")
  ddf1 <- reshape2::melt(ddf[, c(1:4)])
  p <- mybarplot(ddf1, xax = 'variable', yax = 'value', fax = 'variable') +
    geom_bar(stat="identity", position=position_dodge()) +
    facet_wrap(~ set, scales = 'free_y') +
    labs(x = NULL, y = NULL, fill = "Type", subtitle = "Number of cells intersecting")
  pdf(paste0(opt$prefix, "_summary_intersection.pdf"), width = 12, height = 11)
  print(p0)
  print(p)
  graphics.off()

  fname <- list.files(pattern = "_0_table_gex", recursive = TRUE, full.names = TRUE)
  if(!is.null(opt$selected)) fname <- fname[grepl(opt$selected, fname)]
  cat("Summary tables:", length(fname), "\n")
  ddfl <- lapply(fname, function(x){
    y <- read.table(x, sep = "\t")
    y <- y[-nrow(y), -ncol(y)]
    y <- cbind(HT = rownames(y), y)
    ddf <- reshape2::melt(y)
    ddf <- ddf[!(grepl("Doublet|Negative", ddf[, 1]) & grepl("Singlet", ddf[, 2])), ]
    ddf <- ddf[!(grepl("Doublet|Negative|Gex_missed", ddf[, 2]) & ddf$value == 0), ]
    ddf$set <- gsub("_0_table_gex.*", "", basename(x))
    ddf
  })
  ddf <- data.table::rbindlist(ddfl)
  ddf$method <- gsub(".*_", "", ddf$set)
  if(any(ddf$method == "MULTI")) ddf <- ddf[ddf$method == "MULTI", ]
  ddf$set <- gsub("_hto|_multi", "", ddf$set, ignore.case = TRUE)
  ddf$HTN <- gsub("HT([0-9]{1,}).*", "\\1", ddf$HT)
  ddf$HTN <- gsub("^0{1,}", "", ifelse(ddf$HTN == "0", "10", ddf$HTN))
  tvar <- strsplit(as.character(ddf$HT), "-")
  maxln <- max(sapply(tvar, length))
  tvar <- lapply(tvar, function(x){
    if(length(x) < maxln && length(x) == 1) x <- rep(x, length.out = maxln)
    if(length(x) < maxln) x <- rep(paste0(x, collapse = "-"), length.out = maxln)
    as.data.frame(t(x))
  })
  ddf <- cbind(ddf, data.table::rbindlist(tvar))
  if(!is.null(opt$tag_str)) colnames(ddf)[-c(1:5)] <- c("All_Hashtags", translist(opt$tag_str)[[1]])

  cat("Collapsed:", length(colnames(ddf)[-c(1:5)]), "\n")
  for(j in colnames(ddf)[-c(1:5)]){
    cat(" -", j, "\n")
    p <- mybarplot(ddf, xax = j, yax = 'value', fax = 'variable') +
      facet_grid(set ~ method, scales = 'free')
    pdf(paste0(opt$prefix, "_summary_table_collapsed_", j, ".pdf"), width = 7, height = 12)
    print(p)#; print(p + scale_y_log10())
    graphics.off()
  }

  cat("Per method:", length(unique(ddf$method)), "\n")
  for(i in unique(ddf$method)){
    cat(" *", i, "\n")
    p <- mybarplot(ddf[ddf$method == i, ], xax = 'HT', yax = 'value', fax = 'variable') +
      facet_wrap(~ set, scales = 'free_y') +
      labs(x = NULL, y = NULL, fill = "Class", title = casefold(i, upper = TRUE))
    pdf(paste0(opt$prefix, "_summary_table_", i, ".pdf"), width = 12, height = 12)
    print(p)
    print(p + scale_y_log10())
    graphics.off()
  }
}
