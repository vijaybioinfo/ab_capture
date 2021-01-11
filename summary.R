#!/usr/bin/R

########################################
# Antibody-hashtag summary/aggregation #
########################################

# This code will generate a report of the overlap between Gex and
# hashtag libraries as well as the metadata for aggregated libraries

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
    "Eg.: D1-ITU-C0256-6 gives c('donor', 'severity', 'hashtag_name', 'hashtag_id')\n\t\t",
    "or separated by '~', donor~severity~hashtag_name~hashtag_id.\n\t\t",
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
    "or the end ('|~tag'), for suffix. Setting to a random\n\t\t",
    "string turns it off.")
  )
)
optparse <- OptionParser(option_list = optlist)
opt <- parse_args(optparse)

suppressPackageStartupMessages(library(crayon))

cat(cyan("\n************ Vijay Lab - LJI\n"))
cat(cyan("-------------------------------------\n"))
cat(red$bold("------------ Aggregating hashtag information\n"))
str(opt)
### Functions ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R'); rm(.Last)
# Link is loading: dircheck, theObjectSavedIn, remove.factors, translist, get_grouping_cols
mybarplot <- function(mytab, xax, yax, fax){
  ggplot(mytab, aes_string(x = xax, y = yax, fill = fax)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() + theme(
      strip.text.y = element_text(angle = 0),
      axis.text.x = element_text(size=10, angle = 45, hjust = 1, face = "bold")
    ) + labs(x = NULL, y = NULL, fill = "Class") +
    scale_fill_brewer(palette = "Set1")
}

setwd(opt$captures); cat("Working at", getwd(), "\n")

cat("Loading libraries\n") ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(crayon)
}); theme_set(theme_cowplot())

### Getting the subset if necessary ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot_names <- list.files(
  pattern = "_0_annotation.rdata",
  full.names = TRUE, recursive = TRUE
)
libnames <- if(file.exists(opt$selected)){
  clust_annot <- read.csv(opt$selected, stringsAsFactors = FALSE)
  opt$prefix <- paste0(opt$prefix, basename(gsub('.outs.*', '', opt$selected)))
  # substitutions for compatibility between Gex and CITE libraries
  gsub("^[0-9]{1,}_|_Gex", "", as.character(clust_annot[, 1]))
}else{
  dirname(annot_names)
}
opt$selected <- paste0(libnames, collapse = "|")
annot_names <- unlist(sapply( # This preserves the aggregation file order
  X = libnames,
  FUN = function(ssample){
    annot_name <- annot_names[grepl(ssample, annot_names)]
    if(length(annot_name) < 1) annot_name <- "no_file.sad"
    # clust_annot[libnames %in% ssample, ]
    paste0(dircheck(opt$captures), annot_name)
  }
))
cat("Selection pattern:", opt$selected, "\n")
if(opt$prefix == "") opt$prefix <- "all"
cat("ID:", opt$prefix, "\n")

### Summary data.frame ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# You just need to declare 'opt' and 'annot_names'
cat("Summarising annotation:", annot_names, sep = "\n")
captured_htl <- lapply(annot_names, function(x){
  if(!file.exists(x)) return("NULL")
  y <- readfile(x)
  rownames(y) <- gsub("\\-.*", paste0("-", which(annot_names %in% x)), rownames(y))
  y
})
tvar <- !sapply(sapply(captured_htl, nrow), is.null); tvar
captured_htl <- captured_htl[tvar]
cat("Suffixes in barcodes:\n");
tvar <- data.frame(
  t(sapply(captured_htl, function(x){
    gsub(".*\\-(.*)", "-\\1", head(rownames(x), 5))
  }))
)
colnames(tvar) <- paste0("Barcode_", 1:5)
print(tvar)

## Operations
captured_htdf0 <- data.frame(
  data.table::rbindlist(captured_htl, fill = TRUE),
  row.names = unlist(lapply(captured_htl, rownames))
)
table(captured_htdf0$HT_ID)
captured_htdf0$HT_ID.global <- ifelse(
  captured_htdf0$HT_ID %in% c("Doublet", "Negative"), captured_htdf0$HT_ID, "Singlet"
)
tvar <- captured_htdf0$HT_ID.global != captured_htdf0$MULTI_classification.global
if(1){
  cat("Global fold changes:\n")
  print(sapply(unique(captured_htdf0$HT_ID.global), function(x){
    summary(captured_htdf0[captured_htdf0$HT_ID.global == x, "HT_FoldChange"])
  }))
  cat(
    "Discrepancy between our final classification and MULTI:",
    sum(tvar), "/", nrow(captured_htdf0)
  )
  print(table(
    captured_htdf0$HT_ID.global,
    captured_htdf0$MULTI_classification.global
  ))
}

if(!any(grepl(opt$separator, as.character(captured_htdf0$HT_ID)))){
  cat("Changing separator\n"); opt$separator <- "_"
}
cat("Separator for metadata columns: ", opt$separator, ".\n", sep = "")
tvar <- strsplit(as.character(captured_htdf0$HT_ID), opt$separator)
maxln <- max(sapply(tvar, length))
cat("Creating donor metadata with", maxln, "columns.\n")
meta_donor <- data.table::rbindlist(lapply(tvar, function(x){
  if(length(x) < maxln && length(x) == 1) x <- rep(x, length.out = maxln)
  if(length(x) < maxln) x <- rep(paste0(x, collapse = "-"), length.out = maxln)
  as.data.frame(t(x))
}))
meta_donor <- remove.factors(
  data.frame(meta_donor, row.names = rownames(captured_htdf0))
)
if(!is.null(opt$tag_str)) colnames(meta_donor) <- translist(opt$tag_str)[[1]]
headtail(meta_donor)
cat("Checking created columns\n")
lapply(meta_donor, table, useNA = 'always')
rnname <- ifelse(grepl("~", opt$metadata), gsub(".*~", "", opt$metadata), 1)
opt$metadata <- gsub("~.*", "", opt$metadata)
if(isTRUE(file.exists(opt$metadata))){
  cat(
    "Adding given metadata:", opt$metadata,
    "\n- Using row names:", rnname, "\n"
  )
  extra_meta_donor <- read.csv(
    opt$metadata, stringsAsFactors = FALSE, row.names = rnname
  )
  extra_meta_donor[extra_meta_donor == ""] <- NA
  extra_meta_donor_e <- extra_meta_donor[meta_donor$donor, ]
  rownames(extra_meta_donor_e) <- rownames(meta_donor)
  sapply(extra_meta_donor_e, table, useNA = 'always')
  meta_donor <- joindf(meta_donor, extra_meta_donor_e)
  headtail(meta_donor)
}
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
    cols <- paste0(gsub("~.*", ".", opt$meta_vars), colnames(captured_htdf)[tvar])
    colnames(captured_htdf)[tvar] <- cols
  }else{
    cat("As suffix\n")
    cols <- paste0(colnames(captured_htdf)[tvar], gsub(".*~", ".", opt$meta_vars))
    colnames(captured_htdf)[tvar] <- cols
  }
}
if(!is.null(opt$replace)){
  for(i in translist(opt$replace)){
    colnames(captured_htdf) <- gsub(i[[1]], i[[2]], colnames(captured_htdf))
  }
}
str(captured_htdf)
saveRDS(captured_htdf, file = paste0(opt$prefix, ".rds"))
write.csv(captured_htdf, file = paste0(opt$prefix, ".csv"), quote = FALSE)
list.files(pattern = "rds")
cat(red("Written to:"), paste0(opt$prefix, c(".rds", ".csv")), "\n\n")

cat(yellow$bold("Plotting report\n")) #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("- Doublet rates\n") #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(is.character(captured_htdf$in_gex)){
  captured_htdf$in_gex = captured_htdf$in_gex == "inGex"
}
tvar <- grep("HT_ID.global|origlib", colnames(captured_htdf))
libxclass <- table(
  captured_htdf[captured_htdf$in_gex, tvar]
)
libxclass <- round(prop.table(libxclass, margin = 1) * 100, 2)
rownames(libxclass) <- paste(rownames(libxclass), "-", libxclass[, 'Doublet'])
ddf1 <- reshape2::melt(libxclass)
colnames(ddf1) <- c("Library", "variable", "value")
ddf1$variable <- factor(ddf1$variable, levels = c("Negative", "Singlet", "Doublet"))

p_double_rate <- ggplot(data = ddf1, aes(x = Library, y = value, fill = variable)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  theme_minimal() + labs(y = "Doublet rate")
pdf(paste0(opt$prefix, "_0_doublet_rate.pdf"), width = 12, height = 11)
print(p_double_rate);
graphics.off()

cat("- Intersections\n") ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnames <- list.files(
  pattern = "intersection.csv",
  recursive = TRUE, full.names = TRUE
)
fnames <- fnames[grepl(opt$selected, fnames)]
ddfl <- lapply(fnames, function(x){
  y <- read.csv(x, row.names = 1, stringsAsFactors = FALSE)
  y$Library <- basename(dirname(x))
  y
})
ddf <- data.table::rbindlist(ddfl)
ddf$Library <- gsub("_0", "", ddf$Library)
# % that the intersection represents in the total Gex
ddf$Recovered_HT <- ddf$Intersection / ddf$Gex
ddf$Missing_HT <- 1 - ddf$Recovered_HT
ddf$Library <- paste(ddf$Library, "-", round(ddf$Recovered_HT * 100, 1))

ddf1 <- reshape2::melt(ddf[, -c(1:3)])
p0 <- ggplot(data = ddf1, aes(x = Library, y = value, fill = variable)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  theme_minimal() + labs(subtitle = "Percentage of Gex cells with HT information")
ddf1 <- reshape2::melt(ddf[, c(1:4)])
ddf1$variable <- factor(x = ddf1$variable, c("Gex", "Intersection", "HT"))
p <- mybarplot(ddf1, xax = 'variable', yax = 'value', fax = 'variable') +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Library, scales = 'free_y') +
  scale_fill_brewer(palette = "Paired") +
  labs(x = NULL, y = NULL, fill = "Type", subtitle = "Number of cells intersecting")
pdf(paste0(opt$prefix, "_0_intersection.pdf"), width = 12, height = 11)
print(p0); print(p)
graphics.off()
tvar <- file.remove(list.files(pattern = "_summary_intersection"))

cat("- UMI levels\n") #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
captured_htdf$tmp <- factor(
  x = ifelse(!captured_htdf$in_gex, "No Gex", captured_htdf$HT_classification.global),
  levels = c("Negative", "No Gex", "Singlet", "Doublet")
)
p <- ggplot(
    data = captured_htdf, aes(x = tmp, y = log2(nCount_HTO + 1), fill = tmp)
  ) + geom_violin(width = 0.8, alpha = 0.7, scale = 'width') +
  geom_boxplot(
    width = 0.1, fill = "white", alpha = 0.25, outlier.shape = NA, color = "black"
  ) + scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(x = "Category", y = expression("Log"[2]*"(Total hashtags' UMI + 1)")) +
  theme(legend.position = "none")
pdf(paste0(opt$prefix, "_0_vlnplot_per_class.pdf"), width = 14, height = 14);
print(p); print(p + facet_wrap(~ origlib));
graphics.off()

fnames <- list.files(pattern = "_path", recursive = TRUE, full.names = TRUE)
cat("Total files:", length(fnames), "\n")
fnames <- fnames[grepl(opt$selected, fnames)]
cat("- Checking distributions on", fnames, sep = "\n")
cat("  Filtering (>", opt$min_count,"UMI)\n")
summdfl <- lapply(fnames, function(x){
  cat(basename(dirname(x)), "\n")
  htos_count <- Seurat::Read10X(readLines(x))
  filtered_cb <- matrixStats::colMaxs(as.matrix(htos_count)) > opt$min_count
  cat(
    "  ", sum(filtered_cb), "of", ncol(htos_count),
    "who's top feature has more than", opt$min_count, "counts\n"
  )
  cat("  Range:", paste0(range(htos_count), collapse = " to "), "\n");
  ddf <- data.frame(t(as.matrix(htos_count)))
  ddf$Total <- rowSums(ddf)
  ddf$Library <- basename(dirname(x))
  ddf[filtered_cb, ]
})
summdfl <- summdfl[sapply(summdfl, nrow) > 1]
ddf <- reshape2::melt(data.table::rbindlist(summdfl, fill = TRUE))
ddf <- ddf[!is.na(ddf$value), ]
ddf$value[ddf$value == 0] <- 1

cat("- Hashtag per library\n")
p <- ggplot(ddf, aes(x = value, fill = variable)) +
  geom_vline(xintercept = unique(c(opt$min_count, 500))) +
  facet_wrap(~ Library, ncol = ifelse(length(unique(ddf$Library)) > 10, 2, 1)) +
  theme_minimal() +
  labs(fill = "Hashtag", x = "UMI") + scale_x_log10()
pdf(paste0(opt$prefix, "_summary_distribution_htxlib.pdf"), width = 12, height = 15)
print(p + geom_density(alpha = .5) + labs(title = "Not scaled, go to the second page"))
print(p + geom_density(aes(y = ..scaled..), alpha = .5) + labs(title = "Scaled"))
graphics.off()

cat("- Library per hashtag\n")
p <- ggplot(ddf, aes(x = value, fill = Library)) +
  geom_vline(xintercept = unique(c(opt$min_count, 500))) +
  facet_wrap(~ variable, ncol = 2) +
  theme_minimal() +
  labs(fill = "Hashtag", x = "UMI") + scale_x_log10()
pdf(paste0(opt$prefix, "_summary_distribution_libxht.pdf"), width = 12, height = 15)
print(p + geom_density(alpha = .5) + labs(title = "Not scaled, go to the second page"))
print(p + geom_density(aes(y = ..scaled..), alpha = .5) + labs(title = "Scaled"))
graphics.off()
tvar <- file.remove(list.files(pattern = "1pre"))

# Per category generated by the hashtags
fnames <- list.files(pattern = "_0_table_gex", recursive = TRUE, full.names = TRUE)
fnames <- fnames[grepl(opt$selected, fnames)]
cat("-", length(fnames), "libraries with HT data\n")
ddfl <- lapply(fnames, function(x){
  y <- read.table(x, sep = "\t")
  y <- y[-nrow(y), -ncol(y)]
  y <- cbind(HT = rownames(y), y)
  ddf <- reshape2::melt(y)
  ddf <- ddf[!(grepl("Doublet|Negative", ddf[, 1]) & grepl("Singlet", ddf[, 2])), ]
  ddf <- ddf[!(grepl("Doublet|Negative|Gex_missed", ddf[, 2]) & ddf$value == 0), ]
  ddf$Library <- paste0(basename(dirname(x)), gsub("step_|_0_table_gex.*", "", basename(x)))
  ddf
})
ddf <- data.table::rbindlist(ddfl)
ddf$method <- gsub(".*method", "", ddf$Library)
ddf$Library <- gsub("method.*", "", ddf$Library, ignore.case = TRUE)
if(any(ddf$method == "MULTI")){
  print(table(ddf$Library, ddf$method))
  cat("One probably failed, so taking only MULTI will make it consistent\n")
  ddf <- ddf[ddf$method == "MULTI", ]
}
tvar <- strsplit(as.character(ddf$HT), "-")
maxln <- max(sapply(tvar, length))
tvar <- lapply(tvar, function(x){
  if(length(x) < maxln && length(x) == 1) x <- rep(x, length.out = maxln)
  if(length(x) < maxln) x <- rep(paste0(x, collapse = "-"), length.out = maxln)
  as.data.frame(t(x))
})
ddf <- cbind(ddf, data.table::rbindlist(tvar))
highlighted_vars <- -c(2:5)
if(!is.null(opt$tag_str)){
  colnames(ddf)[highlighted_vars] <- c("All_Hashtags", translist(opt$tag_str)[[1]])
}
tvar <- as.character(ddf$variable) == "Gex_missed"
ddf$variable <- ifelse(tvar, "No Gex", as.character(ddf$variable))
tvar <- c("Doublet", "Singlet", "No Gex", "Negative")
ddf$variable <- factor(ddf$variable, tvar)

cat("Collapsed:", length(colnames(ddf)[highlighted_vars]), "\n")
for(j in colnames(ddf)[highlighted_vars]){
  cat(" -", j, "\n")
  p <- mybarplot(ddf, xax = j, yax = 'value', fax = 'variable') +
    facet_grid(Library ~ method, scales = 'free')
  fname <- paste0(opt$prefix, "_summary_table_collapsed_", j, ".pdf")
  pdf(fname, width = 10, height = 12)
  print(p)#; print(p + scale_y_log10())
  graphics.off()
}

# cat("Per method:", length(unique(ddf$method)), "\n")
# for(i in unique(ddf$method)){
#   cat(" *", i, "\n")
#   p <- mybarplot(
#     ddf[ddf$method == i, ], xax = 'All_Hashtags', yax = 'value', fax = 'variable'
#   ) +
#     facet_wrap(~ Library, scales = 'free_y') +
#     labs(x = NULL, y = NULL, fill = "Class", title = casefold(i, upper = TRUE))
#   pdf(paste0(opt$prefix, "_summary_table_", i, ".pdf"), width = 12, height = 12)
#   print(p)
#   print(p + scale_y_log10())
#   graphics.off()
# }
