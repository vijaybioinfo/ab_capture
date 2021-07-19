#!/usr/bin/R

#############################################
# Antibody-hashtag summary/aggregation jobs #
#############################################

library(optparse)

optlist <- list(
  make_option(
    opt_str = c("-y", "--yaml"), type = "character", default = "config.yaml",
    help = "Configuration file."
  ),
  make_option(
    opt_str = c("-s", "--submit"), type = "logical", default = FALSE,
    help = "Submit jobs."
  ),
  make_option(
    opt_str = c("-p", "--pipeline"), type = "character",
    help = "Pipeline location."
  ),
  make_option(
    opt_str = c("-v", "--verbose"), type = "logical", default = TRUE,
    help = "Verbose."
  )
)
optparse <- OptionParser(option_list = optlist)
defargs <- parse_args(optparse)

## Functions ## ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(yaml)
  library(crayon)
})
dircheck <- function(x) ifelse(grepl("\\/$", x), x, paste0(x, "/"))
running_jobs <- function(){
  system("qstat -fu ${USER} | grep -E 'Job_Name|Job Id|job_state' | sed 's/Id: /Id_/g; s/ = /: /g; s/.herman.*/:/g' > ~/.tmp")
  jobs_yaml = yaml::read_yaml("~/.tmp")
  jobs_yaml <- jobs_yaml[sapply(jobs_yaml, function(x) x[['job_state']] ) != "C"]
  jobs_df <- data.frame(
    id = gsub(".*Id_", "", names(jobs_yaml)),
    Name = sapply(jobs_yaml, function(x) x[["Job_Name"]] ),
    stringsAsFactors = FALSE
  ); rownames(jobs_df) <- NULL
  jobs_df
}

## Reading files ## ------------------------------------------------------------
config_file = read_yaml(defargs$yaml)
username <- system("echo ${USER}", intern = TRUE)

output_dir <- paste0(dircheck(config_file$output_dir), config_file$project_id)
setwd(output_dir)
if(defargs$verbose){
  cat(cyan("\n************ Vijay Lab - LJI\n"))
  cat(cyan("-------------------------------------\n"))
  cat(red$bold("------------ Summary/Aggregations\n"))
  str(config_file[!names(config_file) %in% "job"])
  cat("Working at:", getwd(), "\n")
  system("ls -loh")
  cat("\n")
}

template = "https://raw.githubusercontent.com/vijaybioinfo/cellranger_wrappeR/main/routine_template.sh"
template_pbs_con <- file(description = template, open = "r")
template_pbs <- readLines(con = template_pbs_con)
close(template_pbs_con)
if(is.null(defargs$pipeline))
  defargs$pipeline <- dirname(gsub(".*=(.*)", "\\1", grep("--file=", base::commandArgs(), value = TRUE)[1]))
submit = isTRUE(defargs$submit) || isTRUE(grepl("^f|force", config_file$job$submit))

aggregations <- unlist(config_file$aggregation$source)
if(defargs$verbose) cat("Aggregations", length(aggregations), "\n\n")
if(defargs$verbose) cat("--------------------------------------\n")
for(aggr in aggregations){
  my_aggregations <- if(dir.exists(aggr)){
    my_aggregations <- list.files(path = aggr, full.names = TRUE)
    my_aggregations <- paste0(gsub(".outs.*", "", my_aggregations), "/outs/aggregation.csv")
    tvar <- file.exists(my_aggregations); if(tvar) my_aggregations[tvar] else basename(aggr)
  }else{ aggr }
  for(aggr_i in my_aggregations){
    aggregation_name <- basename(gsub("\\|", "_", gsub(".outs.*|.csv", "", aggr_i)))
    if(defargs$verbose) cat(aggregation_name)
    if(file.exists(paste0(aggregation_name, ".rds")) && !submit){
      if(defargs$verbose) cat(green(" - done\n")); next
    }

    sample_patterns <- if(file.exists(aggr_i)) read.csv(aggr_i, stringsAsFactors = FALSE)[, 1] else aggr_i
    # substitutions for compatibility between Gex and CITE libraries
    sample_patterns <- gsub("^[0-9]{1,}_|_Gex", "", sample_patterns)
    sample_patterns <- paste0(paste0("DemuxHT_", sample_patterns), collapse = "|")
    samples <- list.files(path = "scripts", pattern = "sh$")
    selected_samples <- gsub("DemuxHT_|.sh$", "", samples[grepl(sample_patterns, samples)])
    if(defargs$verbose) cat("\n +", length(selected_samples), "samples\n")

    # Parameters
    routine_pbs_fname = "aggrHT"
    params <- paste0(
      config_file$exec_r,
      " ", defargs$pipeline, "/summary.R",
      " --captures=", getwd(),
      " --min_count=", config_file$demux$max_count_min,
      " --selected='", aggr_i, "'",
      " --tag_str=", config_file$tag_str,
      " ", config_file$aggregation$args
    )

    pbs <- gsub("\\{username\\}", username, template_pbs)
    pbs <- gsub("\\{sampleid\\}", aggregation_name, pbs)
    pbs <- gsub("\\{routine_pbs\\}", routine_pbs_fname, pbs)
    pbs <- gsub("\\{outpath\\}...", output_dir, pbs)
    pbs <- gsub("\\{routine_params\\}", params, pbs)
    for(i in names(config_file$job)){
      job_parm <- config_file$job[[i]]
      job_parm <- if(routine_pbs_fname %in% names(job_parm)) job_parm[[routine_pbs_fname]] else job_parm[[1]]
      pbs <- gsub(paste0("\\{", i, "\\}"), job_parm, pbs)
    }

    running <- try(running_jobs(), silent = TRUE)
    if(class(running) == "try-error") running <- list(Name = "none")
    if(any(grepl(paste0(routine_pbs_fname, "_", aggregation_name, "$"), running$Name))){
      if(defargs$verbose) cat(" - running\n"); next
    }

    pbs_file <- paste0(getwd(), "/scripts/", routine_pbs_fname, "_", aggregation_name, ".sh")
    writeLines(text = pbs, con = pbs_file)

    if(isTRUE(config_file$job$submit) || submit){
      depend <- if(isTRUE(config_file$job$depend %in% running$id)) paste0("-W depend=afterok:", config_file$job$depend)
      depend_routine <- paste0(c(depend, running[grepl(sample_patterns, running$Name), ]$id), collapse = ":")
      if(is.null(depend) && depend_routine != "") depend_routine <- paste0("-W depend=afterok:", depend_routine)
      pbs_command <- paste("qsub", depend_routine, pbs_file)
      if(defargs$verbose) cat("\n", pbs_command, "\n"); system(pbs_command)
      void <- suppressWarnings(file.remove(gsub("sh", "out.txt", pbs_file)))
    }
    if(defargs$verbose) cat("\n")
  }
}
if(defargs$verbose) cat("--------------------------------------\n")
if(defargs$verbose) cat("PBS files at:", paste0(getwd(), "/scripts/"), "\n")
