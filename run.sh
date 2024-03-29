#!/bin/bash

###############################
# Hashtag data demultiplexing #
###############################

# This script creates jobs to run the Hashtag demultiplexing anaylisis for a
# given pair of expression and hashtag data sets

# set -euo pipefail

function usage () {
    cat >&2 <<EOF

USAGE: $0 [-y] [options]
  -y <config file> : Path to the YAML config file. Required.
  -v Verbose.
  -h Print the usage info.

EOF
}
# initial : makes this loop silent and now requires '?)'
# ${opt} is each option and ${OPTARG} its the argumet (if a colon is there ${opt}:)
VERBOSE=FALSE
while getopts ":y:vh" opt; do
  case ${opt} in
    y) CONFIG_FILE=${OPTARG};;
    v) VERBOSE=TRUE;;
    h) usage; exit 1;;
    \?) echo "No -${OPTARG} argument found."; usage; exit 1;;
  esac
done
if [[ ${OPTIND} -eq 1 ]] ; then
    usage; exit 1
fi

#### Parameters #### -----------------------------------------------------------
function read_yaml(){
  sed 's/#.*//g' ${1} | grep ${2}: | sed 's/.*:[^:\/\/]//; s/\"//g'
}
OUTPUT_DIR="$(read_yaml ${CONFIG_FILE} output_dir)"
PROJECT_ID="$(read_yaml ${CONFIG_FILE} project_id)"
PROJECT_ID=${PROJECT_ID%_}
MAX_COUNT_MIN="$(read_yaml ${CONFIG_FILE} max_count_min)"
OUTPUT_DIR=${OUTPUT_DIR%/}/${PROJECT_ID}
FOLD_CHANGE="$(read_yaml ${CONFIG_FILE} fold_change)"
ABODIES="$(read_yaml ${CONFIG_FILE} subset_tags)"
SSHEET=${OUTPUT_DIR}/library.csv

## Job coordination
SUBMIT="$(read_yaml ${CONFIG_FILE} submit)"
DEPEND="$(read_yaml ${CONFIG_FILE} dependency)"
WALLTIME="$(read_yaml ${CONFIG_FILE} walltime)"
if [[ "${WALLTIME}" == "0" ]]; then WALLTIME=01:00:00; fi
MEM="$(read_yaml ${CONFIG_FILE} mem)"
NODES="$(read_yaml ${CONFIG_FILE} nodes)"
PPN="$(read_yaml ${CONFIG_FILE} ppn)"
CIDS="NONE"

### Tool variables
EXEC_R="$(read_yaml ${CONFIG_FILE} exec_r)"
PIPELINE_DIR="$(dirname "${0}")"

echo ' '
echo -e "\033[0;36m**** Vijay Lab - LJI 2020\033[0m"

echo -e "\033[0;36m------------------------------- PRESENTING PARAMETERS -------------------------------\033[0m"
echo "Paths to samples: ${SSHEET}"
echo "Output path: ${OUTPUT_DIR}"
echo "Minimum count for top feature count: ${MAX_COUNT_MIN}"
echo "Ratio with the second highest hashtag: ${FOLD_CHANGE}"
echo "Specified antibodies pattern: ${ABODIES}"
echo "Job system requirements"
echo "Rscript: ${EXEC_R}"
echo "Time: ${WALLTIME}"
echo "Memory: ${MEM}"
echo "Nodes: ${NODES}"
echo "Sumitting jobs? ${SUBMIT}"
echo "Processors: ${PPN}"
echo "Submit: ${SUBMIT}"
echo -e "\033[0;36m------------------------------- --------------------- -------------------------------\033[0m"
echo; echo
echo "------------------------ Initialising Hashtag demultiplexing ------------------------"
# Clean and create directory
if [ ! -d ${OUTPUT_DIR} ]; then mkdir --parents $OUTPUT_DIR; fi
if [ ! -d ${OUTPUT_DIR}/scripts ]; then mkdir ${OUTPUT_DIR}/scripts; fi
cd ${OUTPUT_DIR}
echo "Working at: `pwd`"

echo "Fetching samples"
if [[ "${VERBOSE}" == "TRUE" ]]; then
  sh ${PIPELINE_DIR}/create_library_csv.sh -y ${CONFIG_FILE} -v
else
  sh ${PIPELINE_DIR}/create_library_csv.sh -y ${CONFIG_FILE}
fi

EDATA=(`tail -n +2 ${SSHEET} | cut -d, -f1`)
CAPTURE=(`tail -n +2 ${SSHEET} | cut -d, -f2`)
TVAR=`wc -l ${SSHEET} | sed 's/ .*//'`
NSAMPLES=(`seq 0 $((TVAR-2))`)
if [ `head -n 1 ${SSHEET} | grep -o "," | wc -l` -ge 2 ]; then
  echo "Sample names given."
  LNAME=(`tail -n +2 ${SSHEET} | cut -d, -f3`)
else
  LNAME=${NSAMPLES[@]}
fi

echo -e "\033[0;33mAnalysing ${#NSAMPLES[@]} samples\033[0m"
for IT in ${NSAMPLES[@]}; do
  echo -e "Sample: \033[0;32m${LNAME[IT]}\033[0m"
  JOBFILE=${OUTPUT_DIR}/scripts/DemuxHT_${LNAME[IT]}

  if [[ -s ${OUTPUT_DIR}/${LNAME[IT]}/step_0_annotation.rdata ]] && \
     [[ "$(echo "${SUBMIT}" | grep -E "force|f" | wc -l)" == "0" ]] # if it's not force
  then
    echo -e "\033[0;31mResults already present\033[0m"; continue
  fi

  echo "Job file: ${JOBFILE}.sh"

  wget -q https://raw.githubusercontent.com/vijaybioinfo/cellranger_wrappeR/main/routine_template.sh -O ${JOBFILE}.sh

  sed -i 's|{username}|'${USER}'|g' ${JOBFILE}.sh
  sed -i 's|{sampleid}|'${LNAME[IT]}'|g' ${JOBFILE}.sh
  sed -i 's|\/\.\.||g' ${JOBFILE}.sh
  sed -i 's|{routine_pbs}|DemuxHT|' ${JOBFILE}.sh
  sed -i 's|{outpath}|'${OUTPUT_DIR}'|g' ${JOBFILE}.sh
  echo "Pushing critical line..."
  sed -i 's|{routine_params}|'${EXEC_R}' '${PIPELINE_DIR}'/demux_seurat.R --edata='${EDATA[IT]}' --capture='${CAPTURE[IT]}' --outdir='${OUTPUT_DIR}' --min_count='${MAX_COUNT_MIN}' --ratio_second='${FOLD_CHANGE}' --prefix='${LNAME[IT]}' --abodies="'${ABODIES}'"|g' ${JOBFILE}.sh
  sed -i 's|{outpath}|'${OUTPUT_DIR}'|g' ${JOBFILE}.sh

  sed -i 's|{walltime}|'${WALLTIME}'|g' ${JOBFILE}.sh
  sed -i 's|{nodes}|'${NODES}'|g' ${JOBFILE}.sh
  sed -i 's|{ppn}|'${PPN}'|g' ${JOBFILE}.sh
  sed -i 's|{mem}|'${MEM}'|g' ${JOBFILE}.sh
  # ls -loh scripts

  if [[ "$(echo "${SUBMIT}" | grep -E "FALSE|no|n" | wc -l)" != "0" ]]; then
    echo "Check it out"; continue # do not submit
  fi
  rm --force ${JOBFILE}.*.txt
  if [[ "${DEPEND}" != "" ]]; then
    DEPEND=$(qsub -W depend=afterok:${DEPEND} ${JOBFILE}.sh)
  else
    CID=$(qsub ${JOBFILE}.sh)
  fi; CID=`echo ${CID} | sed 's/\..*//'`
  if [ "${CIDS}" != 'NONE' ]; then CIDS=${CIDS}:${CID}; else CIDS=${CID}; fi
  echo "Job ID: ${CID}";
  echo
done

# Take the aggregations
# ROUTINE_PARAMS0="${EXEC_R} ${PIPELINE_DIR}/summary.R -c ${OUTPUT_DIR} -m ${MAX_COUNT_MIN}"
if [[ "$(grep "aggregation:" ${CONFIG_FILE} | wc -l)" != "0" ]]; then
  ${EXEC_R} ${PIPELINE_DIR}/summary_jobs.R -y ${CONFIG_FILE} -p "${PIPELINE_DIR}/"
fi
