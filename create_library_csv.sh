#!/bin/bash

###############################
# Hashtag data demultiplexing #
###############################

# This script creates jobs to run the Hashtag demultiplexing anaylisis for a
# given pair of expression and hashtag data sets

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
sed 's/#.*//g' "${CONFIG_FILE}" > tmp.csv; CONFIG_FILE=tmp.csv
SAMPLE_SHEET=$(grep sample_sheet: ${CONFIG_FILE} | sed 's/.*:[^:\/\/]//; s/\"//g')
OUTPUT_DIR=$(grep output_dir: ${CONFIG_FILE} | sed 's/.*:[^:\/\/]//; s/\"//g')
PROJECT_ID=$(grep project_id: ${CONFIG_FILE} | sed 's/.*:[^:\/\/]//; s/\"//g')
PROJECT_ID=${PROJECT_ID%_}
MAX_COUNT_MIN=$(grep max_count_min: ${CONFIG_FILE} | sed 's/.*:[^:\/\/]//; s/\"//g')
OUTPUT_DIR=${OUTPUT_DIR%/}/${PROJECT_ID}_${MAX_COUNT_MIN}th
COUNT_DIR=$(grep count_dir: ${CONFIG_FILE} | sed 's/.*:[^:\/\/]//; s/\"//g')
COUNT_DIR=${COUNT_DIR%/}
GEX_DATA=$(grep gex_data: ${CONFIG_FILE} | sed 's/.*:[^:\/\/]//; s/\"//g')
FBARCODE_DATA=$(grep fbarcode_data: ${CONFIG_FILE} | sed 's/.*:[^:\/\/]//; s/\"//g')

if [[ ! -s "${SAMPLE_SHEET:-none}" ]] && [[ -d "${COUNT_DIR:-none}" ]]; then
  if [[ "${VERBOSE}" == "TRUE" ]]; then echo -e "\033[0;35m---------- Creating library CSV.\033[0m"; fi
  SAMPLE_SHEET=${OUTPUT_DIR}/library.csv
  echo "gex,capture,name" > "${SAMPLE_SHEET}"
  ABNAMES=$(ls -d "${COUNT_DIR}"/*CITE)
  for ABNAME in ${ABNAMES[@]}; do
    if [[ "${VERBOSE}" == "TRUE" ]]; then printf " + $(basename "${ABNAME}")"; fi
    NNAME="$(basename "${ABNAME}" | sed -E 's/_CITE//; s|^[0-9]{1,}_||')"
    GEXNAME="$(dirname "${ABNAME}")/$(ls $(dirname "${ABNAME}") | grep "${NNAME}_Gex")"
    if [[ "${VERBOSE}" == "TRUE" ]]; then printf " > $(basename "${GEXNAME}")\n"; fi
    if [[ ! -d "${GEXNAME}" ]]; then echo "WARNING: $(basename "${GEXNAME}") folder doens't exist"; fi
    if [[ ! -d "${ABNAME}" ]]; then echo "WARNING: $(basename "${ABNAME}") folder doens't exist"; fi
    echo "${GEXNAME}/outs/${GEX_DATA}_feature_bc_matrix,${ABNAME}/outs/${FBARCODE_DATA}_feature_bc_matrix,${NNAME}" >> ${SAMPLE_SHEET}
  done
elif [[ -s "${SAMPLE_SHEET:-none}" ]]; then
  echo -e "\033[0;35m---------- Library CSV provided.\033[0m"
  sed 's/"//g' ${SAMPLE_SHEET} > ${OUTPUT_DIR}/library.csv
else
  echo -e "\033[0;31mNo Sample sheet.\033[0m"; exit 1
fi
if [[ "${VERBOSE}" == "TRUE" ]]; then echo "Saved into: ${OUTPUT_DIR}/library.csv"; fi
rm tmp.csv
