## Feature Barcoding analysis

#### These scripts help you demultiplexing your feature barcoding data when you have hashtags of donors.

All the information goes into the configuration file (YAML format). There is an example (config.yaml) with comments regarding the files' format.

There are two main R scripts:
- _demux_seurat.R_, performes the demultiplexing of hashtags using Seurat's MULTIseqDemux and further making sure we classify doublets correctly and rescuing barcodes labeled as negative.
- _summary.R_, summarises the classified data and gets the barcode names ready to merge it into an aggregated data.

The files you need to prepare are:
1. Sample sheet ('samples' in the YAML file).
	- Table with columns 'gex' (path), 'capture' (path), and name.
	- Can be avoided if your samples have consistent names you can use to select pairs of Ab capture (CITE) and gene expression (Gex) data. Then you just indicate the Cell Ranger output folder of containing CITE and Gex folders in the 'count_info' part of config.yaml.

NOTES:
1. The pipeline relies on the samples having the follwing patterns in their name:
	- **Gex**: to use Cell Ranger's "count" routine.
	- **CITE**: to use Cell Ranger's "count" routine and identify it as Feature Barcode.
2. For now this has to be a two step process.

### Install
Clone this repository (your ~/bin folder is a good place).
```
git clone https://github.com/vijaybioinfo/ab_capture.git
cd ab_capture
```

Make sure your config template is pointing to where you have the pipeline.
This will also add the run.sh script as an alias.
```
sh locate_pipeline.sh
```

### Run the pipeline
After you've added the necessary information to the YAML file you can call the pipeline.
```
ab_capture -y /path/to/project/config.yaml
Rscript /path/to/ab_capture/summary.R -h
Rscript /path/to/ab_capture/summary.R -c /path/to/yout/results/project_name_100th \
  --selected /path/to/aggr_name/outs/aggregation.csv \ # <- optional
  --metadata /path/to/metadata_donor.csv~donor \ # <- optional
  --tag_str treatment~donor~hashtag_n~hashtag_id # check the help
```
