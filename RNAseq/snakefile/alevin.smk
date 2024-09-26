##########################################################################
# Snakemakefile Version:   1.0
# Description:             Snakemake file to run Salmon alevin module
##########################################################################

################## Context ##################
# launch snakemake -s  snakefile_alevin -c(numberofthreads) --config run=absolutepathoftherundirectory without / at the end of the path
# to launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
# every variable defined in the yaml file can be change
# separate multiple variable with a space (ex  --config run=runname var1=0.05 var2=12)
# also use option --configfile another.yaml to replace and merge existing config.yaml file variables

# use -p to display shell commands
# use --lt to display docstrings of rules

# input file = fastq files
# output file = salmon/alvevin count files (quant.sf)
################## Import libraries ##################

########## Note ########################################################################################
# recipe for Salmon index to perform count extraction
# GRch38 assembly and gencode.vxx.transcripts.fa.gz can be found on https://www.gencodegenes.org/human/
# release can be change (actually v41 is used, relase is 07.2022)
# for GRch38 primary assembly v42 : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
# for gencode transcripts : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz
# for mouse https://www.gencodegenes.org/mouse/releases.html
# grep "^>" <(gunzip -c GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
# sed -i.bak -e 's/>//g' decoys.txt
# zcat gencode.vM33.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > gentrome.fa.gz
# salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
# eventually use --keepDuplicates for isoform analysis
# to generate the transcript/gene table
# zgrep "^>" gentrome.fa.gz| cut -d "|" -f 1,6 --output-delimiter=$'\t' - | sed 's/>//g; s/gene_symbol://g; s/"//g' > txp2gene.tsv
########################################################################################################

import os
import glob
import pandas as pd
import json
import csv
import logging
from shutil import copy2
from datetime import datetime
from itertools import product
from collections import defaultdict

################## Configuration file ##################
configfile: "/app/snakefile/alevin_default.yaml"

####################### FUNCTIONS #####################

def populate_dictionary(dictionary, samples_list, extensions_list, files_list, pattern_include=None, pattern_exclude=None, split_index=0):
	for sample in samples_list:
		for ext in extensions_list:
			for file in files_list:
				file_parts = os.path.basename(file).split(".")
				
				if split_index >= len(file_parts):
					continue  # Skip if split_index is out of range
				
				file_base = file_parts[split_index]
				
				if file_base != sample or not os.path.basename(file).endswith(ext):
					continue

				if pattern_exclude and pattern_exclude in file:
					continue

				if pattern_include and pattern_include not in file:
					continue

				dictionary.setdefault(sample, {})[ext] = file


def filter_files(files_list, filter_in=None, filter_out=None):
	return [file_name for file_name in files_list if (not filter_in or filter_in in file_name) and (not filter_out or filter_out not in file_name)]

			
def find_item_in_dict(sample_list, ext_list, dictionary, include_ext, exclude_ext=None):
	""" Function to search in a dictionary for a non-empty file path by iterating through sample_list and ext_list with inclusion and exclusion filters """
	search_result = ""
	
	for sample in sample_list:
		for ext in ext_list:
			try:
				items = dictionary.get(sample, {}).get(ext, [])
				if include_ext in items and (exclude_ext is None or exclude_ext not in items):
					if os.path.exists(items) and os.path.getsize(items) != 0:
						search_result = items
			except KeyError as e:
				print(f"KeyError encountered: {e}")
	
	return search_result

def searchfiles(directory, search_arg, recursive_arg):
	""" Function to search all files in a directory, adding a search arguement append to the directory and a recursive_search options (True/False) """
	return sorted(filter(os.path.isfile, glob.glob(directory + search_arg, recursive=recursive_arg)))

def extractlistfromfiles(file_list, ext_list, sep, position):
	""" Function for creating list from a file list, with a specific extension, a separator and the position of the string we want to extract """
	return list(set(os.path.basename(files).split(sep)[position] for files in file_list if any(files.endswith(ext) for ext in ext_list)))


### END OF FUNCTIONS ###
serviceName = config['serviceName']
runName = os.path.basename(os.path.normpath(config['run']))
date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")
resultDir = f"/app/res/{runName}/{date_time}"
outputDir = config['OUTPUT_DIR'] if config['OUTPUT_DIR'] else config['run']
directories = [resultDir, outputDir]

for directory in directories:
	os.makedirs(directory, exist_ok=True)

# Search files in repository 
files_list = searchfiles(os.path.normpath(config['run']), config['SEARCH_ARGUMENT'],  config['RECURSIVE_SEARCH'])

# Create sample list
sample_list = extractlistfromfiles(files_list, config['PROCESS_FILE'], '.', 0)

# Exclude samples from the exclude_list , case insensitive
sample_list = [sample for sample in sample_list if not any(sample.upper().startswith(exclude.upper()) for exclude in config['EXCLUDE_SAMPLE'])]

# If filter_sample_list variable is not empty, it will force the sample list
if config['FILTER_SAMPLE']:
	sample_list = list(config['FILTER_SAMPLE'])

runDict = defaultdict(dict)
populate_dictionary(runDict, sample_list, config['EXT_INDEX_LIST'], files_list, None, 'validation')

# Log
logfile = f"{resultDir}/{serviceName}.{date_time}.parameters.log"
logging.basicConfig(filename=logfile, level=config['LOG_LEVEL'], format='%(asctime)s %(message)s')

log_items = [
	('Start of the analysis:', date_time),
	('Analysing run:', runName),
	('List of all samples:', sample_list)
]

for item in log_items:
	if isinstance(item[1], list):
		logging.info(f"{item[0]}\n{'\n'.join(map(str, item[1]))}")
	else:
		logging.info(f"{item[0]} {item[1]}")

print(dict(runDict))
################################################## RULES ##################################################
ruleorder: copy_fastq > bcl_convert

rule all:
	input:
		expand(f"{resultDir}/tmp/{{sample}}/alevin/{{sample}}_alevinReport.html", sample=sample_list)


rule help:
	"""
	General help for alevin module
	Launch snakemake -s  alevin.smk -c(numberofthreads) --config DATA_DIR=absolutepathoftherundirectory (default is data) without / at the end of the path
	To launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
	Every variable defined in the yaml file can be change
	Separate multiple variable with a space (ex  --config DATA_DIR=runname transProb=0.05 var1=0.05 var2=12)
	Also use option --configfile another.yaml to replace and merge existing config.yaml file variables
	Use -p to display shell commands
	Use --lt to display docstrings of rules
	Input file = fastq from ScRNAseq (R1 barecode+umi, R2 reads)
	Output file = quantification files
	"""

rule bcl_convert:
	""" Convert BCL files to FASTQ using bcl-convert """
	input:
		bcl_dir = config['BCL_DIRECTORY']
	output:
		fastqR1 = f"{resultDir}/{{sample}}.R1.fastq.gz",
		fastqR2 = f"{resultDir}/{{sample}}.R2.fastq.gz"
	params:
		sample_sheet = config['SAMPLE_SHEET'], 
		output_dir = resultDir,  
		sample_name = lambda wildcards: wildcards.sample
		#config_file = config['BCL_CONVERT_CONFIG']  # Path to bcl-convert config JSON file
	threads: workflow.cores
	shell:
		"""
		bcl-convert --input-dir {input.bcl_dir} \
					--output-dir {params.output_dir} \
					--sample-sheet {params.sample_sheet} \
					--threads {threads}

		R1_fastq=$(find {params.output_dir} -type f -name "{params.sample_name}*_R1_*.fastq.gz" | head -n 1)
		R2_fastq=$(find {params.output_dir} -type f -name "{params.sample_name}*_R2_*.fastq.gz" | head -n 1)
		
		mv $R1_fastq {output.fastqR1}
		mv $R2_fastq {output.fastqR2}
		"""

rule copy_fastq:
	""" Copy input files """
	output:
		fastqR1=temp(f"{resultDir}/{{sample}}.R1.fastq.gz"),
		fastqR2=temp(f"{resultDir}/{{sample}}.R2.fastq.gz")
	params:
		process = config['PROCESS_CMD'],
		download_link1 = lambda wildcards: runDict[wildcards.sample]['.R1.fastq.gz'],
		download_link2 = lambda wildcards: runDict[wildcards.sample]['.R2.fastq.gz']
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link1} {output.fastqR1} && ln -sfn {params.download_link2} {output.fastqR2} || rsync -azvh {params.download_link1} {output.fastqR1} && rsync -azvh {params.download_link2} {output.fastqR2}"

rule alevin:
	input:
		fastqR1=f"{resultDir}/{{sample}}.R1.fastq.gz", # BARCODE+UMI
		fastqR2=f"{resultDir}/{{sample}}.R2.fastq.gz" # READS
	output:
		f"{resultDir}/tmp/{{sample}}.success"
	params:
		index = config['SALMON_INDEX'],
		misc_options = config['MISC_SALMON_OPTIONS'], 
		library_type = config['ALEVIN_LIBRARY'],
		alevindir= f"{resultDir}/tmp/{{sample}}/",
		gene_map = config['ALEVIN_GENE_TABLE']
	threads: workflow.cores
	shell:
		"""
		salmon alevin --dumpFeatures -l {params.library_type} -1 {input.fastqR1} -2 {input.fastqR2} {params.misc_options} -i {params.index} -p {threads} -o {params.alevindir} --tgMap {params.gene_map} && touch {output}
		"""

rule alevinQC:
	input: f"{resultDir}/tmp/{{sample}}/alevin/quants_mat.gz"
	output: f"{resultDir}/tmp/{{sample}}/alevin/{{sample}}_alevinReport.html"
	conda: "rtools"
	shell: """ Rscript /app/scripts/ScRNASeq/AlevinQC.Rscript --input {input} --output {output} --sample {wildcards.sample} """

rule seurat:
	input: f"{resultDir}/tmp/{{sample}}/alevin/quants_mat.gz"
	output: f"{resultDir}/tmp/"
	params:
		prefix=config['PROJECT_PREFIX'],
		min_genes=config['MIN_GENES'],
		max_genes=config['MAX_GENES'],
		max_mito=config['MAX_MITO'],
		min_housekeeping_expr=config['MIN_HOUSEKEEPING_EXPR'],
		remove_ribo=config['REMOVE_RIBO'],
		species=config['SPECIES']
	conda: "rtools"
	shell:
		"""
		Rscript /app/scripts/ScRNASeq/Seurat.Rscript --files {input} --output {output} --prefix {params.prefix} --sample_id {wildcards.sample} --min_genes {params.min_genes} --max_genes {params.max_genes} --max_mito {params.max_mito} --min_housekeeping_expr {params.min_housekeeping_expr} --remove_ribo {params.remove_ribo} --species {params.species}
		"""

onstart:
	shell(f"touch {os.path.join(outputDir, f'{serviceName}Running.txt')}")
	with open(logfile, "a+") as f:
		f.write("\n")
		f.write("Global parameters of the analysis for debug only")
		json.dump(config, f, ensure_ascii=False, indent=2)
		f.write("\n")

onsuccess:
	include = config['INCLUDE_RSYNC']
	shell(f"rm -f {outputDir}/{serviceName}Running.txt")
	shell(f"touch {outputDir}/{serviceName}Complete.txt")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	with open(logfile, "a+") as f:
		f.write(f"End of the analysis : {date_time_end}\n")
	
	# Copy results to the main output directory
	shell("rsync -azvh --include={include} --exclude='*'  {resultDir}/ {outputDir}")

	# Copy individual sample results to their respective directories
	for sample in sample_list:
		shell(f"cp {outputDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/* {outputDir}/{sample}/{serviceName}/ || true")

onerror:
	include_log = config['INCLUDE_LOG_RSYNC']
	shell(f"touch {outputDir}/{serviceName}Failed.txt")
	shell(f"rm -f {outputDir}/{serviceName}Running.txt")
	shell("rsync -azvh --include={include_log} --exclude='*' {resultDir}/ {outputDir}")