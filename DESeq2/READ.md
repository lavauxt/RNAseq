# To build and launch the docker
docker compose -f deseq2.docker-compose.yml run --name=deseq2-cli deseq2-cli

# Or
docker build -f deseq2_dockerfile --tag deseq2:1.0 .
docker run -itd --name=deseq2 --entrypoint="/bin/bash" -v /home1/data/STARK/databases:/STARK/databases:ro -v /home1/BAS/lavauxt/input:/STARK/input  -v /home1/BAS/lavauxt/output:/STARK/output deseq2:1.0
docker exec -it deseq2 bash

# To launch the analysis
Rscript Count.R --input /path/to/counts --output /path/to/result/folder  --sampletable /path/to/sample_table.csv --compare condition --level2compare nameofthelevelwithtreatment --baselevel nameofthebaselevel

to specify databases --ensembl /path/to/database.tar.gz

# recipe for Salmon index
# GRch38 assembly and gencode.vxx.transcripts.fa.gz can be found on https://www.gencodegenes.org/human/
# release can be change (actually v41 is used)
# for GRch38 primary assembly v42 : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
# for gencode transcripts : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.v41.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode

# Update EnsDb Homo Sapiens
# if needed, install required libraries
BiocManager::install(c("AnnotationHub", "ensembldb", "AnnotationForge"))

library(AnnotationHub)
library(ensembldb)
library(AnnotationForge)
library(devtools)

# start an AnnotationHub instance/connection.
ah <- AnnotationHub()

# query for available Ensembl databases (ex for v107)
EnsDb.Hsapiens.v107 <- query(ah, c("EnsDb", "Homo Sapiens"))
EnsDb.Hsapiens.v107

# Fetch the v107 EnsDb and put it in the cache
EnsDb.Hsapiens.v107 <- EnsDb.Hsapiens.v107[["AH104864"]]

# check (optionnal)
columns(EnsDb.Hsapiens.v107)

# Now copy, and install the database locally.
# By doing so, you can just quickly load the library next time
# see: ?makeEnsembldbPackage
# set working dir to store the file
setwd("./database")

# Copy databse from the cache to working dir
file.copy(AnnotationHub::cache(ah["AH104864"]), "./EnsDb.Hsapiens.sqlite")

# now make it a package. Change name and email accordingly
makeEnsembldbPackage("EnsDb.Hsapiens.sqlite", version="1.0.0", maintainer = "Th LAVAUX <lavauxt@gmail.com>", author = "Th LAVAUX <lavauxt@gmail.com>", destDir=".", license="MIT")
install.packages("./EnsDb.Hsapiens.v107", type = "source", repos = NULL)

# create a source package (TAR.GZ)
build("EnsDb.Hsapiens.v107", binary = TRUE)

# to install that package in R
#install.packages("./EnsDb.Hsapiens.v107.tar.gz", type = "source", repos = NULL)
# load the lib
#library(EnsDb.Hsapiens.v107)
