#!/bin/bash

# based on: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05606-4

# Set default values if parameters are not provided
fragments=${1:-/data/bagel5/db_proteins/bagel5_fragments.faa}
path=${2:-/data/www/users/vicky/bagel5_annotation/classI/Paenibacillus/NZ_CP024795.1}
basename=${3:-NZ_CP024795.1__AOI}
cpu=${4:-6}
mismatch=${5:-1}
kmer=${6:-5}

#pepmatch-preprocess [-h] -p PROTEOME [-n PROTEOME_NAME] -k KMER_SIZE -f PREPROCESS_FORMAT [-P PREPROCESSED_FILES_PATH] [-g GENE_PRIORITY_PROTEOME
cd $path
pepmatch-preprocess -p $basename.all.faa -k $kmer -f pickle
pepmatch-match -q $fragments -p $basename.all.faa -n $cpu -m $mismatch -k $kmer -f csv -o $basename.pepmatch_all

awk -F',' '$2 != ""' $basename.pepmatch_all.csv > $basename.pepmatch.csv
awk -F',' '{print $1"\t"$2"\t"$3"\t"$8"\t"$10"\t"$11}' $basename.pepmatch.csv > $basename.pepmatch.tab

#rm $basename*.pkl

cat $basename.pepmatch.tab

