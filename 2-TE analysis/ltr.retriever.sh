#!/usr/bin/bash

# Help message
usage() {
    echo "Usage: $0 [species] [genome] [genomename]"
    echo
    echo "Parameters:"
    echo "  species     : Species name (arbitrary input)"
    echo "  genome      : Absolute path to genome file"
    echo "  genomename  : Genome file name"
    echo
    echo "Example:"
    echo "  $0 homo_sapiens /path/to/genome.fasta genome"
    echo
    exit 1
}

# Show help if -h or --help is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
fi

# Parameter assignment
species=$1
genome=$2
genomename=$3

# Check required parameters
if [ -z "$species" ] || [ -z "$genome" ] || [ -z "$genomename" ]; then
    echo "Error: Missing required parameters."
    usage
fi

# Create symbolic link to genome
ln -s $genome ./${genomename}

# Build genome index
gt suffixerator -db $genomename -indexname $genomename -tis -suf -lcp -des -ssp -sds -dna > suff.log 2>&1

# Create output directories
mkdir 02_harvest

# Run LTR harvest
gt ltrharvest -index $genomename -minlenltr 100 -maxlenltr 20000 -seqids yes > 02_harvest/${species}.harvest.scn

# Prepare finder directory
mkdir 03_finder
cp $genome 03_finder/
cd 03_finder/

# Run LTR finder
ltr_finder -C -M 0.8 ../${genomename} > ${species}.finder.scn


cd ../

# Prepare retriever directory
mkdir 04_retriever
cd 04_retriever/

# Run LTR retriever with multi-threading
LTR_retriever -genome $genome -inharvest ../02_harvest/${species}.harvest.scn -infinder ../03_finder/${species}.finder.scn -threads 40 > retriever.log 2>&1