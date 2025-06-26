#!/bin/bash

# Download GEO methylation raw data and metadata
echo "Downloading GEO raw IDAT files for GSE42861..."

# Make and move into data folder
mkdir -p data
cd data

# (You will need to manually download the files from GEO if needed)
# Example placeholder (GEO doesn't always provide direct FTP links anymore):
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/suppl/*.tar

echo "Download done. Unpack manually or with: tar -xvf file.tar"

