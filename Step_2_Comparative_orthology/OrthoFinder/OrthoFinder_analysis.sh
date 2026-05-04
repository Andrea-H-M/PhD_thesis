################################################################################################
## Topic: OrthoFinder analysis for TFF filtered dataset                                       ##
## Author: Olga Andrea Hernandez Miranda (Miranda H)                                          ##
## Date: 09/07/2024                                                                           ##
## Description:                                                                               ##
## This script:                                                                               ##
##  (1) Runs OrthoFinder v2.5.5 on the filtered TFF FASTA dataset.                            ##
##  (2) Uses four species: Arabidopsis, Solanum, Vanilla, and Vitis.                          ##
##  (3) Takes PruebaFiltradoTFF/ as the input directory.                                      ##
##  (4) Generates orthogroups and working files in the OrthoFinder results directory.          ##
##  (5) Reports the expected working directory and orthogroup output file based on the log.    ##
################################################################################################

#!/bin/bash

set -e

BASE_DIR="/home/winter/OrthoFinder"
INPUT_DIR="${BASE_DIR}/PruebaFiltradoTFF"

cd "$BASE_DIR" || exit

echo "Starting OrthoFinder analysis..."
echo "Input directory: $INPUT_DIR"

echo "Species used:"
echo "0: Arabidopsis.fa"
echo "1: Solanum.fa"
echo "2: Vanilla.fa"
echo "3: Vitis.fa"

orthofinder.py -f "$INPUT_DIR"

echo "OrthoFinder run completed."

# 🔹 Detectar automáticamente la carpeta de resultados más reciente
RESULTS_DIR=$(ls -td ${INPUT_DIR}/OrthoFinder/Results_* | head -1)
WORKING_DIR="${RESULTS_DIR}/WorkingDirectory"

echo "Working directory:"
echo "$WORKING_DIR"

echo "Orthogroups file:"
echo "${WORKING_DIR}/clusters_OrthoFinder_I1.5.txt_id_pairs.txt"