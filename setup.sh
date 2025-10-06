#!/usr/bin/env bash
set -euo pipefail

VRNA_URL="https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz"

command -v conda >/dev/null || { echo "Conda not found. Please install conda/mamba."; exit 1; }

# Install tools via conda and pip
conda env update -f environment.yml
conda install -y python=3.12.11
conda install -y -c conda-forge -c bioconda snakemake git
conda install -y -c bioconda bowtie2 sam2pairwise mafft
conda install -y -c conda-forge make pkg-config gsl zlib curl
pip install biopython pandas torch joblib einops scikit-learn gget # now includes gget, but only for the dependencies
pip install -q mysql-connector-python biopython
pip install -q --log log git+https://github.com/pachterlab/gget.git@delphy_dev # this build has support for ncbi_virus
tar -xjvf ./target_seqs/original/HUMAN.tar.bz2 -C ./target_seqs/original

# Build & install ViennaRNA into ./src
set -euo pipefail
cd ./src
curl -L "$VRNA_URL" -o ViennaRNA-2.7.0.tar.gz
tar -xzf ViennaRNA-2.7.0.tar.gz
cd ViennaRNA-2.7.0
./configure --prefix="$(pwd)/viennaRNA" --without-perl # perl is broken on my machine, skip it
make
make install
cd ..
rm ViennaRNA-2.7.0.tar.gz

echo
echo "Done installing dependencies."
