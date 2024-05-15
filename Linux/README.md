# Install MicroWorldOmics in Linux
Download from links:https://caiyun.139.com/m/i?1G5C1h8CUoxBx; password:y3NR

## Install MicroWorldOmics_Linux.zip and run
```
# Setup env
conda create -n microworldomics python==3.6
source activate microworldomics
conda install bioconda::aragorn
conda install -c bioconda clinker-py
conda install bioconda::raxml
conda install bioconda::clustalo
conda install bioconda::diamond==0.9.14
conda install bioconda::fastani
conda install bioconda::fasttree
conda install bioconda::iqtree
conda install bioconda::mafft
conda install bioconda::muscle
conda install bioconda::prodigal


# Unzip MicroWorldOmics_Linux.zip and run
unzip MicroWorldOmics_Linux.zip
cd MicroWorldOmics

export QTWEBENGINE_DISABLE_SANDBOX=1
./MicroWorldOmics
```
