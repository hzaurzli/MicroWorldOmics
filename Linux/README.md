# Install MicroWorldOmics in Linux
Download from links:https://caiyun.139.com/m/i?1G5C1h8CUoxBx; password:y3NR

## Install MicroWorldOmics_Linux.zip and run
```
# Setup R env
unzip Rpackages

~/miniconda3/bin/conda create -n R
source ~/miniconda3/bin/activate R
conda install r-base=4.3.2
cp -r /PATH/Rpackages/R/* ~/miniconda3/envs/R/lib/R/library/
cp -r ~/miniconda3/envs/R/lib/R/* /PATH/MicroWorldOmics/Shiny/R/

~/miniconda3/bin/conda create -n R_opt
source ~/miniconda3/bin/activate R_opt
conda install r-base=4.0.3
cp -r /PATH/Rpackages/R_opt/* ~/miniconda3/envs/R_opt/lib/R/library/
cp -r ~/miniconda3/envs/R_opt/lib/R/* /PATH/MicroWorldOmics/Shiny/R_opt/


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

source activate microworldomics
export QTWEBENGINE_DISABLE_SANDBOX=1
python ./MicroWorldOmics.py
```
