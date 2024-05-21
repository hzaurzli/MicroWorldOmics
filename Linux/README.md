# Install MicroWorldOmics in Linux
Download from links:https://caiyun.139.com/m/i?1G5C1h8CUoxBx; password:y3NR

## Install MicroWorldOmics_Linux.zip and run
```
# Setup R env
~/miniconda3/bin/conda create -n R
source ~/miniconda3/bin/activate R
conda install r-base=4.3.2

## R packages install
### By using conda
source ~/miniconda3/bin/activate R
conda install bioconda::bioconductor-biobase
conda install conda-forge::r-devtools
conda install bioconda::bioconductor-phyloseq
conda install bioconda::bioconductor-ggtree
conda install conda-forge::r-nloptr
conda install bioconda::bioconductor-rsamtools
conda install conda-forge::r-phytools
conda install conda-forge::r-lsa
conda install conda-forge::r-ggpubr


### Enter R base
cd ~/miniconda3/envs/R/lib/R/bin
./R

### By using R base
install.packages("BiocManager")
install.packages('shiny')
install.packages('ggplot2')
install.packages('RColorBrewer')
install.packages('r3dmol')
install.packages('colourpicker')
devtools::install_github("xavierdidelot/BactDating")
install.packages('ape')
BiocManager::install('mixOmics')
install.packages('gridExtra')
install.packages('robustbase')
install.packages('biglm')
install.packages('pcaPP')
install.packages('optparse')
install.packages("lme4")
BiocManager::install("Maaslin2")
install.packages('flashClust')
install.packages('adegenet')
BiocManager::install("GenomicAlignments")
install.packages('randomcoloR')
devtools::install_github('gladkia/igvShiny',INSTALL_opts = '--no-lock')
install.packages('htmlwidgets')
install.packages('dplyr')
devtools::install_github('hzaurzli/MicroMCA',INSTALL_opts = '--no-lock')
install.packages('DT')
install.packages('ggrepel')
install.packages('tidyr')
install.packages('shinycssloaders')
install.packages('shinythemes')
install.packages('tidyfst')
install.packages('sna')
install.packages('vegan')
install.packages('shinyFiles')
remotes::install_github("MPBA/r-sparcc")
install.packages('rhierbaps')
devtools::install_github('zdk123/SpiecEasi')
install.packages('shinyBS')
install.packages('metricsgraphics')
install.packages('pls')
devtools::install_github("kevqyzhu/TMscoreAlign")
remotes::install_github("xiaolw95/NetMoss2")
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch")
BiocManager::install("systemfonts")
BiocManager::install("WGCNA")
remotes::install_github("xiaolw95/NetMoss2")
install.packages('igraph')


# Setup R_opt env
~/miniconda3/bin/conda create -n R_opt
source ~/miniconda3/bin/activate R_opt
conda install r-base=4.1.3

## R packages install
### By using conda
conda install conda-forge::r-sf
conda install conda-forge::r-tidyverse
conda install conda-forge::r-ragg
conda install conda-forge::r-sna
conda install conda-forge::r-peptides
conda install bioconda::bioconductor-phyloseq
conda install conda-forge::r-ggpubr


### Enter R base
cd ~/miniconda3/envs/R_opt/lib/R/bin
./R

### By using R base
install.packages('shiny')
install.packages('dplyr')
install.packages('DT')
install.packages('ggrepel')
install.packages('fastmap')
install.packages('shinycssloaders')
install.packages('shinythemes')
install.packages('ape')
install.packages('vegan')
install.packages('remotes')
install.packages('shinyFiles')
install.packages('shinydashboard')
install.packages('waiter')
install.packages('colourpicker')
install.packages('devtools')
require(devtools)
install_version("tidyfst", version = "1.6.5", repos = "http://cran.us.r-project.org")
remotes::install_github("taowenmicro/ggClusterNet")


# Unzip MicroWorldOmics_Linux.zip and run
cp -r ~/miniconda3/envs/R/lib/R/* /.../MicroWorldOmics/Shiny/R/
cp -r ~/miniconda3/envs/R_opt/lib/R/* /.../MicroWorldOmics/Shiny/R_opt/
unzip MicroWorldOmics_Linux.zip
cd MicroWorldOmics

python ./MicroWorldOmics.py
```
