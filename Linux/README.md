# Install MicroWorldOmics in Linux
Download from links:https://caiyun.139.com/m/i?1G5C1h8CUoxBx; password:y3NR

## Install MicroWorldOmics_Linux.zip and run
```
# Setup Python env
~/miniconda3/bin/conda create -n mwo python==3.6
source ~/miniconda3/bin/activate mwo

pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5-tools
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQtWebEngine --trusted-host pypi.tuna.tsinghua.edu.cn
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple biopython --trusted-host pypi.tuna.tsinghua.edu.cn
conda install bioconda::aragorn
conda install numpy
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple matplotlib --trusted-host pypi.tuna.tsinghua.edu.cn
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple pandas --trusted-host pypi.tuna.tsinghua.edu.cn
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5-stubs --trusted-host pypi.tuna.tsinghua.edu.cn
conda install anaconda::paramiko
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple psutil --trusted-host pypi.tuna.tsinghua.edu.cn
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple requests --trusted-host pypi.tuna.tsinghua.edu.cn
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple scikit-learn --trusted-host pypi.tuna.tsinghua.edu.cn
conda install -c bioconda fastani
conda install raxml
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple networkx --trusted-host pypi.tuna.tsinghua.edu.cn
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple markov_clustering --trusted-host pypi.tuna.tsinghua.edu.cn
conda install pytorch torchvision torchaudio cpuonly -c pytorch
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple datasets --trusted-host pypi.tuna.tsinghua.edu.cn
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple transformers==4.11.3 --trusted-host pypi.tuna.tsinghua.edu.cn
conda install -c anaconda keras==2.6.0
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple tensorflow==2.6.0 --trusted-host pypi.tuna.tsinghua.edu.cn
pip uninstall protobuf
pip install protobuf==3.20.0
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple opencv_python==4.1.2.30 --trusted-host pypi.tuna.tsinghua.edu.cn


## Install tools
conda install numpy
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
conda install anaconda::paramiko


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

source ~/miniconda3/bin/activate mwo
export QTWEBENGINE_DISABLE_SANDBOX=1
python ./MicroWorldOmics.py
```
