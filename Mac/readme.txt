# 创建虚拟环境激活
~/miniconda3/bin/conda create -n mwo python==3.6
source ~/miniconda3/bin/activate mwo


# Mac 电脑安装pyqt5及其插件 (若报错,卸载再次安装)
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5 --trusted-host pypi.tuna.tsinghua.edu.cn

pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5-tools --trusted-host pypi.tuna.tsinghua.edu.cn

pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQtWebEngine --trusted-host pypi.tuna.tsinghua.edu.cn


# 检测缺失的依赖库
## 创建一个简单的pyqt5应用,test.py用于检测缺失的依赖库
python test.py


# 安装python依赖模块
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

### 版本一定要对
conda install -c anaconda keras==2.6.0

pip install -i https://pypi.tuna.tsinghua.edu.cn/simple tensorflow==2.6.0 --trusted-host pypi.tuna.tsinghua.edu.cn

### 如果有报错: 'TypeError: Descriptors cannot be created directly.' 则: (optional) 
pip uninstall protobuf
pip install protobuf==3.20.0

## 安装 opencv
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple opencv_python==4.1.2.30 --trusted-host pypi.tuna.tsinghua.edu.cn

# 安装 R
~/miniconda3/bin/conda create -n R
source ~/miniconda3/bin/activate R
conda install -c conda-forge r-base=4.3.2

## 安装 R 包
### 利用conda安装的包
source ~/miniconda3/bin/activate R
conda install bioconda::bioconductor-biobase
conda install conda-forge::r-devtools
conda install bioconda::bioconductor-phyloseq
conda install bioconda::bioconductor-ggtree
conda install conda-forge::r-nloptr
conda install bioconda::bioconductor-rsamtools
conda install conda-forge::r-phytools
conda install conda-forge::r-lsa


### 利用 R base 安装的包
cd /Users/mac/miniconda3/envs/R/lib/R/bin
./R

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
install.packages('igraph')
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

##### 安装 PLSDAbatch 包及其依赖
conda install conda-forge::r-ggpubr
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch")


##### 安装 WGCNA 包和 NetMoss2 包及其依赖
BiocManager::install("systemfonts")
BiocManager::install("WGCNA")
remotes::install_github("xiaolw95/NetMoss2")


##### 安装 sf 包及其依赖
##### 此时 sf 包和其他包可能会产生冲突,因此重新配置环境
~/miniconda3/bin/conda create -n R_opt
source ~/miniconda3/bin/activate R_opt
conda install -c conda-forge r-base=4.0.3

conda install conda-forge::r-sf
conda install conda-forge::r-tidyverse
conda install conda-forge::r-ragg
conda install conda-forge::r-sna
conda install conda-forge::r-peptides


### 利用 R base 安装的包
cd /Users/mac/miniconda3/envs/R_opt/lib/R/bin
./R

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


## 安装 ggClusterNet 及其依赖包
conda install bioconda::bioconductor-phyloseq
conda install conda-forge::r-ggpubr
remotes::install_github("taowenmicro/ggClusterNet")
