# 创建虚拟环境激活
~/miniconda3/bin/conda create -n mwo python==3.6
source ~/miniconda3/bin/activate mwo

#安装pyqt5及其插件 (若报错,卸载再次安装)
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5-tools
pip install PyQtWebEngine


# 检测缺失的依赖库
## 创建一个简单的pyqt5应用,test.py用于检测缺失的依赖库
export QT_DEBUG_PLUGINS=1
python test.py


## 安装依赖库 (optional) 
### 进入pyqt5共享库的文件夹
cd /root/miniconda3/envs/mwo/plugins/platforms/
sudo apt update

### 缺少 libxcb-icccm.so.4
sudo apt-get -y install libxcb-icccm4

### 缺少 libxcb-image.so.0
sudo apt-get -y install libxcb-image0

### 缺少 libxcb-render-util.so.0
sudo apt-get -y install libxcb-render-util0

### 缺少 libxcb-keysyms.so.1
sudo apt-get -y install libxcb-keysyms1

### 缺少 libGL.so.1
sudo apt install libgl1-mesa-glx

### 缺少 libXrender.so.1
sudo apt-get install libxrender1

### 缺少 libXcomposite.so.1
sudo apt-get install libxcomposite-dev

### 缺少 libXcursor.so.1
sudo apt-get install libxcursor-dev

### 缺少 libXdamage.so.1
sudo apt-get install libxdamage1

### 缺少 libXi.so.6
sudo apt-get install libx11-dev libxext-dev libxtst-dev libxrender-dev libxmu-dev libxmuu-dev

### 缺少 libXss.so.1
sudo apt-get install libxss1

### 缺少 libxcb-shape.so.0
sudo apt-get -y install libxcb-shape0

### 缺少 libxcb-xkb.so.1
sudo apt-get -y install libxcb-xkb1

### 缺少 libxkbcommon-x11.so.0
sudo apt -y install libxkbcommon-x11-0

### 缺少 libXrandr.so.2
sudo apt-get install libxrandr2



# 安装python依赖模块
pip install biopython
conda install aragorn
conda install numpy
pip install matplotlib -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install pandas
pip install PyQt5-stubs
pip install paramiko
pip install psutil
pip install requests
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple scikit-learn
conda install -c bioconda fastani
conda install raxml
pip install networkx
pip install markov_clustering
conda install pytorch torchvision torchaudio cpuonly -c pytorch
pip install datasets
pip install transformers
### 版本一定要对
conda install -c anaconda keras==2.6.0
pip install tensorflow==2.6.0
### 如果有报错: 'TypeError: Descriptors cannot be created directly.' 则: (optional) 
pip uninstall protobuf
pip install protobuf==3.20.0

## 安装 opencv
pip install opencv_python==4.1.2.30
## 安装 opencv 可能缺乏的依赖 (optional) 
### 缺少 libEGL.so.1
sudo apt-get -y install libegl1

### 缺少 libOpenGL.so.0
sudo apt-get install libgl1-mesa-dev


# 启动前运行
export QTWEBENGINE_DISABLE_SANDBOX=1


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

### 利用 R base 安装的包
cd /root/miniconda3/envs/R/lib/R/bin
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
BiocManager::install("Maaslin2")
install.packages('flashClust')
install.packages('igraph')
install.packages('adegenet')
BiocManager::install("GenomicAlignments")
install.packages('randomcoloR')
BiocManager::install("igvShiny")
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
remotes::install_github("MPBA/r-sparcc")
install.packages('rhierbaps')
devtools::install_github('zdk123/SpiecEasi')
install.packages('shinyBS')
install.packages('metricsgraphics')
install.packages('pls')
devtools::install_github("kevqyzhu/TMscoreAlign")


##### 安装 PLSDAbatch 包及其依赖
sudo apt install cmake
sudo apt-get install libnlopt-dev
install.packages("nloptr")
install.packages("lme4")
install.packages("ggpubr")
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch")


##### 安装 WGCNA 包及其依赖
sudo apt-get install libfontconfig1-dev
BiocManager::install("systemfonts")
BiocManager::install("WGCNA")


##### 安装 NetMoss2 包及其依赖(用conda先安装pROC和ggforce,然后移动)
source ~/miniconda3/bin/activate R_opt
conda install conda-forge::r-proc
conda install conda-forge::r-ggforce
cp -r /root/miniconda3/envs/R_opt/lib/R/library/pROC/ /root/miniconda3/envs/R/lib/R/library/
cp -r /root/miniconda3/envs/R_opt/lib/R/library/ggforce/ /root/miniconda3/envs/R/lib/R/library/
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
cd /root/miniconda3/envs/R_opt/lib/R/bin
./R
install.packages('shiny')
install.packages('dplyr')
install.packages('DT')
install.packages('ggrepel')
install.packages('fastmap')
install.packages('shinycssloaders')
install.packages('shinythemes')
install.packages('ape')
install.packages('BiocManager')
install.packages('shinyFiles')


## 安装 ggClusterNet 及其依赖包
## 由于在 R 中利用 conda 安装了 phyloseq, 将 R 的 phyloseq, tidyfst, Biostrings, BiocGenerics, S4Vectors, XVector, GenomeInfoDb, RCurl, bitops, GenomeInfoDbData 复制到 R_sf 中
cp -r /root/miniconda3/envs/R/lib/R/library/tidyfst/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/Biostrings/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/phyloseq/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/BiocGenerics/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/S4Vectors/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/XVector/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/GenomeInfoDb/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/RCurl/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/bitops/ /root/miniconda3/envs/R_opt/lib/R/library/
cp -r /root/miniconda3/envs/R/lib/R/library/GenomeInfoDbData/ /root/miniconda3/envs/R_opt/lib/R/library/

remotes::install_github("taowenmicro/ggClusterNet")


