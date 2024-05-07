# 创建虚拟环境激活
~/miniconda3/bin/conda -n mwo python==3.6
source ~/miniconda3/bin/activate mwo

#安装pyqt5及其插件
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PyQt5-tools

#或者(若缺少共享库,一次性安装依赖)
conda install pyqt 
conda install conda-forge::pyqtwebengine


# 检测缺失的依赖库
## 创建一个简单的pyqt5应用,test.py用于检测缺失的依赖库
export QT_DEBUG_PLUGINS=1
python test.py


## 安装依赖库 (optional) 
### 进入pyqt5共享库的文件夹
cd /root/miniconda3/envs/mwo/plugins/platforms/
sudo apt update

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
sudo apt-get install libXScrnSaver


# 安装python依赖模块
pip install biopython
conda install bioconda::aragorn
pip install PyQt5-stubs
pip install biopython
pip install paramiko
pip install psutil
pip install requests
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple opencv-python
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple scikit-learn
conda install -c bioconda fastani
conda install raxml
pip install networkx
pip install markov_clustering


# 启动前运行
export QTWEBENGINE_DISABLE_SANDBOX=1