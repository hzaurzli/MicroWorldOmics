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

pip install opencv_python==4.1.2.30
## 安装 opencv 可能缺乏的依赖 (optional) 
### 缺少 libEGL.so.1
sudo apt-get -y install libegl1

### 缺少 libOpenGL.so.0
sudo apt-get install libgl1-mesa-dev




# 启动前运行
export QTWEBENGINE_DISABLE_SANDBOX=1



