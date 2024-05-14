# Install MicroWorldOmics in Linux
Download from links:https://caiyun.139.com/m/i?1G5C1h8CUoxBx; password:y3NR

## Unzip MicroWorldOmics_Linux.zip and run
```
unzip MicroWorldOmics_Linux.zip
cd MicroWorldOmics

conda create -n mwo
source activate mwo
conda install bioconda::raxml
conda install bioconda::aragorn

export QTWEBENGINE_DISABLE_SANDBOX=1
./MicroWorldOmics
```
