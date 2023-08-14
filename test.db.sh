mkdir -p databases/k2
pip install pyyaml
wget -c https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz -O databases/k2/minikraken2_v2_8GB_201904.tgz
tar xzvf databases/k2/minikraken2_v2_8GB_201904.tgz -C databases/k2
rm -f databases/k2/minikraken2_v2_8GB_201904.tgz