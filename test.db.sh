mkdir -p databases/k2
mkdir -p databases/amrfinderplus
pip install pyyaml
wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz -O databases/krakenStandardDB.tgz
tar xzvf databases/krakenStandardDB.tgz -C databases/k2
rm -f databases/krakenStandardDB.tgz