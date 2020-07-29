

mkdir UCSC_utilities
cd UCSC_utilities
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/ ./
cd ..

pip install mysql-connector-python
pip install nose2 #for testing
pip install unittest