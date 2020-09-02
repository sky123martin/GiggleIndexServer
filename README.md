

mkdir UCSC_utilities
cd UCSC_utilities
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/ ./
cd ..

pip install mysql-connector-python
pip install nose2 #for testing
pip install unittest

conda create --name env
{
git clone https://github.com/sky123martin/GiggleIndexServer
cd GiggleIndexServer
mkdir UCSC_utilities
cd UCSC_utilities
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bigWigToBedGraph ./
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bigBedToBed ./
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bigPslToPsl ./
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/pslToBed ./
cd ..
conda install anaconda
conda create -n env python=3
conda activate env
pip install -r requirements.txt
nose2 --verbose
python setup_indices.py
source deactivate

