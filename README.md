Must have Giggle and bgzip (https://vcf.iobio.io/help.html) setup on local machine.

conda create --name env

git clone https://github.com/sky123martin/GiggleIndexServer

cd GiggleIndexServer

git clone https://github.com/ryanlayer/giggle.git

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

conda install --file requirements.txt

nose2 --verbose

python3 setup_indices.py

source deactivate
