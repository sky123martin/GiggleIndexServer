# Setup
## Setup: Environment

Must have Giggle, pip, conda, and bgzip (https://vcf.iobio.io/help.html) setup on local machine.

Clone directory
```
git clone https://github.com/sky123martin/GiggleIndexServer

cd GiggleIndexServer
```
Clone giggle to use gzip and sorting function
```
git clone https://github.com/ryanlayer/giggle.git
```
Install UCSC utilities for file conversion
```
mkdir UCSC_utilities

cd UCSC_utilities

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bigWigToBedGraph ./

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bigBedToBed ./

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/bigPslToPsl ./

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/pslToBed ./

cd ..
```
Create enviorment
```
conda install anaconda

conda create -n env python=3

conda activate env
```

Install necessary python libararies
```
pip3 install -r requirements.txt 
```
## Setup: Configuration Variables
The file config.py is where you input the genomes you want to index and clustering params.
```
     ###############################
    # USER CONFIGURABLE VARIABLES #
    ###############################
    
    # Max Number of Intervals per Index
    MAX_INTERVALS_PER_INDEX = 1000000

    # Number of proccesses available to server
    AVAILABLE_PROCCESSES = 4

    # Timeout on a file download
    timeout_file_download = 1000
    timeout_file_processing = 100

    # List genomes from UCSC to download and index
    '''
    example:
        CSC_GENOMES = ["hg19", "rn6"]
    '''
    UCSC_GENOMES = []

    # List local genomes to download and index 
    # Format: {<genome name>: [<path to data>, <path metadata file name>],...
    # Example  {"lab data" : ["local/", "local/metadata.csv"]}{"localgenome": ["local/", "", "rn6"]}
    '''
    example:
        LOCAL_GENOMES = [
                         {
                         "project_name": "test",
                         "reference_genome": "hg19",
                         "metadata_path": "",
                         "data_path": "roadmap_sort/"
                         },
                        ...]
    '''
    LOCAL_GENOMES = []


    # Metadata file for local directories
    # Format: .csv format with columns = ["file_name", "short_name", "long_name", "short_info", "long_info"]
    # Example  "metadata.csv"

    # Uncomment UCSC public hubs to index, see https://api.genome.ucsc.edu/list/publicHubs for further info on each hub 
    '''
    example:
        UCSC_HUBS = [
            # "ALFA Hub",
            "Bird assemblies",
            # "Blueprint Hub",
            ...
            ""
            ]
    '''

    UCSC_HUBS = [
                # "ALFA Hub",
                # "Bird assemblies",
                # "Blueprint Hub",
                # "BrainEpigenomeHub",
                ...
                ""
                ]

    
```

## Setup: Test
Errors in config test functions mean that you may have entered invalid config variables 
```
nose2 --verbose
```


# Index Files

```
python3 setup_indices.py
```

# Update Indices
```
python3 update_indices.py
```

# Search Indices

```
python3 query_indices.py -- help

  -h, --help   Show the help
  -i           "returns all currently hosted indices
  -f           "returns all currently hosted files
  -g           "returns all currently hosted genomes
  -s           "returns all currently hosted sources
  --qf         Query a given file given a path. ie: <path> <optional param: source> <optional param: genome>
  --qi         Query a given interval in format 'Chr:#-#' ie. <optional param: source> <optional param: genome>
```

## Database querying functions
Check what indicies are in server:
```
python3 query_indices.py -i
```
Check what files are in server:
```
python3 query_indices.py -f
```
Check what genomes are in server:
```
python3 query_indices.py -g
```
Check what sources are in server:
```
python3 query_indices.py -s
```
## Index querying functions
### Interval

Query an interval on a specfic source
```
python3 query_indices.py --qi 1:10000-20000 rn6
```

Query an interval on a specfic source and genome
```
python3 query_indices.py --qi 1:10000-20000 rn6 UCSC
```
### File

Query a file on a specfic source
```
python3 query_indices.py --qf local/testbed.bed.gz rn6
```

Query a file on a specfic source and genome
```
python3 query_indices.py --qf local/testbed.bed.gz rn6 UCSC
```

## Data
The following that can be used to setup local repositories:

Roadmap Epigenomics: https://s3.amazonaws.com/layerlab/giggle/roadmap/roadmap_sort.tar.gz

Fantom5: https://s3.amazonaws.com/layerlab/giggle/fantom/fantom_sort.tar.gz



### Making Fantom5 a local repo:
Download relavent data into GiggleIndexServer/ directory
```
cd GiggleIndexServer

curl https://s3.amazonaws.com/layerlab/giggle/fantom/fantom_sort.tar.gz

gzip -d < fantom_sort.tar.gz | tar xvf -

rm fantom_sort.tar.gz
```
Add the following entry to the LOCAL_GENOMES list in config.py

```
  {
  "project_name": "fantom5",
  "reference_genome": "hg19",
  "metadata_path": "",
  "data_path": "fantom_sort/"
  }
```