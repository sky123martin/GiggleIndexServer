class config:

    ####################################
    # DEVELOPER CONFIGURABLE VARIABLES #
    ####################################

    # TESTING MODE
    TESTING = False

    # DB file name
    DB = 'Indexing.db'

    # Acceptable file formats. Files not following this format will be ignored
    UCSC_ACCEPTABLE_FILE_FORMATS = ["bw", "bb", "bed", "bigPsl"]
    LOCAL_ACCEPTABLE_FILE_FORMATS = ["bed", "bed.gz"]

    # UCSC Genomes to download
    UCSC_BIG_DATA_LINK = "http://hgdownload.soe.ucsc.edu"
    UCSC_SQL_DB_HOST = "genome-mysql.soe.ucsc.edu"
    UCSC_SQL_DB_USER = "genome"
    UCSC_API = "https://api.genome.ucsc.edu/"

    ###############################
    # USER CONFIGURABLE VARIABLES #
    ###############################

    # Max Number of Intervals per Index
    MAX_INTERVALS_PER_INDEX = 1000000

    # Number of proccesses available to server
    AVAILABLE_PROCCESSES = 4

    # Genome from UCSC to download and index
    UCSC_GENOMES = []#["rn6"]

    # Local genomes to download and index 
    # Format: {<genome name>: [<path to data>, <path metadata file name>],...
    # Example  {"lab data" : ["local/", "local/metadata.csv"]}{"localgenome": ["local/", "", "rn6"]}
    LOCAL_GENOMES = [{"project_name": "test",
                      "reference_genome": "hg19",
                      "metadata_path": "",
                      "data_path": "roadmap_sort/"
                      }]

    # Metadata file
    # Format: .csv format with columns = ["file_name", "short_name", "long_name", "short_info", "long_info"]
    # Example  "metadata.csv"

