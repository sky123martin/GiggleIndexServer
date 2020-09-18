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
    UCSC_BIG_DATA_LINK = ""
    UCSC_SQL_DB_HOST = "genome-mysql.soe.ucsc.edu"
    UCSC_SQL_DB_USER = "genome"
    UCSC_API = "https://api.genome.ucsc.edu/"

    ###############################
    # USER CONFIGURABLE VARIABLES #
    ###############################

    # Max Number of Intervals per Index
    MAX_INTERVALS_PER_INDEX = 10000

    # Number of proccesses available to server
    AVAILABLE_PROCCESSES = 4

    # Genome from UCSC to download and index
    UCSC_GENOMES = []

    # Local genomes to download and index 
    # Format: {<genome name>: <path to data>,...
    # Example  {"lab data" : "local/"}
    LOCAL_GENOMES = {"localgenome": "local"}
