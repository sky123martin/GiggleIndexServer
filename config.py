class config:
    TESTING = True

    # Max Number of Intervals per Index
    MAX_INTERVALS_PER_INDEX = 100000000

    DB = 'Indexing.db'
    # UCSC Genomes to download
    UCSC_GENOMES = ["rn6"]
    UCSC_SQL_DB_HOST = "genome-mysql.soe.ucsc.edu"
    UCSC_SQL_DB_USER = "genome"
    UCSC_API = "https://api.genome.ucsc.edu/"

    # Local genomes to download {genome name: path to data, "rn6" : "local/"}
    LOCAL_GENOMES = {"localgenome": "local"}

