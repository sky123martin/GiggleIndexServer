import os

class config(object):

    # Max Number of Intervals per Index
    MAX_INTERVALS_PER_INDEX = 100_000_0

    # Genomes to download
    UCSC_GENOMES = ["rn6"]
    UCSC_SQL_DB_HOST = "genome-mysql.soe.ucsc.edu"
    UCSC_SQL_DB_USER = "genome"
    UCSC_API = "https://api.genome.ucsc.edu"

    # Check if genomes changed or added new files ever _ minutes
    UPDATE_INDEX_INTERVAL = 10
