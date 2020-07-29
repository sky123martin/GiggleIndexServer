import os

class config(object):

    # Max Number of Intervals per Index
    MAX_INTERVALS_PER_INDEX = 100_000_000

    # Genomes to download
    UCSC_GENOMES = ["rn6"]

    # Check if genomes changed or added new files ever _ minutes
    UPDATE_INDEX_INTERVAL = 10
