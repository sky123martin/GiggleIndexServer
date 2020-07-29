import mysql.connector
import pandas as pd
from contextlib import contextmanager
from application import app
import requests 
import os



class CustomException(Exception):
    pass


@contextmanager
def connect_SQL_db(host, user):
    try:
        db = mysql.connector.connect(
                host=host,
                user=user
                )
        yield db
    except mysql.connector.Error as e:
        raise e
    else:
        db.close()

def target_columns(columns, file_type): 
    """ Returns columns needed to create a BED file.
    Based on genome.ucsc.edu/FAQ/FAQformat.html

    Parameters
    ----------
    columns: list
        column names in file

    Returns
    -------
    target_colunms: list
        columns that need to be extracted to make bed file and base shift
    """
    chrom = ""
    chrom_start = ""
    chrom_end = ""
    base_shift = ""

    # determine chromosome
    if "chrom" in columns:
        chrom = "chrom"

    # determine chromosome base pair start position
    if "chromStart" in columns and "chromEnd" in columns:  # (.bed, pgSNP, nPk, bPk, gappedPeak, tagAlign)
        chrom_start = "chromStart"  # Chromosone start position
        chrom_end = "chromEnd"
        base_shift = 0  # First base is 0
    elif "tStart" in columns and "tEnd" in columns:  # (.psl)
        chrom_start = "tStart"  # Alignment start position in target (.psl)
        chrom_end = "tEnd"
        chrom = "tName"
        base_shift = 0  # First base is 0
    elif "start" in columns and "end" in columns:  # (GFF, MAF)
        chrom_start = "start"  # Starting position in feature
        chrom_end = "end"
        if file_type == "GFF":  # MAF starts 0 and GFF starts 1
            base_shift = -1
        else:
            base_shift = 0
    elif "txStart" in columns and "txEnd" in columns:  # (.genePred)
        chrom_start = "txStart"  # Transcription start position
        chrom_end = "txEnd"
        base_shift = 0  # First base is 0
    elif "genoName" in columns and "genoStart" in columns and "genoEnd" in columns: # rmsk
        chrom_start = "genoStart"  # Repeating Elements by RepeatMasker
        chrom_end = "genoEnd"
        chrom = "genoName"
        base_shift = 0  # First base is 0        

    # if accurate column values not found
    if chrom == "" or chrom_start == "" or chrom_end == "" or base_shift == "":
        print("Given file columns does not have a conversion to .bed for file type {}".format(file_type))
        print(columns)
        print([chrom, chrom_start, chrom_end, base_shift])
        return []

    return [chrom, chrom_start, chrom_end, base_shift]


def cluster_data(source, genome, files_info):
    """ Clusters files to be put into a database
        Parameters
    ----------
    files_info: nested dict
        {file_name : {size, download_function, params}, ...}

    Returns
    -------
    indicies: dict with values of arrays
        {index name : {file_name : {size, download_function, params}}}

    """
    clusters = {}
    max_size = app.config["MAX_INTERVALS_PER_INDEX"]
    num_indices = 1
    curr_index = source + "_" + genome + "." + str(num_indices)
    clusters[curr_index] = {"files": {}}
    curr_size = 0
    while len(files_info) != 0:
        for file_name in files_info:
            if curr_size > max_size and len(files_info) != 0:
                # reset size and name
                clusters[curr_index]["index_size"] = curr_size
                curr_size = 0
                num_indices += 1
                curr_index = source + "_" + genome + "." + str(num_indices)
                clusters[curr_index] = {"files": {}}
                break
            else:
                clusters[curr_index]["files"][file_name] = files_info[file_name]
                curr_size += files_info[file_name]["file_size"]
                del files_info[file_name]
                break
    clusters[curr_index]["index_size"] = curr_size
    return clusters


def UCSC_collect_file_info(genome):
    tracks = requests.get(url="https://api.genome.ucsc.edu/list/tracks?genome={};trackLeavesOnly=1".format(genome)).json()[genome]
    files_info = {}
    with connect_SQL_db("genome-mysql.soe.ucsc.edu", "genome") as db:
        for track, info in tracks.items():
            if "table" in info:
                    track = info["table"]
            if "bigDataUrl" in info:
                files_info[track] = {"file_size": info["itemCount"]}
                files_info[track]["download_function"] = UCSC_get_big_data_file
                files_info[track]["download_params"] = [track, info["bigDataUrl"], info["type"]]
            else: 
                columns = list(pd.read_sql("Show columns from {}.{}".format(genome, track), con=db)["Field"])
                if len(columns) > 2: # Some files just don't have any info on them
                    files_info[track] = {"file_size": info["itemCount"]}
                    bed_columns = target_columns(columns, info["type"])
                    files_info[track]["download_function"] = UCSC_download_sql_file
                    files_info[track]["download_params"] = [track, columns, info["type"]]
    return files_info

def UCSC_download_sql_file(index_name, track_name, params):
    print(index_name, track_name)
    # genome = track_info[0]
    # track = track_info[1]
    # file_type = track_info[2]
    # with connect_SQL_db("genome-ysql.soe.ucsc.edu", "genome") as db:
    #     conn = db.cursor()
    #     columns = list(pd.read_sql("Show columns from {}.{}".format(genome, track), con=db)["Field"])
    #     result = target_columns(columns, track_info["type"])
    #     conn.execute("USE {}".format(genome))
    #     conn.execute()


def UCSC_get_big_data_file(index_name, track_name, params):
    print(index_name, track_name)


def setup_UCSC_indices(genomes):
    for genome in genomes:
        files_info = UCSC_collect_file_info(genome)
        clustered_files_info = cluster_data("UCSC", genome, files_info)
        for index, index_info in clustered_files_info.items():
            for file, file_info in index_info["files"].items():
                file_info["download_function"](index, file, file_info["download_params"])
        
def setup_indices():
    """ Sets up all sources and genomes listed in CONFIG
    """
    if len(app.config["UCSC_GENOMES"]) > 0:
        setup_UCSC_indices(app.config["UCSC_GENOMES"])

