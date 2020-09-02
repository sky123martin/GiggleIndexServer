from config import config
import mysql.connector
import pandas as pd
from contextlib import contextmanager
import requests 
import glob
import os
import os.path
import pybedtools
from alive_progress import alive_bar
from datetime import date
import sqlite3

config = config

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

    return {"chrom": chrom, "start": chrom_start, "end": chrom_end, "base shift": base_shift}


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
    with alive_bar(spinner = 'classic') as bar:
        bar.text("Clustering genome {}".format(genome))
        clusters = {}
        max_size = config.MAX_INTERVALS_PER_INDEX
        num_indices = 1
        curr_index = source + "_" + genome + "." + str(num_indices)
        clusters[curr_index] = {"files": {}}
        curr_size = 0
        while len(files_info) != 0:
            for file_name in files_info:
                if curr_size > max_size and len(files_info) != 0:
                    # reset size and name
                    clusters[curr_index]["full"] = True
                    clusters[curr_index]["index_size"] = curr_size
                    curr_size = 0
                    num_indices += 1
                    curr_index = source + "_" + genome + "." + str(num_indices)
                    clusters[curr_index] = {"files": {}}
                    break
                else:
                    clusters[curr_index]["full"] = False
                    clusters[curr_index]["files"][file_name] = files_info[file_name]
                    curr_size += files_info[file_name]["file_size"]
                    del files_info[file_name]
                    break
        clusters[curr_index]["index_size"] = curr_size
        print("Finished clustering genome {}, total of {} indices".format(genome, len(clusters)))
    return clusters


def UCSC_download_sql_file(index, params):
    """ downloads a specific file from UCSC sql db and stores in index folder
    Parameters
    ----------
    index: string
        ex. hg19_1
    params: list
        [track, genome, extract_columns]
        extract_columns =[chrom, start, end, base_shift]
    Returns
    -------
    nothing

    """
    return
    track = params[0]
    genome = params[1]
    chrom = params[2]["chrom"]
    start = params[2]["start"]
    end = params[2]["end"]
    base_shift = params[2]["base shift"]

    with connect_SQL_db(configUCSC_SQL_DB_HOST, "genome") as db:
        df = pd.read_sql("Select {} AS chrom, {} + {} AS start, {} + {} AS end from {}.{}".format(chrom, start, base_shift, end, base_shift, genome, track), con=db)
        pybedtools.BedTool.from_dataframe(df).saveas('data/{}/{}.bed'.format(index, track))


def UCSC_download_bigDataUrl_file(index, params):
    """ downloads a specific file from URL, converts to bed through,
    UCSC Utilities folder and stores in index folder

    Parameters
    ----------
    index: string
        ex. hg19_1
    params: list
        [track, genome, bigDataUrl, file_type]
    Returns
    -------
    nothing

    """

    file_type = params[3]
    bigDataURL = params[2]
    track = params[0]
    # UCSC_utilities/bigBedToBed input_file output_file
    # os.system("" [track, genome, info["bigDataUrl"], info["type"]]
    if "bigBed" in file_type:
        os.system("UCSC_utilities/bigBedToBed http://hgdownload.soe.ucsc.edu/{} data/{}/{}".format(bigDataURL, index, track + ".bed"))
    elif "bigWig" in file_type:
        os.system("UCSC_utilities/bigWigToBedGraph http://hgdownload.soe.ucsc.edu/{} data/{}/{}".format(bigDataURL, index, track + ".bed"))
        print(track)
    elif "bigPsl" in file_type:
        os.system("UCSC_utilities/bigPslToPsl http://hgdownload.soe.ucsc.edu/{} data/{}/{}".format(bigDataURL, index, track + ".psl"))
        os.system("UCSC_utilities/pslToBed data/{}/{} data/{}/{}".format(index, track + ".psl", index, track + ".bed"))
        os.remove("data/{}/{}".format(index, track + ".psl"))
    else:
        print("File {} of type not able to be converted".format(track, file_type))

def UCSC_collect_file_info(genome):
    """ uses UCSC API to gather download info for each track in a specific genome
    Parameters
    ----------
    genome: string
        ex. hg19

    Returns
    -------
    files_info: dict of file information for each track, look st
    """

    tracks = requests.get(url="{}/list/tracks?genome={};trackLeavesOnly=1".format(config.UCSC_API,genome)).json()[genome]
    files_info = {}
    with alive_bar(len(tracks), bar='bubbles', spinner='classic') as bar:
        for track, info in tracks.items():
            bar.text("Collecting file info for {}: {}".format(genome, track))
            bar()
            if "table" in info:  # sometimes track name !=table name
                track = info["table"]
            if "bigDataUrl" in info: # if stored as a big data file
                files_info[track] = {"file_size": info["itemCount"]}
                                    # "download_function": UCSC_download_bigDataUrl_file}
                files_info[track]["download_function"] = UCSC_download_bigDataUrl_file
                files_info[track]["download_params"] = [track, genome, info["bigDataUrl"], info["type"]]
            else:  # if stored in sql db
                with connect_SQL_db(config.UCSC_SQL_DB_HOST, "genome") as db:
                    columns = list(pd.read_sql("Show columns from {}.{}".format(genome, track), con=db)["Field"])
                    if len(columns) > 2:  # Some files just don't have any info on them
                        files_info[track] = {"file_size": info["itemCount"]}
                        extract_columns = target_columns(columns, info["type"])
                        files_info[track]["download_function"] = UCSC_download_sql_file
                        files_info[track]["download_params"] = [track, genome, extract_columns]
                    elif "fileName" in columns:  # bigWig files do not like big data url in API, only in sql table
                        big_data_URL = list(pd.read_sql("Select fileName from {}.{}".format(genome, track), con=db)["fileName"])[0]
                        files_info[track] ={"file_size": info["itemCount"]}
                        files_info[track]["download_function"] = UCSC_download_bigDataUrl_file
                        files_info[track]["download_params"] = [track, genome, big_data_URL, info["type"]]
    return files_info


def setup_index(source, genome, index, index_info):
    with alive_bar(len(index_info["files"])+3, bar = 'blocks', spinner = 'classic') as bar:
        # create index directory
        try:
            os.mkdir("data/" + index)
        except OSError:
            print("Creation of the directory %s failed" % index)
        # add index info to the Index database
        conn.execute("INSERT INTO INDICES (NAME, ITER, DATE, SOURCE, GENOME, FULL, SIZE) VALUES ('{}', {}, {}, '{}', '{}', {}, {})".format(index, 
            index.split(".")[1], date.today(), source, genome, index_info["full"], index_info["index_size"]))
        # iterate through all files belonging to the index and download them
        for filename, file_info in index_info["files"].items():
            bar.text("Creating {}: Downloading {}".format(index, filename))
            bar()
            file_info["download_function"](index_name, file_info["download_params"])
            conn.execute("INSERT INTO FILES (NAME, DATE, SOURCE, GENOME, SIZE, INDEXNAME) \
                VALUES ('{}', {}, '{}', '{}', {}, '{}')".format(filename, date.today(), 
                source, genome, file_info["file_size"], index))
        # FIXME index downloaded files here
        
        filelist = glob.glob(os.path.join("data/"+index, "*.bed"))

        bar.text("Creating {}: indexing files".format(index))
        bar()
        if index_info["full"]:  # delete files if index is full
            bar.text("Completed {}: deleting files".format(index))
            bar()
            filelist = glob.glob(os.path.join("data/"+index, "*"))
            for f in filelist:
                os.remove(f)
            os.rmdir("data/"+index)
        bar.text("Index {} Completed".format(index))
        print("Index {} Completed".format(index))
        
        cursor = conn.execute("SELECT * from INDICES")
        for row in cursor:
            print(row)

                 
def setup_UCSC_indices(genomes):
    for genome in genomes:
        files_info = UCSC_collect_file_info(genome)
        clustered_files_info = cluster_data("UCSC", genome, files_info)
        # Iterate through each index
        for index, index_info in clustered_files_info.items():
            setup_index("UCSC", genome, index, index_info)

def setup_local_indices(genomes):
    for genome in genomes:
        files_info = local_collect_file_info(genome)


def setup_indices():
    """ Sets up all sources and genomes listed in CONFIG
    """
    setup_UCSC_indices(config.UCSC_GENOMES)
    setup_local_indices(config.LOCAL_GENOMES)

if __name__ == "__main__":
    os.system('python models.py') # set up indexing and files database
    global conn
    conn = sqlite3.connect('Indexing.db')
    setup_indices()

