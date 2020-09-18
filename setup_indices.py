from config import config
from tear_down_indices import delete_directory
import mysql.connector
import pandas as pd
from contextlib import contextmanager
import requests 
import glob
import os
import pybedtools
from alive_progress import alive_bar
from datetime import date
import sqlite3
from shutil import copyfile
import subprocess

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

def extract_bed_columns(columns, file_type): 
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


def cluster_data(source, genome, files_info, index_num = 1):
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
        num_indices = index_num
        curr_index = source + "_" + genome + "_" + str(num_indices)
        clusters[curr_index] = {"files": {}}
        curr_size = 0
        while len(files_info) != 0:
            for file_name in files_info:
                if curr_size > max_size and len(files_info) > 0:
                    # reset size and name
                    clusters[curr_index]["full"] = True
                    clusters[curr_index]["index_size"] = curr_size
                    curr_size = 0
                    num_indices += 1
                    curr_index = source + "_" + genome + "_" + str(num_indices)
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
    track = params[0]
    genome = params[1]
    chrom = params[2]["chrom"]
    start = params[2]["start"]
    end = params[2]["end"]
    base_shift = params[2]["base shift"]

    with connect_SQL_db(config.UCSC_SQL_DB_HOST, "genome") as db:
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

    bigDataURL = params[2]
    file_type = bigDataURL.split(".")[-1]
    track = params[0]
    # UCSC_utilities/bigBedToBed input_file output_file
    if "bb" == file_type:
        os.system("UCSC_utilities/bigBedToBed http://hgdownload.soe.ucsc.edu/{} data/{}/{}.bed".format(bigDataURL, index, track))
    elif "bw" == file_type:
        os.system("UCSC_utilities/bigWigToBedGraph http://hgdownload.soe.ucsc.edu/{} data/{}/{}".format(bigDataURL, index, track + ".bed"))
    elif "bigPsl" == file_type:
        os.system("UCSC_utilities/bigPslToPsl http://hgdownload.soe.ucsc.edu/{} data/{}/{}.psl".format(bigDataURL, index, track))
        os.system("UCSC_utilities/pslToBed data/{}/{} data/{}/{}.bed".format(index, track + ".psl", index, track))
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
    # Check API for file information
    tracks = requests.get(url="{}/list/tracks?genome={};trackLeavesOnly=1".format(config.UCSC_API,genome)).json()[genome]
    files_info = {}  # Create empty dict to fill with info
    with alive_bar(len(tracks), bar='bubbles', spinner='classic') as bar:  # Start progress bar
        for track, info in tracks.items():
            # Progress bar
            bar.text("Collecting file info for {}: {}".format(genome, track))
            bar()
            # Sometimes the track name != table name, restore track as table name
            if "table" in info: #FIXME store table name and track name so file can be retrieved
                track = info["table"]
            
            # If bigData URL is included and not empty then file has been compressed for storage
            if "bigDataUrl" in info.keys() and info["itemCount"] > 0 and info["bigDataUrl"].split(".")[-1] in config.UCSC_ACCEPTABLE_FILE_FORMATS: # if stored as a big data file
                files_info[track] = {"file_size": info["itemCount"],
                                     "download_function": UCSC_download_bigDataUrl_file,
                                     "download_params": [track, genome, info["bigDataUrl"], info["type"]]}       
           
            # If no bigDataURL than it must be uncompressed in sql db
            if "bigDataUrl" not in info.keys() and info["itemCount"] > 0:
                with connect_SQL_db(config.UCSC_SQL_DB_HOST, "genome") as db:
                    # Retrieve column name from sql db track table
                    columns = list(pd.read_sql("Show columns from {}.{}".format(genome, track), con=db)["Field"])
                    if len(columns) > 4:  # Some files just don't have any info in them
                        bed_columns = extract_bed_columns(columns, info["type"])
                        if "" not in bed_columns: # else 
                            files_info[track] = {"file_size": info["itemCount"],
                                                 "download_function": UCSC_download_sql_file,
                                                 "download_params": [track, genome, bed_columns]}
                    # Some bigDataURL files are store in sql db
                    elif "fileName" in columns:  # bigWig files do not like big data url in API, only in sql table
                        big_data_URL = list(pd.read_sql("Select fileName from {}.{}".format(genome, track), con=db)["fileName"])[0]
                        if big_data_URL.split(".")[-1] in config.UCSC_ACCEPTABLE_FILE_FORMATS:
                            files_info[track] = {"file_size": info["itemCount"],
                                                "download_function" : UCSC_download_bigDataUrl_file,
                                                "download_params": [track, genome, big_data_URL, info["type"]]}   

    return files_info


def giggle_index(path, dest):
    """ Gzips and indexes a given directory and destination for index
    Parameters
    ----------
    path: path to a directory of bed files
        ex. data/UCSC_hg19_1/
    dest:
        ex. index/UCSC_hg19_1.d
    Returns
    -------
    nothing
    """

    f = open('giggle_log.txt','w') 
    temp_path = path + "_sorted"
    subprocess.check_call("mkdir -p " + temp_path , shell=True)
    cmd_str = 'giggle/scripts/sort_bed \"' + path + '/*.bed*\" '+temp_path +' 4'
    proc = subprocess.check_call(cmd_str,
                                    stdout=f,
                                    shell=True)
    delete_directory(temp_path)

    cmd_str = 'giggle index -s -f -i \"' + path + '/*.bed.gz\" -o ' + dest + '.d'
    proc = subprocess.check_call(cmd_str,
                                    stdout=f,
                                    shell=True)
    delete_directory(path)
    f.close()


def setup_indices(source, genome, index, index_info, conn):
    with alive_bar(len(index_info["files"])+1, bar = 'blocks', spinner = 'classic') as bar:
        # create index directory
        proc = subprocess.check_call("mkdir -p data/"+index, shell=True)
        # add index info to the Index database
        # iterate through all files belonging to the index and download them
        for track_name, file_info in index_info["files"].items():
            bar.text("Creating {}: Downloading {}".format(index, track_name))
            bar()
            file_info["download_function"](index, file_info["download_params"])
            if os.path.isfile("data/"+index + "/" + track_name + ".bed") or os.path.isfile("data/"+index + "/" + track_name + ".bed.gz"):
                conn.execute("INSERT INTO FILES (NAME, DATE, SOURCE, GENOME, SIZE, INDEXNAME) \
                    VALUES ('{}', {}, '{}', '{}', {}, '{}')".format(track_name, date.today(),
                    source, genome, file_info["file_size"], index))
            else:
                print(file_info)
       

        bar.text("Creating {}: indexing files".format(index))
        if len(glob.glob(os.path.join("data/"+index, "*"))) > 0:
            conn.execute("INSERT INTO INDICES (NAME, ITER, DATE, SOURCE, GENOME, FULL, SIZE) VALUES ('{}', {}, {}, '{}', '{}', {}, {})".format(
                index, index.split("_")[-1], date.today(), source, genome, index_info["full"], index_info["index_size"]))

            giggle_index("data/"+index, "indices/"+index)

        bar()

        bar.text("Completed {}: deleting files".format(index))



        bar.text("Index {} Completed".format(index))
        print("Index {} Completed".format(index))



def setup_UCSC_indices(genomes, conn):
    for genome in genomes:
        files_info = UCSC_collect_file_info(genome)
        clustered_files_info = cluster_data("UCSC", genome, files_info)
        # Iterate through each index
        for index, index_info in clustered_files_info.items():
            setup_indices("UCSC", genome, index, index_info, conn)


def LOCAL_download_file(index, params):
    file_name = params[1]
    path = params[0]
    copyfile(path+"/"+file_name, "data/"+index+"/"+file_name)


def local_collect_file_info(path):
    file_list = glob.glob(os.path.join(path, "*.bed*"))
    # file_list += glob.glob(os.path.join(path, "*.bed.gz"))
    files_info = {}
    for file_name in file_list:
        file_name = file_name.split("/")[-1]
        file_type = file_name.split(".")[1]
        track_name = file_name.split(".")[0]
        interval_count = subprocess.check_output("wc -l {}/{}".format(path, file_name), shell=True, text=True)
        files_info[track_name] = {"file_size": int(interval_count.split()[0]),
                                 "download_function": LOCAL_download_file,
                                 "download_params": [path, file_name, file_type]}
    return files_info

def setup_local_indices(genomes, conn):
    for genome, path in genomes.items():
        files_info = local_collect_file_info(path)
        clustered_files_info = cluster_data("local", genome, files_info)
        for index, index_info in clustered_files_info.items():
            setup_indices("local", genome, index, index_info, conn)


if __name__ == "__main__":
    os.system('python3 tear_down_indices.py')  # delete previous indices
    os.system('python3 models.py')  # set up indexing and files database
    conn = sqlite3.connect('Indexing.db')  # make connection to database

    # make directories for data
    proc = subprocess.check_call("mkdir -p data", shell=True)
    proc = subprocess.check_call("mkdir -p indices", shell=True)

    #setup indices
    setup_UCSC_indices(config.UCSC_GENOMES, conn)
    setup_local_indices(config.LOCAL_GENOMES, conn)

    conn.commit()
    conn.close()
