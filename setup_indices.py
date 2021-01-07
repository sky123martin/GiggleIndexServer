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
import re
import numpy as np

config = config

###########################
# GENERAL SETUP FUNCTIONS #
###########################

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
    cmd_str = "mkdir -p " + temp_path

    proc = subprocess.check_call(cmd_str,
                                    stdout=f,
                                    shell=True)


    cmd_str = 'giggle/scripts/sort_bed \"' + path + '/*.bed*\" ' + temp_path + ' 4'
    proc = subprocess.check_call(cmd_str,
                                    stdout=f,
                                    shell=True)

    delete_directory(path)
    cmd_str = 'mv ' + temp_path + ' ' + path
    proc = subprocess.check_call(cmd_str,
                                    stdout=f,
                                    shell=True)

    cmd_str = 'giggle index -s -f -i \"' + path + '/*.bed.gz\" -o ' + dest + '.d'
    proc = subprocess.check_call(cmd_str,
                                    stdout=f,
                                    shell=True)

    delete_directory(temp_path)
    f.close()


def setup_indices(source, genome, index, index_info, metadata, conn):
    with alive_bar(len(index_info["files"])+1, bar = 'blocks', spinner = 'classic') as bar:
        # create index directory
        proc = subprocess.check_call("mkdir -p data/"+index, shell=True)
        # add index info to the Index database
        # iterate through all files belonging to the index and download them
        for track_name, file_info in index_info["files"].items():
            file_metadata = metadata[metadata["file_name"] == track_name].iloc[0]
            bar.text("Creating {}: Downloading {}".format(index, track_name))
            bar()
            file_info["download_function"](index, file_info["download_params"])
            if os.path.isfile("data/"+index + "/" + track_name + ".bed") or os.path.isfile("data/"+index + "/" + track_name + ".bed.gz"):
                conn.execute("INSERT INTO FILES (NAME, DATE, SOURCE, GENOME, SIZE, INDEXNAME, SHORTNAME, LONGNAME, SHORTINFO, LONGINFO) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                        (
                                        track_name,
                                        date.today(),
                                        source,
                                        genome,
                                        file_info["file_size"],
                                        index,
                                        file_metadata["short_name"],
                                        file_metadata["long_name"],
                                        str(file_metadata["short_info"]),
                                        str(file_metadata["long_info"]),
                                        ))
            else:
                print("file not found:", file_info)

        bar.text("Creating {}: indexing files".format(index))
        if len(glob.glob(os.path.join("data/"+index, "*"))) > 0:
            conn.execute("INSERT INTO INDICES (NAME, ITER, DATE, SOURCE, GENOME, FULL, SIZE) VALUES (?, ?, ?, ?, ?, ?, ?)", 
                (index, index.split("_")[-1], date.today(), source, genome, index_info["full"], index_info["index_size"],))

            giggle_index("data/"+index, "indices/"+index)

        bar()

        bar.text("Completed {}: deleting files".format(index))



        bar.text("Index {} Completed".format(index))
        print("Index {} Completed".format(index))

#########################
# UCSC GENOME FUNCTIONS #
#########################


def setup_UCSC_GENOMES(genomes, conn):
    for genome in genomes:
        # collect file information
        files_info = UCSC_collect_file_info(genome, "")
        # collect metadata for files
        metadata = UCSC_metadata(genome)
        #metadata.to_csv("metadata.csv")
        # cluster files based on hyperparam in .config
        clustered_files_info = cluster_data("UCSC", genome, files_info)
        # iterate through each cluster then download and index cluster files

        for index, index_info in clustered_files_info.items():
            # extract files for the index
            files = list(index_info["files"].keys())
            # grab metadata pertaining to those files
            relevant_metadata = metadata[metadata['file_name'].isin(files)]
            index_metadata = pd.merge(relevant_metadata, 
                                    pd.DataFrame([[i] for i in files], columns=['file_name']),
                                    on ="file_name",
                                    how ="right")
            index_metadata.fillna("")
            # setup this index
            setup_indices("UCSC", genome, index, index_info, index_metadata, conn)


def UCSC_collect_file_info(genome, HUB_EXT):
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
    request_url = "{}/list/tracks?{}genome={};trackLeavesOnly=1".format(config.UCSC_API, HUB_EXT, genome)
    print("REQUEST URL", request_url)
    tracks = requests.get(url=request_url).json()[genome]

    files_info = {}  # Create empty dict to fill with info
    with alive_bar(len(tracks), bar='bubbles', spinner='classic') as bar:  # Start progress bar
        for track, info in tracks.items():
            # Progress bar
            bar.text("Collecting file info for {}: {}".format(genome, track))
            bar()
            # Sometimes the track name != table name, restore track as table name
            track = info["table"] if "table" in info else track

            # If bigData URL is included and not empty then file has been compressed for storage
            if "bigDataUrl" in info.keys() and info["itemCount"] > 0 and info["bigDataUrl"].split(".")[-1] in config.UCSC_ACCEPTABLE_FILE_FORMATS:  # if stored as a big data file
                files_info[track] = {"file_size": info["itemCount"],
                                     "download_function": UCSC_download_bigDataUrl_file,
                                     "download_params": [track, genome, info["bigDataUrl"], info["type"]]
                                    }   

            # If no bigDataURL than it must be uncompressed in sql db
            if "bigDataUrl" not in info.keys() and info["itemCount"] > 0:
                with connect_SQL_db(config.UCSC_SQL_DB_HOST, "genome") as db:
                    # Retrieve column name from sql db track table
                    columns = list(pd.read_sql("Show columns from {}.{}".format(genome, track), con=db)["Field"])
                    if len(columns) > 4:  # Some files just don't have any info in them
                        bed_columns = extract_bed_columns(columns, info["type"])
                        if "" not in bed_columns:  # else 
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
        os.system("UCSC_utilities/bigBedToBed " + config.UCSC_BIG_DATA_LINK + "{} data/{}/{}.bed".format(bigDataURL, index, track))
    elif "bw" == file_type:
        os.system("UCSC_utilities/bigWigToBedGraph " + config.UCSC_BIG_DATA_LINK + "{} data/{}/{}".format(bigDataURL, index, track + ".bed"))
    elif "bigPsl" == file_type:
        os.system("UCSC_utilities/bigPslToPsl " + config.UCSC_BIG_DATA_LINK + "{} data/{}/{}.psl".format(bigDataURL, index, track))
        os.system("UCSC_utilities/pslToBed data/{}/{} data/{}/{}.bed".format(index, track + ".psl", index, track))
        os.remove("data/{}/{}".format(index, track + ".psl"))
    else:
        print("File {} of type not able to be converted".format(track, file_type))


def UCSC_metadata(genome):
    with connect_SQL_db(config.UCSC_SQL_DB_HOST, "genome") as db:
        query = "SELECT tableName AS file_name, shortLabel AS short_name, longLabel as long_name, html as long_info from {}.trackDb order by tableName".format(genome)
        df_metadata = pd.read_sql(query, con=db)

    df_metadata = df_metadata.fillna("")
    df_metadata["short_info"] = ""

    for _, row in df_metadata.iterrows():
        row["long_info"] = row["long_info"].decode('latin-1') 
        if row["long_info"] != "":
            row["short_info"] = UCSC_metadata_extract_description(row["long_info"])

    return df_metadata


def UCSC_metadata_extract_description(html):
    # Format of descriptions <H3>Description</H3> <P> ... <P> <h2>Description</h2> <p>
    if "Description" in html:
        result = re.search('^[ \t]*<[hH]{1}[0-9]{1}>[ \t]*Description[ \t]*<\/[hH]{1}[0-9]{1}>(.*?)<[pP]{1}>[ \t]*(.*?)<\/[pP]{1}>', html, flags = re.DOTALL)
        if result != None:
            htmldescr = result.group(0)
            index = htmldescr.find("<P>")
            if index == -1:
                index = htmldescr.find("<p>")
            return htmldescr[index:len(htmldescr)]
    return ""

##########################
# LOCAL GENOME FUNCTIONS #
##########################

def setup_LOCAL(projects, conn):
    for project in projects: #genome, paths in genomes.items():
        # collect file information
        files_info = local_collect_file_info(project["data_path"])
        # collect metadata for files
        metadata = local_metadata(project["metadata_path"])
        # cluster files based on hyperparam in .config
        clustered_files_info = cluster_data(project["project_name"], project["reference_genome"], files_info)
        # iterate through each cluster then download and index cluster files

        for index, index_info in clustered_files_info.items():
            files = list(index_info["files"].keys())
            index_metadata = pd.merge(metadata, 
                                    pd.DataFrame([[i] for i in files], columns = ['file_name']),
                                    on = "file_name",
                                    how = "right")
            index_metadata.fillna("")
            
            setup_indices(project["project_name"], project["reference_genome"], index, index_info, index_metadata, conn)


def local_collect_file_info(path):
    file_list = glob.glob(os.path.join(path, "*.bed*"))

    files_info = {}
    for file_name in file_list:
        file_name = file_name.split("/")[-1]
        file_type = file_name.split(".")[1]
        track_name = clean_file_name(file_name)
        interval_count = subprocess.check_output("wc -l {}/{}".format(path, file_name), shell=True, text=True)
        files_info[track_name] = {"file_size": int(interval_count.split()[0]),
                                  "download_function": LOCAL_download_file,
                                  "download_params": [path, file_name, file_type]
                                  }
    return files_info

def clean_file_name(file_name):
    return file_name.split("/")[-1].split(".")[0]

def local_metadata(path):
    necesary_columns = ["file_name", "short_name", "long_name", "short_info", "long_info"]
    metadata_exists = False 

    if path != "":
        metadata_exists = True 
        metadata = pd.read_csv(path)
        # error handling if metadata doesn't meet format
        columns = list(metadata.columns)
        for col in necesary_columns:
            if col not in columns:
                metadata_exists = False 
                print("Metadata file {} does not contain column {}. Please check the LOCAL_METADATA formatting guidelines in Config.py. Indexing will continue without metadata.".format(config.LOCAL_METADATA, col))
        metadata["file_name"] = metadata["file_name"].apply(clean_file_name)

    # if metadata doesn't follow standards then replace metadata with empty dataframe
    if metadata_exists == False:
        metadata = pd.DataFrame(columns = necesary_columns)

    return metadata



def LOCAL_download_file(index, params):
    file_name = params[1]
    path = params[0]
    copyfile(path+"/"+file_name, "data/"+index+"/"+file_name)


if __name__ == "__main__":
    os.system('python3 tear_down_indices.py')  # delete previous indices
    os.system('python3 models.py')  # setup indexing and files database
    conn = sqlite3.connect('Indexing.db')  # make connection to database

    # make directories for data and indicies 
    proc = subprocess.check_call("mkdir -p data", shell=True)
    proc = subprocess.check_call("mkdir -p indices", shell=True)

    # setup data sources
    setup_UCSC_GENOMES(config.UCSC_GENOMES, conn)
    setup_LOCAL(config.LOCAL_GENOMES, conn)

    conn.commit()
    conn.close()