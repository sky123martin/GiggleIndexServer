from config import config
import mysql.connector
import pandas as pd
from contextlib import contextmanager
import requests 
import glob
import os
import pybedtools
from multiprocessing.pool import Pool
from functools import partial
from alive_progress import alive_bar
from datetime import date
import sqlite3
import subprocess
import re
import numpy as np
from math import isnan
import math
import shutil


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
        print("Given file columns, {}, does not have a conversion to .bed for file type {}".format(columns, file_type))
        return "bed mapping not found"

    return {chrom: "chrom", chrom_start: "start", chrom_end: "end"}, base_shift

def file_download_handler(path, file_info):
    try:
        file_info["download_function"](path, file_info["download_params"])
    except Exception as e:
        print(e)

def download_all_files(folder, files_info):
    proc = subprocess.check_output("mkdir -p data/{}".format(folder), shell=True)
    with Pool(config.AVAILABLE_PROCCESSES) as p:  # to check multiprocessing.cpu_count()
        output = p.map(partial(file_download_handler, folder), files_info.values())


def giggle_sort(path):
    sorted_path = path + "_sorted"
    temp_path = path + "_temp"

    # make sorted directory
    os.makedirs(sorted_path, exist_ok=True)

    # split into 500 files groups and then sort bed
    batch_size = 500
    file_list = os.listdir(path)
    while len(file_list) > batch_size:
        os.makedirs(temp_path, exist_ok=True)
        # sort batch of files
        for f in file_list[:batch_size]:
            shutil.move(path+"/"+f, temp_path)

        cmd_str = 'giggle/scripts/sort_bed \"{}/*.bed*\" {} {}'.format(temp_path, sorted_path, config.AVAILABLE_PROCCESSES)
        proc = subprocess.check_output(cmd_str,
                                    shell=True,
                                    timeout=2*60)
        proc = subprocess.check_output("rm -R -f {}".format(temp_path), shell=True)
        proc = subprocess.check_output("mkdir -p {}".format(temp_path), shell=True, timeout=config.timeout_file_processing)
        # remove unsorted files from queue
        file_list = file_list[batch_size:]

    # sort remaining files using giggle sort
    cmd_str = 'giggle/scripts/sort_bed \"{}/*.bed*\" {} {}'.format(path, sorted_path, config.AVAILABLE_PROCCESSES)
    print(cmd_str)
    proc = subprocess.check_output(cmd_str,
                                    shell=True,
                                    timeout=2*60)
    
    proc = subprocess.check_output("rm -R -f {}".format(path), shell=True)
    proc = subprocess.check_output("rm -R -f {}".format(temp_path), shell=True)

def giggle_move_index(path, params):
    try:
        index_name = params[0]
        files = params[1]
        print("Making indices")
        index_path = "data/{}".format(index_name)
        os.makedirs(index_path, exist_ok=True)

        for f in files:
            shutil.move("{}/{}.bed.gz".format(path, f),index_path)
            # proc = subprocess.check_output("mv {}/{}.bed.gz {}/".format(path, f, index_path), shell=True)

        cmd_str = 'giggle index -i \"{}/*.bed.gz\" -o indices/{}.d  -f -s'.format(index_path, index_name)
        print(cmd_str)
        proc = subprocess.check_output(cmd_str, timeout=config.timeout_file_processing, shell=True)

        proc = subprocess.check_output("rm -R -f {}".format(index_path), shell=True)
        proc = subprocess.check_output("rm -R -f {}/".format(path), shell=True)
    except Exception as e:
        print(e)
        return "error"    

def download_cluster_index(index_base, source, project, genome, files_info, metadata, conn, current_index_num=1):
    # download all files
    all_files_folder = index_base
    download_all_files(all_files_folder, files_info)
    files_downloaded = os.listdir("data/" + all_files_folder)

    # sort all files
    if len(files_downloaded) > 0:
        giggle_sort("data/{}".format(all_files_folder))
        all_files_folder = "data/{}_sorted".format(all_files_folder)
        files_downloaded = [i.replace(".bed.gz", "") for i in os.listdir(all_files_folder)]
        
        metadata = metadata[metadata["file_name"].isin(files_downloaded)]
        metadata = pd.merge(metadata, 
                            pd.DataFrame([[i] for i in files_downloaded], columns=['file_name']),
                            on ="file_name",
                            how ="right")
        metadata["file_size"] = metadata["file_name"].apply(lambda n: files_info[n]["file_size"] if n in files_info.keys() else 0)

        indices = []
        while metadata.shape[0] > 0:
            current_index_name = "{}_{}".format(index_base, current_index_num)
            current_index_size = 0
            current_index_full = False
            current_index_files = []
            for index, row in metadata.iterrows():
                current_index_size = current_index_size + row["file_size"]
                current_index_files.append(row["file_name"])
                metadata.drop([index], inplace=True)
                # insert file into DB
                conn.execute("INSERT INTO FILES (FILEID, SIZE, GENOME, PROJECTID, INDEXID, SHORTNAME, LONGNAME, SHORTINFO, LONGINFO) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (row["file_name"], row["file_size"], genome, index_base, current_index_name, row["short_name"], row["long_name"], str(row["short_info"]), str(row["long_info"]),))
                
                if current_index_size > config.MAX_INTERVALS_PER_INDEX:
                    current_index_full = True
                    break
            current_index_num = current_index_num + 1
            indices.append([current_index_name, current_index_files])

            conn.execute("INSERT INTO INDICES (INDEXID, ITER, DATE, DSOURCE, PROJECTID, GENOME, FULL, SIZE) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", 
                (current_index_name, int(current_index_name.split("_")[-1]), date.today(), source, index_base, genome, (current_index_full), current_index_size,))

        with Pool(config.AVAILABLE_PROCCESSES) as p:  # to check multiprocessing.cpu_count()
            output = p.map(partial(giggle_move_index, all_files_folder), indices)

#########################
# UCSC GENOME FUNCTIONS #
#########################

def collect_ucsc_files(genome, max_file_size = config.max_file_size, HUB_EXT = ""):
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
    try:
        tracks = requests.get(url=request_url).json()[genome]

        files_info = {}  # Create empty dict to fill with info
        with alive_bar(len(tracks), bar='bubbles', spinner='classic') as bar:  # Start progress bar
            for track, info in tracks.items():
                # Progress bar
                bar.text("Collecting file info for {}: {}".format(genome, track))
                bar()
                # Sometimes the track name != table name, restore track as tablename
                track = info["table"] if "table" in info else track

                # If bigData URL is included and not empty then file has been compressed for storage
                if "bigDataUrl" in info.keys() and info["itemCount"] > 0 and info["bigDataUrl"].split(".")[-1] in config.UCSC_ACCEPTABLE_FILE_FORMATS and max_file_size>info["itemCount"]:  # if stored as a big data file
                    files_info[track] = {"file_size": info["itemCount"],
                                        "download_function": download_linked_file,
                                        "download_params": {"track": track, 
                                                            "download_location": info["bigDataUrl"] if HUB_EXT else (config.UCSC_BIG_DATA_LINK + info["bigDataUrl"]), 
                                                            "file_type": info["bigDataUrl"].replace(".gz","").split(".")[-1]
                                                            }}

                # If no bigDataURL than it must be uncompressed in sql db
                if "bigDataUrl" not in info.keys() and info["itemCount"] > 0 and max_file_size > info["itemCount"]:
                    files_info[track] = {"file_size": info["itemCount"],
                                         "download_function": download_UCSC_sql_file,
                                         "download_params": {"track": track, 
                                                             "genome": genome,
                                                             "file_type": info["type"],
                                                             "hub_ext": HUB_EXT
                                                             }}
                   
        return files_info
    except Exception as e:
        print(e)
        return "error"


def download_UCSC_sql_file(path, params):
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
    track = params["track"]
    genome = params["genome"]
    file_type = params["file_type"]
    HUB_EXT = params["hub_ext"]

    with connect_SQL_db(config.UCSC_SQL_DB_HOST, "genome") as db:
        df = pd.read_sql("Select * from {}.{}".format(genome, track), con=db)
        columns = df.columns
        if "fileName" in columns: # bigWig files do not like big data url in API, only in sql table
            big_data_URL = df["fileName"].iloc[0]
            if big_data_URL.split(".")[-1] in config.UCSC_ACCEPTABLE_FILE_FORMATS:
                download_linked_file(path,{"track":track, "download_location":big_data_URL if HUB_EXT!="" else (config.UCSC_BIG_DATA_LINK+big_data_URL), "file_type":big_data_URL.split(".")[-1]})
        elif len(columns)>3: # some tables are empty
            column_mapping, baseshift = extract_bed_columns(columns, file_type)
            if column_mapping!="bed mapping not found":
                df.rename(columns=column_mapping, inplace=True)
                df = df[["chrom","start","end"]]
                pybedtools.BedTool.from_dataframe(df).saveas('data/{}/{}.bed'.format(path, track))


def download_linked_file(path, params):
    """ downloads a specific file from URL, converts to bed through,
    UCSC Utilities folder and stores in index folder
    Parameters
    ----------
    index: string
        ex. hg19_1
    params: list
        [track, genome, location_path, file_type]
    Returns
    -------
    nothing
    """

    try:
        # download file from
        if 'local' in path:
            shutil.copyfile(params["download_location"], "data/{}/{}.{}".format(path, params["track"], params["file_type"]))
        else:
            cmd_str = "curl {} --output data/{}/{}.{}".format(params["download_location"], path, params["track"], params["file_type"])
            subprocess.check_output(cmd_str, shell=True, timeout = config.timeout_file_download)
        params["download_location"] = "data/{}/{}.{}".format(path, params["track"], params["file_type"])

        try:
            if "bb" == params["file_type"]:
                cmd_str  = "UCSC_utilities/bigBedToBed {} data/{}/{}.bed".format(params["download_location"], path, params["track"])
                subprocess.check_output(cmd_str, shell=True, timeout = config.timeout_file_download)
                os.remove(params["download_location"])
            elif "bw" == params["file_type"]:
                cmd_str  = "UCSC_utilities/bigWigToBedGraph {} data/{}/{}.bed".format(params["download_location"], path, params["track"])
                subprocess.check_output(cmd_str, shell=True, timeout = config.timeout_file_download)
                os.remove(params["download_location"])
            elif "bigPsl" == params["file_type"]:
                cmd_str = "UCSC_utilities/bigPslToPsl {} data/{}/{}.psl".format(params["download_location"], path, params["track"])
                subprocess.check_output(cmd_str, shell=True, timeout = config.timeout_file_download)
                cmd_str = "UCSC_utilities/pslToBed data/{}/{}.psl data/{}/{}.bed".format(path, params["track"], path, params["track"])
                subprocess.check_output(cmd_str, shell=True)
                os.remove(params["download_location"])
                os.remove("data/{}/{}.psl".format(path, params["track"]))
            elif "bed" == params["file_type"]:
                pass
            else:
                print("File {} of type not able to be converted".format(params["track"], params["file_type"]))

        except Exception as e:
            # some files are labeled as .bb but store html with an alternate url
            if int(subprocess.check_output("wc -l {}".format(params["download_location"]), shell=True, text=True).split()[0]) < 15:
                f = open(params["download_location"], "r")
                url_found = False
                for l in f:
                    if "https://" in l and params["file_type"] in l:
                        print("LINE", l )
                        start = l.find("https://")
                        stop = l.find(params["file_type"])
                        if start != -1 and stop !=-1:
                            url_found = True
                            url = l[start:stop+len(params["file_type"])]
                            params["download_location"] = url
                            download_linked_file(path, params)
                            break
                if not url_found:
                    os.remove(params["download_location"])

            else:
                os.remove(params["download_location"])
    except Exception as e:
        print("Exception thrown for track", params["track"], e)


def UCSC_metadata(genome):
    with connect_SQL_db(config.UCSC_SQL_DB_HOST, "genome") as db:
        query = "SELECT tableName AS file_name, shortLabel AS short_name, longLabel as long_name, html as long_info from {}.trackDb order by tableName".format(genome)
        df_metadata = pd.read_sql(query, con=db)
    
    df_metadata = df_metadata.fillna("")
    df_metadata["short_info"] = ""

    for _, row in df_metadata.iterrows():
        row["long_info"] = row["long_info"].decode('latin-1') 
        if row["long_info"] != "":
            row["short_info"] = HTML_extract_description(row["long_info"])
    return df_metadata


def HTML_extract_description(html):
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

######################
# UCSC HUB FUNCTIONS #
######################
def add_genome(genome, science_name, common_name, description, conn):
    # add genome if it has yet to be added
    cursor = conn.execute("SELECT GENOME from GENOMES WHERE GENOMES.GENOME = '{}'".format(genome))
    df = pd.DataFrame(cursor.fetchall())
    if df.empty:
        conn.execute("INSERT INTO GENOMES (GENOME, SCINAME, COMMONNAME, DESCRIPTION) VALUES (?, ?, ?, ?)", (genome, science_name, common_name, description,))


def setup_ucscGenomes(genomes, conn):
    url = "{}/list/ucscGenomes".format(config.UCSC_API)
    genome_info = requests.get(url = url).json()["ucscGenomes"]
    for genome in genomes:
        # add genome to GENOMES
        add_genome(genome, genome_info[genome]["scientificName"], genome_info[genome]["organism"], genome_info[genome]["description"], conn)

        # add UCSC genome info to PROJECTS
        project_id = "ucscGenomes_" + genome.replace(" ","-")
        conn.execute("INSERT INTO PROJECTS (PROJECTID, OSOURCE,  DSOURCE, SHORTNAME, LONGNAME, INFO) VALUES (?, ?, ?, ?, ?, ?)", 
                (project_id, genome_info[genome]["sourceName"], "ucscGenomes", genome, genome_info[genome]["description"], genome_info[genome]["htmlPath"],))
        # collect file information
        files_info = collect_ucsc_files(genome, max_file_size = config.max_setup_file_size)
        if isinstance(files_info, str):  # 404 not found on API
            print("Unable to collect file info from UCSC for {} ref genome".format(genome))
            continue
        # collect metadata for files
        metadata = UCSC_metadata(genome)
        # only setup files that have metadata, take intersection of metadata & files_info
        files_info = {key: files_info[key] for key in list(set(files_info.keys()) & set(metadata['file_name'].to_numpy()))}
        
        download_cluster_index(project_id, "ucscGenomes", genome.replace(" ","-"), genome, files_info, metadata, conn)

def setup_ucscHubs(hub_names, conn):
    hubs = collect_ucscHubs(hub_names)
    for hub in hubs:
        url = "{}/list/hubGenomes?hubUrl={}".format(config.UCSC_API, hub["hub_url"])
        hub_genomes = requests.get(url = url).json()["genomes"]
        for genome in hub["genomes"]:
            if genome not in hub_genomes.keys():
                continue
            # add genome to GENOMES
            scientific_name = hub_genomes[genome]["scientificName"] if "scientificName" in hub_genomes[genome].keys() else ""
            organism = hub_genomes[genome]["organism"] if "organism" in hub_genomes[genome].keys() else ""
            description = hub_genomes[genome]["description"] if "description" in hub_genomes[genome].keys() else ""
            add_genome(genome, scientific_name, organism, description, conn)

            # add hub info to PROJECTS
            project_id = "ucscHubs_{}_{}".format(hub["hub_short_label"].replace(" ","-").replace("_","-"), genome)
            conn.execute("INSERT INTO PROJECTS (PROJECTID, OSOURCE,  DSOURCE, SHORTNAME, LONGNAME, INFO) VALUES (?, ?, ?, ?, ?, ?)", 
                (project_id, hub["hub_url"], "ucscHubs", hub["hub_short_label"], hub["hub_long_label"], hub["descriptionUrl"],))
            # collect file information
            files_info = collect_ucsc_files(genome, max_file_size = config.max_setup_file_size, HUB_EXT="hubUrl={};".format(hub["hub_url"]))
            if isinstance(files_info, str):
                continue
            # iterate through each cluster then download and index cluster files
            metadata = UCSC_hubs_metadata(hub, genome)
            download_cluster_index(project_id, "ucscHubs", hub["hub_short_label"], genome, files_info, metadata, conn)

def collect_ucscHubs(hub_names):
    hubs = requests.get(url = config.UCSC_API + "/list/publicHubs").json()["publicHubs"]
    hubs_info = []
    for hub in hubs:
        if hub["shortLabel"] in hub_names:#hub["dbList"].split(","):
            hubs_info.append({
                "hub_url": hub["hubUrl"],
                "hub_short_label": hub["shortLabel"],
                "hub_long_label": hub["longLabel"],
                "descriptionUrl": hub["descriptionUrl"],
                "genomes": hub["dbList"].split(",")
            })

    return hubs_info


def UCSC_hubs_metadata(hub, genome):
    url = "{}/list/tracks?hubUrl={};genome={};trackLeavesOnly=1".format(config.UCSC_API, hub["hub_url"], genome)
    track_metadata = requests.get(url = url).json()[genome]
    hub_info_string = "This file is from the UCSC public hub <a href=\"{}\">{} : {}</a>".format(hub["descriptionUrl"], hub["hub_short_label"], hub["hub_long_label"])
    # retrieve descriptions through html key
    for track, info in track_metadata.items():
        if "html" in info.keys():
            url = hub["hub_url"][:-7] + info["html"][5:-(len(track)+1)] + "/" + info["html"]
            response = requests.get(url=url)
            track_metadata[track]["short_info"] = HTML_extract_description(str(response.content) )+ "\\n" + hub_info_string 
            track_metadata[track]["long_info"] = response.content
        else:
            track_metadata[track]["short_info"] = hub_info_string 
            track_metadata[track]["long_info"] = ""


    metadata_df = pd.DataFrame(track_metadata).transpose()
    metadata_df.reset_index(level=0, inplace=True)
    metadata_df = metadata_df.rename(columns={"index":"file_name", "shortLabel":"short_name", "longLabel":"long_name"})

    return metadata_df[["file_name", "short_name", "long_name", "short_info", "long_info"]]

##########################
# LOCAL GENOME FUNCTIONS #
##########################

def setup_local(projects, conn):
    for project in projects: #genome, paths in genomes.items():
        # add genome to GENOMES
        add_genome(project["reference_genome"], "", "", "", conn)
        # add UCSC genome info to PROJECTS
        project_id = "local_{}_{}".format(project["project_name"].replace(" ","-"), project["reference_genome"].replace(" ","-"))
        conn.execute("INSERT INTO PROJECTS (PROJECTID, OSOURCE,  DSOURCE, SHORTNAME, LONGNAME, INFO) VALUES (?, ?, ?, ?, ?, ?)", 
                (project_id, project["reference_genome"], "local", "", "","",))
        # collect file information
        files_info = local_collect_file_info(project["data_path"], max_file_size=config.max_setup_file_size)
        # collect metadata for files
        metadata = local_metadata(project["metadata_path"])
        # download files, cluster, then index
        download_cluster_index(project_id, "local", project["project_name"], project["reference_genome"], files_info, metadata, conn)

def local_collect_file_info(path, max_file_size=config.max_file_size):
    file_list = glob.glob(os.path.join(path, "*.bed*"))
    files_info = {}
    for file_name in file_list:
        # record info
        file_name = file_name.split("/")[-1]
        file_type = file_name.replace(".gz","").split(".")[-1]
        track_name = file_name.split(".")[0]
        # count number of intervals per file
        interval_count = int(subprocess.check_output("wc -l {}/{}".format(path, file_name), shell=True, text=True).split()[0])
        # store file info
        if interval_count < max_file_size:
            files_info[track_name] = {"file_size": interval_count,
                                      "download_function": download_linked_file,
                                      "download_params": {"track": track_name, 
                                                          "download_location": "{}{}".format(path, file_name),
                                                          "file_type": file_type
                                                        }  
                                        }
    return files_info

def clean_file_name(file_name):
    return file_name.split("/")[-1].split(".")[0]

def local_metadata(path):
    necesary_columns = ["file_name", "short_name", "long_name", "short_info", "long_info"]
    metadata_exists = False 
    # check if metadata file exists and it is correct format
    if not os.path.exists(path) and ".csv" in path:
        metadata = pd.read_csv(path)
        # error handling if metadata doesn't meet format
        columns = list(metadata.columns)
        if sorted(columns) == sorted(necesary_columns):
            metadata_exists = True 
            metadata["file_name"] = metadata["file_name"].apply(clean_file_name)
        else:
            metadata_exists = False 
            print("Metadata file {} does not contain required columns. Please check the LOCAL_METADATA formatting guidelines in Config.py.".format(path))
    else:
        print("Matadata path not provided or incorrect")
        
    # if metadata doesn't follow standards then replace metadata with empty dataframe
    if metadata_exists == False:
        metadata = pd.DataFrame(columns = necesary_columns)
        print("Indexing will continue without metadata.")

    return metadata


if __name__ == "__main__":
    if os.path.exists("data"):
        shutil.rmtree("data")
    if os.path.exists("indices"):
        shutil.rmtree("indices")
    if os.path.exists("outputs"):
        shutil.rmtree("outputs")
    if os.path.exists("Indexing.db"):
        os.remove("Indexing.db")
    os.system('python3 models.py')  # setup indexing and files database
    conn = sqlite3.connect('Indexing.db')  # make connection to database
    
    # make directories for data and indicies 
    proc = subprocess.check_output("mkdir data/", shell=True)
    proc = subprocess.check_output("mkdir indices/", shell=True)
    proc = subprocess.check_output("mkdir outputs/", shell=True)

    # setup data sources
    setup_ucscGenomes(config.UCSC_GENOMES, conn)
    setup_ucscHubs(config.UCSC_HUBS[:-1], conn)
    setup_local(config.LOCAL_GENOMES, conn)

    conn.commit()
    conn.close()
    proc = subprocess.check_output("python3 query_indices.py -g outputs/genomes.csv",
                                    stderr=None,
                                    shell=True)