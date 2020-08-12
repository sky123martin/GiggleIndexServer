from application.models import Files, Index
import mysql.connector
import pandas as pd
from contextlib import contextmanager
from application import app, db
import requests 
import glob, os, os.path
import pybedtools
from alive_progress import alive_bar

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
        bar("Clustering genome {}".format(genome))
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


def UCSC_collect_file_info(genome):
    tracks = requests.get(url="{}/list/tracks?genome={};trackLeavesOnly=1".format(app.config["UCSC_API"],genome)).json()[genome]
    files_info = {}
    with alive_bar(len(tracks), bar = 'bubbles', spinner = 'classic') as bar:
        for track, info in tracks.items():
            bar("Collecting file info for {}: {}".format(genome, track))
            if "table" in info:
                    track = info["table"]
            if "bigDataUrl" in info: # if stored as a big data file
                files_info[track] = {"file_size": info["itemCount"]}
                files_info[track]["download_function"] = UCSC_get_big_data_file
                files_info[track]["download_params"] = [track, genome, info["bigDataUrl"], info["type"]]
            else: # if stored in sql db
                with connect_SQL_db(app.config["UCSC_SQL_DB_HOST"], "genome") as db:
                    columns = list(pd.read_sql("Show columns from {}.{}".format(genome, track), con=db)["Field"])
                    if len(columns) > 2: # Some files just don't have any info on them
                        files_info[track] = {"file_size": info["itemCount"]}
                        extract_columns = target_columns(columns, info["type"])
                        files_info[track]["download_function"] = UCSC_download_sql_file
                        files_info[track]["download_params"] = [track, genome, extract_columns]
    return files_info


def UCSC_download_sql_file(index, params):
    """ downloads a specific file and puts in index folder
    Parameters
    ----------
    index: string
        ex. hg19_1
    params: list
        [track, genome, extract_columns]
        extract_columns = [chrom, start, end, base_shift]
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

    with connect_SQL_db(app.config["UCSC_SQL_DB_HOST"], "genome") as db:
        df = pd.read_sql("Select {} AS chrom, {} + {} AS start, {} + {} AS end from {}.{}".format(chrom, start, base_shift, end, base_shift, genome, track), con=db)
        pybedtools.BedTool.from_dataframe(df).saveas('data/{}/{}.bed'.format(index, track))


def UCSC_get_big_data_file(index, params):
    pass
    # print(index_name, track_name)

def setup_index(source, genome, index, index_info):
    with alive_bar(len(index_info["files"])+3, bar = 'blocks', spinner = 'classic') as bar:
        # create index directory
        try:
            os.mkdir("data/" + index)
        except OSError:
            print("Creation of the directory %s failed" % index)
        # add index info to the Index database
        db.session.add(Index(id = index, source = source, genome = genome, size = index_info["index_size"], full = index_info["full"]))
        # iterate through all files belonging to the index and download them
        for filename, file_info in index_info["files"].items():
            bar("Creating {}: Downloading {}".format(index, filename))
            file_info["download_function"](index, file_info["download_params"])
            db.session.add(Files(tablename=filename, source="UCSC", genome=genome, size=file_info["file_size"], index_id=index))
        # FIXME index downloaded files here
        bar("Creating {}: indexing files".format(index))
        if index_info["full"]:  # delete files if index is full
            bar("Completed {}: deleting files".format(index))
            filelist = glob.glob(os.path.join("data/"+index, "*.bed"))
            for f in filelist:
                os.remove(f)
            os.rmdir("data/"+index)
        bar("Index {} Completed".format(index))
        print("Index {} Completed".format(index))
        db.session.commit()

                 
def setup_UCSC_indices(genomes):
    for genome in genomes:
        files_info = UCSC_collect_file_info(genome)
        clustered_files_info = cluster_data("UCSC", genome, files_info)
        # Iterate through each index
        for index, index_info in clustered_files_info.items():
            setup_index("UCSC", genome, index, index_info)

        print(Files.query.all())
        print(Index.query.all())


def setup_indices():
    """ Sets up all sources and genomes listed in CONFIG
    """
    if len(app.config["UCSC_GENOMES"]) > 0:
        setup_UCSC_indices(app.config["UCSC_GENOMES"])

