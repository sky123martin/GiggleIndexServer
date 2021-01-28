import os.path
from os import path
import subprocess
from config import config
import sqlite3
from sys import argv
from clize import run
from multiprocessing.pool import Pool
import multiprocessing 
from functools import partial
import pandas as pd

config = config
conn = sqlite3.connect(config.DB)

def GENOMES():
    """"returns all currently hosted genomes"""
    cursor = conn.execute("SELECT DISTINCT GENOME from INDICES")
    df = pd.DataFrame(cursor.fetchall())
    df.columns = [i[0] for i in cursor.description]
    df.to_csv("output.csv", index=False)


def SOURCES():
    """"returns all currently hosted sources"""
    cursor = conn.execute("SELECT DISTINCT SOURCE from INDICES")
    df = pd.DataFrame(cursor.fetchall())
    df.columns = [i[0] for i in cursor.description]
    df.to_csv("output.csv", index=False)


def INDICES():
    """"returns all currently hosted indices"""
    cursor = conn.execute("SELECT * from INDICES")
    df = pd.DataFrame(cursor.fetchall())
    df.columns = [i[0] for i in cursor.description]
    df.to_csv("output.csv", index=False)


def FILES(genome="", output_path=""):
    """"returns currently hosted files information, <optional param: genome> <optional param: output_path>"""
    if genome == "":
        cursor = conn.execute("SELECT * from FILES")
    else:
        validate_genome(genome)
        cursor = conn.execute("SELECT * from FILES WHERE GENOME='{}'".format(genome))

    df = pd.DataFrame(cursor.fetchall())
    df.columns = [i[0].replace(" ", "") for i in cursor.description]
    
    if output_path == "":
        df.to_csv("output.csv", index=False)
    else:
        df.to_csv(output_path, index=False)


def validate_source(source):
    out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = '{}'".format(source))
    source_valid = False
    for i in out:
        source_valid = True
        break

    if not source_valid:
        print("ERROR: Source {} not found in system, use -s to find valid sources".format(source))
        raise


def validate_genome(genome):
    out = conn.execute("SELECT NAME from INDICES WHERE GENOME = '{}'".format(genome))
    genome_valid = False
    for i in out:
        genome_valid = True
        break

    if not genome_valid:
        print("ERROR: Genome {} not found in system , use -g to find valid genomes".format(genome))
        raise


def validate_interval(interval):
    if len(interval.split("-")) != 2 or len(interval.split(":")) != 2:
        print("ERROR: Interval not in correct format Chr:#-# ex 1:200457776-200457776")
        raise
    try:
        chr = interval.split(":")[0]
        LB = interval.split(":")[1].split("-")[0]
        UP = interval.split(":")[1].split("-")[1]
    except:
        print("ERROR: Interval not in correct format Chr:#-# ex 1:200457776-200457776")
        raise


def query_interval_index(interval, index):
    index_path = "indices/" + index + ".d"
    cmd_str = 'giggle search -i ' + index_path + ' -r ' + interval
    search_out = subprocess.check_output(cmd_str, shell=True)
    cleaned_out = ""
    for x in str(search_out).split('#')[2:]:
        formatted_string = x[:-2].replace('\\t', ", ").replace('data/', "").replace("\\", "")
        cleaned_out += formatted_string.replace("size:","").replace("overlaps:","") + "\n"

    return cleaned_out


def QUERY_INTERVAL(interval, genome, output_path=""):
    """Query a given interval in format 'Chr:#-#' <genome> <optional param: output_path>"""
    validate_interval(interval)
    validate_genome(genome)

    out = conn.execute("SELECT NAME from INDICES WHERE GENOME = '{}'".format(genome))

    indices = list([i[0] for i in out])  # reformat sql results into list

    # Multiproccesing used to query indices
    with Pool(config.AVAILABLE_PROCCESSES) as p:  # to check multiprocessing.cpu_count()
        output = p.map(partial(query_interval_index, interval), indices)

    with open("output.txt" if output_path == "" else output_path, "w") as f:
    # Format output
        f.write("file,size,overlaps\n")
        for result in output:
            f.write(result)


def validate_query_file(path_):
    # validate if path exists
    if not path.exists(path_):
        print("ERROR: Path given, {} , does not exist".format(path_))
        raise ImportError
    
    # validate if path leads to a file
    if not path.isfile(path_):
        print("ERROR: Path given, {} , is not a file".format(path_))
        raise ImportError
    
    # GIGGLE accepts query files in either bed.gz or vcf.gz
    p = path_.split(".")
    if not ((p[-2] == "bed" or p[-2] == "vcf") and p[-1] == "gz"):
         print("ERROR: File given, {} , is not in format bed.gz or vcf.gz".format(path_.split("/")[-1]))
         raise TypeError


def query_file_index(path, index):
    index_path = "indices/" + index + ".d"
    cmd_str = 'giggle search -i ' + index_path + ' -q ' + path + " -s"
    search_out = str(subprocess.check_output(cmd_str, shell=True))
    cleaned_out = ""
    # delete header words in query file
    rows = search_out.replace("b'#","").replace("","").split("\\n")
    # headers are in rows[0]
    cleaned_rows = [i.replace("\\t",",")+"\n" for i in rows]
    cleaned_rows = [i.replace(",\n","\n") for i in cleaned_rows]
    return cleaned_rows


def QUERY_FILE(path, genome, output_path=""):
    """Query a given file given a path. ie: <path> <genome> <optional param: output_path>"""
    
    validate_query_file(path)
    validate_genome(genome)

    out = conn.execute("SELECT NAME from INDICES WHERE GENOME = '{}'".format(genome))

    indices = list([i[0] for i in out])

    with Pool(config.AVAILABLE_PROCCESSES) as p: # to check multiprocessing.cpu_count()
        output = p.map(partial(query_file_index, path), indices)

    with open("output.txt" if output_path == "" else output_path, "w") as f:
        # Format output
        f.write(output[0][0])
        for result in output:
            for row in result[1:-1]:
                # print(row)
                f.write(row)


if __name__ == "__main__":
    run(FILES, alt={
                'i': INDICES,
                'f': FILES,
                'g': GENOMES,
                's': SOURCES,
                'qf': QUERY_FILE,
                'qi': QUERY_INTERVAL,
                })