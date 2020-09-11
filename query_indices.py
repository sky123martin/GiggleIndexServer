import os.path
from os import path
import subprocess
from config import config
import sqlite3
from sys import argv
from clize import run

config = config
conn = sqlite3.connect(config.DB)

def GENOMES():
    """"returns all currently hosted genomes"""
    out = conn.execute("SELECT DISTINCT GENOME from INDICES")
    for row in out:
        print(row)


def SOURCES():
    """"returns all currently hosted sources"""
    out = conn.execute("SELECT DISTINCT SOURCE from INDICES")
    for row in out:
        print(row)


def INDICES():
    """"returns all currently hosted indices"""
    out = conn.execute("SELECT * from INDICES")
    for row in out:
        print(row)


def FILES():
    """"returns all currently hosted files"""
    out = conn.execute("SELECT * from FILES")
    for row in out:
        print(row)


def validate_source(source):
    out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = '{}'".format(source))
    source_valid = False
    for i in out:
        print("HHH")
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

def QUERY_INTERVAL(interval, source="", genome=""):
    """Query a given interval in format 'Chr:#-#' ex. 1:200457776-200457776"""
    validate_interval(interval)

    if source != "":
        validate_source(source)

    if genome != "":
        validate_genome(genome)

    if source != "" and genome != "":
        out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = '{}' AND GENOME = '{}'".format(source, genome))
    elif source != "":
        out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = '{}'".format(source))
    else:
        out = conn.execute("SELECT NAME from INDICES")

    for index in out:
        # TODO Multithread and concat results
        print("RESULTS FOR",index[0])
        index_path = "indices/" + index[0] + ".d"
        cmd_str = 'giggle search -i ' + index_path + ' -r ' + interval
        output = subprocess.check_output(cmd_str, shell=True)
        print(output)
    
    print("QUERY_INTERVAL", interval, genome)

def validate_query_file(path_):
    # validate if path exists
    if not path.exists(path_):
        print("ERROR: Path given, {} , does not exist".format(path_))
        raise
    
    # validate if path leads to a file
    if not path.isfile(path_):
        print("ERROR: Path given, {} , is not a file".format(path_))
        raise
    
    # GIGGLE accepts query files in either bed.gz or vcf.gz
    p = path_.split(".")
    if not ((p[-2] == "bed" or p[-2] == "vcf") and p[-1] == "gz"):
         print("ERROR: File given, {} , is not in format bed.gz or vcf.gz".format(path_.split("/")[-1]))
         raise


def QUERY_FILE(path, source="", genome=""):
    """Query a given file given a path"""
    
    validate_query_file(path)

    if source != "":
        validate_source(source)

    if genome != "":
        validate_genome(genome)

    if source != "" and genome != "":
        out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = '{}' AND GENOME = '{}'".format(source, genome))
    elif source != "":
        out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = '{}'".format(source))
    else:
        out = conn.execute("SELECT NAME from INDICES")

    for index in out:
        # TODO Multithread and concat results
        print("RESULTS FOR",index[0])
        index_path = "indices/" + index[0] + ".d"
        cmd_str = 'giggle search -i ' + index_path + ' -q ' + path + ' -s'
        output = subprocess.check_output(cmd_str, shell=True)
        print(output)


if __name__ == "__main__":
    run(FILES, alt={
    'i': INDICES,
    'f': FILES,
    'g': GENOMES,
    's': SOURCES,
    'qf': QUERY_FILE,
    'qi': QUERY_INTERVAL,
    })