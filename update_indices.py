from config import config
from setup_indices import local_collect_file_info, UCSC_collect_file_info, cluster_data, setup_indices
import sqlite3
import glob
import os

config = config

def update_UCSC(genome):
    files_info = UCSC_collect_file_info(config.LOCAL_GENOMES[genome])
    files_out = conn.execute("SELECT NAME, INDEXNAME from FILES WHERE SOURCE = 'UCSC' AND GENOME = '{}'".format(genome))
    re_indexed_files = {}
    
    # Determine if there is an unfull index that could be deleted and recompiled with new files
    indices_out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = 'UCSC' AND GENOME = '{}' and FULL = False ".format(genome))
    
    open_index = ""
    for row in indices_out:
        open_index = row[0]
    
    # Determine which files are new
    num_new_files = len(files_info)
    for row in files_out:
        file_name = row[0]
        index_name = row[1]
        if file_name in files_info.keys():
            num_new_files -= 1
            if index_name != open_index:
                del files_info[file_name]

    if open_index != "" and num_new_files > 0:
        # delete stored index
        try:
            index_list = glob.glob(os.path.join("indices/"+open_index+".d", "*"))
            for f in index_list :
                os.remove(f)
            os.rmdir("indices/"+open_index+".d")
        except OSError:
            print("folder {} not found".format("indices/"+open_index))
        # delete current sql db logs of files that are being re-indexing
        conn.execute("DELETE FROM FILES WHERE INDEXNAME = '{}'".format(open_index))
        conn.execute("DELETE FROM INDICES WHERE NAME = '{}'".format(open_index))

    print("HEKDJDJ", num_new_files)

    if num_new_files > 0:
        # Re cluster and index files
        clustered_files_info = cluster_data("local", genome, files_info, open_index.split(".")[-1])
        for index, index_info in clustered_files_info.items():
                setup_indices("local", genome, index, index_info, conn)



def update_LOCAL(genome):
    files_info = local_collect_file_info(config.LOCAL_GENOMES[genome])
    files_out = conn.execute("SELECT NAME, INDEXNAME from FILES WHERE SOURCE = 'local' AND GENOME = '{}'".format(genome))
    re_indexed_files = {}
    
    # Determine if there is an unfull index that could be deleted and recompiled with new files
    indices_out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = 'local' AND GENOME = '{}' and FULL = False ".format(genome))
    
    open_index = ""
    for row in indices_out:
        open_index = row[0]
    
    # Determine which files are new
    num_new_files = len(files_info)
    for row in files_out:
        file_name = row[0]
        index_name = row[1]
        if file_name in files_info.keys():
            num_new_files -= 1
            if index_name != open_index:
                del files_info[file_name]

    if open_index != "" and num_new_files > 0:
        # delete stored index
        try:
            index_list = glob.glob(os.path.join("indices/"+open_index+".d", "*"))
            for f in index_list :
                os.remove(f)
            os.rmdir("indices/"+open_index+".d")
        except OSError:
            print("folder {} not found".format("indices/"+open_index))
        # delete current sql db logs of files that are being re-indexing
        conn.execute("DELETE FROM FILES WHERE INDEXNAME = '{}'".format(open_index))
        conn.execute("DELETE FROM INDICES WHERE NAME = '{}'".format(open_index))

    print("HEKDJDJ", num_new_files)

    if num_new_files > 0:
        # Re cluster and index files
        clustered_files_info = cluster_data("local", genome, files_info, open_index.split(".")[-1])
        for index, index_info in clustered_files_info.items():
                setup_indices("local", genome, index, index_info, conn)


def update(source, genome):
    print("HELLO")
    if source == "UCSC":
        update_UCSC(genome)
    elif source == "local":
        print("HELLO")
        update_LOCAL(genome)


if __name__ == "__main__":
    global conn
    conn = sqlite3.connect(config.DB)
    out = conn.execute("SELECT DISTINCT SOURCE, GENOME from FILES")
    for row in out:
        source = row[0]
        genome = row[1]
        update(source, genome)
        print("Fff")

        conn.commit()
    conn.close()