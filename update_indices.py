from config import config
import pandas as pd
from setup_indices import local_collect_file_info, UCSC_collect_file_info, cluster_data, setup_indices, local_metadata, UCSC_metadata, UCSC_hubs_metadata, UCSC_collect_hubs
import sqlite3
import glob
import os

config = config

def update_UCSC(genome):
    newly_scraped_files = UCSC_collect_file_info(genome)
    indexed_files = conn.execute("SELECT NAME, INDEXNAME from FILES WHERE SOURCE = 'UCSC' AND PROJECT = 'UCSC Genomes' AND GENOME = '{}'".format(genome))
    
    # Determine if there is an unfull index that could be deleted and recompiled with new files
    indices_out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = 'UCSC' AND PROJECT = 'UCSC Genomes' AND GENOME = '{}' and FULL = False ".format(genome))
    
    open_index = ""
    for row in indices_out:
        open_index = row[0]
    
    # Determine which files are new
    num_new_files = len(newly_scraped_files)
    for row in indexed_files:
        file_name = row[0]
        index_name = row[1]
        if file_name in newly_scraped_files.keys():
            num_new_files -= 1
            if index_name != open_index:
                del newly_scraped_files[file_name]

    if num_new_files > 0:
        if open_index != "":
            # delete stored index
            try:
                index_list = glob.glob(os.path.join("indices/"+open_index+".d", "*"))
                for f in index_list :
                    os.remove(f)
                os.rmdir("indices/" + open_index + ".d")
            except OSError:
                print("folder {} not found".format("indices/"+open_index))
            # delete current sql db logs of files that are being re-indexing
            conn.execute("DELETE FROM FILES WHERE INDEXNAME = '{}' AND  SOURCE = 'UCSC' AND PROJECT = 'UCSC Genomes' AND GENOME = '{}'".format(open_index, genome))
            conn.execute("DELETE FROM INDICES WHERE NAME = '{}' AND  SOURCE = 'UCSC' AND PROJECT = 'UCSC Genomes' AND GENOME = '{}' ".format(open_index, genome))

        # Re cluster and index files
        clustered_files_info = cluster_data("local", genome, newly_scraped_files, open_index.split("_")[-1])
        # collect metadata for files
        metadata = UCSC_metadata(genome)

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
            setup_indices("UCSC", "UCSC Genomes", genome, index, index_info, index_metadata, conn)

def update_UCSC_HUBS(project, genome):
    
    # collect hub info
    hubs = requests.get(url = config.UCSC_API + "/list/publicHubs").json()["publicHubs"]
    for hub in hubs:
        if hub["shortLabel"] == project:
            break

    indexed_files = conn.execute("SELECT NAME, INDEXNAME from FILES WHERE SOURCE = 'UCSC' AND PROJECT = '{}' AND GENOME = '{}'".format(project, genome))
    
    # Determine if there is an unfull index that could be deleted and recompiled with new files
    indices_out = conn.execute("SELECT NAME, HUBEXT from INDICES WHERE SOURCE = 'UCSC' AND PROJECT = '{}' AND GENOME = '{}' AND FULL = False ".format(project, genome))
    
    open_index = ""
    for row in indices_out:
        open_index = row[0]
    
    newly_scraped_files = UCSC_collect_file_info(genome, "hubUrl={};".format(hub_ext)

    # Determine which files are new
    num_new_files = len(newly_scraped_files)
    for row in indexed_files:
        file_name = row[0]
        index_name = row[1]
        if file_name in newly_scraped_files.keys():
            num_new_files -= 1
            if index_name != open_index:
                del newly_scraped_files[file_name]

    if num_new_files > 0:
        if open_index != "":
            # delete stored index
            try:
                index_list = glob.glob(os.path.join("indices/"+open_index+".d", "*"))
                for f in index_list:
                    os.remove(f)
                os.rmdir("indices/" + open_index + ".d")
            except OSError:
                print("folder {} not found".format("indices/"+open_index))
            # delete current sql db logs of files that are being re-indexing
            conn.execute("DELETE FROM FILES WHERE INDEXNAME = '{}' AND SOURCE = 'UCSC' AND PROJECT = '{}' AND GENOME = '{}'".format(open_index, project, genome))
            conn.execute("DELETE FROM INDICES WHERE NAME = '{}' AND SOURCE = 'UCSC' AND PROJECT = '{}' AND GENOME = '{}'".format(open_index, project, genome))

            # cluster files based on hyperparam in .config
            clustered_files_info = cluster_data(hub["hub_short_label"].replace(" ",""), genome, files_info)
            # iterate through each cluster then download and index cluster files
            metadata = UCSC_hubs_metadata(hub, genome)

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
                setup_indices("UCSC", hub["hub_short_label"], genome, index, index_info, metadata, conn)


def update_LOCAL(project_name, genome):
    for project in config.LOCAL_GENOMES:
        if project["project_name"] == project_name:
            break
    files_info = local_collect_file_info(project["data_path"])
    files_out = conn.execute("SELECT NAME, INDEXNAME from FILES WHERE SOURCE = 'LOCAL' AND PROJECT = '{}' AND GENOME = '{}'".format(project, genome))
    
    # Determine if there is an unfull index that could be deleted and recompiled with new files
    indices_out = conn.execute("SELECT NAME from INDICES WHERE SOURCE = 'LOCAL' AND PROJECT = '{}' AND GENOME = '{}' AND FULL = False ".format(project, genome))
    
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

    if num_new_files > 0:
        if open_index != "":
            # delete stored index
            try:
                index_list = glob.glob(os.path.join("indices/"+open_index+".d", "*"))
                for f in index_list :
                    os.remove(f)
                os.rmdir("indices/"+open_index+".d")
            except OSError:
                print("folder {} not found".format("indices/"+open_index))
            # delete current sql db logs of files that are being re-indexing
            conn.execute("DELETE FROM FILES WHERE INDEXNAME = '{}' AND SOURCE = 'LOCAL' AND PROJECT = '{}' AND GENOME = '{}'".format(open_index, project, genome))
            conn.execute("DELETE FROM INDICES WHERE NAME = '{}' AND SOURCE = 'LOCAL' AND PROJECT = '{}' AND GENOME = '{}'".format(open_index, project, genome))

        # Re cluster and index files
        # collect metadata for files
        metadata = local_metadata(project["metadata_path"])
        # cluster files based on hyperparam in .config
        clustered_files_info = cluster_data(project["project_name"], project["reference_genome"], files_info, open_index.split("_")[-1])
        # iterate through each cluster then download and index cluster files

        for index, index_info in clustered_files_info.items():
            files = list(index_info["files"].keys())
            index_metadata = pd.merge(metadata, 
                                    pd.DataFrame([[i] for i in files], columns = ['file_name']),
                                    on = "file_name",
                                    how = "right")
            index_metadata.fillna("")
            
            setup_indices("LOCAL", project["project_name"], project["reference_genome"], index, index_info, index_metadata, conn)


def update(source, project, genome):
    if source == "UCSC" and project == "UCSC Genomes": # UCSC Genome
        update_UCSC(genome)
    elif source == "UCSC" and project != "UCSC Genomes": # UCSC hub
        update_UCSC_HUBS(project, genome)
    elif source == "local": # local project
        update_LOCAL(genome)
    else:
        print("No compatible update function {}, {}, {}".format(source, project, genome))


if __name__ == "__main__":
    global conn
    conn = sqlite3.connect(config.DB)
    out = conn.execute("SELECT DISTINCT SOURCE, PROJECT, GENOME from FILES")
    for row in out:
        source = row[0]
        project = row[1]
        genome = row[2]
        update(source, project, genome)

        conn.commit()
    conn.close()