from config import config
import pandas as pd
import subprocess
from setup_indices import *
import sqlite3
import glob
import os
import shutil

config = config

def update_ucscGenomes(genomes, conn):
    print(genomes)
    # check if project has already been setup, else setup
    setup_projects = [i[0] for i in conn.execute("select DISTINCT PROJECTID from PROJECTS  where DSOURCE='ucscGenomes'")]
    config_projects = ["ucscGenomes_" + g.replace(" ","-") for g in genomes]
    print("HERE", config_projects, setup_projects)
    for genome in genomes:
        project_id = "ucscGenomes_" + genome.replace(" ","-")
        if project_id in setup_projects: # project has already been setup, check for new files
            # scrap file info for genome
            scrapped_files = collect_ucsc_files(genome)
            if isinstance(scrapped_files, str):  # 404 not found on API
                print("Unable to collect file info from UCSC for {} ref genome".format(genome))
                continue
            # collect files that have been indexed
            setup_files = [i[0] for i in conn.execute("select FILEID from FILES where PROJECTID='{}'".format(project_id))]
            # scan for new files
            new_files = {}
            for f in scrapped_files.keys():
                if f not in setup_files:
                    new_files[f] = scrapped_files[f]
            # if no new files are found then move to next genome
            if len(new_files) == 0:
                continue
            # if new files then find an unfilled index
            unfull_index_iter = [i[0] for i in conn.execute("select MIN(ITER) from INDICES where PROJECTID='{}' and FULL=0".format(project_id))][0]
            unfull_files = [i[0] for i in conn.execute("select FILEID from FILES as f left join INDICES as i on f.INDEXID=i.INDEXID where f.PROJECTID='{}' and i.FULL=0".format(project_id))]
            # create list of files being indexed
            for f in unfull_files:
                new_files[f] = scrapped_files[f]
            # collect metadata for files
            metadata = UCSC_metadata(genome)
            # only setup files that have metadata, take intersection of metadata & files_info
            new_files = {key: new_files[key] for key in list(set(new_files.keys()) & set(metadata['file_name'].to_numpy()))}
            if len(new_files)==0:
                continue
            # delete SQL records for indicies being reindexed
            conn.execute("DELETE from INDICES where PROJECTID='{}' and FULL=0".format(project_id))
            conn.execute("DELETE from FILES where INDEXID IN (select INDEXID from INDICES where PROJECTID='{}' and FULL=0)".format(project_id))
            download_cluster_index(project_id, "ucscGenomes", genome.replace(" ","-"), genome, new_files, metadata, conn, current_index_num=unfull_index_iter)

        elif project_id not in setup_projects: # project has not been setup
            setup_ucscGenomes([genome], conn)

    # if project is no longer in config then delete
    for setup_id in setup_projects:
        genome = setup_id.split("_",1)[-1]
        if setup_id not in config_projects:
            # delete indices
            files = [i[0] for i in conn.execute("select INDEXID from INDICES where PROJECTID='{}'".format(setup_id))]
            for f in files:
                shutil.rmtree("indices/{}.d".format(f))

            # delete sql metadata
            conn.execute("DELETE FROM PROJECTS WHERE PROJECTID = '{}'".format(setup_id))
            conn.execute("DELETE FROM INDICES WHERE PROJECTID = '{}'".format(setup_id))
            conn.execute("DELETE FROM FILES WHERE PROJECTID = '{}'".format(setup_id))

            # delete genome if there is no other project with that genome
            cursor = conn.execute("SELECT GENOME from GENOMES WHERE GENOMES.GENOME = '{}' LIMIT 1".format(genome))
            df = pd.DataFrame(cursor.fetchall())
            if df.empty:
                conn.execute("DELETE FROM GENOMES WHERE GENOME = '{}'".format(genome))

def update_ucscHubs(hub_names, conn):
    # check if project has already been setup, else setup
    setup_projects = [i[0] for i in conn.execute("select DISTINCT PROJECTID from PROJECTS where DSOURCE='ucscHubs'")]
    config_projects = []
    print(hub_names, setup_projects)
    hubs = collect_ucscHubs(hub_names)
    for hub in hubs:
        for genome in hub["genomes"]:
            project_id = "ucscHubs_{}_{}".format(hub["hub_short_label"].replace(" ","-").replace("_","-"), genome)
            print(project_id)
            config_projects.append(project_id)
            if genome == "":
                continue
            # add geno
            if project_id in setup_projects: # project has already been setup, check for new files
                # scrap file info for genome
                scrapped_files = collect_ucsc_files(genome, HUB_EXT="hubUrl={};".format(hub["hub_url"]))
                if isinstance(scrapped_files, str):  # 404 not found on API
                    print("Unable to collect file info from ucscHub {}".format(project_id))
                    continue
                # collect files that have been indexed
                setup_files = [i[0] for i in conn.execute("select FILEID from FILES where PROJECTID='{}'".format(project_id))]
                # scan for new files
                new_files = {}
                for f in scrapped_files.keys():
                    if f not in setup_files:
                        new_files[f] = scrapped_files[f]
                # if no new files are found then move to next genome
                if len(new_files) == 0:
                    print("no files")
                    continue
                # if new files then find an unfilled index
                unfull_index_iter = [i[0] for i in conn.execute("select MIN(ITER) from INDICES where PROJECTID='{}' and FULL=0".format(project_id))][0]
                unfull_files = [i[0] for i in conn.execute("select FILEID from FILES as f left join INDICES as i on f.INDEXID=i.INDEXID where f.PROJECTID='{}' and i.FULL=0".format(project_id))]
                # delete SQL records for indicies being reindexed
                conn.execute("DELETE from INDICES where PROJECTID='{}' and FULL=0".format(project_id))
                conn.execute("DELETE from FILES where INDEXID IN (select INDEXID from INDICES where PROJECTID='{}' and FULL=0)".format(project_id))
                # create list of files being indexed
                for f in unfull_files:
                    new_files[f] = scrapped_files[f]
                print("HERE", new_files)

                # collect metadata for files
                metadata = UCSC_hubs_metadata(hub, genome)

                download_cluster_index(project_id, "ucscHubs", hub["hub_short_label"], genome, new_files, metadata, conn)

            elif project_id not in setup_projects: # project has not been setup
                print("NEW INDEX")
                setup_ucscHubs([hub["hub_short_label"]], conn)
                break

    # if project is no longer in config then delete
    for setup_id in setup_projects:
        genome = setup_id.split("_",1)[-1]
        if setup_id not in config_projects:
            # delete indices
            files = [i[0] for i in conn.execute("select INDEXID from INDICES where PROJECTID='{}'".format(setup_id))]
            for f in files:
                shutil.rmtree("indices/{}.d".format(f))

            # delete sql metadata
            conn.execute("DELETE FROM PROJECTS WHERE PROJECTID = '{}'".format(setup_id))
            conn.execute("DELETE FROM INDICES WHERE PROJECTID = '{}'".format(setup_id))
            conn.execute("DELETE FROM FILES WHERE PROJECTID = '{}'".format(setup_id))

            # delete genome if there is no other project with that genome
            cursor = conn.execute("SELECT GENOME from GENOMES WHERE GENOMES.GENOME = '{}' LIMIT 1".format(genome))
            df = pd.DataFrame(cursor.fetchall())
            if df.empty:
                conn.execute("DELETE FROM GENOMES WHERE GENOME = '{}'".format(genome))


def update_local(projects, conn):
    # check if project has already been setup, else setup
    setup_projects = [i[0] for i in conn.execute("select DISTINCT PROJECTID from PROJECTS where DSOURCE='local'")]
    config_projects = ["local_{}_{}".format(p["project_name"].replace(" ","-"), p["reference_genome"].replace(" ","-")) for p in projects]

    for project in projects:
        project_id = "local_{}_{}".format(project["project_name"].replace(" ","-"), project["reference_genome"].replace(" ","-"))
        if project_id in setup_projects: # project has already been setup, check for new files
            # scrap file info for genome
            scrapped_files = local_collect_file_info(project["data_path"])
            if isinstance(scrapped_files, str):  # 404 not found on API
                print("Unable to collect file info from local for {} project".format(project["project_name"]))
                continue
            # collect files that have been indexed
            setup_files = [i[0] for i in conn.execute("select FILEID from FILES where PROJECTID='{}'".format(project_id))]
            # scan for new files
            new_files = {}
            for f in scrapped_files.keys():
                if f not in setup_files:
                    new_files[f] = scrapped_files[f]
            # if no new files are found then move to next genome
            if len(new_files) == 0:
                continue
            # if new files then find an unfilled index
            unfull_index_iter = [i[0] for i in conn.execute("select MIN(ITER) from INDICES where PROJECTID='{}' and FULL=0".format(project_id))][0]
            unfull_files = [i[0] for i in conn.execute("select FILEID from FILES as f left join INDICES as i on f.INDEXID=i.INDEXID where f.PROJECTID='{}' and i.FULL=0".format(project_id))]
            # create list of files being indexed
            for f in unfull_files:
                new_files[f] = scrapped_files[f]
            # collect metadata for files
            metadata = UCSC_metadata(genome)
            # only setup files that have metadata, take intersection of metadata & files_info
            new_files = {key: new_files[key] for key in list(set(new_files.keys()) & set(metadata['file_name'].to_numpy()))}
            if len(new_files)==0:
                continue
            # delete SQL records for indicies being reindexed
            conn.execute("DELETE from INDICES where PROJECTID='{}' and FULL=0".format(project_id))
            conn.execute("DELETE from FILES where INDEXID IN (select INDEXID from INDICES where PROJECTID='{}' and FULL=0)".format(project_id))
            download_cluster_index(project_id, "local", project["project_name"], project["reference_genome"], new_files, metadata, conn, current_index_num=unfull_index_iter)

        elif project_id not in setup_projects: # project has not been setup
            setup_local([project], conn)

    # if project is no longer in config then delete
    for setup_id in setup_projects:
        genome = setup_id.split("_",1)[-1]
        if setup_id not in config_projects:
            # delete indices
            files = [i[0] for i in conn.execute("select INDEXID from INDICES where PROJECTID='{}'".format(setup_id))]
            for f in files:
                shutil.rmtree("indices/{}.d".format(f))

            # delete sql metadata
            conn.execute("DELETE FROM PROJECTS WHERE PROJECTID = '{}'".format(setup_id))
            conn.execute("DELETE FROM INDICES WHERE PROJECTID = '{}'".format(setup_id))
            conn.execute("DELETE FROM FILES WHERE PROJECTID = '{}'".format(setup_id))

            # delete genome if there is no other project with that genome
            cursor = conn.execute("SELECT GENOME from GENOMES WHERE GENOMES.GENOME = '{}' LIMIT 1".format(genome))
            df = pd.DataFrame(cursor.fetchall())
            if df.empty:
                conn.execute("DELETE FROM GENOMES WHERE GENOME = '{}'".format(genome))


if __name__ == "__main__":
    conn = sqlite3.connect(config.DB)

    update_ucscGenomes(config.UCSC_GENOMES, conn)
    update_ucscHubs(config.UCSC_HUBS[:-1], conn)
    update_local(config.LOCAL_GENOMES, conn)

    # delete any projects that are maintained but no longer listed in the config
    conn.commit()
    conn.close()

    proc = subprocess.check_output("python3 query_indices.py -g outputs/genomes.csv",
                                    stderr=None,
                                    shell=True)