# project/test_utility.py

from application import app, utility

from contextlib import contextmanager
import pandas as pd
import unittest
import mysql.connector
import requests
import random


class BasicTests(unittest.TestCase):

    ######################
    # setup and teardown #
    ######################

    # executed prior to each test
    def setUp(self):
        # app.config['TESTING'] = True
        # app.config['DEBUG'] = False
        # self.app = app.test_client()
        pass

    # executed after each test
    def tearDown(self):
        pass

    #############################
    # UCSC API Structure Tests #
    #############################
    # If any of these tests break try to check and see if UCSC has
    # modified API structure.
    def test_UCSC_API_genomes(self):
        # testing if genome directory exists
        r = requests.get('https://api.genome.ucsc.edu/list/ucscGenomes')
        self.assertEqual(r.status_code, requests.codes.ok)
        r = requests.get('https://api.genome.ucsc.edu/list/fakeURL')
        self.assertEqual(r.status_code, 400)

    def test_UCSC_API_track_structure(self):
        # testing if track structure holds
        genome = "rn6"
        r = requests.get('https://api.genome.ucsc.edu/list/tracks?genome={};trackLeavesOnly=1'.format("rn6"))
        self.assertEqual(r.status_code, requests.codes.ok)
        rn6_leaves = r.json()["rn6"]
        r = requests.get('https://api.genome.ucsc.edu/list/tracks?genome={}'.format("rn6"))
        rn6_full = r.json()["rn6"]
        self.assertTrue(len(rn6_leaves) > len(rn6_full))

    def test_UCSC_API_all_tracks(self):
        # testing that /genomes listed can link to /tracks
        r = requests.get('https://api.genome.ucsc.edu/list/ucscGenomes')
        self.assertEqual(r.status_code, requests.codes.ok)
        all_genomes = r.json()["ucscGenomes"]
        # Only chosing acouple to check because hg19 delay, normal delay
        for key in random.choices(list(all_genomes.keys()), k=3):
            track = requests.head('https://api.genome.ucsc.edu/list/tracks?genome={}'.format(key), timeout=45.00)
            self.assertEqual(track.status_code, requests.codes.ok)

    #############################
    # UCSC MySQL Database Tests #
    #############################
    def test_false_SQL_connection(self):
        # test that false connection throws connection error
        with self.assertRaises(mysql.connector.Error):
            host = "fake"
            user = "password"
            with utility.connect_SQL_db(host, user):
                pass

    def test_UCSC_SQL_connection(self):
        # test that current UCSC connection works
        host = "genome-mysql.soe.ucsc.edu"
        user = "genome"
        with utility.connect_SQL_db(host, user) as db:
            mycursor = db.cursor()
            mycursor.execute("SHOW DATABASES")
            self.assertTrue(mycursor)

    def test_UCSC_SQL_databases(self):
        # test thatincorrect query throws
        host = "genome-mysql.soe.ucsc.edu"
        user = "genome"
        with utility.connect_SQL_db(host, user) as db:
            mycursor = db.cursor()
            mycursor.execute("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA\
                              WHERE SCHEMA_NAME = 'hg19'")
        self.assertTrue(mycursor)

    ######################
    # Data Scapper Tests #
    ######################
    def test_target_columns(self):
        # Local test of target_columns
        file_type = "chain ailMel1"
        columns = ['bin', 'score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStrand', 'qStart', 'qEnd', 'id']
        expected_result = ['tName', 'tStart', 'tEnd', 0]
        result = utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "netAlign panTro5 chainPanTro5"
        columns = ['bin', 'level', 'tName', 'tStart', 'tEnd', 'strand', 'qName', 'qStart', 'qEnd', 'chainId', 'ali', 'score', 'qOver', 'qFar', 'qDup', 'type', 'tN', 'qN', 'tR', 'qR', 'tNewR', 'qNewR', 'tOldR', 'qOldR', 'tTrf', 'qTrf']
        expected_result = ['tName', 'tStart', 'tEnd', 0]
        result = utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "bed 4 +"
        columns = ['bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'period', 'copyNum', 'consensusSize', 'perMatch', 'perIndel', 'score', 'A', 'C', 'G', 'T', 'entropy', 'sequence']
        expected_result = ['chrom', 'chromStart', 'chromEnd', 0]
        result = utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "bed 4"
        columns = ['chrom', 'chromStart', 'chromEnd', 'name']
        expected_result = ['chrom', 'chromStart', 'chromEnd', 0]
        result = utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "psl est"
        columns = ['bin', 'matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts']
        expected_result = ['tName', 'tStart', 'tEnd', 0]
        result = utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "genePred genscanPep"
        columns = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']
        expected_result = ['chrom', 'txStart', 'txEnd', 0]
        result =  utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "rmsk"
        columns = ['bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName', 'genoStart', 'genoEnd', 'genoLeft', 'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'repLeft', 'id']
        expected_result = ['genoName', 'genoStart', 'genoEnd', 0]
        result =  utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "GFF"
        columns = ["start", "end", "chrom"]
        expected_result = ["chrom", "start", "end", -1]
        result =  utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

        file_type = "MAF"
        columns = ["start", "end", "chrom"]
        expected_result = ["chrom", "start", "end", 0]
        result =  utility.target_columns(columns, file_type)
        self.assertEqual(result, expected_result)

    def test_UCSC_target_columns(self):
        # Integration test of both target_columns and UCSC servers
        # Iterate through a genome in UCSC and make sure target_columns runs
        with utility.connect_SQL_db("genome-mysql.soe.ucsc.edu", "genome") as db:
            test_genomes = ["rn6"]
            for genome in test_genomes:  # iterate through all genomes being checked
                tracks_info = requests.get(url="https://api.genome.ucsc.edu/list/tracks?genome={};trackLeavesOnly=1".format(genome)).json()[genome]
                # unzip track details and names
                for track, track_info in tracks_info.items():
                    # sometimes table lists the track name instead of the key of the dict
                    if "table" in track_info:
                        track = track_info["table"]
                    # bigDataUrl then the file needs to be downloaded/parsed
                    if "bigDataUrl" not in track_info:
                        columns = list(pd.read_sql("Show columns from {}.{}".format(genome, track), con=db)["Field"])
                        if len(columns) > 2: # Some files just don't have any info on them
                            result = utility.target_columns(columns, track_info["type"])
                            self.assertTrue([result])
    ####################
    # Config Variables #
    ####################
    def test_CONFIG_VAR_MAX_INTERVAL_PER_INDEX(self):
        self.assertGreater(app.config["MAX_INTERVALS_PER_INDEX"], 0)
        genome = "rn6"
        tracks_info = requests.get(url="https://api.genome.ucsc.edu/list/tracks?genome={};trackLeavesOnly=1".format(genome)).json()[genome]
        # unzip track details and names
        for track, track_info in tracks_info.items():
            self.assertGreaterEqual(app.config["MAX_INTERVALS_PER_INDEX"], track_info["itemCount"])

    def test_cluster_data(self):
        testing_size = app.config["MAX_INTERVALS_PER_INDEX"]/2
        files_info = {"bannana": {
                    "file_size": testing_size, 
                    "download_function": "ba", 
                    "download_params": [1,1,1]
                    },
                 "mellon": {
                    "file_size": testing_size * 1.1, 
                    "download_function": "", 
                    "download_params": [2,2,2]
                    },

                 "apple": {
                    "file_size": testing_size, 
                    "download_function": "", 
                    "download_params": [3,3,3]
                }}
        genome = "fruits"
        source = "kitchen"
        clusters = utility.cluster_data(source, genome, files_info)
        expected_clusters = {source + "_" + genome+".1": {
                                "files": {
                                    "bannana": {
                                        "file_size": testing_size, 
                                        "download_function": "ba", 
                                        "download_params": [1,1,1]
                                        },
                                    "mellon": {
                                        "file_size": testing_size * 1.1, 
                                        "download_function": "", 
                                        "download_params": [2,2,2]
                                    }},
                                "index_size": testing_size + 1.1*testing_size},
                             source + "_" + genome+".2": {
                                 "files": {
                                    "apple": {
                                        "file_size": testing_size, 
                                        "download_function": "", 
                                        "download_params": [3,3,3]
                                    }},
                                "index_size": testing_size}
                          }                     
        self.assertTrue(expected_clusters == clusters)


    ##################
    # Indexing Tests #
    ##################

    ##################
    # Updating Tests #
    ##################

    #################
    # Routing Tests #
    #################
    # def test_main_page(self):
    #     response = app.post('/', follow_redirects=True)
    #     self.assertEqual(response.status_code, 200)

if __name__ == "__main__":
    unittest.main()
